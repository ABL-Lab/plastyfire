"""fit_v4 generalised to N pathway models (PROTOCOLS.md). fit_v4.py is not modified; with --joint and no --extra this
file gives the same objective and the same chi2 as fit_v4.py --joint.

--extra NAME:GROUPS:DIRS:BASIS[:GEOM[:PAIRS]] (repeatable) adds one pathway model, scored with the SAME parameter
vector (one uniform set; pathway differences only via each synapse's own traces and basis). GROUPS / DIRS are comma
lists (load_targets groups, extracted dirs); BASIS the pathway's basis_results_edges_* dir (relative to the
plastyfire root or absolute). GEOM: optional pair-geometry csv (pregid, postgid, all_protocols and boolean *_distal
location columns -> loc map, pairs = all_protocols rows). PAIRS: optional file with a comma list of pairs (intersected
with GEOM pairs). --joint is shorthand for --extra l23:paired_l23l5:<--l23-dirs>:<L23_BASIS>:<GEOM>.
The first model (l5, --dirs/--groups) supplies the rules_row peaks. Each model gets its own targets, so pathway groups
that share protocol ids are separate models. --drop-targets entries "target|cond" apply to every model,
"NAME/target|cond" to model NAME only. Outputs: <save>.csv (l5), <save>_<NAME>.csv, json chi2_<NAME>,
n_targets_<NAME>, chi2_total, models. With --maxiter 0 the seed's per-model chi2 are compared with its json.

Original fit_v4 notes:
fit_v3 with the v4 post-rule candidates (rho_v4 / gpu_v4_rho; RHO_REDESIGN.md). fit_v3.py is not modified.

New free filters (add them to --free-filters): rho_gamma (A, box 0-1), tau_fast (B, 5-300 ms, log; needs --fast),
rho_sigma (C, 0.005-0.3, log). --set '{"rho_gamma": 0.0}' fixes an option for every candidate (goes into filters).
--seed-set '{"tau_fast": 278.3}' fills a missing option in the seed fits (a v3 json has none).
--joint: also score paired_l23l5 on the L2/3 -> L5 dirs (own basis dir, loc selections) and minimise the sum
(diagnostic only: can one uniform parameter set fit both pathways).
--vgate: (D) pot also needs G > 1/2, G = own unweighted spine VDCC events (cev) through tau_E1 (rho_v4d.py).
--check-v3: before the DE, print chi2 of the seed fits under gpu_v3 and under v4 (must agree when v4 is off).
--fit-gamma: the rho rates gamma_d (20-250), gamma_p (100-600) join the DE vector (one value each, uniform); seeds
without them start at GAMMA_D_GB, GAMMA_P_GB. With --check-v3 the seed chi2 with the rates kernel at the GB values
must equal the constant-rate kernel (asserted to 1e-9).
ROUND2 C1 / C2 (scan_vgate_amp.py): --set {"vamp_mode": 1 or 2, "theta_V": x} fixes the VDCC-amplitude gate on
potentiation (V = own -ica_VDCC through tau_E1, / i_scale; pot needs V > theta_V; C1 a failed crossing depresses, C2 it
is neutral); --free-filters theta_V frees it (0-50, seeds without it start at 0 = A0). vamp_mode 0 = off (current kernel).
--fix-params '{"a01": 1.17, "A_NO": 0}': as fit_v4.py (fixed values leave the DE vector; a00..a11 and A_mglu / A_NO in
linear units fill their slots, 0 -> log10 -inf -> exactly 0; a filter is dropped from the free list and set in the
filters). Default {} = unchanged. With --maxiter 0 the free slots come from the first seed fit.
"""
import argparse, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); V2 = os.path.dirname(HERE); ROOT = os.path.dirname(V2)
sys.path.insert(0, V2); sys.path.insert(0, HERE)
import fit_v2 as F                    # noqa: E402
import fit_v3 as F3                   # noqa: E402  (rules_row)
import batch_v2                       # noqa: E402
from batch_v2 import BatchV2          # noqa: E402
from targets import load_targets      # noqa: E402
import model_v2 as MV                 # noqa: E402
import rho_v4                         # noqa: E402

F.FILTER_BOX.update(rho_gamma=(0.0, 1.0), tau_fast=(5.0, 300.0), rho_sigma=(0.005, 0.3),   # runtime only
                    gamma_d=(20.0, 250.0), gamma_p=(100.0, 600.0),    # widened 2026-10-01 per GAMMA_LIT.md (was 50-200 / 150-300)
                    theta_V=(0.0, 50.0))                              # ROUND2 C1 / C2 VDCC-amplitude gate
F.LOG_FILTERS |= {"tau_fast", "rho_sigma"}
for k, v in {**rho_v4.V4_DEFAULTS, **rho_v4.RATE_DEFAULTS, **rho_v4.VAMP_DEFAULTS}.items():   # so fit_v2.pack finds them for old seeds
    MV.DEFAULTS.setdefault(k, v)

GEOM = os.environ.get("GEOM_L23") or os.path.join(ROOT, "ebner", "pair_geometry_L23PC_L5TTPC.csv")   # GEOM_L23: alternative geometry csv
L5_BASIS = batch_v2.BASIS_DIR
L23_BASIS = os.environ.get("L23_BASIS_DIR") or os.path.join(ROOT, "basis_results_edges_ebner_l23l5_delta_rs")   # L23_BASIS_DIR: other emodel (run_new_emodel_chain.sh)
_M = []            # [(model, name)]; _M[0] is L5 (rules use its peaks)
_CFG = None
_POP = {}
_FIX = None        # --fix-params: (keep indices into the full vector, {index: value}, full length); None = no fixing


def expand(X):
    """reduced DE vector(s) (n_free,) or (n_free, S) -> full vector(s) with the fixed slots filled (fit_v4.expand)."""
    if _FIX is None:
        return X
    keep, fixv, n = _FIX
    Xf = np.empty((n,) + X.shape[1:])
    Xf[keep] = X
    for i, v in fixv.items():
        Xf[i] = v
    return Xf


def build(dirs, groups, conds, fil, pairs, basis, fast, loc=None, vgate=False, rates=False, drop=()):
    import gpu_v4_rho
    batch_v2.BASIS_DIR = basis
    T = {k: v for k, v in load_targets(tuple(groups.split(","))).items() if k[1] in conds and f"{k[0]}|{k[1]}" not in drop}
    protos = sorted({k[0].split("@")[0] for k in T})
    sig = ("vdcc",) if fil.get("pre_drive") else ("shaft_cai",)
    B = BatchV2(dirs.split(","), protocols=protos, pairs=pairs, fast=False, signals=sig)
    have = {r["proto"] for r in B.recs}
    T = {k: v for k, v in T.items() if k[0].split("@")[0] in have}
    peak = np.concatenate([r["peak"] for r in B.recs])
    B._stack()
    G = gpu_v4_rho.GPUModelV4(B, T, fil, fast=fast, loc=loc, vgate=vgate, rates=rates)
    for r in B.recs:
        for k in ("effcai", "vdcc", "shaft_cai", "cacr", "t", "cev", "cev_lo", "cev_hi", "vev", "_c4", "_v3", "arr"):
            r.pop(k, None)
    B._E = B._H = None
    import gc; gc.collect()
    print(f"{groups}: {len(T)} targets, {len(B.recs)} records", flush=True)
    return G, T, peak


def objective_gpu(X):
    X = expand(np.asarray(X))
    xs = X.T if X.ndim == 2 else X[None]
    out = np.zeros(len(xs))
    AP = [F.unpack(x, _CFG) for x in xs]
    Ps_all = [{**_CFG["filters"], **P} for _, P in AP]
    G0 = _M[0][0]
    td, tp = G0.thetas_all([a for a, _ in AP], Ps_all)
    idx = []
    for i, (a, P) in enumerate(AP):
        pen = F3.rules_row(td[i], tp[i], _CFG["peak"])
        if not pen:
            tz, tn = P.get("tau_Z", MV.DEFAULTS["tau_Z"]), P.get("tau_NO", MV.DEFAULTS["tau_NO"])
            if tz <= tn:
                pen = 1e3 * (1.0 + (tn - tz) / tn)
        if pen:
            out[i] = pen
        else:
            idx.append(i)
    if idx:
        t0 = time.time()
        Ps = [Ps_all[i] for i in idx]; A = [AP[i][0] for i in idx]
        c = G0.evaluate_th(td[idx], tp[idx], Ps)
        for G, _ in _M[1:]:
            c = c + G.evaluate(A, Ps)
        out[idx] = np.where(np.isfinite(c), c, 1e5)
        _POP["sec"] = time.time() - t0; _POP["n"] = len(idx)
    return out if X.ndim == 2 else out[0]


def table(G, T, a, P):
    import pandas as pd
    c, preds = G.chi2(*G.rho_dpre([a], [P]), return_pred=True, Ps=[P])
    return float(c[0]), pd.DataFrame([dict(target=k[0], condition=k[1], target_mean=T[k][0], target_sem=T[k][1],
                                           pred=float(v[0]), z=(float(v[0]) - T[k][0]) / T[k][1]) for k, v in preds.items()])


def main():
    global _CFG
    ap = argparse.ArgumentParser()
    ap.add_argument("--dirs", required=True)
    ap.add_argument("--groups", default="paired_l5,sjostrom07")
    ap.add_argument("--filters", default="{}")
    ap.add_argument("--set", default="{}", help="v4 options fixed for every candidate (merged into filters)")
    ap.add_argument("--free-filters", default="")
    ap.add_argument("--tie", default="{}")
    ap.add_argument("--conditions", default="control,mglu_block,post_nmdar,nmdar_block,no_block")
    ap.add_argument("--a-lo", type=float, default=0.0); ap.add_argument("--a-hi", type=float, default=5.0)
    ap.add_argument("--x0", default=None)
    ap.add_argument("--maxiter", type=int, default=150)
    ap.add_argument("--popsize", type=int, default=8)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--pairs", default=None)
    ap.add_argument("--save", required=True)
    ap.add_argument("--seed-fits", default=None)
    ap.add_argument("--seed-set", default="{}")
    ap.add_argument("--resume", default=None)
    ap.add_argument("--fast", action="store_true", help="compile the tau_fast kernel (B)")
    ap.add_argument("--vgate", action="store_true", help="(D) potentiation needs a recent own spine VDCC event")
    ap.add_argument("--joint", action="store_true")
    ap.add_argument("--l23-dirs", default=os.path.join(V2, "extracted", "ebner_l23l5_delta-prefire-vseg-rs"))
    ap.add_argument("--check-v3", action="store_true")
    ap.add_argument("--drop-targets", default="", help="target|condition,... left out of the objective (validation only)")
    ap.add_argument("--fit-gamma", action="store_true", help="free gamma_d, gamma_p (rho rates, uniform)")
    ap.add_argument("--extra", action="append", default=[], help="NAME:GROUPS:DIRS:BASIS[:GEOM[:PAIRS]], repeatable")
    ap.add_argument("--fix-params", default="{}", help="json {name: value}: removed from the DE vector, held fixed")
    args = ap.parse_args()
    global _FIX

    conds = set(args.conditions.split(","))
    fil = {**json.loads(args.filters), **json.loads(args.set)}
    if args.vgate:
        fil["vgate"] = 1
    fix = {k: float(v) for k, v in json.loads(args.fix_params).items()}
    head = F.NAMES_A + ["A_mglu", "A_NO"]
    free = [k for k in args.free_filters.split(",") if k]
    if args.fit_gamma:
        free += [k for k in ("gamma_d", "gamma_p") if k not in free]
    bad = [k for k in fix if k not in head and k not in free]
    assert not bad, f"--fix-params: {bad} not free in this fit (a00..a11, A_mglu, A_NO or a free filter)"
    fil.update({k: v for k, v in fix.items() if k not in head})
    free = [k for k in free if k not in fix]
    rates = args.fit_gamma or "gamma_d" in fil or "gamma_p" in fil
    specs = []                                   # (name, groups, dirs, basis, geom, pairs_file)
    if args.joint:
        specs.append(("l23", "paired_l23l5", args.l23_dirs, L23_BASIS, GEOM, ""))
    for s in args.extra:
        f = s.split(":")
        if len(f) < 4 or len(f) > 6:
            raise SystemExit(f"--extra {s!r}: need NAME:GROUPS:DIRS:BASIS[:GEOM[:PAIRS]]")
        f += [""] * (6 - len(f))
        if not os.path.isabs(f[3]) and not os.path.isdir(f[3]):
            f[3] = os.path.join(ROOT, f[3])
        specs.append(tuple(f))
    names = ["l5"] + [s[0] for s in specs]
    if len(set(names)) != len(names):
        raise SystemExit(f"model names must be unique: {names}")
    drops = [t for t in args.drop_targets.split(",") if t]

    def drop_for(name):
        return set(t.split("/", 1)[1] if "/" in t else t for t in drops if "/" not in t or t.split("/", 1)[0] == name)

    G5, T5, peak = build(args.dirs, args.groups, conds, fil, set(args.pairs.split(",")) if args.pairs else None,
                         L5_BASIS, args.fast, vgate=args.vgate, rates=rates, drop=drop_for("l5"))
    _M.append((G5, "l5"))
    TM = {"l5": T5}
    for name, groups, dirs, basis, geom, pf in specs:
        loc = pairs = None
        if geom:
            import pandas as pd
            g = pd.read_csv(geom); g["pair"] = g.pregid.astype(str) + "-" + g.postgid.astype(str); g = g[g.all_protocols]
            lc = [c for c in g.columns if c.endswith("_distal")]       # pair_geometry_L23PC_L5TTPC: letzkus_distal, sh_distal
            loc = {p: {c: bool(v) for c, v in zip(lc, row)} for p, row in zip(g.pair, g[lc].itertuples(index=False))}
            pairs = set(g.pair)
        if pf:
            pp = set(open(pf).read().strip().split(","))
            pairs = pp if pairs is None else pairs & pp
        G, T, _ = build(dirs, groups, conds, fil, pairs, basis, args.fast, loc=loc, vgate=args.vgate,
                        drop=drop_for(name), rates=rates)
        _M.append((G, name)); TM[name] = T
        print(f"model {name}: groups {groups}, basis {basis}, geom {geom or '-'}, {len(T)} targets", flush=True)
    _CFG = dict(model="v2", rho="full", filters=fil, free_filters=free, filter_steps=8, loc=None, peak=peak,
                tie=json.loads(args.tie), continuous=True)
    F._CFG = _CFG

    bounds = [(args.a_lo, args.a_hi)] * 4 + [(-4.0, 0.0), (-1.0, 4.0)] + [(0.0, 1.0)] * len(free)
    lo, hi = np.array(bounds).T
    with np.errstate(divide="ignore"):       # A_mglu / A_NO = 0 -> log10 slot -inf -> 10 ** -inf = 0 exactly
        fixv = {i: (fix[k] if i < 4 else float(np.log10(fix[k]))) for i, k in enumerate(head) if k in fix}
    keep = [i for i in range(len(bounds)) if i not in fixv]
    if fixv:
        _FIX = (keep, fixv, len(bounds))
    if fix:
        print(f"fixed ({len(fix)}): {fix}; free ({len(keep)}): {[head[i] for i in keep if i < 6] + free}", flush=True)
    seeds = []
    for f in (args.seed_fits.split(",") if args.seed_fits else []):
        fj = json.load(open(f)); pre = {**json.loads(args.seed_set), **fj["pre"]}
        x = np.clip(F.pack(fj["a"], pre, _CFG), lo, hi)
        seeds.append((f, fj, expand(x[keep])))      # full vector; fixed slots hold the --fix-params values
    bounds = [bounds[i] for i in keep]; lo, hi = lo[keep], hi[keep]
    if args.check_v3:
        for f, fj, x in seeds:
            a, P = F.unpack(x, _CFG); P = {**fil, **P}
            c4 = float(G5.evaluate([a], [P])[0])
            print(f"check: {os.path.basename(f)} json chi2 {fj['chi2']:.6f}  v4 {c4:.6f}  (options {rho_v4.opts(P)}, "
                  f"rates {rho_v4.rates(P)}, vamp {rho_v4.vamp(P)})", flush=True)
            if rates:      # rates kernel at the exact GB values vs the constant-rate kernel (gpu_v3 when v4 is off)
                import gpu_v3, gpu_v4_rho
                Pg = {**P, **rho_v4.RATE_DEFAULTS}
                cg = float(G5.evaluate([a], [Pg])[0])
                kr = G5._kern; G5.rates = False
                G5._kern = (gpu_v4_rho._make_kernel_fast(G5.t_mode, G5.no_mode, G5.use_s, G5.fast, G5.vgate, False, G5.vamp)
                            if (G5.fast or G5.vgate or G5.vamp) else gpu_v3._make_kernel(G5.t_mode, G5.no_mode, G5.use_s))
                c0 = float(G5.evaluate([a], [Pg])[0])
                G5.rates = True; G5._kern = kr
                print(f"check gamma: rates kernel at GB {cg:.12f}  constant-rate kernel {c0:.12f}  diff {abs(cg - c0):.1e}",
                      flush=True)
                assert abs(cg - c0) < 1e-9, (cg, c0)
    init = "latinhypercube"
    if (args.x0 or seeds) and keep:
        rng = np.random.default_rng(args.seed)
        pop = lo + rng.random((args.popsize * len(bounds), len(bounds))) * (hi - lo)
        if args.x0:
            a0 = json.load(open(args.x0)); a0 = a0.get("params", a0)
            for j, i in enumerate(keep):
                if i < 4:
                    pop[0, j] = a0[F.NAMES_A[i]]
        for i, (f, fj, x) in enumerate(seeds):
            pop[1 + i] = x[keep]; print(f"seed member {1 + i}: {f} (chi2 {fj['chi2']:.2f})", flush=True)
        init = pop
    if args.resume:
        init = np.load(args.resume)["x"]; print(f"resume from {args.resume}: {len(init)} members", flush=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.save)), exist_ok=True)
    from scipy.optimize import differential_evolution
    tg = [time.time()]

    def checkpoint(intermediate_result):
        ir = intermediate_result
        np.savez(args.save + "_ckpt.npz", x=ir.population, f=ir.population_energies, best=ir.x, nit=ir.nit)
        now = time.time()
        print(f"  ckpt nit {ir.nit}: best {ir.fun:.3f}, objective call {_POP.get('sec', 0):.2f}s "
              f"({_POP.get('n', 0)} admissible), generation {now - tg[0]:.2f}s", flush=True)
        tg[0] = now
    t0 = time.time()
    if args.maxiter > 0:
        r = differential_evolution(objective_gpu, bounds, maxiter=args.maxiter, popsize=args.popsize, init=init,
                                   polish=False, seed=args.seed, disp=True, vectorized=True, updating="deferred",
                                   tol=1e-6, callback=checkpoint)
        xb, fun, nfev = expand(r.x), float(r.fun), int(r.nfev)
    else:
        assert seeds or not keep, "--maxiter 0 needs a seed fit unless every parameter is fixed"
        xb, fun, nfev = (seeds[0][2] if seeds else expand(np.zeros(0))), float("nan"), 0
    a, P = F.unpack(xb, _CFG); Pf = {**fil, **P}
    chi2, tab = table(G5, T5, a, Pf)
    tab.to_csv(args.save + ".csv", index=False)
    out = dict(model="v2", a=a, pre=P, chi2=chi2, de_fun=fun, nfev=nfev, n_targets=len(tab), minutes=(time.time() - t0) / 60,
               args=vars(args), bap_gate=batch_v2.BAP_GATE, objective="gpu_v4_rho", v4=dict(zip(("rho_gamma", "tau_fast",
               "rho_sigma"), rho_v4.opts(Pf)), vgate=int(bool(args.vgate)), vamp_mode=rho_v4.vamp(Pf)[0], theta_V=rho_v4.vamp(Pf)[1]), n_free=len(bounds),
               fixed=fix,
               gamma_d=rho_v4.rates(Pf)[0], gamma_p=rho_v4.rates(Pf)[1], fit_gamma=int(bool(args.fit_gamma)))
    print(tab.round(3).to_string(index=False))
    print(f"chi2 L5 {chi2:.2f} over {len(tab)} targets (gamma_d {out['gamma_d']:.3f}, gamma_p {out['gamma_p']:.3f}, "
          f"vamp {rho_v4.vamp(Pf)})")
    tot = chi2
    for G, name in _M[1:]:
        cm, tm = table(G, TM[name], a, Pf); tm.to_csv(f"{args.save}_{name}.csv", index=False)
        out[f"chi2_{name}"] = cm; out[f"n_targets_{name}"] = len(tm); tot += cm
        print(tm.round(3).to_string(index=False)); print(f"chi2 {name} {cm:.2f} over {len(tm)} targets")
    out["chi2_total"] = tot
    out["models"] = [dict(name="l5", groups=args.groups, dirs=args.dirs, basis=L5_BASIS)] + [
        dict(name=n, groups=g, dirs=d, basis=b, geom=ge, pairs=p) for n, g, d, b, ge, p in specs]
    print(f"chi2 total {tot:.6f} over {len(_M)} models")
    if args.maxiter <= 0 and seeds:          # reproduction check against the seed json
        fj = seeds[0][1]; ok = True
        for _, name in _M:
            k = "chi2" if name == "l5" else f"chi2_{name}"
            if k in fj:
                d = abs(out[k] - fj[k]); ok &= d < 1e-9
                print(f"repro {name}: json {fj[k]:.12f}  now {out[k]:.12f}  diff {d:.1e}", flush=True)
        print("REPRO OK (1e-9)" if ok else "REPRO DIFF", flush=True)
    json.dump(out, open(args.save + ".json", "w"), indent=1)
    print(f"-> {args.save}.json")


if __name__ == "__main__":
    main()
