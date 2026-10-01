"""v5 fitter (V5_DESIGN.md): fit_v4n.py (N pathway models, --extra) with the v5 shared-pool rule of gpu_v5_rho.
fit_v4n.py / fit_v4.py are not modified.

v5 = Chindemi 2022 + one own-spine VDCC-Ca pool V (tau_E1, i_scale) with one threshold theta_V that both gates
potentiation and lets an own pre spike drive eCB-LTD (d -> d_min). No NO arm, rho_gamma 1. Set by the filters
(run_fit_v5.sh): v5_mode 1, vamp_mode 1, rho_gamma 1, A_eCB, dpre_min, tau_E1, i_scale. v5_mode 2 = v5b ((1 - b)-weighted
eCB trigger pool, threshold theta_V); + --free-filters theta_eCB = v5c (separate eCB threshold, log box 0.05-50);
v5_mode 3 = v5-E2 (v5b with only VDCC excursions that start in an unprimed pool, V <= theta_V, feeding the trigger).
DE vector: a00, a01, a10, a11 + --free-filters (theta_V) + gamma_d, gamma_p (--fit-gamma). The v4 A_mglu / A_NO slots
are not in the vector: they are held at 0 exactly (log slot -inf, fit_v2.unpack gives 10**-inf = 0.0), so the v4
eCB / NO chain is off in every v4 evaluation made here.

--check-v4: before the DE, at the first seed (a's, gammas, theta_V of the seed; rho_gamma 1, A_mglu = A_NO = 0), for
every pathway model assert (1e-9 in chi2)
  M1: v5 kernel (A_eCB 0)               == v4 vamp kernel (vamp_mode 1, gpu_v4_rho)
  M0: v5 kernel (A_eCB 0, theta_V 0)    == v4 rates kernel (vamp_mode 0) == v4 vamp kernel at theta_V 0
and print the v4 M0 chi2 next to results/v4_LM0chk.json (ladder check 22135018) when that file exists.
With --maxiter 0 the output is the seed scored with the full v5 rule (the filters).

--extra NAME:GROUPS:DIRS:BASIS[:GEOM[:PAIRS]] (repeatable) and --joint as in fit_v4n.py (one uniform parameter set;
pathway differences only via each synapse's own traces and basis). --drop-targets "target|cond" (all models) or
"NAME/target|cond". Outputs: <save>.csv (l5), <save>_<NAME>.csv, json chi2_<NAME>, n_targets_<NAME>, chi2_total.
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
import gpu_v5_rho                     # noqa: E402

F.FILTER_BOX.update(rho_gamma=(0.0, 1.0), tau_fast=(5.0, 300.0), rho_sigma=(0.005, 0.3),   # runtime only
                    gamma_d=(20.0, 250.0), gamma_p=(100.0, 600.0),    # GAMMA_LIT.md, as fit_v4n
                    theta_V=(0.0, 50.0),                              # shared pool threshold (flagged fitted)
                    theta_eCB=(0.05, 50.0))                           # v5c separate eCB threshold on W (log box)
F.LOG_FILTERS |= {"tau_fast", "rho_sigma", "theta_eCB"}
for k, v in {**rho_v4.V4_DEFAULTS, **rho_v4.RATE_DEFAULTS, **rho_v4.VAMP_DEFAULTS, **gpu_v5_rho.V5_DEFAULTS}.items():
    MV.DEFAULTS.setdefault(k, v)

GEOM = os.environ.get("GEOM_L23") or os.path.join(ROOT, "ebner", "pair_geometry_L23PC_L5TTPC.csv")   # GEOM_L23: alternative geometry csv
L5_BASIS = batch_v2.BASIS_DIR
L23_BASIS = os.environ.get("L23_BASIS_DIR") or os.path.join(ROOT, "basis_results_edges_ebner_l23l5_delta_rs")   # L23_BASIS_DIR: other emodel (run_new_emodel_chain.sh)
LM0CHK = os.path.join(HERE, "results", "v4_LM0chk.json")
_M = []            # [(model, name)]; _M[0] is L5 (rules use its peaks)
_CFG = None
_POP = {}


def expand(x):
    """DE vector (no A_mglu / A_NO) -> fit_v2 vector with both log slots at -inf (amplitudes exactly 0)."""
    return np.concatenate([x[:4], [-np.inf, -np.inf], x[4:]])


def unpack5(x):
    a, P = F.unpack(expand(np.asarray(x, float)), _CFG)
    return a, P


def build(dirs, groups, conds, fil, pairs, basis, fast, loc=None, vgate=False, rates=True, drop=()):
    batch_v2.BASIS_DIR = basis
    T = {k: v for k, v in load_targets(tuple(groups.split(","))).items() if k[1] in conds and f"{k[0]}|{k[1]}" not in drop}
    protos = sorted({k[0].split("@")[0] for k in T})
    sig = ("vdcc",) if fil.get("pre_drive") else ("shaft_cai",)
    B = BatchV2(dirs.split(","), protocols=protos, pairs=pairs, fast=False, signals=sig)
    have = {r["proto"] for r in B.recs}
    T = {k: v for k, v in T.items() if k[0].split("@")[0] in have}
    peak = np.concatenate([r["peak"] for r in B.recs])
    B._stack()
    G = gpu_v5_rho.GPUModelV5(B, T, fil, fast=fast, loc=loc, vgate=vgate, rates=rates)
    for r in B.recs:
        for k in ("effcai", "vdcc", "shaft_cai", "cacr", "t", "cev", "cev_lo", "cev_hi", "vev", "_c4", "_v3", "arr"):
            r.pop(k, None)
    B._E = B._H = None
    import gc; gc.collect()
    print(f"{groups}: {len(T)} targets, {len(B.recs)} records (v5 kernel {G.v5})", flush=True)
    return G, T, peak


def objective_gpu(X):
    xs = X.T if X.ndim == 2 else X[None]
    out = np.zeros(len(xs))
    AP = [unpack5(x) for x in xs]
    Ps_all = [{**_CFG["filters"], **P} for _, P in AP]
    G0 = _M[0][0]
    td, tp = G0.thetas_all([a for a, _ in AP], Ps_all)
    idx = []
    for i, (a, P) in enumerate(AP):
        pen = F3.rules_row(td[i], tp[i], _CFG["peak"])
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


WEIGHTS = {}       # env WEIGHTS / --weights: {"ltd": w, "ltp": w, "target|cond": w, "NAME/target|cond": w}; empty = unweighted


def target_weight(name, k, mean):
    """z^2 weight of one target: ltd (exp mean < 1.0) / ltp class weight, overridden by an explicit key."""
    w = WEIGHTS.get("ltd", 1.0) if mean < 1.0 else WEIGHTS.get("ltp", 1.0)
    for key in (f"{k[0]}|{k[1]}", f"{name}/{k[0]}|{k[1]}"):
        w = WEIGHTS.get(key, w)
    return float(w)


def install_weights(G, T, name):
    """Instance-level chi2 wrapper: weighted sum of z^2 for the DE objective (class chi2 untouched; only called when
    WEIGHTS is set, so the default path is bit-identical)."""
    orig = type(G).chi2; W = {k: target_weight(name, k, T[k][0]) for k in T}

    def chi2w(rho, d, return_pred=False, Ps=None):
        _, preds = orig(G, rho, d, True, Ps)
        out = np.zeros(rho.shape[0])
        for k, v in preds.items():
            z = (v - T[k][0]) / T[k][1]
            out += W[k] * np.where(np.isfinite(z), z ** 2, 0.0)
        return (out, preds) if return_pred else out
    G.chi2 = chi2w


def table(G, T, a, P, name="l5"):
    """Unweighted chi2 (class chi2, also when a weighting wrapper is installed) + per-target table."""
    import pandas as pd
    c, preds = type(G).chi2(G, *G.rho_dpre([a], [P]), return_pred=True, Ps=[P])
    df = pd.DataFrame([dict(target=k[0], condition=k[1], target_mean=T[k][0], target_sem=T[k][1],
                            pred=float(v[0]), z=(float(v[0]) - T[k][0]) / T[k][1]) for k, v in preds.items()])
    if WEIGHTS:
        df["ltd"] = df.target_mean < 1.0
        df["w"] = [target_weight(name, (t, c_), m) for t, c_, m in zip(df.target, df.condition, df.target_mean)]
    return float(c[0]), df


def check_v4(seed, fil):
    """M1 / M0 limits of the v5 kernel against the v4 kernels at the seed (see module doc). Raises on a mismatch."""
    import gpu_v4_rho
    f, fj, x = seed
    a, P = unpack5(x); P = {**fil, **P}
    assert P["A_mglu"] == 0.0 and P["A_NO"] == 0.0, (P["A_mglu"], P["A_NO"])
    P1 = {**P, "A_eCB": 0.0}                      # M1: C1 gate on the pool, no eCB, no NO
    P0 = {**P1, "theta_V": 0.0}                   # M0: pure Chindemi (fitted rates)
    tot = {"m1_v4": 0.0, "m1_v5": 0.0, "m0_v4": 0.0, "m0_v5": 0.0, "v5": 0.0}
    for G, name in _M:
        assert G.v5 and G.vamp == 1 and G.rates, (name, G.v5, G.vamp, G.rates)
        kv = G._kern
        kr = gpu_v4_rho._make_kernel_fast(G.t_mode, G.no_mode, G.use_s, False, False, True, 0)   # v4 rates, vamp 0

        def ev(v5, Pq, kern=None):
            G.v5 = v5
            if kern is not None:
                G._kern = kern
            try:
                return float(G.evaluate([a], [Pq])[0])
            finally:
                G._kern = kv; G.v5 = True
        c = dict(m1_v4=ev(False, P1), m1_v5=ev(True, P1), m0_v4=ev(False, P0, kr), m0_v4g=ev(False, P0),
                 m0_v5=ev(True, P0), v5=ev(True, P))
        for k in tot:
            tot[k] += c[k]
        print(f"check-v4 {name}: M1 v4 {c['m1_v4']:.12f}  v5 {c['m1_v5']:.12f}  diff {abs(c['m1_v4'] - c['m1_v5']):.1e} | "
              f"M0 v4 {c['m0_v4']:.12f}  v4-vamp(theta_V 0) {c['m0_v4g']:.12f}  v5 {c['m0_v5']:.12f}  "
              f"diff {abs(c['m0_v4'] - c['m0_v5']):.1e} | v5 full rule at seed {c['v5']:.6f}", flush=True)
        assert abs(c["m1_v4"] - c["m1_v5"]) < 1e-9, (name, "M1", c)
        assert abs(c["m0_v4"] - c["m0_v5"]) < 1e-9 and abs(c["m0_v4"] - c["m0_v4g"]) < 1e-9, (name, "M0", c)
    print(f"check-v4 total over {len(_M)} models: M1 v4 {tot['m1_v4']:.9f} v5 {tot['m1_v5']:.9f}; "
          f"M0 v4 {tot['m0_v4']:.9f} v5 {tot['m0_v5']:.9f}; v5 full {tot['v5']:.6f}  (seed {os.path.basename(f)}: "
          f"json chi2_total {fj.get('chi2_total', fj['chi2'] + fj.get('chi2_l23', 0.0)):.6f})", flush=True)
    if os.path.exists(LM0CHK):
        lj = json.load(open(LM0CHK))
        ref = lj.get("chi2", np.nan) + lj.get("chi2_l23", 0.0)
        print(f"check-v4 vs ladder {os.path.basename(LM0CHK)}: json total {ref:.9f}  v4 M0 here {tot['m0_v4']:.9f}  "
              f"diff {abs(ref - tot['m0_v4']):.1e} (equal only if the ladder used these a's / gammas)", flush=True)
    else:
        print(f"check-v4: {LM0CHK} not there yet (ladder check pending), in-job v4 kernels used", flush=True)
    print("CHECK-V4 OK (1e-9)", flush=True)
    return tot


def pd_read(f):
    import pandas as pd
    return pd.read_csv(f)


def main():
    global _CFG
    ap = argparse.ArgumentParser()
    ap.add_argument("--dirs", required=True)
    ap.add_argument("--groups", default="paired_l5,sjostrom07")
    ap.add_argument("--filters", default="{}")
    ap.add_argument("--set", default="{}", help="options fixed for every candidate (merged into filters)")
    ap.add_argument("--free-filters", default="theta_V")
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
    ap.add_argument("--joint", action="store_true")
    ap.add_argument("--l23-dirs", default=os.path.join(V2, "extracted", "ebner_l23l5_delta-prefire-vseg-rs"))
    ap.add_argument("--check-v4", action="store_true")
    ap.add_argument("--drop-targets", default="", help="target|condition,... left out of the objective (validation only)")
    ap.add_argument("--fit-gamma", action="store_true", help="free gamma_d, gamma_p (rho rates, uniform)")
    ap.add_argument("--weights", default=os.environ.get("WEIGHTS", ""),
                    help='json {"ltd": w, "ltp": w, "target|cond": w}: z^2 weights in the DE objective only (env WEIGHTS)')
    ap.add_argument("--extra", action="append", default=[], help="NAME:GROUPS:DIRS:BASIS[:GEOM[:PAIRS]], repeatable")
    args = ap.parse_args()

    conds = set(args.conditions.split(","))
    fil = {**json.loads(args.filters), **json.loads(args.set)}
    free = [k for k in args.free_filters.split(",") if k]
    if args.fit_gamma:
        free += [k for k in ("gamma_d", "gamma_p") if k not in free]
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
                         L5_BASIS, False, drop=drop_for("l5"))
    _M.append((G5, "l5"))
    TM = {"l5": T5}
    for name, groups, dirs, basis, geom, pf in specs:
        loc = pairs = None
        if geom:
            import pandas as pd
            g = pd.read_csv(geom); g["pair"] = g.pregid.astype(str) + "-" + g.postgid.astype(str); g = g[g.all_protocols]
            lc = [c for c in g.columns if c.endswith("_distal")]
            loc = {p: {c: bool(v) for c, v in zip(lc, row)} for p, row in zip(g.pair, g[lc].itertuples(index=False))}
            pairs = set(g.pair)
        if pf:
            pp = set(open(pf).read().strip().split(","))
            pairs = pp if pairs is None else pairs & pp
        G, T, _ = build(dirs, groups, conds, fil, pairs, basis, False, loc=loc, drop=drop_for(name))
        _M.append((G, name)); TM[name] = T
        print(f"model {name}: groups {groups}, basis {basis}, geom {geom or '-'}, {len(T)} targets", flush=True)
    if args.weights:
        WEIGHTS.update(json.loads(args.weights))
        for G, name in _M:
            install_weights(G, TM[name], name)
        for G, name in _M:
            ws = [target_weight(name, k, TM[name][k][0]) for k in TM[name]]
            print(f"weights {name}: {WEIGHTS}; {sum(w != 1.0 for w in ws)} of {len(ws)} targets weighted", flush=True)
    _CFG = dict(model="v2", rho="full", filters=fil, free_filters=free, filter_steps=8, loc=None, peak=peak,
                tie={}, continuous=True)
    F._CFG = _CFG

    bounds = [(args.a_lo, args.a_hi)] * 4 + [(0.0, 1.0)] * len(free)
    lo, hi = np.array(bounds).T
    print(f"v5 free ({len(bounds)}): {F.NAMES_A + free}; fixed filters {fil}", flush=True)
    seeds = []
    for f in (args.seed_fits.split(",") if args.seed_fits else []):
        fj = json.load(open(f)); pre = {**json.loads(args.seed_set), **fj["pre"]}
        with np.errstate(divide="ignore"):
            xf = F.pack(fj["a"], pre, _CFG)
        seeds.append((f, fj, np.clip(np.delete(xf, [4, 5]), lo, hi)))
    if args.check_v4:
        assert seeds, "--check-v4 needs --seed-fits"
        check_v4(seeds[0], fil)
    init = "latinhypercube"
    if args.x0 or seeds:
        rng = np.random.default_rng(args.seed)
        pop = lo + rng.random((args.popsize * len(bounds), len(bounds))) * (hi - lo)
        if args.x0:
            a0 = json.load(open(args.x0)); a0 = a0.get("params", a0); pop[0, :4] = [a0[k] for k in F.NAMES_A]
        for i, (f, fj, x) in enumerate(seeds):
            pop[1 + i] = x; print(f"seed member {1 + i}: {f} (chi2 {fj['chi2']:.2f})", flush=True)
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
        xb, fun, nfev = r.x, float(r.fun), int(r.nfev)
    else:
        xb, fun, nfev = seeds[0][2], float("nan"), 0
    a, P = unpack5(xb); Pf = {**fil, **P}
    chi2, tab = table(G5, T5, a, Pf, "l5")
    tab.to_csv(args.save + ".csv", index=False)
    out = dict(model="v2", a=a, pre=P, chi2=chi2, de_fun=fun, nfev=nfev, n_targets=len(tab), minutes=(time.time() - t0) / 60,
               args=vars(args), bap_gate=batch_v2.BAP_GATE, objective="gpu_v5_rho", n_free=len(bounds),
               v5=dict(v5_mode=int(Pf.get("v5_mode", 0)), theta_V=float(Pf["theta_V"]), A_eCB=float(Pf.get("A_eCB", 0.0)),
                       theta_eCB=float(Pf.get("theta_eCB", 0.0)), tau_d_NMDA=float(Pf.get("tau_d_NMDA", MV.DEFAULTS["tau_d_NMDA"])),
                       dpre_min=float(Pf["dpre_min"]), tau_E1=float(Pf["tau_E1"]), rho_gamma=rho_v4.opts(Pf)[0]),
               gamma_d=rho_v4.rates(Pf)[0], gamma_p=rho_v4.rates(Pf)[1], fit_gamma=int(bool(args.fit_gamma)))
    print(tab.round(3).to_string(index=False))
    print(f"chi2 L5 {chi2:.2f} over {len(tab)} targets (gamma_d {out['gamma_d']:.3f}, gamma_p {out['gamma_p']:.3f}, "
          f"v5 {out['v5']})")
    tot = chi2
    for G, name in _M[1:]:
        cm, tm = table(G, TM[name], a, Pf, name); tm.to_csv(f"{args.save}_{name}.csv", index=False)
        out[f"chi2_{name}"] = cm; out[f"n_targets_{name}"] = len(tm); tot += cm
        print(tm.round(3).to_string(index=False)); print(f"chi2 {name} {cm:.2f} over {len(tm)} targets")
    out["chi2_total"] = tot
    out["models"] = [dict(name="l5", groups=args.groups, dirs=args.dirs, basis=L5_BASIS)] + [
        dict(name=n, groups=g, dirs=d, basis=b, geom=ge, pairs=p) for n, g, d, b, ge, p in specs]
    k = len(bounds); n = sum(out.get(f"n_targets_{m}", len(tab)) if m != "l5" else len(tab) for _, m in _M)
    out["aic"] = tot + 2 * k; out["bic"] = tot + k * np.log(n)
    if WEIGHTS:      # unweighted chi2 above; split LTD / LTP and weighted objective value of the final point
        tabs = [(tab, "l5")] + [(pd_read(f"{args.save}_{m}.csv"), m) for _, m in _M[1:]]
        z2 = lambda t, msk: float((t.z[msk] ** 2).sum())
        out["weights"] = WEIGHTS
        out["chi2_ltd"] = sum(z2(t, t.ltd) for t, _ in tabs); out["chi2_ltp"] = sum(z2(t, ~t.ltd) for t, _ in tabs)
        out["n_ltd"] = int(sum(t.ltd.sum() for t, _ in tabs)); out["n_ltp"] = int(sum((~t.ltd).sum() for t, _ in tabs))
        out["chi2_weighted"] = float(sum((t.w * t.z ** 2).sum() for t, _ in tabs))
        print(f"UNWEIGHTED chi2 total {tot:.6f} = LTD {out['chi2_ltd']:.3f} ({out['n_ltd']}) + LTP {out['chi2_ltp']:.3f} "
              f"({out['n_ltp']}); weighted objective {out['chi2_weighted']:.3f} (de_fun {fun:.3f})", flush=True)
    print(f"chi2 total {tot:.6f} over {len(_M)} models, {n} targets, k {k}: AIC {out['aic']:.2f}, BIC {out['bic']:.2f}")
    if args.maxiter <= 0 and seeds and seeds[0][1].get("objective") == "gpu_v5_rho":    # repro of a v5 json
        fj = seeds[0][1]; ok = True
        for _, name in _M:
            kk = "chi2" if name == "l5" else f"chi2_{name}"
            if kk in fj:
                d = abs(out[kk] - fj[kk]); ok &= d < 1e-9
                print(f"repro {name}: json {fj[kk]:.12f}  now {out[kk]:.12f}  diff {d:.1e}", flush=True)
        print("REPRO OK (1e-9)" if ok else "REPRO DIFF", flush=True)
    json.dump(out, open(args.save + ".json", "w"), indent=1)
    print(f"-> {args.save}.json")


if __name__ == "__main__":
    main()
