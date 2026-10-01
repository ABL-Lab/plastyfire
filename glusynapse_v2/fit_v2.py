"""Analytical fit on all protocols (tsk T11: --model post, the cooker rule; T17: --model v2).

    post  x = [a00, a01, a10, a11]                   apical tied to basal, gammas fixed (cooker)
    v2    x = [a00, a01, a10, a11, log10 A_mglu, log10 A_NO] + the filter parameters in --free-filters

Filter parameters change the dpre features (one pass over shaft_cai), so each free filter parameter is
quantised to --filter-steps levels in its box and the features are cached per worker by value.
Objective: sum over targets.py targets of ((pred - target)/SEM)^2, plus fit.py's rules 1-3.

    python glusynapse_v2/fit_v2.py --dirs glusynapse_v2/extracted/markram_delta-cooker \
        --groups markram --model post --workers 30 --maxiter 100 --save glusynapse_v2/results/fit_post_markram
"""
import argparse, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import batch_v2                        # noqa: E402
from batch_v2 import BatchV2          # noqa: E402
from targets import load_targets      # noqa: E402
import model_v2 as MV                 # noqa: E402

NAMES_A = ["a00", "a01", "a10", "a11"]
MIN_ACTIVE, MIN_POT = 0.30, 0.15      # fit.py rules 2 and 3
FILTER_BOX = dict(theta_T=(5e-5, 2e-3), tau_T=(2.0, 500.0), theta_NO=(5e-4, 8e-3),
                  tau_NO=(2.0, 50.0), tau_Z=(5.0, 100.0), ca_scale=(1e-4, 1e-2), theta_Z=(0.0, 3.0),
                  # pre_drive 1 (spine VDCC influx, nA): per-synapse AP peaks span 1e-8..2.5e-3 (median 3e-5,
                  # ebner_preview), so the boxes cover the whole range
                  theta_Ti=(1e-8, 1e-3), theta_NOi=(1e-8, 1e-3), i_scale=(1e-7, 1e-3),
                  # no_drive 1 (effcai, mM): rise after a pre spike is 0.03-0.04 for single pairings and
                  # 0.06-0.09 for bursts / 50 Hz (logs_no_drive.py)
                  theta_NOe=(0.005, 0.5), e_scale=(1e-3, 0.1),
                  # no_drive 2 (spine Ca above rest, mM): AP peaks 0.7-3e-3 (logs_cacr_drive.py); theta_N is the
                  # threshold on N (charge above theta_NOc / c_scale, ~2-20 for c_scale 1e-3)
                  theta_NOc=(1e-5, 3e-3), c_scale=(1e-4, 1e-2), theta_N=(0.0, 100.0),   # N per bAP ~14 at i_scale 1e-5
                  # t_drive 1 (T25, eCB): effcai / c_post, 1 at a single-bAP peak; dpre_min is the LTD floor
                  # (L5-L5 tLTD saturates near 0.65-0.73, Sjostrom 2001/2003)
                  # effcai / c_post before a pre spike: 0.07 at +10 ms, 0.65 at -10 ms (0.1 Hz), 30-400 inside trains
                  theta_Te=(0.03, 300.0), dpre_min=(-0.8, -0.05),
                  # t_drive 2 (post-AP trace S, +1 per AP, tau_E1 100 ms): theta_Te on S (0.5-1.5: single AP vs
                  # burst), theta_Tg the gate threshold on T (ms-scale: T ~ integral of pos(S - theta_Te))
                  theta_Tg=(0.0, 20.0))
# time constants are log-spaced too (v2.2 preview put tau_T on its old 5 ms bound; 2 ms is the 0.25 ms-grid limit)
LOG_FILTERS = {"theta_T", "theta_NO", "ca_scale", "theta_Ti", "theta_NOi", "i_scale", "tau_T", "tau_NO", "tau_Z",
               "theta_NOe", "e_scale", "theta_NOc", "c_scale", "theta_Te"}
_B = _T = _CFG = None


def unpack(x, cfg):
    a = dict(zip(NAMES_A, x[:4]))
    a.update(a20=a["a00"], a21=a["a01"], a30=a["a10"], a31=a["a11"])
    P = dict(cfg["filters"])
    if cfg["model"] == "v2":
        P["A_mglu"], P["A_NO"] = 10 ** x[4], 10 ** x[5]
        for k, xi in zip(cfg["free_filters"], x[6:]):
            lo, hi = FILTER_BOX[k]; n = cfg["filter_steps"]
            q = xi if cfg.get("continuous") else np.round(xi * (n - 1)) / (n - 1)   # [0, 1] -> (quantised) level
            P[k] = float(lo * (hi / lo) ** q if k in LOG_FILTERS else lo + q * (hi - lo))
    return a, tie(P, cfg)


def pack(a, P, cfg):
    """inverse of unpack (continuous): a fit json's a and pre -> DE vector (to seed a population)."""
    x = [a[k] for k in NAMES_A]
    P = {**MV.DEFAULTS, **cfg["filters"], **P}     # older fits lack newer filters (theta_N ...)
    if cfg["model"] == "v2":
        x += [np.log10(P["A_mglu"]), np.log10(P["A_NO"])]
        for k in cfg["free_filters"]:
            lo, hi = FILTER_BOX[k]
            x.append(np.log(P[k] / lo) / np.log(hi / lo) if k in LOG_FILTERS else (P[k] - lo) / (hi - lo))
    return np.array(x, float)


def tie(P, cfg):
    """--tie {"theta_NOi": "theta_Ti"}: copy one filter parameter onto another (shared threshold)."""
    for k, src in cfg.get("tie", {}).items():
        if src in P:
            P[k] = P[src]
    return P


def rules(a):
    td = np.concatenate([_B.thetas(r, a)[0] for r in _B.recs])
    tp = np.concatenate([_B.thetas(r, a)[1] for r in _B.recs])
    viol = np.maximum(td - tp, 0.0)
    if viol.max() > 0:
        return 1e4 * (1.0 + viol.mean() / max(np.median(tp), 1e-9))
    frac = float((td < _CFG["peak"]).mean())
    if frac < MIN_ACTIVE:
        return 1e3 * (1.0 + (MIN_ACTIVE - frac))
    fpot = float((tp < _CFG["peak"]).mean())
    if fpot < MIN_POT:
        return 1e3 * (1.0 + (MIN_POT - fpot))
    return 0.0


def objective(x):
    a, P = unpack(x, _CFG)
    pen = rules(a)
    if pen:
        return pen
    # RULE 4 (v2): the NO pathway is causal (pre before the Ca event) only if the presynaptic trace
    # outlasts the NO signal; with tau_Z < tau_NO post-then-pre potentiates (tests/test_units.py).
    if _CFG["model"] == "v2":
        tz, tn = P.get("tau_Z", MV.DEFAULTS["tau_Z"]), P.get("tau_NO", MV.DEFAULTS["tau_NO"])
        if tz <= tn:
            return 1e3 * (1.0 + (tn - tz) / tn)
    chi2, _ = _B.objective(a, P, _T, mode=_CFG["rho"], loc=_CFG["loc"])
    return chi2 if np.isfinite(chi2) else 1e5


def penalty(a, P):
    """Rules 1-3 (rules) and rule 4 (tau_Z > tau_NO) as in objective(); 0 if the candidate is admissible."""
    pen = rules(a)
    if pen:
        return pen
    if _CFG["model"] == "v2":
        tz, tn = P.get("tau_Z", MV.DEFAULTS["tau_Z"]), P.get("tau_NO", MV.DEFAULTS["tau_NO"])
        if tz <= tn:
            return 1e3 * (1.0 + (tn - tz) / tn)
    return 0.0


_G = None
_POP = {}          # timing of the last GPU objective call (for the checkpoint line)


def objective_gpu(X):
    """Vectorised objective for DE (vectorized=True): X (D, S) -> (S,). Admissible candidates go to the GPU
    in one batch."""
    xs = X.T if X.ndim == 2 else X[None]
    out = np.zeros(len(xs)); A, Ps, idx = [], [], []
    for i, x in enumerate(xs):
        a, P = unpack(x, _CFG)
        pen = penalty(a, P)
        if pen:
            out[i] = pen
        else:
            A.append(a); Ps.append({**_CFG["filters"], **P}); idx.append(i)
    if idx:
        t0 = time.time()
        c = _G.evaluate(A, Ps)
        out[idx] = np.where(np.isfinite(c), c, 1e5)
        _POP["sec"] = time.time() - t0
    return out if X.ndim == 2 else out[0]


def main():
    global _B, _T, _CFG
    ap = argparse.ArgumentParser()
    ap.add_argument("--dirs", required=True, help="comma-separated npz dirs")
    ap.add_argument("--groups", default="markram,nevian,ebner")
    ap.add_argument("--model", choices=["post", "v2"], default="post")
    ap.add_argument("--rho", choices=["fast", "full"], default="full",
                    help="full: exact ODE (default; validate_v2 shows fast mean-field is off by up to 0.15 in ratio)")
    ap.add_argument("--filters", default="{}", help="json: fixed filter parameters (model_v2.DEFAULTS keys)")
    ap.add_argument("--free-filters", default="", help="comma list from " + ",".join(FILTER_BOX))
    ap.add_argument("--filter-steps", type=int, default=8)
    ap.add_argument("--tie", default="{}", help='json {"dst": "src"}: dst filter = src filter, e.g. {"theta_NOi": "theta_Ti"}')
    ap.add_argument("--conditions", default="control,mglu_block,post_nmdar,nmdar_block")
    ap.add_argument("--loc-csv", default=None, help="csv pair,loc (Letzkus proximal/distal)")
    ap.add_argument("--a-lo", type=float, default=0.0); ap.add_argument("--a-hi", type=float, default=5.0)
    ap.add_argument("--x0", default=None, help="json of a-params to seed the population (e.g. fit_results/delta-cooker.json)")
    ap.add_argument("--workers", type=int, default=30)
    ap.add_argument("--maxiter", type=int, default=100)
    ap.add_argument("--popsize", type=int, default=8)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--pairs", default=None)
    ap.add_argument("--save", required=True, help="output prefix (.json, .csv)")
    ap.add_argument("--backend", choices=["cpu", "jax"], default="cpu",
                    help="jax: GPU objective (jax_v2), continuous filters, whole population per call")
    ap.add_argument("--seed-fits", default=None, help="comma list of fit jsons (same model and free filters) put "
                    "into the initial population after --x0, e.g. optima from subset runs")
    ap.add_argument("--resume", default=None, help="jax only: <save>_ckpt.npz of an earlier run: start DE from its "
                    "population (same bounds and free filters)")
    ap.add_argument("--chunk-gb", type=float, default=None,
                    help="jax only: stream traces to the GPU in chunks of this size (record sets larger than GPU "
                         "memory, e.g. all 120 pairs); host traces are freed and the table comes from the GPU")
    args = ap.parse_args()

    conds = set(args.conditions.split(","))
    _T = {k: v for k, v in load_targets(tuple(args.groups.split(","))).items() if k[1] in conds}
    protos = sorted({k[0].split("@")[0] for k in _T})
    if not _T:
        sys.exit(f"no targets for groups {args.groups!r} and conditions {sorted(conds)}")
    fil = json.loads(args.filters)
    sig = () if args.model == "post" else (("vdcc",) if fil.get("pre_drive") else ("shaft_cai",))
    if fil.get("no_drive") == 2 and args.backend == "cpu":
        sig += ("cacr",)                       # the GPU scan recovers spine Ca from effcai itself
    _B = BatchV2(args.dirs.split(","), protocols=protos, pairs=set(args.pairs.split(",")) if args.pairs else None,
                 fast=args.rho == "fast", signals=sig)
    loc = None
    if args.loc_csv:
        import pandas as pd
        loc = dict(pd.read_csv(args.loc_csv)[["pair", "loc"]].values)
    have = {r["proto"] for r in _B.recs}
    _T = {k: v for k, v in _T.items() if k[0].split("@")[0] in have}
    print(f"targets with records: {len(_T)}   protocols: {sorted(have)}", flush=True)
    peak = np.concatenate([r["peak"] for r in _B.recs])
    if args.rho == "full":
        _B._stack()          # build before forking so workers share it copy-on-write
    free = [k for k in args.free_filters.split(",") if k]
    _CFG = dict(model=args.model, rho=args.rho, filters=json.loads(args.filters), free_filters=free,
                filter_steps=args.filter_steps, loc=loc, peak=peak, tie=json.loads(args.tie))

    if args.backend == "jax":
        global _G
        import jax_v2
        _CFG["continuous"] = True
        t0 = time.time()
        _G = jax_v2.GPUModel(_B, _T, _CFG["filters"], chunk_gb=args.chunk_gb)
        if args.chunk_gb:                        # traces now live in _G's chunks only
            for r in _B.recs:
                for k in ("effcai", "vdcc", "shaft_cai", "cacr", "t"):
                    r.pop(k, None)
            _B._E = _B._H = None; import gc; gc.collect()
        print(f"GPUModel: {_G.L} lanes, buckets (steps, lanes) {[(b['T'], len(b['idx'])) for b in _G.buckets]}, "
              f"padded/real work {_G.pad_work:.2f}, built in {(time.time()-t0)/60:.1f} min", flush=True)
    elif args.model == "v2" and free:
        import itertools
        levels = {k: [unpack(np.r_[np.zeros(6), [q if kk == k else 0 for kk in free]], _CFG)[1][k]
                      for q in np.linspace(0, 1, args.filter_steps)] for k in free}
        grid = [dict(zip(free, v)) for v in itertools.product(*levels.values())]
        t0 = time.time()
        n = _B.precompute([tie({**_CFG["filters"], **g}, _CFG) for g in grid], workers=args.workers)
        print(f"precomputed {n} feature passes for {len(grid)} filter sets in {(time.time()-t0)/60:.1f} min",
              flush=True)

    bounds = [(args.a_lo, args.a_hi)] * 4
    if args.model == "v2":
        bounds += [(-4.0, 0.0), (-1.0, 4.0)] + [(0.0, 1.0)] * len(free)
    init = "latinhypercube"
    if args.x0:
        a0 = json.load(open(args.x0)); a0 = a0.get("params", a0)
        rng = np.random.default_rng(args.seed)
        n = args.popsize * len(bounds)
        lo, hi = np.array(bounds).T
        pop = lo + rng.random((n, len(bounds))) * (hi - lo)
        pop[0, :4] = [a0[k] for k in NAMES_A]
        for i, f in enumerate(args.seed_fits.split(",") if args.seed_fits else []):
            fj = json.load(open(f)); x = pack(fj["a"], fj["pre"], _CFG)
            assert np.allclose(unpack(x, _CFG)[1]["A_NO"], fj["pre"]["A_NO"]), f
            pop[1 + i] = np.clip(x, lo, hi)
            print(f"seed member {1 + i}: {f} (chi2 {fj['chi2']:.2f})", flush=True)
        init = pop

    os.makedirs(os.path.dirname(os.path.abspath(args.save)), exist_ok=True)
    if args.resume:
        init = np.load(args.resume)["x"]
        print(f"resume from {args.resume}: {len(init)} members", flush=True)
    from scipy.optimize import differential_evolution

    def checkpoint(intermediate_result):
        """per generation: the DE population (for --resume) and the objective call time"""
        ir = intermediate_result
        np.savez(args.save + "_ckpt.npz", x=ir.population, f=ir.population_energies, best=ir.x, nit=ir.nit)
        print(f"  ckpt nit {ir.nit}: best {ir.fun:.3f}, objective call {_POP.get('sec', 0):.0f}s", flush=True)
    from multiprocessing import Pool
    t0 = time.time()
    if args.backend == "jax":
        r = differential_evolution(objective_gpu, bounds, maxiter=args.maxiter, popsize=args.popsize, init=init,
                                   polish=False, seed=args.seed, disp=True, vectorized=True,
                                   updating="deferred", tol=1e-6, callback=checkpoint)
    else:
        with Pool(args.workers) as pool:          # fork: workers share _B, _T, _CFG
            r = differential_evolution(objective, bounds, maxiter=args.maxiter, popsize=args.popsize, init=init,
                                       polish=False, seed=args.seed, disp=True, workers=pool.map,
                                       updating="deferred", tol=1e-6)
    a, P = unpack(r.x, _CFG)
    if args.chunk_gb:
        import pandas as pd
        c, preds = _G.chi2(*_G.rho_dpre([a], [{**_CFG["filters"], **P}]), return_pred=True)
        chi2 = float(c[0])
        table = pd.DataFrame([dict(target=k[0], condition=k[1], target_mean=_T[k][0], target_sem=_T[k][1],
                                   pred=float(v[0]), z=(float(v[0]) - _T[k][0]) / _T[k][1]) for k, v in preds.items()])
    else:
        chi2, table = _B.objective(a, {**_CFG["filters"], **P} if args.model == "v2" else {}, _T, mode="full" if args.rho == "full" else "fast", loc=loc)
    table.to_csv(args.save + ".csv", index=False)
    json.dump(dict(model=args.model, a=a, pre=P, chi2=chi2, de_fun=float(r.fun), nfev=int(r.nfev),
                   n_targets=len(table), minutes=(time.time() - t0) / 60, args=vars(args),
                   bap_gate=batch_v2.BAP_GATE),
              open(args.save + ".json", "w"), indent=1)
    print(table.round(3).to_string(index=False))
    print(f"chi2 {chi2:.2f} over {len(table)} targets -> {args.save}.json")


if __name__ == "__main__":
    main()
