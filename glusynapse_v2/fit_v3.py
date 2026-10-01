"""fit_v2.py --backend jax with the v3 GPU objective (gpu_v3.GPUModelV3). Same DE, bounds, rules, seeding,
checkpoints (--resume takes fit_v2 checkpoints and vice versa) and output files; see FIT_SPEED.md.

Differences from fit_v2 (speed only): the objective is gpu_v3 (one CUDA kernel, traces resident on the GPU);
rules 1-3 use td, tp computed once for the whole population (same values as fit_v2.rules); host traces are
freed after the GPU copy, so the final table comes from the GPU model (as fit_v2 --chunk-gb). --chunk-gb is
accepted and ignored.
"""
import argparse, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fit_v2 as F                    # noqa: E402  (unpack, pack, tie, FILTER_BOX, rule constants)
import batch_v2                       # noqa: E402
from batch_v2 import BatchV2          # noqa: E402
from targets import load_targets      # noqa: E402
import model_v2 as MV                 # noqa: E402

_G = _CFG = None
_POP = {}


def rules_row(td, tp, peak):
    """fit_v2.rules on one candidate's concatenated td, tp (same operations, same values)."""
    viol = np.maximum(td - tp, 0.0)
    if viol.max() > 0:
        return 1e4 * (1.0 + viol.mean() / max(np.median(tp), 1e-9))
    frac = float((td < peak).mean())
    if frac < F.MIN_ACTIVE:
        return 1e3 * (1.0 + (F.MIN_ACTIVE - frac))
    fpot = float((tp < peak).mean())
    if fpot < F.MIN_POT:
        return 1e3 * (1.0 + (F.MIN_POT - fpot))
    return 0.0


def objective_gpu(X):
    """As fit_v2.objective_gpu: X (D, S) -> (S,); admissible candidates go to the GPU in one batch."""
    xs = X.T if X.ndim == 2 else X[None]
    out = np.zeros(len(xs))
    AP = [F.unpack(x, _CFG) for x in xs]
    td, tp = _G.thetas_all([a for a, _ in AP])
    Ps, idx = [], []
    for i, (a, P) in enumerate(AP):
        pen = rules_row(td[i], tp[i], _CFG["peak"])
        if not pen:
            tz, tn = P.get("tau_Z", MV.DEFAULTS["tau_Z"]), P.get("tau_NO", MV.DEFAULTS["tau_NO"])
            if tz <= tn:
                pen = 1e3 * (1.0 + (tn - tz) / tn)
        if pen:
            out[i] = pen
        else:
            Ps.append({**_CFG["filters"], **P}); idx.append(i)
    if idx:
        t0 = time.time()
        c = _G.evaluate_th(td[idx], tp[idx], Ps)
        out[idx] = np.where(np.isfinite(c), c, 1e5)
        _POP["sec"] = time.time() - t0; _POP["n"] = len(idx)
    return out if X.ndim == 2 else out[0]


def main():
    global _G, _CFG
    ap = argparse.ArgumentParser()
    ap.add_argument("--dirs", required=True)
    ap.add_argument("--groups", default="markram,nevian,ebner")
    ap.add_argument("--model", choices=["v2"], default="v2")
    ap.add_argument("--filters", default="{}")
    ap.add_argument("--free-filters", default="")
    ap.add_argument("--tie", default="{}")
    ap.add_argument("--conditions", default="control,mglu_block,post_nmdar,nmdar_block")
    ap.add_argument("--a-lo", type=float, default=0.0); ap.add_argument("--a-hi", type=float, default=5.0)
    ap.add_argument("--x0", default=None)
    ap.add_argument("--maxiter", type=int, default=100)
    ap.add_argument("--popsize", type=int, default=8)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--pairs", default=None)
    ap.add_argument("--save", required=True)
    ap.add_argument("--seed-fits", default=None)
    ap.add_argument("--resume", default=None)
    ap.add_argument("--chunk-gb", type=float, default=None, help="ignored (v2 streaming option)")
    ap.add_argument("--backend", default="jax", help="ignored")
    args = ap.parse_args()

    conds = set(args.conditions.split(","))
    T = {k: v for k, v in load_targets(tuple(args.groups.split(","))).items() if k[1] in conds}
    protos = sorted({k[0].split("@")[0] for k in T})
    if not T:
        sys.exit(f"no targets for groups {args.groups!r} and conditions {sorted(conds)}")
    fil = json.loads(args.filters)
    sig = ("vdcc",) if fil.get("pre_drive") else ("shaft_cai",)
    B = BatchV2(args.dirs.split(","), protocols=protos, pairs=set(args.pairs.split(",")) if args.pairs else None,
                fast=False, signals=sig)
    have = {r["proto"] for r in B.recs}
    T = {k: v for k, v in T.items() if k[0].split("@")[0] in have}
    print(f"targets with records: {len(T)}   protocols: {sorted(have)}", flush=True)
    peak = np.concatenate([r["peak"] for r in B.recs])
    B._stack()
    free = [k for k in args.free_filters.split(",") if k]
    _CFG = dict(model="v2", rho="full", filters=fil, free_filters=free, filter_steps=8, loc=None, peak=peak,
                tie=json.loads(args.tie), continuous=True)
    F._B, F._T, F._CFG = B, T, _CFG           # fit_v2.unpack / pack read _CFG only through their argument

    import gpu_v3
    _G = gpu_v3.GPUModelV3(B, T, fil)
    for r in B.recs:                           # traces now live on the GPU only
        for k in ("effcai", "vdcc", "shaft_cai", "cacr", "t", "cev", "cev_lo", "cev_hi", "vev", "_c4", "_v3", "arr"):
            r.pop(k, None)
    B._E = B._H = None; import gc; gc.collect()

    bounds = [(args.a_lo, args.a_hi)] * 4 + [(-4.0, 0.0), (-1.0, 4.0)] + [(0.0, 1.0)] * len(free)
    init = "latinhypercube"
    if args.x0:
        a0 = json.load(open(args.x0)); a0 = a0.get("params", a0)
        rng = np.random.default_rng(args.seed)
        n = args.popsize * len(bounds)
        lo, hi = np.array(bounds).T
        pop = lo + rng.random((n, len(bounds))) * (hi - lo)
        pop[0, :4] = [a0[k] for k in F.NAMES_A]
        for i, f in enumerate(args.seed_fits.split(",") if args.seed_fits else []):
            fj = json.load(open(f)); x = F.pack(fj["a"], fj["pre"], _CFG)
            assert np.allclose(F.unpack(x, _CFG)[1]["A_NO"], fj["pre"]["A_NO"]), f
            pop[1 + i] = np.clip(x, lo, hi)
            print(f"seed member {1 + i}: {f} (chi2 {fj['chi2']:.2f})", flush=True)
        init = pop
    os.makedirs(os.path.dirname(os.path.abspath(args.save)), exist_ok=True)
    if args.resume:
        init = np.load(args.resume)["x"]
        print(f"resume from {args.resume}: {len(init)} members", flush=True)
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
    r = differential_evolution(objective_gpu, bounds, maxiter=args.maxiter, popsize=args.popsize, init=init,
                               polish=False, seed=args.seed, disp=True, vectorized=True,
                               updating="deferred", tol=1e-6, callback=checkpoint)
    a, P = F.unpack(r.x, _CFG)
    import pandas as pd
    c, preds = _G.chi2(*_G.rho_dpre([a], [{**_CFG["filters"], **P}]), return_pred=True)
    chi2 = float(c[0])
    table = pd.DataFrame([dict(target=k[0], condition=k[1], target_mean=T[k][0], target_sem=T[k][1],
                               pred=float(v[0]), z=(float(v[0]) - T[k][0]) / T[k][1]) for k, v in preds.items()])
    table.to_csv(args.save + ".csv", index=False)
    json.dump(dict(model="v2", a=a, pre=P, chi2=chi2, de_fun=float(r.fun), nfev=int(r.nfev),
                   n_targets=len(table), minutes=(time.time() - t0) / 60, args=vars(args),
                   bap_gate=batch_v2.BAP_GATE, objective="gpu_v3"),
              open(args.save + ".json", "w"), indent=1)
    print(table.round(3).to_string(index=False))
    print(f"chi2 {chi2:.2f} over {len(table)} targets -> {args.save}.json")


if __name__ == "__main__":
    main()
