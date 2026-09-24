"""Fit a-params (and optionally gamma_d/gamma_p) to Markram 1997 in-vitro STDP.

Parallelism
-----------
The Batch holds ~800 MB (effcai + its transpose). Sending that to 60 workers
would be 48 GB, so it is a MODULE-LEVEL global built before the pool forks:
on Linux fork shares it copy-on-write and the workers only read it. The
objective is therefore a top-level function, not a closure -- a closure over
`b` would be pickled and defeat the whole thing.

DE evaluates its whole population per generation with updating='deferred',
which is what makes `--workers` actually useful. CMA-ES uses ask/tell with an
explicit pool.

Search choice: the EPSP stage binarises rho at 0.5, so the objective is
piecewise-constant -- no gradient method applies. CMA-ES adapts step size and is
the primary; --de is an independent check (agreement between them is evidence a
minimum is real).

    python fit.py --time-eval --fit-gamma
    python fit.py --fit-gamma --de --workers 60 --maxiter 50 --verbose
    python fit.py --fit-gamma --workers 60 --maxiter 80          # CMA-ES

Methods (the "cpre" method is the 2-a-param one)
------------------------------------------------
default    theta = a*c_pre + b*c_post, 4 a-params (apical tied to basal).
--cpre-only  theta_d = a00*c_pre, theta_p = a10*c_pre. The c_post coefficients
           are pinned to 0, so the search is 2 a-params (+ gamma_d/gamma_p with
           --fit-gamma). Motivated by 22.1% of synapses sitting at the resting
           c_post floor, where the c_post term carries no information anyway.
           Note theta_p > theta_d then reduces exactly to a10 > a00.
"""
import argparse, multiprocessing as mp, os, time
import numpy as np
from constants import INVITRO, TIED_IC4, GAMMA_D_GB, GAMMA_P_GB
from batch import Batch

FIT_DT  = [-10, 5, 10]
NAMES_A = ["a00", "a01", "a10", "a11"]
NAMES_A_CPRE = ["a00", "a10"]      # --cpre-only: a01 = a11 = 0
MIN_ACTIVE = 0.30      # min fraction of synapses with theta_d below peak effcai
# Physiological bounds on the plasticity rates (user-supplied). Earlier fits
# used [10,400] / [10,800] and railed gamma_p at 483 then 739 -- outside range,
# so those results are invalid. gamma_p in particular is tight: the ic4 value
# (199.77) and chindemi (216.18) sit inside, the mod-file default (450) does not.
GAMMA_D_BOUNDS = (50.0, 200.0)
GAMMA_P_BOUNDS = (150.0, 300.0)

MIN_POT    = 0.15      # min fraction with theta_p below peak effcai, i.e. able
                       # to potentiate at all. Without this the search finds
                       # "depress everything": theta_p above the calcium so pot
                       # never fires, which fits -10 while killing +5/+10.
# Outside the STDP window there should be little plasticity. Anchoring the tails
# at ~1.0 blocks solutions that depress every dt uniformly. Weight is loose
# because these are not measured points, just the physiological expectation.
TAIL_DT     = [-50, -30, 30, 50]
TAIL_TARGET = 1.0
TAIL_SEM    = 0.15

_B   = None            # Batch, shared with workers via fork
_CFG = {"fit_gamma": False, "verbose": False,
        "anchor_tails": False, "sign": False, "cpre_only": False,
        "fit_dt": FIT_DT}


def _init_worker(cfg, src, limit, stride=1):
    """Pool initializer.

    Under 'fork' the child already inherited _B and this is a no-op guard.
    Under 'spawn' the child re-imports this module as __mp_main__, so the
    __main__ block never ran and _B is None -- rebuild it there. We force fork
    below precisely to avoid 60 workers each re-reading 381 MB.
    """
    global _B, _CFG
    _CFG.update(cfg)
    if _B is None:
        _B = Batch(limit=limit, extracted=src, stride=stride)


def unpack(x, fit_gamma, cpre_only=False):
    """apical tied to basal (a20=a00, a21=a01, a30=a10, a31=a11).

    cpre_only: x is [a00, a10(, gd, gp)] and the c_post coefficients are pinned
    to zero, so theta_d = a00*c_pre and theta_p = a10*c_pre.
    """
    if cpre_only:
        a00, a10 = x[:2]
        a01 = a11 = 0.0
        na = 2
    else:
        a00, a01, a10, a11 = x[:4]
        na = 4
    a = dict(a00=a00, a01=a01, a10=a10, a11=a11,
             a20=a00, a21=a01, a30=a10, a31=a11)
    gd, gp = (x[na], x[na + 1]) if fit_gamma else (GAMMA_D_GB, GAMMA_P_GB)
    return a, gd, gp


def objective(x):
    """Top-level so it pickles cheaply; reads the forked global _B."""
    b = _B
    a, gd, gp = unpack(x, _CFG["fit_gamma"], _CFG["cpre_only"])
    td, tp = b.thetas(a)

    # RULE 1: theta_p > theta_d for EVERY synapse. Coefficient tests (a10>a00)
    # are neither necessary nor sufficient -- theta = a*c_pre + b*c_post, so the
    # ordering depends on each synapse's own (c_pre, c_post). Graded so the
    # optimiser is pushed back rather than meeting an invisible cliff.
    viol = np.maximum(td - tp, 0.0)
    if viol.max() > 0:
        return 1e4 * (1.0 + viol.mean() / max(np.median(tp), 1e-9))

    # RULE 2: thresholds must sit inside the calcium that actually occurs.
    # theta_d above peak effcai -> nothing crosses, rho frozen, every dt = 1.0.
    # That flat basin is where an unconstrained search parks itself.
    frac = float((td < b.peak).mean())
    if frac < MIN_ACTIVE:
        return 1e3 * (1.0 + (MIN_ACTIVE - frac))

    # RULE 3: potentiation must be reachable for a real fraction of synapses.
    fpot = float((tp < b.peak).mean())
    if fpot < MIN_POT:
        return 1e3 * (1.0 + (MIN_POT - fpot))

    summ, _ = b.curve(b.rho_full(td, tp, gd, gp))
    m = dict(zip(summ.dt, summ["mean"]))
    err = 0.0
    for dt in _CFG["fit_dt"]:
        if dt not in m:
            return 1e3
        tgt, sem = INVITRO[dt]
        err += ((m[dt] - tgt) / sem) ** 2
    if _CFG["anchor_tails"]:
        for dt in TAIL_DT:
            if dt in m:
                err += ((m[dt] - TAIL_TARGET) / TAIL_SEM) ** 2
    if _CFG["sign"]:
        # the qualitative shape: depression below 1 at -10, potentiation above
        # 1 at +5/+10. Penalise the wrong side even when the magnitude is close.
        err += 50.0 * max(0.0, m[-10] - 1.0) ** 2 / 0.01
        for dt in (5, 10):
            err += 50.0 * max(0.0, 1.0 - m[dt]) ** 2 / 0.01
    if _CFG["verbose"]:
        print(f"  {np.array2string(x, precision=3)} -> "
              f"{'  '.join(f'{d:+d}:{m[d]:.3f}' for d in _CFG['fit_dt'])}  err={err:.2f}",
              flush=True)
    return err


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", type=int, default=0, help="limit records (0=all)")
    ap.add_argument("--stride", type=int, default=1,
                    help="decimate effcai on load. With extracted_d0p25 (0.25 ms): "
                         "4 -> 1 ms (fitting), 2 -> 0.5 ms (validation). See dt_scan.py")
    ap.add_argument("--workers", type=int, default=60)
    ap.add_argument("--maxiter", type=int, default=60)
    ap.add_argument("--sigma", type=float, default=0.3)
    ap.add_argument("--a-lo", type=float, default=1.0, help="lower bound, all a-params")
    ap.add_argument("--a-hi", type=float, default=5.0, help="upper bound, all a-params")
    ap.add_argument("--gd-lo", type=float, default=GAMMA_D_BOUNDS[0])
    ap.add_argument("--gd-hi", type=float, default=GAMMA_D_BOUNDS[1])
    ap.add_argument("--gp-lo", type=float, default=GAMMA_P_BOUNDS[0])
    ap.add_argument("--gp-hi", type=float, default=GAMMA_P_BOUNDS[1])
    ap.add_argument("--de", action="store_true")
    # DE convergence: scipy stops when
    #     std(pop_energies) <= atol + tol * |mean(pop_energies)|
    # Defaults are scipy's own (tol=0.01, atol=0). tol is the useful lever:
    # lower it to keep the population evolving. RAISING atol stops it SOONER.
    # Note the objective is piecewise-constant (rho is binarised at 0.5), so a
    # stalled population can be a genuine plateau rather than a loose tolerance
    # -- a smaller tol then just burns evaluations without improving f(x).
    # seed and popsize were hardcoded (1 and 8). The objective is
    # piecewise-constant, so DE routinely collapses the population onto a
    # plateau and exits early regardless of tol -- f(x) repeating unchanged for
    # many steps is that signature. A different seed explores a different
    # trajectory; a larger popsize covers the box better. Sweeping seeds is the
    # cheapest way to tell a real floor from a stagnation.
    ap.add_argument("--seed", type=int, default=1,
                    help="DE/CMA-ES seed (was hardcoded 1). Sweep it to test "
                         "whether a plateau is a real optimum.")
    ap.add_argument("--popsize", type=int, default=8,
                    help="DE popsize multiplier (was hardcoded 8); "
                         "individuals = popsize * n_free_params.")
    ap.add_argument("--tol", type=float, default=0.01,
                    help="DE relative convergence tolerance (scipy default "
                         "0.01). Lower => runs longer.")
    ap.add_argument("--atol", type=float, default=0.0,
                    help="DE absolute convergence tolerance (scipy default 0). "
                         "Higher => stops sooner.")
    ap.add_argument("--fit-gamma", action="store_true")
    ap.add_argument("--time-eval", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--anchor-tails", action="store_true",
                    help="also require dt=+-30,+-50 to sit near 1.0")
    ap.add_argument("--fit-dt", type=str, default="",
                    help="comma-separated dt values to score, e.g. '-10,10'. "
                         "Default (empty) keeps the standard -10,5,10. Only the "
                         "listed dt enter the objective; --anchor-tails is "
                         "independent of this.")
    ap.add_argument("--cpre-only", action="store_true",
                    help="fit theta from c_pre alone: a01=a11=0, so only "
                         "a00/a10 (+ gammas) are free")
    ap.add_argument("--sign", action="store_true",
                    help="penalise wrong-sign plasticity (-10 must be <1, +5/+10 >1)")
    args = ap.parse_args()

    # keep BLAS single-threaded: parallelism is across candidates, and nested
    # threads would oversubscribe 60 cores badly
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(v, "1")

    try:
        mp.set_start_method("fork")      # so workers inherit _B copy-on-write
    except RuntimeError:
        pass
    print(f"multiprocessing start method: {mp.get_start_method()}")

    _B = Batch(limit=args.pairs or None, stride=args.stride)
    _CFG["fit_gamma"] = args.fit_gamma
    _CFG["verbose"]      = args.verbose
    _CFG["anchor_tails"] = args.anchor_tails
    _CFG["sign"]         = args.sign
    _CFG["cpre_only"]    = args.cpre_only
    if args.fit_dt:
        _CFG["fit_dt"] = [int(v) for v in args.fit_dt.split(",") if v.strip()]
        missing = [d for d in _CFG["fit_dt"] if d not in INVITRO]
        if missing:
            raise SystemExit(f"--fit-dt {missing} has no in-vitro target; "
                             f"available: {sorted(INVITRO)}")
        print(f"scoring only dt = {_CFG['fit_dt']} "
              f"(default is {FIT_DT})")
    NAMES = NAMES_A_CPRE if args.cpre_only else NAMES_A
    NA    = len(NAMES)
    if args.cpre_only:
        print("method: cpre-only (theta_d = a00*c_pre, theta_p = a10*c_pre; "
              "a01 = a11 = 0)")
    print(f"{_B.N} synapses / {len(_B.meta)} records / T={_B.T} / "
          f"dt={_B.dt_ms:.3f} ms (stride {args.stride})")

    x0 = [TIED_IC4[k] for k in NAMES] + ([GAMMA_D_GB, GAMMA_P_GB] if args.fit_gamma else [])
    x0 = np.array(x0)
    td0, tp0 = _B.thetas(unpack(x0, args.fit_gamma, args.cpre_only)[0])
    print(f"seed: theta_d med {np.median(td0):.4f}  theta_p med {np.median(tp0):.4f}  "
          f"peak effcai med {np.median(_B.peak):.4f}")
    print(f"  synapses with theta_d < peak effcai: {100*(td0 < _B.peak).mean():.0f}%  "
          f"(need >={100*MIN_ACTIVE:.0f}%)")
    print(f"  synapses with theta_p < peak effcai: {100*(tp0 < _B.peak).mean():.0f}%  "
          f"(need >={100*MIN_POT:.0f}%, else no potentiation is possible)")

    t0 = time.time(); e0 = objective(x0); dt1 = time.time() - t0
    n_par = (NA + (2 if args.fit_gamma else 0)) * args.popsize
    print(f"seed err {e0:.2f}   ({dt1:.1f} s / eval, {args.workers} workers)")
    if args.time_eval:
        print(f"  DE: {args.maxiter} gens x {n_par} pop / {args.workers} cores "
              f"~ {args.maxiter*np.ceil(n_par/args.workers)*dt1/60:.0f} min")
        raise SystemExit

    ng = [args.gd_lo, args.gp_lo] if args.fit_gamma else []
    xg = [args.gd_hi, args.gp_hi] if args.fit_gamma else []
    lo = np.array([args.a_lo]*NA + ng)
    hi = np.array([args.a_hi]*NA + xg)
    sig0 = np.array([0.25*(args.a_hi-args.a_lo)]*NA
                    + ([0.25*(args.gd_hi-args.gd_lo),
                        0.25*(args.gp_hi-args.gp_lo)] if args.fit_gamma else []))

    # Feasibility: theta_d = a00*c_pre + a01*c_post must fall below peak effcai
    # for at least MIN_ACTIVE of synapses, else nothing crosses and the whole
    # box is the flat no-plasticity basin.
    cp_med, cq_med = np.median(_B.c_pre), np.median(_B.c_post)
    print(f"\n  a-params bounded to [{args.a_lo}, {args.a_hi}]")
    print(f"  median c_pre {cp_med:.4f}  c_post {cq_med:.5f}  peak effcai {np.median(_B.peak):.4f}")
    q_med = 0.0 if args.cpre_only else cq_med       # c_post term is pinned off
    for lab, av in (("at a_lo", args.a_lo), ("at a_hi", args.a_hi)):
        th = av*cp_med + av*q_med
        print(f"    theta {lab:8s} ~ {th:.4f}   "
              f"{'BELOW peak (active)' if th < np.median(_B.peak) else 'ABOVE peak (inert)'}")
    cq_term = 0.0 if args.cpre_only else args.a_lo*_B.c_post
    fa = (( args.a_lo*_B.c_pre + cq_term) < _B.peak).mean()
    print(f"  best case (both a at lower bound): {100*fa:.0f}% of synapses active"
          f"   -- need >={100*MIN_ACTIVE:.0f}%")
    if args.fit_gamma:
        print(f"  gamma_d in [{args.gd_lo}, {args.gd_hi}]   "
              f"gamma_p in [{args.gp_lo}, {args.gp_hi}]")
    if fa < MIN_ACTIVE:
        print("  WARNING: even the lower bound leaves too few synapses active; "
              "no feasible point exists in this box.")

    if args.de:
        from scipy.optimize import differential_evolution
        from multiprocessing import Pool
        cfg = dict(_CFG)
        with Pool(args.workers, initializer=_init_worker,
                  initargs=(cfg, _B.src, args.pairs or None, args.stride)) as pool:
            r = differential_evolution(objective, list(zip(lo, hi)),
                                       maxiter=args.maxiter, popsize=args.popsize,
                                       polish=False,
                                       seed=args.seed, disp=True,
                                       tol=args.tol, atol=args.atol,
                                       workers=pool.map, updating="deferred")
        xb, eb, nev = r.x, r.fun, r.nfev
    else:
        try:
            import cma
        except ImportError:
            raise SystemExit("cma not installed:  pip install --user cma   (or --de)")
        from multiprocessing import Pool
        cfg = dict(_CFG)
        es = cma.CMAEvolutionStrategy(x0, args.sigma,
                                      {"bounds": [lo.tolist(), hi.tolist()],
                                       "CMA_stds": (sig0/args.sigma).tolist(),
                                       "maxiter": args.maxiter, "seed": args.seed})
        nev = 0
        with Pool(args.workers, initializer=_init_worker,
                  initargs=(cfg, _B.src, args.pairs or None, args.stride)) as pool:
            while not es.stop():
                X = es.ask()
                es.tell(X, pool.map(objective, X))
                nev += len(X); es.disp()
        xb, eb = es.result.xbest, es.result.fbest

    a, gd, gp = unpack(xb, args.fit_gamma, args.cpre_only)
    print("\n=== best ===")
    for k in NAMES: print(f"  {k} = {a[k]:.6f}")
    if args.cpre_only:
        print("  a01 = a11 = 0.000000   (cpre-only method)")
    print("  (apical tied: a20=a00, a21=a01, a30=a10, a31=a11)")
    print(f"  gamma_d = {gd:.4f}   gamma_p = {gp:.4f}"
          f"{'' if args.fit_gamma else '   (fixed)'}")
    print(f"  weighted err {eb:.3f}   evaluations {nev}")
    if args.fit_gamma:
        for nm, val, blo, bhi in (("gamma_d", gd, args.gd_lo, args.gd_hi),
                                  ("gamma_p", gp, args.gp_lo, args.gp_hi)):
            frac = (val - blo) / (bhi - blo)
            if frac > 0.95 or frac < 0.05:
                print(f"  WARNING: {nm} = {val:.3f} is railed at "
                      f"{100*frac:.0f}% of [{blo}, {bhi}] -- the objective wants "
                      f"to leave the physiological range, so the fit is "
                      f"compensating for something else.")
    td, tp = _B.thetas(a)
    summ, _ = _B.curve(_B.rho_full(td, tp, gd, gp))
    print(f"\n  {'dt':>5} {'fitted':>9} {'in vitro':>9}")
    for dt in sorted(summ.dt):
        tgt = f"{INVITRO[dt][0]:.4f}" if dt in INVITRO else "-"
        print(f"  {dt:>5} {summ[summ.dt==dt]['mean'].values[0]:>9.4f} {tgt:>9}")
