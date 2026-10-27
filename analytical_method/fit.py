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
import argparse, glob, json, multiprocessing as mp, os, time
import numpy as np
from constants import INVITRO, TIED_IC4, GAMMA_D_GB, GAMMA_P_GB, TAU_EFFCA_GB
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
        "anchor_tails": False, "sign": False, "cpre_only": False, "untie": False,
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


NAMES_A8 = ["a00", "a01", "a10", "a11", "a20", "a21", "a30", "a31"]   # --untie-apical


def unpack(x, fit_gamma, cpre_only=False, untie=False):
    """apical tied to basal (a20=a00, a21=a01, a30=a10, a31=a11).

    cpre_only: x is [a00, a10(, gd, gp)] and the c_post coefficients are pinned
    to zero, so theta_d = a00*c_pre and theta_p = a10*c_pre.
    untie (OGK-03): x is the 8 a-params in NAMES_A8 order (apical free).
    """
    if untie:
        a = dict(zip(NAMES_A8, (float(v) for v in x[:8])))
        gd, gp = (x[8], x[9]) if fit_gamma else (GAMMA_D_GB, GAMMA_P_GB)
        return a, gd, gp
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


# ---------------------------------------------------------------------------
# OGK-03: protocol_id-keyed objective (behind --targets; the --fit-dt path above
# is untouched). One Batch per protocol (own T, npz set, basis dir); the
# objective is the sum over protocols of ((pred - target)/sem)^2. The three
# validity rules (theta_p>theta_d, MIN_ACTIVE, MIN_POT) are evaluated on the
# POOL of all synapses of all batches, which is exactly what the single
# --fit-dt Batch does, so the two paths agree to the last digit.
# ---------------------------------------------------------------------------
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROTO_CSV = os.path.join(ROOT, "biodata", "paired_recordings.csv")
_BS = []               # [dict(id, batch, tgt, sem, dt, tail)] -- shared with workers via fork
_POOL = {}             # pooled c_pre / c_post / peak for the startup printout


def load_targets(ids):
    """protocol_id -> dict(proto='<f>Hz_<dt>ms' npz tag, dt, tgt, sem, pre, post)."""
    import csv
    rows = {r["protocol_id"]: r for r in csv.DictReader(open(PROTO_CSV))}
    out = {}
    for pid in ids:
        if pid not in rows:
            raise SystemExit(f"--targets: {pid} not in {PROTO_CSV}")
        r = rows[pid]
        f, d = float(r["frequency_train"]), int(float(r["dt_train"]))
        out[pid] = dict(proto=f"{f:g}Hz_{d}ms", dt=d, tgt=float(r["mean_epsp_ratio"]),
                        sem=float(r["sem_epsp_ratio"]), pre=r["pre_mtype"], post=r["post_mtype"])
    return out


def build_proto_batches(ids, proto_map, anchor_tails, stride, limit):
    """proto_map: {protocol_id: {extracted, basis, suffix, proto}} (all keys optional;
    defaults = $ANALYTICAL_EXTRACTED / $ANALYTICAL_BASIS_DIR, suffix '', proto from the CSV)."""
    tg = load_targets(ids)
    d_ext = os.environ.get("ANALYTICAL_EXTRACTED")
    d_bas = os.environ.get("ANALYTICAL_BASIS_DIR")
    bs = []
    for pid in ids:
        pm = proto_map.get(pid, {})
        ext = pm.get("extracted", d_ext); bas = pm.get("basis", d_bas)
        if not ext: raise SystemExit(f"{pid}: no extracted dir (proto-map or ANALYTICAL_EXTRACTED)")
        proto = pm.get("proto", tg[pid]["proto"])
        b = Batch(limit=limit, extracted=ext, stride=stride, basis_dir=bas,
                  pattern=f"*__{proto}{pm.get('suffix', '')}.npz")
        bs.append(dict(id=pid, batch=b, tgt=tg[pid]["tgt"], sem=tg[pid]["sem"], w=float(pm.get("weight", 1.0)),   # opt-in per-protocol weight (proto-map "weight"; default 1 = unchanged)
                       dt=b.meta[0]["dt"], tail=False, proto=proto, extracted=ext, basis=bas))
    if anchor_tails:      # same semantics as the old path: a tail counts only if its npz exist
        for dt in TAIL_DT:
            proto = f"10Hz_{dt}ms"
            if not d_ext or not glob.glob(os.path.join(d_ext, f"*__{proto}.npz")):
                print(f"anchor-tails: no npz for {proto} in {d_ext} -> tail skipped", flush=True)
                continue
            b = Batch(limit=limit, extracted=d_ext, stride=stride, basis_dir=d_bas,
                      pattern=f"*__{proto}.npz")
            bs.append(dict(id=f"tail{dt:+d}", batch=b, tgt=TAIL_TARGET, sem=TAIL_SEM,
                           dt=dt, tail=True, proto=proto, extracted=d_ext, basis=d_bas))
    return bs


def load_apical_ids(path):
    """Apical synapse global ids: .json (list, or {gid: bool}), .npy/.npz ('syn'), or text (one id per line)."""
    if path.endswith(".json"):
        j = json.load(open(path))
        ids = [int(k) for k, v in j.items() if v] if isinstance(j, dict) else [int(v) for v in j]
    elif path.endswith((".npy", ".npz")):
        z = np.load(path)
        ids = (z["syn"] if path.endswith(".npz") else z).astype(np.int64).tolist()
    else:
        ids = [int(s) for s in open(path).read().split()]
    return set(ids)


class PooledBatch:
    """Read-only stand-in for the startup printout: the synapses of all protocol batches pooled."""
    def __init__(self, bs):
        B = [p["batch"] for p in bs]
        self.N = sum(b.N for b in B)
        self.meta = [m for b in B for m in b.meta]
        self.T = "/".join(str(b.T) for b in B)
        self.dt_ms = B[0].dt_ms
        self.c_pre = np.concatenate([b.c_pre for b in B])
        self.c_post = np.concatenate([b.c_post for b in B])
        self.peak = np.concatenate([b.peak for b in B])
        self.src = "; ".join(sorted({p["extracted"] for p in bs}))
        self._B = B

    def thetas(self, a):
        t = [b.thetas(a) for b in self._B]
        return np.concatenate([x[0] for x in t]), np.concatenate([x[1] for x in t])


def objective_proto(x):
    """Sum over protocols of ((pred-target)/sem)^2, rules on the pooled synapses."""
    a, gd, gp = unpack(x, _CFG["fit_gamma"], _CFG["cpre_only"], _CFG["untie"])
    tds, tps = [], []
    for p in _BS:
        td, tp = p["batch"].thetas(a); tds.append(td); tps.append(tp)
    td = np.concatenate(tds); tp = np.concatenate(tps)
    viol = np.maximum(td - tp, 0.0)
    if viol.max() > 0:
        return 1e4 * (1.0 + viol.mean() / max(np.median(tp), 1e-9))
    frac = float((td < _POOL["peak"]).mean())
    if frac < MIN_ACTIVE:
        return 1e3 * (1.0 + (MIN_ACTIVE - frac))
    fpot = float((tp < _POOL["peak"]).mean())
    if fpot < MIN_POT:
        return 1e3 * (1.0 + (MIN_POT - fpot))
    err = 0.0
    preds = {}
    for p, d_, t_ in zip(_BS, tds, tps):
        b = p["batch"]
        summ, _ = b.curve(b.rho_full(d_, t_, gd, gp))
        m = dict(zip(summ.dt, summ["mean"]))
        if p["dt"] not in m:
            if p["tail"]: continue          # same as the old path: a missing tail is skipped
            return 1e3
        preds[p["id"]] = m[p["dt"]]
        err += p.get("w", 1.0) * ((m[p["dt"]] - p["tgt"]) / p["sem"]) ** 2
    if _CFG["verbose"]:
        print(f"  {np.array2string(x, precision=3)} -> "
              f"{'  '.join(f'{k}:{v:.3f}' for k, v in preds.items())}  err={err:.2f}", flush=True)
    return err


def proto_predictions(x):
    """{protocol_id: predicted ratio} at x (final report / --eval-only)."""
    a, gd, gp = unpack(x, _CFG["fit_gamma"], _CFG["cpre_only"], _CFG["untie"])
    out = {}
    for p in _BS:
        b = p["batch"]; td, tp = b.thetas(a)
        summ, _ = b.curve(b.rho_full(td, tp, gd, gp))
        m = dict(zip(summ.dt, summ["mean"]))
        out[p["id"]] = float(m[p["dt"]]) if p["dt"] in m else None
    return out


def seeds_from_json(paths, names, fit_gamma, lo, hi):
    """Each fit json (its 'params' block) -> one x vector in the active layout, clipped to bounds."""
    out = []
    for pth in paths:
        pr = json.load(open(pth))["params"]
        v = [pr[k] for k in names]
        if fit_gamma:
            v += [pr["gamma_d_GB_GluSynapse"], pr["gamma_p_GB_GluSynapse"]]
        v = np.array(v, dtype=float)
        c = np.clip(v, lo, hi)
        if (c != v).any():
            print(f"  --x0 {pth}: clipped to bounds at {[i for i in range(len(v)) if c[i] != v[i]]}")
        out.append(c)
    return out


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
    ap.add_argument("--save-json", default=None,
                    help="write the best params (all ten, apical written out, plus "
                         "tau_effca) here; run_de_fit2_pool.py --params-json reads it")
    # ---- OGK-03 (all new; none of these touch the --fit-dt path) ----
    ap.add_argument("--targets", default="",
                    help="comma-separated protocol_id from biodata/paired_recordings.csv, "
                         "e.g. mrk97_07,mrk97_08. One Batch per protocol, objective = sum of "
                         "((pred-target)/sem)^2. Replaces --fit-dt (do not combine).")
    ap.add_argument("--proto-map", default="",
                    help="JSON file or inline JSON: {protocol_id: {extracted, basis, suffix, proto}}. "
                         "extracted/basis default to $ANALYTICAL_EXTRACTED/$ANALYTICAL_BASIS_DIR; "
                         "suffix is appended to the npz tag ('<f>Hz_<dt>ms'+suffix+'.npz'); "
                         "proto overrides the tag. Lets two pathways use different dirs.")
    ap.add_argument("--x0", action="append", default=[],
                    help="fit json whose 'params' seed the DE population (repeatable; "
                         "CMA-ES starts from the first). Clipped to the bounds.")
    ap.add_argument("--untie-apical", action="store_true",
                    help="8 a-params (a20,a21,a30,a31 free). Needs --apical-syns.")
    ap.add_argument("--apical-syns", default="",
                    help="apical synapse global ids (.json list or {gid: bool}, .npz 'syn', "
                         "or text one id/line); all other synapses are basal. The extracted "
                         "npz carry no location, so this is required with --untie-apical.")
    ap.add_argument("--eval-only", action="store_true",
                    help="evaluate the objective once at x0 (first --x0, else the default "
                         "seed), print full-precision err + predictions, exit. No DE.")
    ap.add_argument("--dump-eval", default="",
                    help="with --eval-only --targets: write <PREFIX>.csv (one row per pair x protocol: offline "
                         "ratio, n_syn, mean c_pre/c_post/peak, rho0/rho_pred/rho_obs means, n apical) and "
                         "<PREFIX>_syn.npz (per-synapse arrays). OGK-13 diagnosis; default behaviour unchanged.")
    args = ap.parse_args()
    proto_mode = bool(args.targets)
    if proto_mode and (args.fit_dt or args.sign):
        raise SystemExit("--targets replaces --fit-dt and does not support --sign")
    if args.untie_apical and (args.cpre_only or not args.apical_syns):
        raise SystemExit("--untie-apical needs --apical-syns and excludes --cpre-only")

    # keep BLAS single-threaded: parallelism is across candidates, and nested
    # threads would oversubscribe 60 cores badly
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(v, "1")

    try:
        mp.set_start_method("fork")      # so workers inherit _B copy-on-write
    except RuntimeError:
        pass
    print(f"multiprocessing start method: {mp.get_start_method()}")

    if proto_mode:
        ids = [s.strip() for s in args.targets.split(",") if s.strip()]
        pmap = {}
        if args.proto_map:
            pmap = json.loads(args.proto_map if args.proto_map.lstrip().startswith("{")
                              else open(args.proto_map).read())
        _BS[:] = build_proto_batches(ids, pmap, args.anchor_tails, args.stride, args.pairs or None)
        apical = load_apical_ids(args.apical_syns) if args.apical_syns else None
        for p in _BS:
            if apical is not None:
                p["batch"].is_apical = np.isin(p["batch"].syn, np.fromiter(apical, np.int64))
            print(f"protocol {p['id']:9s} {p['proto']:12s} target {p['tgt']:.4f} +- {p['sem']:.4f}  "
                  f"{p['batch'].N} syn / {len(p['batch'].meta)} rec / T={p['batch'].T} / "
                  f"dt={p['batch'].dt_ms:.3f} ms  src={p['extracted']}  basis={p['basis']}")
        _B = PooledBatch(_BS)
        _POOL["peak"] = _B.peak
    else:
        _B = Batch(limit=args.pairs or None, stride=args.stride)
        if args.untie_apical:
            _B.is_apical = np.isin(_B.syn, np.fromiter(load_apical_ids(args.apical_syns), np.int64))
    OBJ = objective_proto if proto_mode else objective
    _CFG["untie"] = args.untie_apical
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
    NAMES = NAMES_A_CPRE if args.cpre_only else (NAMES_A8 if args.untie_apical else NAMES_A)
    NA    = len(NAMES)
    if args.cpre_only:
        print("method: cpre-only (theta_d = a00*c_pre, theta_p = a10*c_pre; "
              "a01 = a11 = 0)")
    print(f"{_B.N} synapses / {len(_B.meta)} records / T={_B.T} / "
          f"dt={_B.dt_ms:.3f} ms (stride {args.stride})")

    x0 = [TIED_IC4[k] for k in NAMES] + ([GAMMA_D_GB, GAMMA_P_GB] if args.fit_gamma else [])
    x0 = np.array(x0)
    # bounds (also built again below for the optimiser; needed here to clip --x0 seeds)
    _lo = np.array([args.a_lo]*NA + ([args.gd_lo, args.gp_lo] if args.fit_gamma else []))
    _hi = np.array([args.a_hi]*NA + ([args.gd_hi, args.gp_hi] if args.fit_gamma else []))
    seeds = seeds_from_json(args.x0, NAMES, args.fit_gamma, _lo, _hi) if args.x0 else []
    if seeds:
        x0 = seeds[0]
        print(f"--x0: {len(seeds)} seed(s) injected into the initial population; x0 = {args.x0[0]}")
    td0, tp0 = _B.thetas(unpack(x0, args.fit_gamma, args.cpre_only, args.untie_apical)[0])
    print(f"seed: theta_d med {np.median(td0):.4f}  theta_p med {np.median(tp0):.4f}  "
          f"peak effcai med {np.median(_B.peak):.4f}")
    print(f"  synapses with theta_d < peak effcai: {100*(td0 < _B.peak).mean():.0f}%  "
          f"(need >={100*MIN_ACTIVE:.0f}%)")
    print(f"  synapses with theta_p < peak effcai: {100*(tp0 < _B.peak).mean():.0f}%  "
          f"(need >={100*MIN_POT:.0f}%, else no potentiation is possible)")

    t0 = time.time(); e0 = OBJ(x0); dt1 = time.time() - t0
    n_par = (NA + (2 if args.fit_gamma else 0)) * args.popsize
    print(f"seed err {e0:.2f}   ({dt1:.1f} s / eval, {args.workers} workers)")
    if args.eval_only:
        print(f"EVAL err {e0!r}")
        if proto_mode:
            for k, v in proto_predictions(x0).items():
                print(f"EVAL pred {k} {v!r}")
            if args.dump_eval:
                import pandas as pd
                a_, gd_, gp_ = unpack(x0, args.fit_gamma, args.cpre_only, args.untie_apical)
                rows, arr = [], {}
                for p in _BS:
                    b = p["batch"]; td, tp = b.thetas(a_)
                    rho = b.rho_full(td, tp, gd_, gp_)
                    _, df = b.curve(rho); rat = dict(zip(df.pair, df.ratio))
                    off = 0
                    for m in b.meta:
                        sl = slice(off, off + m["n"]); off += m["n"]
                        ap_ = b.is_apical[sl] if b.is_apical is not None else np.zeros(m["n"], bool)
                        rows.append(dict(protocol=p["id"], pair=m["pair"], ratio_off=rat.get(m["pair"], np.nan),
                                         n_syn=m["n"], n_apical=int(ap_.sum()),
                                         c_pre=b.c_pre[sl].mean(), c_post=b.c_post[sl].mean(), peak=b.peak[sl].mean(),
                                         theta_d=td[sl].mean(), theta_p=tp[sl].mean(),
                                         rho0=b.rho0[sl].mean(), rho_pred=rho[sl].mean(), rho_obs=b.rho_obs[sl].mean()))
                        pre = f"{p['id']}__{m['pair']}__"
                        for k_, v_ in (("syn", b.syn[sl]), ("c_pre", b.c_pre[sl]), ("c_post", b.c_post[sl]),
                                       ("peak", b.peak[sl]), ("rho0", b.rho0[sl]), ("rho_pred", rho[sl]),
                                       ("rho_obs", b.rho_obs[sl]), ("apical", ap_)):
                            arr[pre + k_] = np.asarray(v_)
                pd.DataFrame(rows).to_csv(args.dump_eval + ".csv", index=False)
                np.savez_compressed(args.dump_eval + "_syn.npz", **arr)
                print(f"EVAL dumped {len(rows)} records -> {args.dump_eval}.csv / _syn.npz")
        else:
            td, tp = _B.thetas(unpack(x0, args.fit_gamma, args.cpre_only, args.untie_apical)[0])
            gdd, gpp = unpack(x0, args.fit_gamma, args.cpre_only, args.untie_apical)[1:]
            summ, _ = _B.curve(_B.rho_full(td, tp, gdd, gpp))
            for dtv, mv in zip(summ.dt, summ["mean"]):
                print(f"EVAL pred dt{int(dtv):+d} {float(mv)!r}")
        raise SystemExit
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
        de_init = {}
        if seeds:
            # population = Latin hypercube over the box, first rows replaced by the seeds.
            # (Only with --x0; without it scipy's own init is used, as before.)
            from scipy.stats import qmc
            S = max(args.popsize * len(lo), len(seeds))
            pop = lo + qmc.LatinHypercube(d=len(lo), seed=args.seed).random(S) * (hi - lo)
            pop[:len(seeds)] = np.array(seeds)
            de_init = {"init": pop}
        cfg = dict(_CFG)
        with Pool(args.workers, initializer=_init_worker,
                  initargs=(cfg, _B.src, args.pairs or None, args.stride)) as pool:
            r = differential_evolution(OBJ, list(zip(lo, hi)),
                                       maxiter=args.maxiter, popsize=args.popsize,
                                       polish=False, **de_init,
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
                es.tell(X, pool.map(OBJ, X))
                nev += len(X); es.disp()
        xb, eb = es.result.xbest, es.result.fbest

    a, gd, gp = unpack(xb, args.fit_gamma, args.cpre_only, args.untie_apical)
    print("\n=== best ===")
    for k in NAMES: print(f"  {k} = {a[k]:.6f}")
    if args.cpre_only:
        print("  a01 = a11 = 0.000000   (cpre-only method)")
    print("  (apical free: a20,a21,a30,a31 fitted)" if args.untie_apical else
          "  (apical tied: a20=a00, a21=a01, a30=a10, a31=a11)")
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
    if proto_mode:
        pred = proto_predictions(xb)
        print(f"\n  {'protocol':>9} {'proto':>12} {'fitted':>9} {'target':>9} {'sem':>7}")
        for p in _BS:
            v = pred[p["id"]]
            print(f"  {p['id']:>9} {p['proto']:>12} {('%.4f' % v) if v is not None else 'n/a':>9} "
                  f"{p['tgt']:>9.4f} {p['sem']:>7.4f}")
        summ = None
    else:
        td, tp = _B.thetas(a)
        summ, _ = _B.curve(_B.rho_full(td, tp, gd, gp))
        print(f"\n  {'dt':>5} {'fitted':>9} {'in vitro':>9}")
        for dt in sorted(summ.dt):
            tgt = f"{INVITRO[dt][0]:.4f}" if dt in INVITRO else "-"
            print(f"  {dt:>5} {summ[summ.dt==dt]['mean'].values[0]:>9.4f} {tgt:>9}")

    if args.save_json:
        # apical written out explicitly: simulator_edges reads a missing a-key as 1.0
        params = {"gamma_d_GB_GluSynapse": float(gd), "gamma_p_GB_GluSynapse": float(gp),
                  **{k: float(a[k]) for k in ("a00", "a01", "a10", "a11")},
                  "a20": float(a["a20"]), "a21": float(a["a21"]),
                  "a30": float(a["a30"]), "a31": float(a["a31"]),
                  "tau_effca_GB_GluSynapse": TAU_EFFCA_GB}
        extra = {}
        if proto_mode:      # protocol_id-keyed curve; offline_curve keys are protocol ids here
            curve = {k: v for k, v in pred.items()}
            extra = {"protocols": {p["id"]: dict(proto=p["proto"], target=p["tgt"], sem=p["sem"],
                                                 pred=pred[p["id"]], extracted=p["extracted"],
                                                 basis=p["basis"]) for p in _BS}}
        else:
            curve = {int(dt): float(summ[summ.dt == dt]["mean"].values[0]) for dt in summ.dt}
        with open(args.save_json, "w") as f:
            json.dump({"params": params, "weighted_err": float(eb), "evaluations": int(nev),
                       "offline_curve": curve, **extra,
                       "extracted": _B.src, "argv": os.sys.argv}, f, indent=2)
        print(f"\n  params -> {args.save_json}")
