"""Stage 1 (expensive, run once): traces -> compact per-(pair,protocol) npz.

The 343 MB simulation_traces.pkl files hold cai_CR at 0.025 ms. effcai_GB is a
leaky integral of it (tau=278 ms), so it is smooth on a millisecond scale and
decimating to DECIM_MS loses nothing that matters for time-above-threshold.
That turns ~235 GB of traces into a few hundred MB that the optimiser can hold.

Nothing here depends on the a-params — that is the whole point. Run once, then
every candidate parameter set is evaluated from these arrays.

    python extract.py --workers 12
"""
import argparse, glob, os, pickle, sys
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from scipy.signal import lfilter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from constants import MIN_CA_CR, TAU_EFFCA_GB, PROTOCOLS

ROOT     = "/lustre06/project/6077694/dhuruva/plastyfire"
SIMS     = os.path.join(ROOT, "refitting_results/fitting/n100/seed19091997/"
                              "L5TTPC_L5TTPC_STDP/simulations")
CACHE    = os.path.join(ROOT, "cpre_cpost_cache/ion_channels_tau278.pkl")
HERE     = os.path.dirname(os.path.abspath(__file__))
DECIM_MS = 2.0   # default; override with --decim


def out_dir_for(decim):
    """Each decimation gets its own directory so resolutions never mix."""
    return os.path.join(HERE, "extracted" if decim == 2.0
                        else f"extracted_d{decim:g}".replace(".", "p"))


def effcai_from_cai(cai, dt_ms):
    """Exact exponential integrator for  effcai' = -effcai/tau + (cai - min_ca)."""
    a = np.exp(-dt_ms / TAU_EFFCA_GB)
    b = TAU_EFFCA_GB * (1.0 - a)
    u = np.asarray(cai, dtype=np.float64) - MIN_CA_CR
    return lfilter([0.0, b], [1.0, -a], u)


def _one(job):
    pair, proto, trace_f, res_f, cpre, cpost, decim, out = job
    try:
        tr = pickle.load(open(trace_f, "rb"))
        rs = pickle.load(open(res_f, "rb"))
        t  = tr["t"]
        dt = float(np.median(np.diff(t[:1000])))
        step = max(1, int(round(decim / dt)))

        rho0 = dict(zip(rs["global_ids"], rs["initial_rho"]))
        fin  = dict(zip(rs["global_ids"], rs["final_rho"]))

        sids, eff, r0, rf, cp, cq = [], [], [], [], [], []
        for sid, cai in tr["cai_CR"].items():
            sid = int(sid)
            if sid not in cpre or sid not in rho0:
                continue
            sids.append(sid)
            eff.append(effcai_from_cai(cai, dt)[::step].astype(np.float32))
            r0.append(rho0[sid]); rf.append(fin[sid])
            cp.append(cpre[sid]); cq.append(cpost[sid])
        if not sids:
            return pair, proto, 0, "no overlapping synapses"

        os.makedirs(out, exist_ok=True)
        np.savez_compressed(
            os.path.join(out, f"{pair}__{proto}.npz"),
            syn=np.array(sids, dtype=np.int64),
            effcai=np.stack(eff),                       # (n_syn, n_t) float32
            dt_ms=np.float64(dt * step),
            rho0=np.array(r0), rho_obs=np.array(rf),
            c_pre=np.array(cp), c_post=np.array(cq),
        )
        return pair, proto, len(sids), None
    except Exception as e:
        return pair, proto, 0, str(e)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--limit", type=int, default=0, help="first N pairs only")
    ap.add_argument("--param-hash", default="cdf3a1e1db98")
    ap.add_argument("--decim", type=float, default=DECIM_MS,
                    help="decimation in ms (default 2.0). 0.25 gives 8x finer "
                         "threshold-crossing resolution at 8x the storage.")
    args = ap.parse_args()

    cache = pickle.load(open(CACHE, "rb"))
    cpre, cpost = {}, {}
    for v in cache.values():
        for s, x in v["c_pre"].items():  cpre[int(s)]  = x
        for s, x in v["c_post"].items(): cpost[int(s)] = x

    out = out_dir_for(args.decim)
    jobs = []
    pairs = sorted(p for p in os.listdir(SIMS) if "-" in p)
    if args.limit:
        pairs = pairs[:args.limit]
    for pair in pairs:
        for proto in PROTOCOLS:
            tf = glob.glob(f"{SIMS}/{pair}/{proto}/*/simulation_traces.pkl")
            rf = f"{SIMS}/{pair}/{proto}/simulation_edges_{args.param_hash}.pkl"
            if tf and os.path.isfile(rf):
                jobs.append((pair, proto, tf[0], rf, cpre, cpost, args.decim, out))

    print(f"jobs: {len(jobs)}  workers: {args.workers}  decim: {args.decim} ms")
    print(f"out : {out}")
    ok = bad = 0
    with ProcessPoolExecutor(args.workers) as ex:
        for f in as_completed([ex.submit(_one, j) for j in jobs]):
            pair, proto, n, err = f.result()
            if err: bad += 1;  print(f"  FAIL {pair}/{proto}: {err}", flush=True)
            else:   ok  += 1
            if (ok + bad) % 50 == 0:
                print(f"  [{ok+bad}/{len(jobs)}] ok={ok} fail={bad}", flush=True)
    print(f"done: {ok} written, {bad} failed -> {out}")


if __name__ == "__main__":
    main()
