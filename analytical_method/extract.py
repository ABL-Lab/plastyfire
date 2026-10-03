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

ROOT     = "/project/rrg-emuller/dhuruva/plastyfire"
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


def window_index(n, dt, events, fine_step, coarse_step, pre_ms, post_ms):
    """Sample indices of a variable-rate grid over n raw samples: every `fine_step`-th sample
    within [e - pre_ms, e + post_ms] of any event time e (pre/post spikes), every `coarse_step`-th
    elsewhere. The first and last raw samples are always kept."""
    fine = np.zeros(n, bool)
    for e in np.asarray(events, dtype=np.float64):
        i0 = max(0, int(np.floor((e - pre_ms) / dt))); i1 = min(n, int(np.ceil((e + post_ms) / dt)) + 1)
        fine[i0:i1] = True
    k = np.arange(n)
    keep = np.where(fine, k % fine_step == 0, k % coarse_step == 0)
    keep[0] = keep[-1] = True
    return np.flatnonzero(keep).astype(np.int64)


def protocols_from_args(protocols, index_csv, default):
    """Protocol ids to run: comma list `protocols`, else the protocol_id column of a simwriter
    index csv (e.g. index_Ebner2019_L5TTPC_L5TTPC.csv), else `default` (the 7 Markram 10Hz_* ids)."""
    if protocols:
        return [p.strip() for p in protocols.split(",") if p.strip()]
    if index_csv:
        import csv
        with open(index_csv) as f:
            return list(dict.fromkeys(row["protocol_id"] for row in csv.DictReader(f)))
    return list(default)


def _one(job):
    pair, proto, trace_f, res_f, cpre, cpost, decim, out, win = job
    try:
        tr = pickle.load(open(trace_f, "rb"))
        rs = pickle.load(open(res_f, "rb"))
        t  = tr["t"]
        dt = float(np.median(np.diff(t[:1000])))
        step = max(1, int(round(decim / dt)))
        if win:  # variable-rate grid: fine around pre/post spikes, coarse in the inter-pairing gaps
            # float32 t makes the median step 0.02499962: round it, or t_ms drifts 0.6 ms per 42 s
            dt = round(dt, 5)
            pre_t = np.asarray(tr.get("prespikes", rs.get("prespikes", [])), dtype=np.float64)
            post_t = np.asarray(rs.get("postspikes", []), dtype=np.float64)
            keep = window_index(len(t), dt, np.concatenate([pre_t, post_t]), step,
                                max(step, int(round(win["coarse"] / dt))), win["pre"], win["post"])
        else:
            keep = slice(None, None, step)

        rho0 = dict(zip(rs["global_ids"], rs["initial_rho"]))
        fin  = dict(zip(rs["global_ids"], rs["final_rho"]))

        sids, eff, r0, rf, cp, cq = [], [], [], [], [], []
        for sid, cai in tr["cai_CR"].items():
            sid = int(sid)
            if sid not in cpre or sid not in rho0:
                continue
            sids.append(sid)
            eff.append(effcai_from_cai(cai, dt)[keep].astype(np.float32))
            r0.append(rho0[sid]); rf.append(fin[sid])
            cp.append(cpre[sid]); cq.append(cpost[sid])
        if not sids:
            return pair, proto, 0, "no overlapping synapses"

        os.makedirs(out, exist_ok=True)
        extra = {}
        if win:
            # t_ms: sample times (non-uniform). Integrate rho with h[k] = t_ms[k+1] - t_ms[k]
            # (left point); dt_ms is the fine step, coarse_dt_ms the gap step.
            extra = dict(t_ms=keep * dt, coarse_dt_ms=np.float64(dt * max(step, int(round(win["coarse"] / dt)))),
                         win_pre_ms=np.float64(win["pre"]), win_post_ms=np.float64(win["post"]),
                         pre_t_ms=pre_t, post_t_ms=post_t, raw_dt_ms=np.float64(dt))
        np.savez_compressed(
            os.path.join(out, f"{pair}__{proto}.npz"),
            syn=np.array(sids, dtype=np.int64),
            effcai=np.stack(eff),                       # (n_syn, n_t) float32
            dt_ms=np.float64(dt * step),
            rho0=np.array(r0), rho_obs=np.array(rf),
            c_pre=np.array(cp), c_post=np.array(cq),
            **extra,
        )
        return pair, proto, len(sids), None
    except Exception as e:
        return pair, proto, 0, str(e)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--limit", type=int, default=0, help="first N pairs only")
    ap.add_argument("--param-hash", default="cdf3a1e1db98")
    ap.add_argument("--sims", default=None, help="..._STDP/simulations dir (default: n100 set)")
    ap.add_argument("--cache", default=CACHE, help="cpre/cpost cache pkl")
    ap.add_argument("--out", default=None, help="output dir (default: per-decim dir here)")
    ap.add_argument("--protocols", default=None,
                    help="comma-separated protocol ids; default: constants.PROTOCOLS (10Hz_*)")
    ap.add_argument("--index-csv", default=None,
                    help="simwriter index csv: extract every protocol_id in it")
    ap.add_argument("--decim", type=float, default=DECIM_MS,
                    help="decimation in ms (default 2.0). 0.25 gives 8x finer "
                         "threshold-crossing resolution at 8x the storage.")
    ap.add_argument("--window", action="store_true",
                    help="variable-rate grid for long protocols: --decim within [-win-pre, +win-post] ms "
                         "of every pre/post spike, --coarse ms elsewhere; adds t_ms (sample times) to the npz")
    ap.add_argument("--win-pre", type=float, default=50.0, help="ms before each spike kept fine (default 50)")
    ap.add_argument("--win-post", type=float, default=1000.0,
                    help="ms after each spike kept fine (default 1000 = 3.6 tau_effca)")
    ap.add_argument("--coarse", type=float, default=5.0, help="gap sampling in ms (default 5)")
    ap.add_argument("--skip-existing", action="store_true", help="skip jobs whose npz exists")
    ap.add_argument("--pairs", default=None, help="comma-separated pair dirs (<pre>-<post>); default: all")
    args = ap.parse_args()
    win = dict(pre=args.win_pre, post=args.win_post, coarse=args.coarse) if args.window else None

    SIMS = args.sims or globals()["SIMS"]
    cache = pickle.load(open(args.cache, "rb"))
    cpre, cpost = {}, {}
    for v in cache.values():
        for s, x in v["c_pre"].items():  cpre[int(s)]  = x
        for s, x in v["c_post"].items(): cpost[int(s)] = x

    out = args.out or out_dir_for(args.decim)
    jobs = []
    pairs = sorted(p for p in os.listdir(SIMS) if "-" in p)
    if args.pairs:
        pairs = [p for p in pairs if p in set(args.pairs.split(","))]
    if args.limit:
        pairs = pairs[:args.limit]
    protocols = protocols_from_args(args.protocols, args.index_csv, PROTOCOLS)
    for pair in pairs:
        for proto in protocols:
            # traces of THIS prefire run: per-hash subdir, or the legacy shared one
            tf = [f for f in (f"{SIMS}/{pair}/{proto}/bluecellulab_results_{args.param_hash}/simulation_traces.pkl",
                              f"{SIMS}/{pair}/{proto}/bluecellulab_results_optimizer/simulation_traces.pkl")
                  if os.path.isfile(f)][:1]
            rf = f"{SIMS}/{pair}/{proto}/simulation_edges_{args.param_hash}.pkl"
            if tf and os.path.isfile(rf):
                if args.skip_existing and os.path.isfile(os.path.join(out, f"{pair}__{proto}.npz")):
                    continue
                jobs.append((pair, proto, tf[0], rf, cpre, cpost, args.decim, out, win))

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
