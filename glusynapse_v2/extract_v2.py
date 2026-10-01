"""Traces -> compact npz for the GluSynapse_v2 offline model (tsk T9).

Same as analytical_method/extract.py, plus what the presynaptic pathways need:
  shaft_cai  (section Ca from CaDynamics, drives the pre-LTD trace T and the NO trace N)
  prespikes  (presynaptic spike times: mGluR steps and the Z trace)
  vdcc       (-ica_VDCC, nA, spine VDCC influx: the pre_drive 1 signal). An AP gives a ~1 ms pulse, so
             it is stored as the mean over each grid step [keep[k], keep[k+1]) (charge-preserving),
             not point-sampled like the slow traces.
effcai is re-derived from the SAME trace file, so effcai, shaft_cai and rho_obs always belong
to one BCL run (the Markram 10Hz_* traces on disk are from the delta-cooker run: the shared
bluecellulab_results_optimizer dir was overwritten after delta-prefire).

    python glusynapse_v2/extract_v2.py --param-hash delta-cooker --workers 30 \
        --out glusynapse_v2/extracted/markram_delta-cooker
"""
import argparse, os, pickle, sys
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

ROOT = "/project/rrg-emuller/dhuruva/plastyfire"
sys.path.insert(0, os.path.join(ROOT, "analytical_method"))
from extract import effcai_from_cai, protocols_from_args, window_index   # noqa: E402
import model_v2 as MV                                         # noqa: E402  (v_events, t_drive 3 constants)
from constants import PROTOCOLS                               # noqa: E402

SIMS  = os.path.join(ROOT, "refitting_results/fitting/n120/seed20262009/"
                           "Sabrina_L5TTPC_L5TTPC_STDP/simulations")
CACHE = os.path.join(ROOT, "cpre_cpost_cache/sabrina_n120_delta.pkl")


def step_mean(x, keep, ends):
    """Mean of x over [keep[k], ends[k]) for each grid sample k (float32)."""
    c = np.concatenate([[0.0], np.cumsum(x)])
    return ((c[ends] - c[keep]) / (ends - keep)).astype(np.float32)


def _one(job):
    pair, proto, trace_f, res_f, cpre, cpost, decim, out, win, k_ca = job
    try:
        tr = pickle.load(open(trace_f, "rb"))
        rs = pickle.load(open(res_f, "rb"))
        t  = np.asarray(tr["t"])
        dt = float(np.median(np.diff(t[:1000])))
        step = max(1, int(round(decim / dt)))
        pre_t = np.asarray(tr.get("prespikes", rs.get("prespikes", [])), dtype=np.float64)
        post_t = np.asarray(rs.get("postspikes", []), dtype=np.float64)
        n_raw = len(next(iter(tr["cai_CR"].values())))
        if win:   # same grid as analytical_method/extract.py --window
            keep = window_index(n_raw, dt, np.concatenate([pre_t, post_t]), step,
                                max(step, int(round(win["coarse"] / dt))), win["pre"], win["post"])
        else:
            keep = np.arange(0, n_raw, step)
        # sample times as index * NEURON dt. The median of float32 diffs gives 0.02499962, which drifts
        # 0.6 ms by 42 s and would misplace spike arrivals; effcai keeps dt as extract.py does.
        t_ms = float(t[0]) + keep * round(dt, 5)
        rho0 = dict(zip(rs["global_ids"], rs["initial_rho"]))
        fin  = dict(zip(rs["global_ids"], rs["final_rho"]))
        has_shaft = "shaft_cai" in tr
        has_ica = "ica_VDCC" in tr
        has_v = "v_seg" in tr
        vevs = []
        cevs = {k: [] for k in ("cev", "cev_lo", "cev_hi")}
        ends = np.append(keep[1:], n_raw)
        sids, eff, sh, vd, r0, rf, cp, cq = [], [], [], [], [], [], [], []
        for sid, cai in tr["cai_CR"].items():
            sid = int(sid)
            if sid not in cpre or sid not in rho0:
                continue
            sids.append(sid)
            eff.append(effcai_from_cai(cai, dt)[keep].astype(np.float32))
            if has_shaft:
                sh.append(np.asarray(tr["shaft_cai"][sid], dtype=np.float32)[keep])
            if has_ica:
                vd.append(step_mean(-np.asarray(tr["ica_VDCC"][sid], dtype=np.float64), keep, ends))
            if has_v:   # t_drive 3: crossing times from the FULL-resolution trace; the raw v is not stored
                vevs.append(MV.v_events(tr["v_seg"][sid], round(dt, 5), t0=float(t[0])))
            if has_ica:   # t_drive 4: crossings of -ica_VDCC through K (and K/2, 2K) at the recorded resolution
                q = -np.asarray(tr["ica_VDCC"][sid], dtype=np.float64)
                for k, m in (("cev", 1.0), ("cev_lo", 0.5), ("cev_hi", 2.0)):
                    cevs[k].append(MV.ca_events(q, round(dt, 5), t0=float(t[0]), K=k_ca * m))
            r0.append(rho0[sid]); rf.append(fin[sid]); cp.append(cpre[sid]); cq.append(cpost[sid])
        if not sids:
            return pair, proto, 0, "no overlapping synapses"
        if not has_shaft:
            return pair, proto, 0, "no shaft_cai in traces"
        extra = {}
        if has_v:
            nm = max([len(e) for e in vevs] + [1]); vev = np.full((len(vevs), nm), np.nan)
            for i, e in enumerate(vevs):
                vev[i, :len(e)] = e
            extra["vev"] = vev
        if has_ica:
            for k, L in cevs.items():
                nm = max([len(e) for e in L] + [1]); a = np.full((len(L), nm), np.nan)
                for i, e in enumerate(L):
                    a[i, :len(e)] = e
                extra[k] = a
        os.makedirs(out, exist_ok=True)
        np.savez_compressed(
            os.path.join(out, f"{pair}__{proto}.npz"),
            syn=np.array(sids, dtype=np.int64),
            effcai=np.stack(eff), shaft_cai=np.stack(sh), **({"vdcc": np.stack(vd)} if has_ica else {}), **extra,
            dt_ms=np.float64(dt * step), t0_ms=np.float64(t[0]), t_ms=t_ms, raw_dt_ms=np.float64(round(dt, 5)),
            K_ca=np.float64(k_ca), windowed=np.bool_(bool(win)), prespikes=pre_t, postspikes=post_t,
            rho0=np.array(r0), rho_obs=np.array(rf),
            c_pre=np.array(cp), c_post=np.array(cq))
        return pair, proto, len(sids), None
    except Exception as e:
        return pair, proto, 0, f"{type(e).__name__}: {e}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=30)
    ap.add_argument("--param-hash", required=True,
                    help="names simulation_edges_<hash>.pkl and bluecellulab_results_<hash>/")
    ap.add_argument("--sims", default=SIMS)
    ap.add_argument("--cache", default=CACHE)
    ap.add_argument("--out", required=True)
    ap.add_argument("--protocols", default=None)
    ap.add_argument("--index-csv", default=None)
    ap.add_argument("--decim", type=float, default=0.25)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--pairs", default=None, help="comma-separated pair dirs; default all")
    ap.add_argument("--window", action="store_true", help="windowed grid, as extract.py --window")
    ap.add_argument("--win-pre", type=float, default=50.0)
    ap.add_argument("--win-post", type=float, default=1000.0)
    ap.add_argument("--coarse", type=float, default=5.0)
    ap.add_argument("--skip-existing", action="store_true")
    ap.add_argument("--k-ca", type=float, default=None,
                    help="t_drive 4 influx threshold K (nA) for cev (cev_lo K/2, cev_hi 2K); default model_v2.DEFAULTS['K_ca'] "
                         "(ljp 0). Scales with the VDCC driving force: use the ljp-25 value for ljp-25 traces.")
    args = ap.parse_args()
    k_ca = MV.DEFAULTS["K_ca"] if args.k_ca is None else args.k_ca
    win = dict(pre=args.win_pre, post=args.win_post, coarse=args.coarse) if args.window else None

    cache = pickle.load(open(args.cache, "rb"))
    cpre, cpost = {}, {}
    for v in cache.values():
        for s, x in v["c_pre"].items():  cpre[int(s)] = x
        for s, x in v["c_post"].items(): cpost[int(s)] = x

    protocols = protocols_from_args(args.protocols, args.index_csv, PROTOCOLS)
    pairs = sorted(p for p in os.listdir(args.sims) if "-" in p)
    if args.pairs:
        pairs = [p for p in pairs if p in set(args.pairs.split(","))]
    pairs = pairs[: args.limit or None]
    jobs = []
    for pair in pairs:
        for proto in protocols:
            d = f"{args.sims}/{pair}/{proto}"
            tf = [f for f in (f"{d}/bluecellulab_results_{args.param_hash}/simulation_traces.pkl",
                              f"{d}/bluecellulab_results_optimizer/simulation_traces.pkl")
                  if os.path.isfile(f)][:1]
            rf = f"{d}/simulation_edges_{args.param_hash}.pkl"
            if args.skip_existing and os.path.isfile(os.path.join(args.out, f"{pair}__{proto}.npz")):
                continue
            if tf and os.path.isfile(rf):
                jobs.append((pair, proto, tf[0], rf, cpre, cpost, args.decim, args.out, win, k_ca))
    print(f"jobs {len(jobs)}  workers {args.workers}  decim {args.decim} ms  -> {args.out}", flush=True)
    ok = bad = 0
    with ProcessPoolExecutor(args.workers) as ex:
        for f in as_completed([ex.submit(_one, j) for j in jobs]):
            pair, proto, n, err = f.result()
            if err: bad += 1; print(f"  FAIL {pair}/{proto}: {err}", flush=True)
            else:   ok += 1
    print(f"done: {ok} written, {bad} failed", flush=True)


if __name__ == "__main__":
    main()
