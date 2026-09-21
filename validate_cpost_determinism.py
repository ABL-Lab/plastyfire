"""
Is c_post reproducible at all?

validate_cpre_cpost_recompute.py found c_pre reproduces exactly (68/68 within 1%)
while c_post scatters (32/68 within 10%, one pair off by 206x). Two very different
explanations:

  (a) the cache is stale -- the emodel changed since it was built
  (b) c_post is not deterministic -- recomputing the SAME pair on the SAME emodel
      gives a different answer each time

They need different fixes, and (a) cannot be diagnosed while (b) is unresolved. This
recomputes each pair TWICE in one job, same emodel, same workdir, and compares
run-to-run. Any spread here is (b).

Also records the chosen stimulus (amp, width) and the resulting spike count per run,
because c_post is calcium from a *single* backpropagating AP -- if the two runs pick
different amplitudes, or fire different numbers of spikes, that is the mechanism.
"""

import argparse
import multiprocessing as mp
import os
import pickle
import random
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RESULTS_DIR = os.path.join(
    HERE, "refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations")

FIT_PARAMS = {
    "tau_effca_GB_GluSynapse": 278.3177658387,
    "gamma_d_GB_GluSynapse": 101.5, "gamma_p_GB_GluSynapse": 216.2,
    "a00": 1.002, "a01": 1.954, "a10": 1.159, "a11": 2.483,
    "a20": 1.127, "a21": 2.456, "a30": 5.236, "a31": 1.782,
}


def _worker(arg):
    pair, rep, workdir, cfg = arg
    import plastyfire.simulator_edges as sim_mod
    try:
        r = sim_mod.compute_cpre_cpost_for_workdir(
            workdir, fit_params=dict(FIT_PARAMS), circuit_config=cfg)
        if r is None or "c_pre" not in r:
            return pair, rep, None, (r or {}).get("error", "None")
        return pair, rep, {"c_pre": r["c_pre"], "c_post": r["c_post"],
                           "n_spikes": r.get("n_post_spikes"),
                           "amp": r.get("pulse_amp")}, None
    except Exception as e:  # noqa: BLE001
        return pair, rep, None, str(e)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cache", default="cpre_cpost_cache/ion_channels_tau278.pkl")
    ap.add_argument("--circuit-config",
                    default="data/dhuruva_modified_ion_channels_circuit_config.json")
    ap.add_argument("--pairs", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--protocol", default="10Hz_-10ms")
    args = ap.parse_args()

    cache = pickle.load(open(os.path.join(HERE, args.cache), "rb"))
    random.seed(args.seed)
    chosen = sorted(random.sample(sorted(cache.keys()), min(args.pairs, len(cache))))

    jobs = []
    for p in chosen:
        wd = os.path.join(RESULTS_DIR, "%d-%d" % p, args.protocol)
        if os.path.isfile(os.path.join(wd, "prefire_simulation_config.json")):
            for rep in (1, 2):
                jobs.append((p, rep, wd, os.path.join(HERE, args.circuit_config)))

    print("recomputing %d pairs x 2 replicates = %d runs\n" % (len(jobs) // 2, len(jobs)))

    got = {}
    ctx = mp.get_context("spawn")
    with ProcessPoolExecutor(max_workers=min(args.workers, len(jobs)),
                             mp_context=ctx) as ex:
        futs = [ex.submit(_worker, j) for j in jobs]
        for f in as_completed(futs):
            pair, rep, res, err = f.result()
            got[(pair, rep)] = (res, err)

    print("%-16s %6s %6s %9s %8s %12s %12s %9s" % (
        "pair", "amp1", "amp2", "spikes", "n", "cpost_run1", "cpost_run2", "ratio"))
    print("-" * 92)
    all_pre, all_post = [], []
    for p in chosen:
        r1, e1 = got.get((p, 1), (None, "missing"))
        r2, e2 = got.get((p, 2), (None, "missing"))
        if r1 is None or r2 is None:
            print("%-16s FAILED  r1=%s r2=%s" % ("%d-%d" % p, e1, e2))
            continue
        gids = sorted(set(r1["c_post"]) & set(r2["c_post"]))
        q1 = np.array([r1["c_post"][g] for g in gids])
        q2 = np.array([r2["c_post"][g] for g in gids])
        p1 = np.array([r1["c_pre"][g] for g in gids])
        p2 = np.array([r2["c_pre"][g] for g in gids])
        all_pre.extend((p2 / p1).tolist())
        all_post.extend((q2 / q1).tolist())
        print("%-16s %6.3f %6.3f %4s/%-4s %8d %12.4e %12.4e %9.3f" % (
            "%d-%d" % p, r1["amp"] or np.nan, r2["amp"] or np.nan,
            r1["n_spikes"], r2["n_spikes"], len(gids),
            np.median(q1), np.median(q2), np.median(q2 / q1)))

    def summarize(name, arr):
        a = np.array([x for x in arr if np.isfinite(x)])
        if not a.size:
            print("  %-22s (no data)" % name)
            return
        print("  %-22s n=%3d  median %.4f  identical: %d/%d  within 1%%: %d/%d" % (
            name, a.size, np.median(a),
            int((np.abs(a - 1) < 1e-9).sum()), a.size,
            int((np.abs(a - 1) <= 0.01).sum()), a.size))

    print("\n--- run2 / run1 (same emodel, same workdir) ---")
    summarize("c_pre", all_pre)
    summarize("c_post", all_post)
    b = np.array([x for x in all_post if np.isfinite(x)])
    if b.size:
        print("\nVERDICT: c_post is %s" % (
            "DETERMINISTIC -- scatter vs the cache means the cache is stale"
            if (np.abs(b - 1) < 1e-9).all() else
            "NOT deterministic -- the cache is one draw, not a wrong answer; the "
            "threshold search does not converge to the same stimulus each run"))
    multi = [(p, got[(p, r)][0]["n_spikes"]) for p in chosen for r in (1, 2)
             if got.get((p, r), (None,))[0] and (got[(p, r)][0]["n_spikes"] or 0) > 1]
    if multi:
        print("\nWARNING: %d run(s) measured c_post with MORE THAN ONE spike:" % len(multi))
        for p, n in multi:
            print("  %d-%d: %d spikes -- c_post is not single-AP calcium here" % (p[0], p[1], n))


if __name__ == "__main__":
    main()
