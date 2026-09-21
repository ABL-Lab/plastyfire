"""
Validate a c_pre/c_post cache by RECOMPUTING it and diffing, per synapse.

The cache is keyed on (pre_gid, post_gid) only -- no emodel, no nseg, no tau. So a
cache silently survives any change to the circuit it was built on. The only way to
know whether a cache still describes the current emodel is to recompute some of it
and compare.

This runs the same entry point that built the cache
(simulator_edges.compute_cpre_cpost_for_workdir -> _find_cpre_cpost) on N random
pairs with the given --circuit-config, then reports per-synapse ratios.

Interpreting the output
-----------------------
  ratio ~1.00 everywhere      -> cache matches the current emodel
  ratio systematically != 1   -> cache is stale w.r.t. this circuit config; every
                                 downstream theta was computed from the wrong calcium
  a few synapses off          -> check whether they are apical (bAP-sensitive) --
                                 c_post is far more emodel-sensitive than c_pre

Note c_post depends on the single-AP threshold search, which is itself emodel
dependent; a pair that now fails to fire is reported rather than silently skipped.

Usage:
    python validate_cpre_cpost_recompute.py \
        --cache cpre_cpost_cache/ion_channels_tau278.pkl \
        --circuit-config data/dhuruva_modified_ion_channels_circuit_config.json \
        --pairs 10 --seed 0 --workers 10
"""

import argparse
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import pickle
import random
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

RESULTS_DIR = os.path.join(
    HERE, "refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations")
EDGES_H5 = os.path.join(HERE, "data/dhuruva_modified_edges.h5")
EDGE_POP = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"
APICAL = 3

# Same a-param preset the caches were built with: tau_effca is the only value that
# affects c_pre/c_post; the a-params are present only to trigger the code path.
FIT_PARAMS = {
    "tau_effca_GB_GluSynapse": 278.3177658387,
    "gamma_d_GB_GluSynapse": 101.5,
    "gamma_p_GB_GluSynapse": 216.2,
    "a00": 1.002, "a01": 1.954, "a10": 1.159, "a11": 2.483,
    "a20": 1.127, "a21": 2.456, "a30": 5.236, "a31": 1.782,
}


def _worker(arg):
    pair, workdir, circuit_config = arg
    import plastyfire.simulator_edges as sim_mod
    try:
        res = sim_mod.compute_cpre_cpost_for_workdir(
            workdir, fit_params=dict(FIT_PARAMS), circuit_config=circuit_config)
        if res is None or "c_pre" not in res:
            return pair, None, (res or {}).get("error", "returned None")
        return pair, {"c_pre": res["c_pre"], "c_post": res["c_post"]}, None
    except Exception as e:  # noqa: BLE001
        return pair, None, str(e)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cache", required=True)
    ap.add_argument("--circuit-config", required=True)
    ap.add_argument("--pairs", type=int, default=10)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--protocol", default="10Hz_-10ms")
    args = ap.parse_args()

    cache = pickle.load(open(os.path.join(HERE, args.cache), "rb"))
    with h5py.File(EDGES_H5, "r") as f:
        sec = f["edges/%s/0/afferent_section_type" % EDGE_POP][:]

    keys = sorted(cache.keys())
    random.seed(args.seed)
    chosen = sorted(random.sample(keys, min(args.pairs, len(keys))))

    jobs = []
    for (pre, post) in chosen:
        wd = os.path.join(RESULTS_DIR, "%d-%d" % (pre, post), args.protocol)
        if os.path.isfile(os.path.join(wd, "prefire_simulation_config.json")):
            jobs.append(((pre, post), wd, os.path.join(HERE, args.circuit_config)))

    print("cache          : %s  (%d pairs)" % (args.cache, len(cache)))
    print("circuit config : %s" % args.circuit_config)
    print("recomputing    : %d pairs (seed=%d)\n" % (len(jobs), args.seed))

    # ProcessPoolExecutor, not mp.Pool: compute_cpre_cpost_for_workdir forks a
    # grandchild for the NEURON run, and mp.Pool workers are daemonic ("daemonic
    # processes are not allowed to have children"). Executor workers are not.
    # spawn keeps each worker free of inherited NEURON/HOC state.
    ctx = mp.get_context("spawn")
    results = []
    with ProcessPoolExecutor(max_workers=min(args.workers, len(jobs)),
                             mp_context=ctx) as ex:
        futs = [ex.submit(_worker, j) for j in jobs]
        for f in as_completed(futs):
            results.append(f.result())

    print("%-18s %-7s %5s %11s %11s %8s %11s %11s %8s" % (
        "pair", "loc", "n", "cpre_cache", "cpre_new", "ratio", "cpost_cache",
        "cpost_new", "ratio"))
    print("-" * 104)

    pre_ratios, post_ratios, post_ratios_ap, post_ratios_ba = [], [], [], []
    failed = []
    for pair, got, err in results:
        if got is None:
            failed.append((pair, err))
            continue
        old = cache[pair]
        for loc_name, want_ap in (("apical", True), ("basal", False)):
            gids = [g for g in sorted(old["c_pre"]) if (sec[g] == APICAL) == want_ap]
            if not gids:
                continue
            op = np.array([old["c_pre"][g] for g in gids])
            npv = np.array([got["c_pre"].get(g, np.nan) for g in gids])
            oq = np.array([old["c_post"][g] for g in gids])
            nq = np.array([got["c_post"].get(g, np.nan) for g in gids])
            rp = np.nanmedian(npv / op)
            rq = np.nanmedian(nq / oq)
            pre_ratios.extend((npv / op).tolist())
            post_ratios.extend((nq / oq).tolist())
            (post_ratios_ap if want_ap else post_ratios_ba).extend((nq / oq).tolist())
            print("%-18s %-7s %5d %11.4e %11.4e %8.3f %11.4e %11.4e %8.3f" % (
                "%d-%d" % pair, loc_name, len(gids),
                np.median(op), np.nanmedian(npv), rp,
                np.median(oq), np.nanmedian(nq), rq))

    def summarize(name, arr):
        a = np.array([x for x in arr if np.isfinite(x)])
        if not a.size:
            print("  %-26s (no data)" % name)
            return
        print("  %-26s n=%3d  median %.3f  within +/-10%%: %d/%d  within +/-1%%: %d/%d" % (
            name, a.size, np.median(a),
            int((np.abs(a - 1) <= 0.10).sum()), a.size,
            int((np.abs(a - 1) <= 0.01).sum()), a.size))

    print("\n--- recomputed / cached ---")
    summarize("c_pre  (all)", pre_ratios)
    summarize("c_post (all)", post_ratios)
    summarize("c_post (apical)", post_ratios_ap)
    summarize("c_post (basal)", post_ratios_ba)

    a = np.array([x for x in pre_ratios if np.isfinite(x)])
    b = np.array([x for x in post_ratios if np.isfinite(x)])
    tol = 0.01
    if a.size and b.size:
        ok = (np.abs(a - 1) <= tol).all() and (np.abs(b - 1) <= tol).all()
        print("\nVERDICT: %s" % (
            "cache MATCHES the current emodel (all synapses within 1%)" if ok else
            "cache DOES NOT match the current emodel -- recompute before trusting "
            "any theta derived from it"))

    if failed:
        print("\nFailed to recompute (%d):" % len(failed))
        for pair, err in failed:
            print("  %d-%d: %s" % (pair[0], pair[1], str(err)[:110]))


if __name__ == "__main__":
    main()
