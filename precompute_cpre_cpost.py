"""
Precompute c_pre and c_post for all unique pairs in a results directory and
save them to a cache pickle.

The cache is keyed by (pre_gid, post_gid) and stores:
    {(pre_gid, post_gid): {"c_pre": {syn_gid: float}, "c_post": {syn_gid: float}}}

c_pre/c_post depend only on the pair and the HOC params (tau_effca, gamma_d,
gamma_p) — NOT on the protocol (freq/dt). So one cache entry covers ALL
protocols for a given pair.

Usage:
    python precompute_cpre_cpost.py --params fitting2 --results-dir <root>
    python precompute_cpre_cpost.py --params fitting2 \\
        --results-dir /lustre06/project/6077694/dhuruva/plastyfire/working_new_refitting_results \\
        --output /project/ctb-emuller/dhuruva/plastyfire/cpre_cpost_cache/fitting2.pkl

Then pass to pairrunner_edges.py:
    python pairrunner_edges.py --params fitting2 \\
        --cpre-cpost-cache /project/ctb-emuller/dhuruva/plastyfire/cpre_cpost_cache/fitting2.pkl
"""

import argparse
import glob
import json
import logging
import multiprocessing
import os
import pickle
import tempfile
import time

import bluecellulab
from bluepysnap import Simulation

import plastyfire.simulator_edges as sim_mod
from plastyfire.simulator_edges import (
    _find_cpre_cpost,
    MECHANISMS_PATH,
    EDGES_H5_DEFAULT,
    EDGE_POP_DEFAULT,
    NODE_POP_DEFAULT,
)

# ── Parameter presets (keep in sync with submit_edges_sims.py) ────────────────
PARAM_PRESETS = {
    "chindemi": {
        "tau_effca_GB_GluSynapse": 278.3177658387,
        "gamma_d_GB_GluSynapse":   101.5387594661,
        "gamma_p_GB_GluSynapse":   216.1841700668,
    },
    # Same tau as chindemi but includes a-params so the cache-building path is triggered.
    # Use this to pre-populate cpre_cpost_cache.pkl for the edges GA optimizer
    # (which always runs with FITTED_TAU=278.32 and a-params from the GA).
    "chindemi_aparams": {
        "tau_effca_GB_GluSynapse": 278.3177658387,
        "gamma_d_GB_GluSynapse":   101.5,
        "gamma_p_GB_GluSynapse":   216.2,
        "a00": 1.002, "a01": 1.954,
        "a10": 1.159, "a11": 2.483,
        "a20": 1.127, "a21": 2.456,
        "a30": 5.236, "a31": 1.782,
    },
    "fitted": {
        "tau_effca_GB_GluSynapse": 200.4818,
        "gamma_d_GB_GluSynapse":    60.6964,
        "gamma_p_GB_GluSynapse":   178.2431,
        "a00": 0.7967, "a01": 1.1547,
        "a10": 1.1382, "a11": 1.9027,
        "a20": 2.2637, "a21": 7.0244,
        "a30": 3.5947, "a31": 6.2696,
    },
    "fitting2": {
        "tau_effca_GB_GluSynapse": 308.0155,
        "gamma_d_GB_GluSynapse":    76.1153,
        "gamma_p_GB_GluSynapse":   181.0415,
        "a00": 1.0378, "a01": 1.2593,
        "a10": 1.2523, "a11": 2.2936,
        "a20": 1.9727, "a21": 8.3983,
        "a30": 6.0229, "a31": 1.4087,
    },
    "ga_best": {
        # Best-fit parameters from GA optimizer (DEAP, gen 28, fitness=0.3942)
        "tau_effca_GB_GluSynapse": 278.3177658387,
        "gamma_d_GB_GluSynapse":    94.698202,
        "gamma_p_GB_GluSynapse":   211.404335,
        "a00": 1.000406, "a01": 1.954000,
        "a10": 1.160326, "a11": 2.590750,
        "a20": 1.164831, "a21": 2.489699,
        "a30": 3.670398, "a31": 1.566718,
    },
}

log = logging.getLogger(__name__)


def find_workdirs(results_dir):
    """Find all valid simulation workdirs under results_dir."""
    pattern = os.path.join(
        results_dir, "fitting", "*", "seed*", "*_STDP",
        "simulations", "*-*", "*Hz_*ms",
    )
    dirs = sorted(glob.glob(pattern))
    return [d for d in dirs
            if os.path.isfile(os.path.join(d, "prefire_simulation_config.json"))
            and os.path.isfile(os.path.join(d, "prefire_prespikes.h5"))]


def get_pre_post(workdir):
    """Return (pre_gid, post_gid) by reading node_sets from the workdir config."""
    sim_config = os.path.join(workdir, "prefire_simulation_config.json")
    with open(sim_config) as f:
        cfg = json.load(f)
    node_sets_local = os.path.join(workdir, "node_sets.json")
    if os.path.exists(node_sets_local):
        cfg["node_sets_file"] = node_sets_local
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".json", dir=tempfile.gettempdir(), delete=False)
    json.dump(cfg, tmp)
    tmp.close()
    try:
        snap = Simulation(tmp.name)
        pre_gid  = snap.node_sets.content["precell"]["node_id"][0]
        post_gid = snap.node_sets.content["postcell"]["node_id"][0]
    finally:
        os.unlink(tmp.name)
    return int(pre_gid), int(post_gid)


def _worker(conn, workdir, fit_params, node_pop, edge_pop, fixhp):
    """Subprocess: build tmp config, run mini-sims, return c_pre/c_post."""
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    _logger = logging.getLogger(__name__)
    bluecellulab.neuron.load_mechanisms(MECHANISMS_PATH)

    try:
        sim_config = os.path.join(workdir, "prefire_simulation_config.json")
        with open(sim_config) as f:
            cfg = json.load(f)
        node_sets_local = os.path.join(workdir, "node_sets.json")
        if os.path.exists(node_sets_local):
            cfg["node_sets_file"] = node_sets_local
        out_dir_local = os.path.join(workdir, "out")
        os.makedirs(out_dir_local, exist_ok=True)
        cfg.setdefault("output", {})["output_dir"] = out_dir_local
        cfg.pop("reports", None)
        cfg.get("conditions", {}).get("mechanisms", {}).pop("GluSynapse", None)
        tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".json", dir=tempfile.gettempdir(), delete=False)
        json.dump(cfg, tmp)
        tmp.close()
        sim_config_eff = tmp.name

        # Determine pre/post gids from the already-loaded cfg
        snap_cfg = json.load(open(sim_config))
        if os.path.exists(node_sets_local):
            snap_cfg["node_sets_file"] = node_sets_local
        snap_tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".json", dir=tempfile.gettempdir(), delete=False)
        json.dump(snap_cfg, snap_tmp)
        snap_tmp.close()
        try:
            snap = Simulation(snap_tmp.name)
            pre_gid  = int(snap.node_sets.content["precell"]["node_id"][0])
            post_gid = int(snap.node_sets.content["postcell"]["node_id"][0])
        finally:
            os.unlink(snap_tmp.name)

        c_pre, c_post, _ = _find_cpre_cpost(
            sim_config_eff, cfg, fit_params,
            pre_gid, post_gid, node_pop, edge_pop, fixhp, _logger,
            sim_config_orig=sim_config,
        )
        conn.send({"c_pre": c_pre, "c_post": c_post, "pre_gid": pre_gid, "post_gid": post_gid})
    except Exception as e:
        _logger.error("Worker failed: %s", e, exc_info=True)
        conn.send(None)
    finally:
        conn.close()
        try:
            os.unlink(sim_config_eff)
        except Exception:
            pass


def compute_pair(workdir, fit_params, node_pop, edge_pop, fixhp, circuit_config=None,
                 extracellular_calcium=None):
    """Run mini-sims in a subprocess; return the full result dict.

    Delegates to simulator_edges.compute_cpre_cpost_for_workdir so the
    subprocess / pipe logic lives in one place.
    """
    result = sim_mod.compute_cpre_cpost_for_workdir(
        workdir, fit_params=fit_params, node_pop=node_pop,
        edge_pop=edge_pop, fixhp=fixhp, circuit_config=circuit_config,
        extracellular_calcium=extracellular_calcium,
    )
    if result is None:
        raise RuntimeError(f"Mini-sim subprocess returned None for {workdir}")
    if "error" in result and "c_pre" not in result:
        # Subprocess caught an exception; re-raise with the original message so
        # the caller can classify it (e.g. no-spike vs. crash).
        raise RuntimeError(result["error"])
    return result


def _pool_worker(args):
    """Top-level picklable worker for multiprocessing.Pool."""
    workdir, fit_params, node_pop, edge_pop, fixhp, circuit_config, extracellular_calcium = args
    try:
        res = compute_pair(
            workdir, fit_params, node_pop, edge_pop, fixhp, circuit_config=circuit_config,
            extracellular_calcium=extracellular_calcium,
        )
        return (res["pre_gid"], res["post_gid"], res["c_pre"], res["c_post"],
                {"n_post_spikes": res.get("n_post_spikes"),
                 "n_syn": res.get("n_syn"),
                 "n_syn_at_floor": res.get("n_syn_at_floor")}, None)
    except Exception as e:
        # Extract gids from workdir name as fallback
        pair_str = os.path.basename(os.path.dirname(workdir))
        parts = pair_str.split("-")
        pre = int(parts[0]) if len(parts) == 2 else -1
        post = int(parts[1]) if len(parts) == 2 else -1
        # A threshold-search failure ("could not fire post cell ...") is the
        # explicit no-spike case; keep it separable from other crashes.
        _msg = str(e)
        # NB: this covers "never found a single-AP amplitude", which is not the
        # same as "cell is unexcitable" — the cell may fire 2+ APs at every
        # amplitude tried. Reported as "no usable single-AP stimulus".
        no_spike = ("could not fire post cell" in _msg
                    or "did not spike" in _msg
                    or "c_post invalid" in _msg)
        return pre, post, None, None, {"no_spike": no_spike}, str(e)


class _NoDaemonProcess(multiprocessing.Process):
    """Pool worker that can itself spawn child processes (NEURON subprocess)."""
    @property
    def daemon(self):
        return False
    @daemon.setter
    def daemon(self, value):
        pass


class _NoDaemonPool(multiprocessing.pool.Pool):
    def Process(self, *args, **kwds):
        proc = super().Process(*args, **kwds)
        proc.__class__ = _NoDaemonProcess
        return proc


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--params", required=True, choices=list(PARAM_PRESETS.keys()),
                        help="Parameter preset to use for mini-sims")
    parser.add_argument("--results-dir", required=True,
                        help="Root results directory containing fitting/ tree")
    parser.add_argument("--output", default=None,
                        help="Output cache pkl path (default: <results-dir>/cpre_cpost_<params>.pkl)")
    parser.add_argument("--edges-h5", default=EDGES_H5_DEFAULT)
    parser.add_argument("--node-pop", default=NODE_POP_DEFAULT)
    parser.add_argument("--edge-pop", default=EDGE_POP_DEFAULT)
    parser.add_argument("--no-fixhp", action="store_true",
                        help="Disable SK_E2 removal (fixhp=False)")
    parser.add_argument("--workers", type=int, default=1,
                        help="Number of parallel workers (default: 1 — serial, avoids OOM). "
                             "Increase with caution: each bluecellulab instance uses ~4-8 GB RAM.")
    parser.add_argument("--circuit-config", default=None,
                        help="Override the 'network' field in each prefire_simulation_config.json "
                             "(e.g. for new ion channels)")
    parser.add_argument("--extracellular-calcium", type=float, default=None,
                        help="Extracellular calcium (mM) for the mini-sims. c_pre/c_post are "
                             "calcium-dependent (cao_CR sets Pf_NMDA and Eca_syn in "
                             "GluSynapse.mod), so a cache is only valid at the Ca it was "
                             "built at. Default: whatever the workdir config carries (2.0).")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")

    fit_params = dict(PARAM_PRESETS[args.params])
    fixhp = not args.no_fixhp
    output = args.output or os.path.join(args.results_dir, f"cpre_cpost_{args.params}.pkl")

    # Find all workdirs and collect unique (pre_gid, post_gid) → one representative workdir
    log.info("Scanning workdirs under: %s", args.results_dir)
    workdirs = find_workdirs(args.results_dir)
    log.info("Found %d workdirs total", len(workdirs))

    # Load existing cache if resuming
    cache = {}
    if os.path.exists(output):
        with open(output, "rb") as f:
            cache = pickle.load(f)
        log.info("Loaded existing cache with %d entries from %s", len(cache), output)

    # Map unique pairs to one representative workdir (first protocol found)
    pair_to_workdir = {}
    for wd in workdirs:
        pair_str = os.path.basename(os.path.dirname(wd))
        parts = pair_str.split("-")
        if len(parts) == 2:
            try:
                pre, post = int(parts[0]), int(parts[1])
                if (pre, post) not in pair_to_workdir:
                    pair_to_workdir[(pre, post)] = wd
            except ValueError:
                pass

    log.info("Unique pairs: %d", len(pair_to_workdir))
    todo = [((pre, post), wd) for (pre, post), wd in pair_to_workdir.items()
            if (pre, post) not in cache]
    log.info("Pairs needing computation: %d (already cached: %d)",
             len(todo), len(cache))

    if not todo:
        log.info("All pairs already cached.")
        return

    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)

    worker_args = [
        (wd, fit_params, args.node_pop, args.edge_pop, fixhp, args.circuit_config,
         args.extracellular_calcium)
        for (_, wd) in todo
    ]

    n_workers = min(args.workers, len(todo))
    log.info("Running %d pairs with %d parallel workers", len(todo), n_workers)

    done = failed = 0
    # Spike bookkeeping. A pair is only a usable c_post measurement if the post
    # cell fired AND the resulting bAP deposited calcium at ≥1 synapse; a spike
    # with every synapse at the resting floor measures nothing (see
    # simulator_edges.C_POST_FLOOR_M).
    n_spiked = n_no_spike = n_all_floor = 0
    syn_total = syn_at_floor = 0
    no_spike_pairs, all_floor_pairs = [], []

    with _NoDaemonPool(processes=n_workers) as pool:
        for pre_gid, post_gid, c_pre, c_post, stats, err in pool.imap_unordered(_pool_worker, worker_args):
            if err is not None:
                log.error("FAILED pair (%d, %d): %s", pre_gid, post_gid, err)
                failed += 1
                if (stats or {}).get("no_spike"):
                    n_no_spike += 1
                    no_spike_pairs.append((pre_gid, post_gid))
            else:
                cache[(pre_gid, post_gid)] = {"c_pre": c_pre, "c_post": c_post}
                done += 1
                stats = stats or {}
                n_syn = stats.get("n_syn") or len(c_post)
                n_floor = stats.get("n_syn_at_floor") or 0
                syn_total    += n_syn
                syn_at_floor += n_floor
                if (stats.get("n_post_spikes") or 0) > 0:
                    n_spiked += 1
                if n_syn and n_floor == n_syn:
                    n_all_floor += 1
                    all_floor_pairs.append((pre_gid, post_gid))
                log.info("[%d/%d] Cached pair (%d, %d) — c_pre syns=%d, "
                         "post spikes=%s, c_post at floor=%d/%d",
                         done, len(todo), pre_gid, post_gid, len(c_pre),
                         stats.get("n_post_spikes"), n_floor, n_syn)
                # Save incrementally after each completed pair
                with open(output, "wb") as f:
                    pickle.dump(cache, f, protocol=-1)

    attempted = len(todo)
    log.info("Done: %d cached, %d failed. Cache → %s", done, failed, output)
    log.info("=" * 70)
    log.info("c_post spike report (%d pairs attempted this run)", attempted)
    log.info("  post cell FIRED               : %d / %d", n_spiked, attempted)
    log.info("  no usable single-AP stimulus  : %d / %d", n_no_spike, attempted)
    if failed - n_no_spike:
        log.info("  failed for other reasons      : %d", failed - n_no_spike)
    log.info("  fired but bAP reached NO synapse (all c_post at floor <%.0e): %d / %d",
             sim_mod.C_POST_FLOOR_M, n_all_floor, attempted)
    log.info("  --> usable c_post measurements: %d / %d", n_spiked - n_all_floor, attempted)
    if syn_total:
        log.info("  synapses at resting floor     : %d / %d (%.1f%%)",
                 syn_at_floor, syn_total, 100.0 * syn_at_floor / syn_total)
    if no_spike_pairs:
        log.info("  no-spike pairs: %s",
                 ", ".join(f"({a},{b})" for a, b in sorted(no_spike_pairs)[:20])
                 + (" …" if len(no_spike_pairs) > 20 else ""))
    if all_floor_pairs:
        log.info("  all-floor pairs: %s",
                 ", ".join(f"({a},{b})" for a, b in sorted(all_floor_pairs)[:20])
                 + (" …" if len(all_floor_pairs) > 20 else ""))
    log.info("=" * 70)


if __name__ == "__main__":
    main()
