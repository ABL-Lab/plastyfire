"""
Prefire connected-pair simulator using dhuruva_modified_edges.h5.

Runs the prefire simulation (pairing protocol only, no C01/C02 test pulses) using
prefire_simulation_config.json + prefire_prespikes.h5.

bluecellulab auto-loads rho0_GB, Use_d/p, gmax_d/p from the circuit.  We only
inject theta_d_GB / theta_p_GB, which are stored under different names in edges.h5.

After simulation, writes bluecellulab_results/rho.h5 in SONATA report format so that
plot_refitting_stdp_basis.py can read it directly via the existing basis pipeline.
"""

import json
import os
import pickle
import tempfile
import time
import logging
import multiprocessing

import h5py
import numpy as np
from libsonata import SpikeReader
from bluepysnap import Simulation
import bluecellulab

from plastyfire.simulator import (
    _get_spikes, _map_syn_idx, _set_global_params, spike_threshold_finder,
)

logger = logging.getLogger(__name__)

EDGES_H5_DEFAULT  = "/project/ctb-emuller/dhuruva/plastyfire/data/dhuruva_modified_edges.h5"
EDGE_POP_DEFAULT  = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"
NODE_POP_DEFAULT  = "S1nonbarrel_neurons"
MECHANISMS_PATH   = "/project/ctb-emuller/dhuruva/DEES_cell_packages/"

# Keys that trigger per-synapse theta override via c_pre/c_post computation
_A_PARAM_KEYS = {"a00", "a01", "a10", "a11", "a20", "a21", "a30", "a31"}

# effcai_GB resting value is ~2-5e-5. A c_post at or below this means no
# backpropagating AP calcium arrived at that synapse — the measurement failed
# rather than returning a genuinely small number. Used only for reporting.
C_POST_FLOOR_M = 1e-4

# effcai_GB never approaches this, so a synapse given this threshold never arms
# its WATCH and stays frozen. Used by PLASTYFIRE_ONLY_LOC (see
# _apply_theta_from_a_params) to study one dendritic population in isolation.
THETA_SILENCE = 1e9

# c_post single-AP threshold search. These mirror thresholdfinder.py::run()
# exactly — same widths, same amplitude range, same number of levels — so the
# c_post stimulus here is the one the reference implementation would have found.
#
# The binary search in spike_threshold_finder assumes the leftmost amplitude
# firing >=nspikes also fires ==nspikes. For nspikes=1 that holds only while the
# pulse is short enough that rheobase produces a single AP rather than a train,
# which is why the widths stop at 5 ms. Do not extend this list upward.
C_POST_PULSE_WIDTHS_MS = [1.5, 3, 5]
C_POST_MIN_AMP_NA      = 0.05
C_POST_MAX_AMP_NA      = 5.
C_POST_AMP_LEVELS      = 100   # -> 0.05 nA grid over [0.05, 5.] nA


# ---------------------------------------------------------------------------
# HDF5 helpers
# ---------------------------------------------------------------------------

def _load_thresholds(edges_h5, edge_pop, global_ids):
    """
    Load theta_d and theta_p from edges.h5 for the given global synapse IDs.
    Returns {global_id: (theta_d, theta_p)}.
    """
    pop_path   = f"edges/{edge_pop}/0"
    sorted_ids = np.array(sorted(set(global_ids)), dtype=np.int64)

    with h5py.File(edges_h5, "r") as f:
        td_path = f"{pop_path}/theta_d"
        tp_path = f"{pop_path}/theta_p"
        if td_path not in f or tp_path not in f:
            raise KeyError(
                f"theta_d / theta_p not found under {pop_path} in {edges_h5}. "
                "Run inject_thresholds.py first."
            )
        td_vals = f[td_path][sorted_ids]
        tp_vals = f[tp_path][sorted_ids]

    return {
        int(gid): (float(td), float(tp))
        for gid, td, tp in zip(sorted_ids, td_vals, tp_vals)
    }


def _apply_thresholds(cell, global_ids, thresholds):
    """Set theta_d_GB / theta_p_GB on each synapse from the precomputed thresholds."""
    for (syn_id, synapse), gid in zip(cell.synapses.items(), global_ids):
        td, tp = thresholds.get(gid, (-1.0, -1.0))
        synapse.hsynapse.theta_d_GB = td if td > 0 else -1.0
        synapse.hsynapse.theta_p_GB = tp if tp > 0 else -1.0


def _write_rho_h5(path, node_pop, global_ids, initial_rho, final_rho):
    """
    Write a minimal SONATA synapse report with initial and final rho values.
    Shape: data (2, n_synapses) — data[0]=initial, data[-1]=final.
    Compatible with plot_refitting_stdp_basis.py's _read_rho_h5().
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    data = np.array([initial_rho, final_rho], dtype=np.float32)
    with h5py.File(path, "w") as f:
        grp     = f.require_group(f"report/{node_pop}")
        grp.create_dataset("data", data=data)
        mapping = grp.require_group("mapping")
        mapping.create_dataset("element_ids", data=np.array(global_ids, dtype=np.int64))
        mapping.create_dataset("time", data=np.array([0.0, 1.0, 1.0], dtype=np.float64))


# ---------------------------------------------------------------------------
# Validation helper
# ---------------------------------------------------------------------------

def validate_edges_loading(cell, global_ids, edges_h5, edge_pop, tol=1e-4):
    """
    Compare what bluecellulab loaded from the circuit against dhuruva_modified_edges.h5.
    Checks rho0_GB, Use_d, Use_p, gmax_d_AMPA, gmax_p_AMPA, theta_d_GB, theta_p_GB.
    Logs a warning per mismatch and returns list of (global_id, field, bcl_val, edges_val).
    """
    pop_path   = f"edges/{edge_pop}/0"
    sorted_ids = np.array(sorted(set(global_ids)), dtype=np.int64)

    field_map = {
        "rho0_GB":    "rho0_GB",
        "Use_d_TM":   "Use_d",
        "Use_p_TM":   "Use_p",
        "gmax_d_AMPA":"gmax_d_AMPA",
        "gmax_p_AMPA":"gmax_p_AMPA",
        "theta_d":    "theta_d_GB",
        "theta_p":    "theta_p_GB",
    }

    edges_vals = {}
    with h5py.File(edges_h5, "r") as f:
        for field in field_map:
            ds_path = f"{pop_path}/{field}"
            if ds_path in f:
                arr = f[ds_path][sorted_ids]
                edges_vals[field] = dict(zip(sorted_ids.tolist(), arr.tolist()))

    mismatches = []
    for (syn_id, synapse), gid in zip(cell.synapses.items(), global_ids):
        h = synapse.hsynapse
        for field, mod_attr in field_map.items():
            if field not in edges_vals:
                continue
            edges_v = edges_vals[field].get(gid)
            if edges_v is None:
                continue
            try:
                bcl_v = float(getattr(h, mod_attr))
            except AttributeError:
                continue
            if abs(bcl_v - edges_v) > tol:
                mismatches.append((gid, mod_attr, bcl_v, edges_v))
                logger.warning(
                    "Mismatch syn %d  %-14s  bluecellulab=%.6f  edges=%.6f  Δ=%.2e",
                    gid, mod_attr, bcl_v, edges_v, abs(bcl_v - edges_v),
                )

    if not mismatches:
        logger.info("validate_edges_loading: all %d synapses match edges.h5", len(global_ids))
    else:
        logger.warning(
            "validate_edges_loading: %d mismatches across %d synapses",
            len(mismatches), len(global_ids),
        )
    return mismatches


# ---------------------------------------------------------------------------
# c_pre / c_post finders (mini-sims, same approach as simulator.py)
# ---------------------------------------------------------------------------

def _read_syn_extra_params(sim_config, pre_gid, post_gid, edge_pop, _logger):
    """Per-synapse GluSynapse params from the SONATA edge file, keyed by global edge id.

    Same source thresholdfinder.py uses in its no-recipe branch (read_sonata_params),
    so the values match what produced the published thresholds.
    """
    from bluepysnap import Simulation
    from plastyfire.thresholdfinder import read_sonata_params

    edges = Simulation(sim_config).circuit.edges[edge_pop]
    syn_df = read_sonata_params(edges, np.array([pre_gid]), post_gid)
    params = syn_df.drop(columns=["@source_node"]).to_dict(orient="index")
    _logger.info("Loaded per-synapse params for %d synapses (volume_CR %.3f-%.3f)",
                 len(params),
                 min(p["volume_CR"] for p in params.values()) if params else float("nan"),
                 max(p["volume_CR"] for p in params.values()) if params else float("nan"))
    return params



def _find_cpre_cpost(sim_config_eff, cfg, fit_params, pre_gid, post_gid,
                     node_pop, edge_pop, fixhp, _logger, sim_config_orig=None):
    """
    Compute c_pre and c_post per synapse by running two short NEURON sims:
      - c_pre : deliver 1 pre spike, record peak effcai_GB per synapse
      - c_post: inject current to fire 1 AP in post cell, record peak effcai_GB

    fit_params (with the new tau_effca / gamma_d / gamma_p) are applied so the
    calcium transients themselves reflect the fitted dynamics.

    Returns:
        c_pre  : dict {global_syn_id (int): float}
        c_post : dict {global_syn_id (int): float}
        df     : DataFrame mapping local_syn_idx → global_syn_id
    """
    from plastyfire.simulator import _map_syn_idx, _set_global_params, _set_local_params

    # Per-synapse parameters from the SONATA edge file. Without these every synapse
    # runs on the mod file's global defaults, and volume_CR in particular (the spine
    # volume that scales calcium concentration) varies 0.10-0.36 across a single
    # connection — so apical c_post stayed at the resting floor ~7e-5 instead of the
    # reference ~0.073. Both simulator.py c_pre_finder/c_post_finder pass these via
    # syn_extra_params; this path did not, which is what that 1000x gap was.
    syn_extra_params = _read_syn_extra_params(sim_config_orig, pre_gid, post_gid,
                                              edge_pop, _logger)

    # ── c_pre: 1 pre spike, no post firing ──────────────────────────────────
    _logger.info("Computing c_pre (1 pre spike sim) …")
    sim_pre = bluecellulab.CircuitSimulation(sim_config_eff, base_seed=cfg["run"]["random_seed"])
    sim_pre.instantiate_gids(
        [(node_pop, post_gid)],
        add_synapses=True, add_minis=False, add_pulse_stimuli=False,
        intersect_pre_gids=[(node_pop, pre_gid)],
        pre_spike_trains={(node_pop, pre_gid): [1000.0]},
    )
    cell_pre = sim_pre.cells[(node_pop, post_gid)]
    if fixhp:
        for sec in cell_pre.somatic + cell_pre.axonal:
            sec.uninsert("SK_E2")
    if fit_params:
        _set_global_params(fit_params)
    # Map local->global first: _set_local_params needs the global id to find each
    # synapse's parameters.
    syn_idx_pre = [syn_id[1] for syn_id in cell_pre.synapses]
    df = _map_syn_idx(sim_config_orig, post_gid, syn_idx_pre, edge_pop)

    rec_effcai_pre = {}
    for syn_id, synapse in cell_pre.synapses.items():
        v = bluecellulab.neuron.h.Vector()
        v.record(synapse.hsynapse._ref_effcai_GB, 1.0)
        rec_effcai_pre[syn_id[1]] = v
        _gid = int(df.loc[df["local_syn_idx"] == syn_id[1]].index[0])
        if syn_extra_params and _gid in syn_extra_params:
            _set_local_params(synapse, fit_params, syn_extra_params[_gid])
        # Force potentiated state to match compute_thresholds.py calibration:
        # edges.h5 thresholds were computed with rho0_GB=1, Use_p=1 so that
        # c_pre reflects max NMDA Ca2+ influx. Using mixed initial state here
        # gives 10-80x smaller c_pre → 10-80x smaller thresholds → wrong rho dynamics.
        # Force fully-potentiated state to match compute_thresholds.py calibration.
        # compute_thresholds.py sets rho0_GB=1, Use_p=1, Use=1, gmax0_AMPA=gmax_p_AMPA
        # so c_pre reflects maximum NMDA Ca2+ influx.  We must do the same here.
        # Exactly the four properties simulator.py::_c_pre_finder_process sets
        # (lines 262-265). Use_GB and gmax_AMPA are deliberately NOT set: they appear
        # only in the induction setup, not in either c_pre reference implementation.
        h = synapse.hsynapse
        h.rho0_GB    = 1.0
        h.Use_p      = 1.0
        h.Use        = h.Use_p
        h.gmax0_AMPA = h.gmax_p_AMPA
        h.theta_d_GB = -1.0
        h.theta_p_GB = -1.0
    sim_pre.run(1500.0, cvode=True)
    # Take max after the spike (t=1000ms), matching compute_thresholds.py's arr[i0:].max()
    c_pre = {
        int(df.loc[df["local_syn_idx"] == sid].index[0]): float(
            max(list(rec_effcai_pre[sid])[1000:])  # recording at dt=1ms → index 1000 = t=1000ms
        )
        for sid in syn_idx_pre
    }
    _logger.info("c_pre: %s", {gid: f"{v:.4e}" for gid, v in c_pre.items()})

    # ── c_post: 1 post AP (stimulus found by threshold search) ──────────────
    _logger.info("Computing c_post (1 post AP sim) …")

    # 1500 ms, matching simulator.py::_c_post_finder_process — the lineage that
    # produced the published thresholds. (compute_thresholds.py uses 3000 ms, but it
    # is the reimplementation, not the reference.)
    C_POST_TSTOP_MS = 1500.0

    # Find the stimulus the way thresholdfinder.py does, instead of scraping
    # amp_start from the workdir config. That scraped value is the *induction train*
    # amplitude (the single_cells pkl records nspikes=5, freq=10), which is not a
    # single-spike threshold: it failed to fire the cell on 4 of 5 pairs tested.
    # amp and width are searched TOGETHER — the finder returns the amp that fires
    # exactly one AP at a given width, so an amp is only valid with the width it
    # was found for.
    stim_amp, stim_width = None, None
    for pulse_width in C_POST_PULSE_WIDTHS_MS:
        simres = spike_threshold_finder(
            sim_config_eff, post_gid, 1, 0.1, pulse_width, 1000.,
            C_POST_MIN_AMP_NA, C_POST_MAX_AMP_NA, C_POST_AMP_LEVELS,
            node_pop, fixhp)
        if simres is not None:
            stim_amp, stim_width = float(simres["amp"]), float(simres["width"])
            _logger.info("c_post stimulus: %.4f nA, %.1f ms (threshold search)",
                         stim_amp, stim_width)
            break
        _logger.info("c_post: no single AP at width=%.1f ms up to %.1f nA — widening pulse",
                     pulse_width, C_POST_MAX_AMP_NA)
    if stim_amp is None:
        raise RuntimeError(
            f"spike_threshold_finder could not fire post cell {post_gid} with "
            f"exactly one AP at any width in {C_POST_PULSE_WIDTHS_MS} ms up to "
            f"{C_POST_MAX_AMP_NA} nA — c_post invalid")

    sim_post = bluecellulab.CircuitSimulation(sim_config_eff, base_seed=cfg["run"]["random_seed"])
    sim_post.instantiate_gids(
        [(node_pop, post_gid)],
        add_synapses=True, add_minis=False, add_pulse_stimuli=False,
        intersect_pre_gids=[(node_pop, pre_gid)],
    )
    cell_post = sim_post.cells[(node_pop, post_gid)]
    if fixhp:
        for sec in cell_post.somatic + cell_post.axonal:
            sec.uninsert("SK_E2")
    # h.TStim + train(), exactly as simulator.py::_c_post_finder_process does — NOT a
    # bare IClamp. spike_threshold_finder searched the amp against a TStim train, so
    # the amp is only valid when delivered the same way.
    tstim = bluecellulab.neuron.h.TStim(0.5, sec=cell_post.soma)
    _stim_duration = (1 - 1) * 1000. / 0.1 + stim_width   # nspikes=1, freq=0.1
    tstim.train(1000., _stim_duration, stim_amp, 0.1, stim_width)
    cell_post.persistent.append(tstim)
    if fit_params:
        _set_global_params(fit_params)
    rec_effcai_post, syn_idx_post = {}, []
    for syn_id, synapse in cell_post.synapses.items():
        syn_idx_post.append(syn_id[1])
        v = bluecellulab.neuron.h.Vector()
        v.record(synapse.hsynapse._ref_effcai_GB, 1.0)
        rec_effcai_post[syn_id[1]] = v
        # Same per-synapse params as c_pre — volume_CR above all, since c_post is a
        # calcium concentration and volume_CR is what it is divided by.
        _gid = int(df.loc[df["local_syn_idx"] == syn_id[1]].index[0])
        if syn_extra_params and _gid in syn_extra_params:
            _set_local_params(synapse, fit_params, syn_extra_params[_gid])
        synapse.hsynapse.theta_d_GB = -1.0
        synapse.hsynapse.theta_p_GB = -1.0
    sim_post.run(C_POST_TSTOP_MS, cvode=True)

    # Did the cell actually fire? c_post is calcium from the BACKPROPAGATING AP, so
    # no spike means no measurement — and it fails silently, returning the resting
    # effcai_GB floor (~3e-5) instead of ~0.07, which looks like a small number
    # rather than a broken run. The stimulus is now searched against these very
    # emodels, so a failure here means the search itself did not transfer (small
    # integration differences between the threshold sim and this one) rather than
    # a stale amplitude.
    _t = np.array(sim_post.get_time())
    _v = np.array(sim_post.get_voltage_trace((node_pop, post_gid)))
    _spikes = _get_spikes(_t, _v)
    _n_post_spikes = int(len(_spikes))
    if _n_post_spikes == 0:
        # c_post_finder does exactly this: bump by 0.05 nA and retry once before
        # giving up, for the borderline case where the threshold sim and this sim
        # land on opposite sides of the spike.
        _logger.warning(
            "c_post: post cell %d did not spike at the searched amp=%.4f nA "
            "(Vm peak %.1f mV) — retrying at +0.05 nA",
            post_gid, stim_amp, float(_v.max()) if _v.size else float("nan"))
        # TStim has no settable .amp — rebuild the train at the bumped amplitude.
        tstim = bluecellulab.neuron.h.TStim(0.5, sec=cell_post.soma)
        tstim.train(1000., _stim_duration, stim_amp + 0.05, 0.1, stim_width)
        cell_post.persistent.append(tstim)
        for _v_rec in rec_effcai_post.values():
            _v_rec.resize(0)
        sim_post.run(C_POST_TSTOP_MS, cvode=True)
        _t = np.array(sim_post.get_time())
        _v = np.array(sim_post.get_voltage_trace((node_pop, post_gid)))
        _n_post_spikes = int(len(_get_spikes(_t, _v)))
        if _n_post_spikes == 0:
            raise RuntimeError(
                f"post cell {post_gid} did not spike at {stim_amp:.4f} nA nor at "
                f"{stim_amp + 0.05:.4f} nA ({stim_width} ms) — c_post invalid")
        stim_amp += 0.05
    _logger.info("c_post: post cell %d fired %d spike(s) at amp=%.4f nA (%.1f ms)",
                 post_gid, _n_post_spikes, stim_amp, stim_width)

    c_post = {
        int(df.loc[df["local_syn_idx"] == sid].index[0]): float(rec_effcai_post[sid].max())
        for sid in syn_idx_post
    }
    _logger.info("c_post: %s", {gid: f"{v:.4e}" for gid, v in c_post.items()})

    # A somatic spike is NOT proof that c_post was measured. c_post is calcium
    # deposited by the BACKPROPAGATING AP, so a cell that fires at the soma but
    # whose bAP dies before reaching the synapses returns the effcai_GB resting
    # floor (~2-5e-5) with n_post_spikes=1 and no error. Count the synapses left
    # at that floor so the caller can distinguish "small c_post" from
    # "no measurement at all".
    _n_syn = len(c_post)
    _n_at_floor = sum(1 for v in c_post.values() if v < C_POST_FLOOR_M)
    if _n_syn and _n_at_floor == _n_syn:
        _logger.warning(
            "c_post: post cell %d fired %d spike(s) but ALL %d synapses are at the "
            "resting floor (<%.0e) — the bAP did not reach any synapse, so c_post "
            "carries no signal for this pair",
            post_gid, _n_post_spikes, _n_syn, C_POST_FLOOR_M)
    elif _n_at_floor:
        _logger.info("c_post: %d/%d synapses at resting floor (<%.0e)",
                     _n_at_floor, _n_syn, C_POST_FLOOR_M)

    # Stashed on df rather than widening the return tuple, so both existing
    # call sites (which unpack 3 values) keep working unchanged.
    df.attrs = getattr(df, "attrs", {})
    df.attrs["n_post_spikes"] = _n_post_spikes
    df.attrs["n_syn"] = _n_syn
    df.attrs["n_syn_at_floor"] = _n_at_floor
    df.attrs["pulse_amp"] = float(stim_amp)
    df.attrs["pulse_width"] = float(stim_width)
    return c_pre, c_post, df


def _apply_theta_from_a_params(cell, global_ids, fit_params, c_pre, c_post,
                                edges_h5, edge_pop, _logger):
    """
    Override theta_d_GB / theta_p_GB per synapse using a00-a31 coefficients:
        basal:  theta_d = a00*c_pre + a01*c_post
                theta_p = a10*c_pre + a11*c_post
        apical: theta_d = a20*c_pre + a21*c_post
                theta_p = a30*c_pre + a31*c_post

    Section type 3 = apical dendrite; all others treated as basal.
    """
    with h5py.File(edges_h5, "r") as f:
        sec_types = f[f"edges/{edge_pop}/0/afferent_section_type"][:]

    # Optional location isolation: PLASTYFIRE_ONLY_LOC=apical|basal silences the
    # other population so one can be studied alone. Silencing is done by pushing
    # both thresholds out of reach (THETA_SILENCE), NOT by setting them to -1:
    # GluSynapse.mod arms `WATCH (effcai_GB > theta_d_GB)`, so a negative theta
    # fires immediately and pins dep_GB = pot_GB = 1 — the opposite of silencing.
    # With both watches unreachable, dep_GB = pot_GB = 0 and rho_GB' collapses to
    # the bistable drift term, so rho stays at its initial attractor.
    only_loc = os.environ.get("PLASTYFIRE_ONLY_LOC", "").strip().lower() or None
    if only_loc not in (None, "apical", "basal"):
        raise ValueError(f"PLASTYFIRE_ONLY_LOC must be apical|basal, got {only_loc!r}")
    if only_loc:
        _logger.info("PLASTYFIRE_ONLY_LOC=%s — silencing the other population "
                     "(theta pushed to %.0f)", only_loc, THETA_SILENCE)
    n_silenced = 0

    for (syn_id, synapse), gid in zip(cell.synapses.items(), global_ids):
        cp = c_pre.get(gid, 0.0)
        cq = c_post.get(gid, 0.0)
        is_apical = bool(sec_types[gid] == 3)

        if only_loc is not None:
            keep = is_apical if only_loc == "apical" else (not is_apical)
            if not keep:
                synapse.hsynapse.theta_d_GB = THETA_SILENCE
                synapse.hsynapse.theta_p_GB = THETA_SILENCE
                n_silenced += 1
                continue

        if is_apical:
            td = fit_params.get("a20", 1.0) * cp + fit_params.get("a21", 1.0) * cq
            tp = fit_params.get("a30", 1.0) * cp + fit_params.get("a31", 1.0) * cq
        else:
            td = fit_params.get("a00", 1.0) * cp + fit_params.get("a01", 1.0) * cq
            tp = fit_params.get("a10", 1.0) * cp + fit_params.get("a11", 1.0) * cq

        synapse.hsynapse.theta_d_GB = td if td > 0 else -1.0
        synapse.hsynapse.theta_p_GB = tp if tp > 0 else -1.0
        _logger.info(
            "  syn gid=%d  loc=%s  cp=%.4e  cq=%.4e  theta_d=%.4e  theta_p=%.4e",
            gid, "apical" if is_apical else "basal", cp, cq, td, tp,
        )

    if only_loc:
        _logger.info("PLASTYFIRE_ONLY_LOC=%s: %d/%d synapses silenced, %d active",
                     only_loc, n_silenced, len(global_ids),
                     len(global_ids) - n_silenced)


# ---------------------------------------------------------------------------
# Subprocess
# ---------------------------------------------------------------------------

def _run_prefire_edges_process(conn, workdir, fit_params, edges_h5,
                                pre_gid, post_gid, t_end,
                                fastforward, node_pop, edge_pop, fixhp,
                                traces_output_dir=None,
                                bluecellulab_output_dir=None,
                                cpre_cpost_cache=None,
                                circuit_config=None,
                                bcl_subdir="bluecellulab_results",
                                extracellular_calcium=None,
                                lean=False,
                                force=False):
    """
    Runs the prefire simulation inside a subprocess.
    Uses prefire_simulation_config.json + prefire_prespikes.h5.
    Injects theta_d_GB / theta_p_GB from edges.h5.
    Writes out/rho.h5 for downstream basis extrapolation.
    """
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(os.path.join(workdir, "edges_prefire.log")),
        ],
    )
    _logger = logging.getLogger(__name__)
    t0 = time.time()

    bluecellulab.neuron.load_mechanisms(MECHANISMS_PATH)

    try:
        sim_config      = os.path.join(workdir, "prefire_simulation_config.json")
        pre_spikes_path = os.path.join(workdir, "prefire_prespikes.h5")

        # bluecellulab does not support 'synapse' type reports; strip them
        # so configure_all_reports doesn't raise NotImplementedError.
        # We collect rho manually, so no information is lost.
        with open(sim_config) as f:
            cfg = json.load(f)
        # Fix hardcoded paths that may point to an old/moved location
        node_sets_local = os.path.join(workdir, "node_sets.json")
        if os.path.exists(node_sets_local):
            cfg["node_sets_file"] = node_sets_local
        if circuit_config:
            cfg["network"] = os.path.abspath(circuit_config)
        # Extracellular calcium has to be set in BOTH places or the two halves of
        # the model disagree: conditions.extracellular_calcium drives bluecellulab's
        # Ca-dependent release-probability scaling, while mechanisms.GluSynapse.cao_CR
        # becomes the cao_CR_GluSynapse HOC global that sets the VDCC/NMDA driving
        # force. Setting only one silently mixes two calcium levels.
        if extracellular_calcium is not None:
            _ca = float(extracellular_calcium)
            cfg.setdefault("conditions", {})["extracellular_calcium"] = _ca
            cfg["conditions"].setdefault("mechanisms", {}).setdefault("GluSynapse", {})["cao_CR"] = _ca
            _logger.info("extracellular calcium overridden -> %.4g mM "
                         "(conditions.extracellular_calcium and GluSynapse.cao_CR)", _ca)
        out_dir_local = os.path.join(workdir, "out")
        os.makedirs(out_dir_local, exist_ok=True)
        cfg.setdefault("output", {})["output_dir"] = out_dir_local
        filtered_reports = {
            k: v for k, v in cfg.get("reports", {}).items()
            if v.get("type") != "synapse"
        }
        if filtered_reports:
            cfg["reports"] = filtered_reports
        else:
            cfg.pop("reports", None)  # bluecellulab rejects empty reports dict
        # Save GluSynapse params before stripping (applied manually after instantiation)
        glusyn_conds = dict(cfg.get("conditions", {}).get("mechanisms", {}).get("GluSynapse", {}))
        # BCL calls set_global_condition_parameters() before mechanisms are compiled,
        # so GluSynapse HOC globals don't exist yet → LookupError for init_depleted/minis_single_vesicle.
        # Strip the whole GluSynapse block; we apply all params manually after instantiate_gids().
        cfg.get("conditions", {}).get("mechanisms", {}).pop("GluSynapse", None)
        tmp = tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", dir=tempfile.gettempdir(), delete=False
        )
        json.dump(cfg, tmp)
        tmp.close()
        sim_config_eff = tmp.name

        synapse_seed = int(cfg.get("run", {}).get("synapse_seed") or 0)

        _logger.info("Loading prefire circuit: %s", sim_config)
        sim        = bluecellulab.CircuitSimulation(sim_config_eff, base_seed=cfg["run"]["random_seed"])
        pre_spikes = SpikeReader(pre_spikes_path)[node_pop].get_dict()["timestamps"]

        # bluecellulab auto-loads rho0_GB, Use_d/p, gmax_d/p from dhuruva_circuit_config.json
        sim.instantiate_gids(
            [(node_pop, post_gid)],
            add_synapses=True, add_minis=False, add_pulse_stimuli=True,
            intersect_pre_gids=[(node_pop, pre_gid)],
            pre_spike_trains={(node_pop, pre_gid): pre_spikes},
        )
        cell = sim.cells[(node_pop, post_gid)]

        if fixhp:
            for sec in cell.somatic + cell.axonal:
                sec.uninsert("SK_E2")

        # Apply ALL GluSynapse global params from conditions.mechanisms.
        # BCL's set_global_condition_parameters only covers cao_CR (from extracellular_calcium),
        # but we stripped that block entirely to avoid a pre-compile LookupError.
        # Now that mechanisms are loaded, set everything explicitly.
        _GLUSYN_HOC_MAP = {
            "tau_effca_GB":         "tau_effca_GB_GluSynapse",
            "gamma_d_GB":           "gamma_d_GB_GluSynapse",
            "gamma_p_GB":           "gamma_p_GB_GluSynapse",
            "init_depleted":        "init_depleted_GluSynapse",
            "cao_CR":               "cao_CR_GluSynapse",
            "minis_single_vesicle": "minis_single_vesicle_GluSynapse",
        }
        for param, hoc_attr in _GLUSYN_HOC_MAP.items():
            if param in glusyn_conds:
                val = float(glusyn_conds[param])
                try:
                    setattr(bluecellulab.neuron.h, hoc_attr, val)
                    _logger.info("Set %s = %g", hoc_attr, val)
                except Exception as e:
                    _logger.warning("Could not set %s: %s", hoc_attr, e)

        # Apply fit_params AFTER glusyn_conds so fitted values override config defaults.
        if fit_params:
            _set_global_params(fit_params)
            _logger.info("fit_params applied (override glusyn_conds defaults)")

        # Map local synapse indices → global circuit IDs
        syn_idx    = [syn_id[1] for syn_id in cell.synapses]
        df         = _map_syn_idx(sim_config, post_gid, syn_idx, edge_pop)
        global_ids = [int(df.loc[df["local_syn_idx"] == i].index[0]) for i in syn_idx]

        # Inject theta_d_GB / theta_p_GB.
        # If fit_params contains a00-a31 coefficients, compute c_pre/c_post via
        # proper mini-sims (like simulator.py's c_pre_finder/c_post_finder) and
        # override thresholds per synapse. Otherwise fall back to edges.h5 values.
        use_a_params = fit_params and bool(_A_PARAM_KEYS & set(fit_params.keys()))
        if use_a_params:
            # Try to load c_pre/c_post from precomputed cache first
            cache_hit = False
            if cpre_cpost_cache and os.path.exists(cpre_cpost_cache):
                with open(cpre_cpost_cache, "rb") as _f:
                    _cache = pickle.load(_f)
                _key = (pre_gid, post_gid)
                if _key in _cache:
                    c_pre  = _cache[_key]["c_pre"]
                    c_post = _cache[_key]["c_post"]
                    cache_hit = True
                    _logger.info("Loaded c_pre/c_post from cache for pair (%d, %d)", pre_gid, post_gid)
            if not cache_hit:
                _logger.info("a-params detected — computing c_pre/c_post via mini-sims")
                c_pre, c_post, _ = _find_cpre_cpost(
                    sim_config_eff, cfg, fit_params,
                    pre_gid, post_gid, node_pop, edge_pop, fixhp, _logger,
                    sim_config_orig=sim_config,
                )
            _logger.info("Overriding theta_d/theta_p with a-param formula")
            _apply_theta_from_a_params(cell, global_ids, fit_params, c_pre, c_post,
                                       edges_h5, edge_pop, _logger)
        else:
            _logger.info("Injecting theta_d/theta_p from edges.h5 for %d synapses", len(global_ids))
            thresholds = _load_thresholds(edges_h5, edge_pop, global_ids)
            _apply_thresholds(cell, global_ids, thresholds)

        # Per-synapse RNG seeding to match Neurodamus Random123 pattern.
        # ND seeds each synapse as (tgid, 100000+local_idx, synapse_seed+200).
        _logger.info("Setting per-synapse RNG seeds (synapse_seed=%d)", synapse_seed)
        for syn_id, synapse in cell.synapses.items():
            local_idx = int(syn_id[1])
            s1 = post_gid
            s2 = 100000 + local_idx
            s3 = 200 + synapse_seed
            synapse.randseed1 = s1
            synapse.randseed2 = s2
            synapse.randseed3 = s3
            synapse.hsynapse.setRNG(s1, s2, s3)
            _logger.info("  syn %s: setRNG(%d, %d, %d)", syn_id, s1, s2, s3)

        # Record synapse variables per synapse for full time-series output.
        # cai_CR    : spine calcium (mM) — the primary signal for JAX CICR fitting
        # shaft_cai : dendritic shaft calcium from the local section (mM)
        # ica_NMDA  : NMDA calcium current (nA) — CICR model input
        # ica_VDCC  : VDCC calcium current (nA) — CICR model input
        # rho_GB    : synaptic weight state variable
        rho_vecs      = {}
        cai_vecs      = {}
        shaft_vecs    = {}
        ica_nmda_vecs = {}
        ica_vdcc_vecs = {}
        for syn_id, synapse in cell.synapses.items():
            h_syn = synapse.hsynapse

            v = bluecellulab.neuron.h.Vector()
            v.record(h_syn._ref_rho_GB)
            rho_vecs[syn_id] = v

            v = bluecellulab.neuron.h.Vector()
            v.record(h_syn._ref_cai_CR)
            cai_vecs[syn_id] = v

            v = bluecellulab.neuron.h.Vector()
            v.record(h_syn._ref_shaft_cai)
            shaft_vecs[syn_id] = v

            v = bluecellulab.neuron.h.Vector()
            v.record(h_syn._ref_ica_NMDA)
            ica_nmda_vecs[syn_id] = v

            v = bluecellulab.neuron.h.Vector()
            v.record(h_syn._ref_ica_VDCC)
            ica_vdcc_vecs[syn_id] = v

        # Collect initial rho from bluecellulab-loaded rho0_GB
        initial_rho = [s.hsynapse.rho0_GB for s in cell.synapses.values()]

        # Run
        if fastforward is not None and fastforward < t_end:
            _logger.info("Fast-forwarding to %.0f ms", fastforward)
            sim.run(fastforward, cvode=True)
            # Snap rho and use as the initial state for recording
            initial_rho = []
            for synapse in cell.synapses.values():
                h       = synapse.hsynapse
                snapped = 1.0 if h.rho_GB >= 0.5 else 0.0
                h.rho_GB = snapped
                if snapped == 1.0:
                    h.Use = h.Use_p;  h.Use_GB = h.Use_p
                    h.gmax_AMPA = h.gmax_p_AMPA; h.gmax0_AMPA = h.gmax_p_AMPA
                else:
                    h.Use = h.Use_d;  h.Use_GB = h.Use_d
                    h.gmax_AMPA = h.gmax_d_AMPA; h.gmax0_AMPA = h.gmax_d_AMPA
                initial_rho.append(snapped)
            bluecellulab.neuron.h.cvode_active(1)
            bluecellulab.neuron.h.continuerun(t_end)
        else:
            if fastforward is not None:
                _logger.warning(
                    "fastforward=%.0f >= tstop=%.0f — ignoring fastforward. "
                    "rho0_GB from edges.h5 is already the equilibrated state.",
                    fastforward, t_end,
                )
            sim.run(t_end, cvode=True)

        _logger.info("Simulation done in %.1f s", time.time() - t0)

        # Collect final rho
        final_rho = [s.hsynapse.rho_GB for s in cell.synapses.values()]

        # Determine output directory for rho.h5 / rho_timeseries.npy.
        # If bluecellulab_output_dir is given, write to:
        #   <bluecellulab_output_dir>/<pair_name>/<protocol>/bluecellulab_results/
        # Otherwise fall back to bluecellulab_results/ inside the workdir.
        pair_name = os.path.basename(os.path.dirname(workdir))
        protocol  = os.path.basename(workdir)
        if bluecellulab_output_dir is not None:
            bcl_out = os.path.join(bluecellulab_output_dir, pair_name, protocol, bcl_subdir)
        else:
            bcl_out = os.path.join(workdir, bcl_subdir)
        os.makedirs(bcl_out, exist_ok=True)

        # Write rho.h5 for plot_refitting_stdp_basis.py
        if not lean:
            rho_h5_path = os.path.join(bcl_out, "rho.h5")
            _write_rho_h5(rho_h5_path, node_pop, global_ids, initial_rho, final_rho)
            _logger.info("Wrote %s", rho_h5_path)

        # Write full rho time-series for BCL vs ND comparison
        t = np.array(sim.get_time())
        # lean skips every per-workdir artefact: rho.h5 (7 kB), rho_timeseries.npy
        # (~3.5 MB) and simulation_traces.pkl (~275 MB). rho.h5 is skipped not for
        # size but because its path is fixed per workdir -- a calcium sweep would
        # otherwise overwrite the 2.0 mM rho.h5 three times over. The returned
        # results dict already carries global_ids/initial_rho/final_rho, which is
        # everything rho.h5 holds.
        if lean:
            _logger.info("lean=True: skipping rho.h5, rho_timeseries.npy and simulation_traces.pkl")
        else:
            rho_ts = np.column_stack([t] + [np.array(rho_vecs[sid]) for sid in cell.synapses])
            ts_path = os.path.join(bcl_out, "rho_timeseries.npy")
            np.save(ts_path, rho_ts)
            np.save(ts_path.replace("rho_timeseries.npy", "rho_global_ids.npy"),
                    np.array(global_ids, dtype=np.int64))
            _logger.info("Wrote %s  shape=%s", ts_path, rho_ts.shape)

            # Write simulation_traces.pkl — JAX-ready calcium traces.
            # Interpolated to a uniform time grid (RECORD_DT ms) so JAX can treat
            # each trace as a fixed-length float32 array without variable step logic.
            # Format matches plastyfitting/cicr_common.py _load_pkl() "new format":
            #   cai_CR    : dict {global_id (int): float32 array (T,)}
            #   shaft_cai : dict {global_id (int): float32 array (T,)}
            #   ica_NMDA  : dict {global_id (int): float32 array (T,)}
            #   ica_VDCC  : dict {global_id (int): float32 array (T,)}
            #   rho_GB    : dict {global_id (int): float32 array (T,)}
            RECORD_DT = 0.025  # ms — matches default NEURON dt; fine enough for CICR dynamics
            t_uniform = np.arange(0.0, t[-1] + RECORD_DT, RECORD_DT, dtype=np.float32)

            def _interp_dict(vec_dict, syn_order, gids, t_raw, t_uni):
                """Interpolate NEURON Vector recordings onto a uniform grid, keyed by global_id."""
                out = {}
                for syn_id, gid in zip(syn_order, gids):
                    raw = np.array(vec_dict[syn_id], dtype=np.float64)
                    out[gid] = np.interp(t_uni, t_raw, raw).astype(np.float32)
                return out

            syn_order = list(cell.synapses.keys())
            traces_pkl = {
                "t":         t_uniform,
                "cai_CR":    _interp_dict(cai_vecs,      syn_order, global_ids, t, t_uniform),
                "shaft_cai": _interp_dict(shaft_vecs,    syn_order, global_ids, t, t_uniform),
                "ica_NMDA":  _interp_dict(ica_nmda_vecs, syn_order, global_ids, t, t_uniform),
                "ica_VDCC":  _interp_dict(ica_vdcc_vecs, syn_order, global_ids, t, t_uniform),
                "rho_GB":    _interp_dict(rho_vecs,      syn_order, global_ids, t, t_uniform),
                "prespikes": np.array(pre_spikes, dtype=np.float32),
                "global_ids": global_ids,
            }

            # Determine output path for simulation_traces.pkl.
            # If traces_output_dir is given, write to:
            #   <traces_output_dir>/<pair_name>/<protocol>/simulation_traces.pkl
            # (the structure plastyfitting/cicr_common.py _load_pkl expects).
            # Otherwise fall back to bluecellulab_results/ inside the workdir.
            if traces_output_dir is not None:
                traces_dir = os.path.join(traces_output_dir, pair_name, protocol)
            else:
                traces_dir = bcl_out

            os.makedirs(traces_dir, exist_ok=True)
            traces_path = os.path.join(traces_dir, "simulation_traces.pkl")
            with open(traces_path, "wb") as _f:
                pickle.dump(traces_pkl, _f, protocol=-1)
            _logger.info(
                "Wrote %s  n_syn=%d  T=%d  dt=%.3f ms",
                traces_path, len(global_ids), len(t_uniform), RECORD_DT,
            )

        v = np.array(sim.get_voltage_trace((node_pop, post_gid)))
        post_spikes = _get_spikes(t, v)

        # Induction spike-count guardrail.
        #
        # The pairing protocol only induces plasticity if the postsynaptic cell
        # actually fires on every repetition. The pulse amplitude comes from the
        # single_cells calibration, which only ever verified that the cell fires
        # nspikes APs *in isolation over ~1.4 s* — not that it still fires them
        # during a 40 s train with synaptic input. Those are different questions,
        # and the second one was never asked: runs were silently completing with
        # zero post spikes and producing no plasticity, which is indistinguishable
        # from "these parameters give no plasticity" once the trace is discarded.
        #
        # Expected count is derived from the config rather than hard-coded, so it
        # follows nspikes / nreps / T if the protocol changes: each "pulse" input
        # fires duration_ms * frequency_Hz / 1000 times.
        n_expected = 0
        for _inp in cfg.get("inputs", {}).values():
            if _inp.get("module") == "pulse":
                n_expected += int(round(_inp["duration"] * _inp["frequency"] / 1000.0))
        n_post = int(len(post_spikes))
        guardrail_forced = False
        if n_expected and n_post != n_expected:
            msg = ("induction delivered %d post spikes, expected %d (%.0f%%) — "
                   "pair %d->%d, Vm peak %.1f mV. The calibrated pulse amplitude no "
                   "longer fires this cell as calibrated during the pairing train, so "
                   "any rho change from this run is suspect."
                   % (n_post, n_expected, 100.0 * n_post / n_expected,
                      pre_gid, post_gid, float(v.max()) if v.size else float("nan")))
            if not force:
                raise RuntimeError(msg)
            # force=True: keep the run, but mark it. `guardrail_forced` lands in the
            # output pkl so downstream consumers can tell a forced result from a
            # clean one — the whole point of the check was that a silently kept
            # bad induction is indistinguishable from "these params give no
            # plasticity".
            guardrail_forced = True
            _logger.warning("FORCED past induction guardrail: %s", msg)
        _logger.info("induction: %d/%d post spikes%s", n_post, n_expected,
                     "  (FORCED)" if guardrail_forced else "")

        results = {
            "n_post_spikes":     n_post,
            "n_post_spikes_exp": n_expected,
            "guardrail_forced":  guardrail_forced,
            "global_ids":  global_ids,
            "initial_rho": initial_rho,
            "final_rho":   final_rho,
        }
        if not lean:
            # t/v/spike trains are ~500 kB per workdir; rho alone is ~1 kB.
            results.update({"t": t, "v": v,
                            "prespikes": pre_spikes, "postspikes": post_spikes})
        conn.send(results)

    except Exception as e:
        _logger.error("Subprocess failed: %s", e, exc_info=True)
        raise
    finally:
        conn.close()
        try:
            os.unlink(sim_config_eff)
        except Exception:
            pass


# ---------------------------------------------------------------------------
# Public API — c_pre/c_post cache builder
# ---------------------------------------------------------------------------

def _cpre_cpost_only_process(conn, workdir, fit_params, node_pop, edge_pop, fixhp,
                              circuit_config=None, extracellular_calcium=None):
    """
    Subprocess target: run ONLY the two mini-sims (c_pre, c_post) and send
    the results back through *conn*.  Does NOT run the full induction.
    """
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(os.path.join(workdir, "cpre_cpost_compute.log")),
        ],
    )
    _logger = logging.getLogger(__name__)
    bluecellulab.neuron.load_mechanisms(MECHANISMS_PATH)

    try:
        sim_config = os.path.join(workdir, "prefire_simulation_config.json")
        with open(sim_config) as f:
            cfg = json.load(f)
        node_sets_local = os.path.abspath(os.path.join(workdir, "node_sets.json"))
        if os.path.exists(node_sets_local):
            cfg["node_sets_file"] = node_sets_local
        if circuit_config is not None:
            cfg["network"] = os.path.abspath(circuit_config)
        # c_pre/c_post are calcium-dependent: cao_CR enters GluSynapse.mod both as
        # Pf_NMDA = (4*cao_CR)/(4*cao_CR + (1/1.38)*120)*0.6 and via nernst() in
        # Eca_syn, so peak effcai_GB -- and hence theta_d/theta_p -- move with Ca.
        # A cache built at one calcium level is INVALID at another.
        # Unlike the induction path, this one pops GluSynapse conditions and never
        # re-applies them as HOC globals; that is fine here because bluecellulab's
        # set_global_condition_parameters sets cao_CR_GluSynapse from
        # conditions.extracellular_calcium alone (neuron_globals.py).
        if extracellular_calcium is not None:
            _ca = float(extracellular_calcium)
            cfg.setdefault("conditions", {})["extracellular_calcium"] = _ca
            _logger.info("cpre/cpost extracellular calcium -> %.4g mM", _ca)
        # Strip synapse reports (bcl raises NotImplementedError for them)
        cfg.pop("reports", None)
        # Save and strip GluSynapse conditions to apply manually post-load
        glusyn_conds = dict(cfg.get("conditions", {}).get("mechanisms", {}).get("GluSynapse", {}))
        cfg.get("conditions", {}).get("mechanisms", {}).pop("GluSynapse", None)
        tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".json", dir=tempfile.gettempdir(), delete=False)
        json.dump(cfg, tmp)
        tmp.close()
        sim_config_eff = tmp.name

        try:
            _logger.info("Loading circuit for mini-sims: %s", sim_config)
            # Use bluepysnap.Simulation to read node_sets (same as runconnectedpair_prefire_from_edges)
            _snap_tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".json", dir=tempfile.gettempdir(), delete=False)
            json.dump(cfg, _snap_tmp)
            _snap_tmp.close()
            try:
                _snap = Simulation(_snap_tmp.name)
                pre_gid  = int(_snap.node_sets.content["precell"]["node_id"][0])
                post_gid = int(_snap.node_sets.content["postcell"]["node_id"][0])
            finally:
                os.unlink(_snap_tmp.name)

            c_pre, c_post, _df = _find_cpre_cpost(
                sim_config_eff, cfg, fit_params,
                pre_gid, post_gid, node_pop, edge_pop, fixhp, _logger,
                sim_config_orig=sim_config,
            )
            _attrs = getattr(_df, "attrs", {}) or {}
            conn.send({"pre_gid": pre_gid, "post_gid": post_gid,
                       "c_pre": c_pre, "c_post": c_post,
                       "n_post_spikes": _attrs.get("n_post_spikes"),
                       "n_syn": _attrs.get("n_syn"),
                       "n_syn_at_floor": _attrs.get("n_syn_at_floor"),
                       "pulse_amp": _attrs.get("pulse_amp")})
        finally:
            try:
                os.unlink(sim_config_eff)
            except OSError:
                pass
    except Exception as exc:
        _logger.exception("cpre_cpost_only_process failed")
        # Send the reason back rather than a bare None: the parent classifies
        # no-spike failures separately from crashes, and cannot tell them apart
        # from the generic "subprocess returned None" wrapper alone.
        try:
            conn.send({"error": str(exc), "error_type": type(exc).__name__})
        except Exception:
            conn.send(None)
    finally:
        conn.close()


def compute_cpre_cpost_for_workdir(workdir, fit_params=None,
                                    node_pop=NODE_POP_DEFAULT,
                                    edge_pop=EDGE_POP_DEFAULT,
                                    fixhp=True,
                                    circuit_config=None,
                                    extracellular_calcium=None):
    """
    Run the two mini-sims for *workdir* in an isolated subprocess and return
    ``{"pre_gid": …, "post_gid": …, "c_pre": {…}, "c_post": {…}}`` or None on failure.

    Running in a subprocess keeps the NEURON/HOC state isolated so that
    multiple pairs can be processed sequentially without state leakage.
    """
    parent_conn, child_conn = multiprocessing.Pipe(duplex=False)
    proc = multiprocessing.Process(
        target=_cpre_cpost_only_process,
        args=(child_conn, workdir, fit_params, node_pop, edge_pop, fixhp),
        kwargs={"circuit_config": circuit_config,
                "extracellular_calcium": extracellular_calcium},
    )
    proc.start()
    child_conn.close()
    try:
        result = parent_conn.recv()
    except EOFError:
        result = None
    finally:
        parent_conn.close()
    proc.join()
    if proc.exitcode != 0:
        raise RuntimeError(
            f"cpre_cpost subprocess failed (exit code {proc.exitcode}) for {workdir}"
        )
    return result


# ---------------------------------------------------------------------------
# Public API — induction simulation
# ---------------------------------------------------------------------------

def runconnectedpair_prefire_from_edges(workdir, fit_params=None,
                                        edges_h5=EDGES_H5_DEFAULT,
                                        fastforward=None,
                                        node_pop=NODE_POP_DEFAULT,
                                        edge_pop=EDGE_POP_DEFAULT,
                                        fixhp=True,
                                        traces_output_dir=None,
                                        bluecellulab_output_dir=None,
                                        cpre_cpost_cache=None,
                                        circuit_config=None,
                                        bcl_subdir="bluecellulab_results",
                                        extracellular_calcium=None,
                                        lean=False,
                                        force=False):
    """
    Run the prefire STDP induction simulation using dhuruva_modified_edges.h5.

    Uses prefire_simulation_config.json (pairing protocol only, no test pulses).
    theta_d_GB / theta_p_GB injected from edges.h5; all other synapse params
    loaded automatically by bluecellulab.

    If bluecellulab_output_dir is given, rho.h5 / rho_timeseries.npy are written to:
        <bluecellulab_output_dir>/<pair_name>/<protocol>/bluecellulab_results/
    Otherwise they are written to bluecellulab_results/ inside the workdir.

    If traces_output_dir is given, simulation_traces.pkl is written to:
        <traces_output_dir>/<pair_name>/<protocol>/simulation_traces.pkl
    so that plastyfitting/cicr_common.py can load it directly.
    Otherwise it falls back to bluecellulab_results/ inside the workdir.

    force=True downgrades the induction spike-count guardrail from a raise to a
    warning: the run completes and the pkl is written even when the postsynaptic
    cell fired the wrong number of APs during the pairing train. The result then
    carries guardrail_forced=True. Use it to recover near-misses (51-53 spikes
    out of 50, where the EPSP summates with the calibrated pulse at short |dt|),
    not to paper over a cell that never fired.
    """
    sim_config = os.path.join(workdir, "prefire_simulation_config.json")
    # Patch hardcoded paths (config may reference old/moved locations)
    with open(sim_config) as f:
        _cfg_snap = json.load(f)
    node_sets_local = os.path.join(workdir, "node_sets.json")
    if os.path.exists(node_sets_local):
        _cfg_snap["node_sets_file"] = node_sets_local
    if circuit_config:
        _cfg_snap["network"] = os.path.abspath(circuit_config)
    _snap_tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".json", dir=tempfile.gettempdir(), delete=False
    )
    json.dump(_cfg_snap, _snap_tmp)
    _snap_tmp.close()
    try:
        snap = Simulation(_snap_tmp.name)
        pre_gid  = snap.node_sets.content["precell"]["node_id"][0]
        post_gid = snap.node_sets.content["postcell"]["node_id"][0]
        t_end    = snap.config["run"]["tstop"]
    finally:
        os.unlink(_snap_tmp.name)

    parent_conn, child_conn = multiprocessing.Pipe(duplex=False)
    proc = multiprocessing.Process(
        target=_run_prefire_edges_process,
        args=(child_conn, workdir, fit_params, edges_h5,
              pre_gid, post_gid, t_end,
              fastforward, node_pop, edge_pop, fixhp,
              traces_output_dir, bluecellulab_output_dir, cpre_cpost_cache),
        kwargs={"circuit_config": circuit_config, "bcl_subdir": bcl_subdir,
                "extracellular_calcium": extracellular_calcium, "lean": lean,
                "force": force},
    )
    proc.start()
    child_conn.close()
    try:
        results = parent_conn.recv()
    except EOFError:
        results = None
    finally:
        parent_conn.close()
    proc.join()
    if proc.exitcode != 0:
        raise RuntimeError(f"Subprocess failed (exit code {proc.exitcode})")
    return results
