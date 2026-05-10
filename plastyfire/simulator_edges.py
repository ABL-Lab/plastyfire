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
    _get_spikes, _map_syn_idx, _set_global_params,
)

logger = logging.getLogger(__name__)

EDGES_H5_DEFAULT  = "/project/ctb-emuller/dhuruva/plastyfire/data/dhuruva_modified_edges.h5"
EDGE_POP_DEFAULT  = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"
NODE_POP_DEFAULT  = "S1nonbarrel_neurons"
MECHANISMS_PATH   = "/project/ctb-emuller/dhuruva/DEES_cell_packages/"


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
# Subprocess
# ---------------------------------------------------------------------------

def _run_prefire_edges_process(conn, workdir, fit_params, edges_h5,
                                pre_gid, post_gid, t_end,
                                fastforward, node_pop, edge_pop, fixhp,
                                traces_output_dir=None):
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
        filtered_reports = {
            k: v for k, v in cfg.get("reports", {}).items()
            if v.get("type") != "synapse"
        }
        cfg["reports"] = filtered_reports
        # Save GluSynapse params before stripping (applied manually after instantiation)
        glusyn_conds = dict(cfg.get("conditions", {}).get("mechanisms", {}).get("GluSynapse", {}))
        # BCL calls set_global_condition_parameters() before mechanisms are compiled,
        # so GluSynapse HOC globals don't exist yet → LookupError for init_depleted/minis_single_vesicle.
        # Strip the whole GluSynapse block; we apply all params manually after instantiate_gids().
        cfg.get("conditions", {}).get("mechanisms", {}).pop("GluSynapse", None)
        tmp = tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", dir=workdir, delete=False
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

        if fit_params:
            _set_global_params(fit_params)

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

        # Map local synapse indices → global circuit IDs
        syn_idx    = [syn_id[1] for syn_id in cell.synapses]
        df         = _map_syn_idx(sim_config, post_gid, syn_idx, edge_pop)
        global_ids = [int(df.loc[df["local_syn_idx"] == i].index[0]) for i in syn_idx]

        # Inject theta_d_GB / theta_p_GB (not auto-loaded; name mismatch in edges.h5)
        _logger.info("Injecting theta_d/theta_p for %d synapses", len(global_ids))
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

        # Write bluecellulab_results/rho.h5 for plot_refitting_stdp_basis.py
        rho_h5_path = os.path.join(workdir, "bluecellulab_results", "rho.h5")
        _write_rho_h5(rho_h5_path, node_pop, global_ids, initial_rho, final_rho)
        _logger.info("Wrote %s", rho_h5_path)

        # Write full rho time-series for BCL vs ND comparison
        t = np.array(sim.get_time())
        rho_ts = np.column_stack([t] + [np.array(rho_vecs[sid]) for sid in cell.synapses])
        ts_path = os.path.join(workdir, "bluecellulab_results", "rho_timeseries.npy")
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
            pair_name = os.path.basename(os.path.dirname(workdir))
            protocol  = os.path.basename(workdir)
            traces_dir = os.path.join(traces_output_dir, pair_name, protocol)
        else:
            traces_dir = os.path.join(workdir, "bluecellulab_results")

        os.makedirs(traces_dir, exist_ok=True)
        traces_path = os.path.join(traces_dir, "simulation_traces.pkl")
        with open(traces_path, "wb") as _f:
            pickle.dump(traces_pkl, _f, protocol=-1)
        _logger.info(
            "Wrote %s  n_syn=%d  T=%d  dt=%.3f ms",
            traces_path, len(global_ids), len(t_uniform), RECORD_DT,
        )

        v = np.array(sim.get_voltage_trace((node_pop, post_gid)))
        results = {
            "t":           t,
            "v":           v,
            "prespikes":   pre_spikes,
            "postspikes":  _get_spikes(t, v),
            "global_ids":  global_ids,
            "initial_rho": initial_rho,
            "final_rho":   final_rho,
        }
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
# Public API
# ---------------------------------------------------------------------------

def runconnectedpair_prefire_from_edges(workdir, fit_params=None,
                                        edges_h5=EDGES_H5_DEFAULT,
                                        fastforward=None,
                                        node_pop=NODE_POP_DEFAULT,
                                        edge_pop=EDGE_POP_DEFAULT,
                                        fixhp=True,
                                        traces_output_dir=None):
    """
    Run the prefire STDP induction simulation using dhuruva_modified_edges.h5.

    Uses prefire_simulation_config.json (pairing protocol only, no test pulses).
    theta_d_GB / theta_p_GB injected from edges.h5; all other synapse params
    loaded automatically by bluecellulab.

    Writes bluecellulab_results/rho.h5 (SONATA format, compatible with plot_refitting_stdp_basis.py)
    and returns a results dict with initial_rho, final_rho, v, t.

    If traces_output_dir is given, simulation_traces.pkl is written to:
        <traces_output_dir>/<pair_name>/<protocol>/simulation_traces.pkl
    so that plastyfitting/cicr_common.py can load it directly.
    Otherwise it falls back to bluecellulab_results/ inside the workdir.
    """
    sim_config = os.path.join(workdir, "prefire_simulation_config.json")
    snap       = Simulation(sim_config)
    pre_gid    = snap.node_sets.content["precell"]["node_id"][0]
    post_gid   = snap.node_sets.content["postcell"]["node_id"][0]
    t_end      = snap.config["run"]["tstop"]

    parent_conn, child_conn = multiprocessing.Pipe(duplex=False)
    proc = multiprocessing.Process(
        target=_run_prefire_edges_process,
        args=(child_conn, workdir, fit_params, edges_h5,
              pre_gid, post_gid, t_end,
              fastforward, node_pop, edge_pop, fixhp,
              traces_output_dir),
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
