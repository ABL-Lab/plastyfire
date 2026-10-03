#!/usr/bin/env python3
"""
Compute EPSP basis results for one pre-post pair using dhuruva_modified_edges.h5.

Reads gmax_d/p_AMPA and Use_d/p directly from bluecellulab (auto-loaded from
dhuruva_circuit_config.json -> dhuruva_modified_edges.h5).  No ParamsGenerator step.

For each rho configuration (all-depressed, each singleton-potentiated, all-potentiated),
runs N_TRIALS test-pulse trials and records the mean EPSP amplitude.

Three EPSP measures per test pulse, each relative to v at the pre spike time:
  mean / std          peak: raw maximum within EPSP_WINDOW_MS (what the analytical fit uses)
  mean_pavg / std_pavg  peak averaged over +-1.5 ms around the peak (Tazerart lab measure)
  mean_area / std_area  area under the EPSP over EPSP_WINDOW_MS (mV ms, Tazerart lab measure)
--waveforms NPZ also saves each config's mean aligned EPSP (trial 0) for plotting.

Usage (from the simulation workdir or with --sim-config):
    python run_basis_pair_edges.py \\
        --pre-gid 180164 --post-gid 197248 \\
        --sim-config /path/to/simulation_config.json \\
        --output-csv /path/to/basis_results_edges/basis_180164_197248.csv
"""

import argparse
import json
import logging
import multiprocessing
import os
import tempfile
import traceback

import numpy as np
import pandas as pd
import bluecellulab
from libsonata import SpikeReader

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(name)s %(levelname)s %(message)s",
)
logger = logging.getLogger(__name__)

NODE_POP        = "S1nonbarrel_neurons"
EDGE_POP        = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"
MECHANISMS_PATH = "/project/rrg-emuller/dhuruva/DEES_cell_packages/"
C01_DURATION_MS = 2.0 * 60.0 * 1000.0   # 2 min in ms
N_EPSP          = 30
EPSP_WINDOW_MS  = 100.0
SPIKE_THR_MV    = -30.0
PEAK_AVG_MS     = 1.5       # +- window around the peak (Tazerart lab)
DT_RESAMPLE     = 0.025     # uniform grid for the cvode trace (ms)
WAVE_T          = np.arange(-5.0, EPSP_WINDOW_MS + 1e-9, 0.1)   # saved waveform grid (ms from pre spike)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _strip_synapse_reports(sim_config_path, circuit_config=None):
    """
    Write a temp copy of the config with:
    - synapse-type reports removed (bluecellulab doesn't support them)
    - conditions.mechanisms removed entirely: bluecellulab tries to set every
      entry (init_depleted_GluSynapse, minis_single_vesicle_GluSynapse, …) as
      a HOC global during __init__, before any .mod files are compiled/loaded.
      We set all synapse parameters manually after instantiate_gids, so nothing
      is lost.
    - network overridden if circuit_config is provided (for new ion channels)
    """
    with open(sim_config_path) as f:
        cfg = json.load(f)
    cfg["reports"] = {
        k: v for k, v in cfg.get("reports", {}).items()
        if v.get("type") != "synapse"
    }
    try:
        cfg["conditions"].pop("mechanisms", None)
    except (KeyError, AttributeError):
        pass
    if circuit_config is not None:
        cfg["network"] = os.path.abspath(circuit_config)
    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".json",
        dir=tempfile.gettempdir(), delete=False,
    )
    json.dump(cfg, tmp)
    tmp.close()
    return tmp.name


def _get_epsp_vector(t, v, spikes):
    """-> (peak, peak_avg, area, waveforms) per pre spike; the cvode trace is resampled to DT_RESAMPLE."""
    n = len(spikes)
    peak, pavg, area = np.zeros(n), np.zeros(n), np.zeros(n)
    waves = np.zeros((n, len(WAVE_T)))
    tw = np.arange(0.0, EPSP_WINDOW_MS + 1e-9, DT_RESAMPLE)
    k = int(round(PEAK_AVG_MS / DT_RESAMPLE))
    for i, t_spike in enumerate(spikes):
        if t_spike + EPSP_WINDOW_MS > t[-1]:
            continue
        vw = np.interp(t_spike + tw, t, v)
        if vw.max() > SPIKE_THR_MV:
            raise RuntimeError(f"Postsynaptic spike at t={t_spike:.0f} ms")
        d = vw - vw[0]
        j = int(np.argmax(d))
        peak[i] = d[j]
        pavg[i] = d[max(0, j - k):j + k + 1].mean()
        area[i] = np.trapz(d, tw)
        waves[i] = np.interp(t_spike + WAVE_T, t, v) - vw[0]
    return peak, pavg, area, waves


def _measure_epsp(t, v, pre_spikes):
    """Mean of each measure over the last N_EPSP spikes within C01_DURATION_MS.
    -> dict(epsp, epsp_pavg, epsp_area, wave)"""
    t = np.asarray(t, dtype=np.float64)
    v = np.asarray(v, dtype=np.float64)
    spikes = np.asarray(pre_spikes, dtype=np.float64)
    spikes = spikes[(spikes <= t[-1]) & (spikes <= C01_DURATION_MS)]
    if len(spikes) == 0:
        return dict(epsp=0.0, epsp_pavg=0.0, epsp_area=0.0, wave=np.zeros(len(WAVE_T)))
    peak, pavg, area, waves = _get_epsp_vector(t, v, spikes)
    return dict(epsp=float(np.mean(peak[-N_EPSP:])), epsp_pavg=float(np.mean(pavg[-N_EPSP:])),
                epsp_area=float(np.mean(area[-N_EPSP:])), wave=waves[-N_EPSP:].mean(0))


def _generate_configs(n):
    """all-0, each singleton-1, all-1."""
    configs = [[0] * n]
    for i in range(n):
        c = [0] * n; c[i] = 1
        configs.append(c)
    configs.append([1] * n)
    return [",".join(map(str, c)) for c in configs]


# ---------------------------------------------------------------------------
# Per-trial worker (top-level so it is picklable)
# ---------------------------------------------------------------------------

def _run_trial(args):
    sim_config_eff, sim_config_orig, pre_gid, post_gid, rho_config_str, trial = args
    try:
        import bluecellulab
        from libsonata import SpikeReader
        import numpy as np

        bluecellulab.neuron.load_mechanisms(MECHANISMS_PATH)
        np.random.seed(trial)
        sim = bluecellulab.CircuitSimulation(sim_config_eff, base_seed=trial)

        workdir    = os.path.dirname(os.path.abspath(sim_config_orig))
        pre_spikes = SpikeReader(os.path.join(workdir, "prespikes.h5"))[NODE_POP].get_dict()["timestamps"]

        sim.instantiate_gids(
            [(NODE_POP, post_gid)],
            add_synapses=True, add_minis=False, add_pulse_stimuli=True,
            intersect_pre_gids=[(NODE_POP, pre_gid)],
            pre_spike_trains={(NODE_POP, pre_gid): pre_spikes},
        )
        cell = sim.cells[(NODE_POP, post_gid)]

        for sec in cell.somatic + cell.axonal:
            sec.uninsert("SK_E2")

        rho_vals = [int(x) for x in rho_config_str.split(",")]
        for (syn_id, synapse), rho in zip(cell.synapses.items(), rho_vals):
            h = synapse.hsynapse
            if rho >= 1:
                h.rho0_GB = 1.0;  h.rho_GB = 1.0
                h.Use = h.Use_p;  h.Use_GB = h.Use_p
                h.gmax_AMPA  = h.gmax_p_AMPA
                h.gmax0_AMPA = h.gmax_p_AMPA
            else:
                h.rho0_GB = 0.0;  h.rho_GB = 0.0
                h.Use = h.Use_d;  h.Use_GB = h.Use_d
                h.gmax_AMPA  = h.gmax_d_AMPA
                h.gmax0_AMPA = h.gmax_d_AMPA
            h.theta_d_GB = -1.0
            h.theta_p_GB = -1.0

        bluecellulab.neuron.h.cvode_active(1)
        sim.run(C01_DURATION_MS, cvode=True)

        t = np.array(sim.get_time())
        v = np.array(sim.get_voltage_trace((NODE_POP, post_gid)))
        m = _measure_epsp(t, v, pre_spikes)
        del sim
        return dict(config=rho_config_str, trial=trial, **m)

    except Exception as e:
        logger.error("Trial %d config %s failed: %s", trial, rho_config_str, e)
        traceback.print_exc()
        return dict(config=rho_config_str, trial=trial, epsp=float("nan"), epsp_pavg=float("nan"),
                    epsp_area=float("nan"), wave=None)


# ---------------------------------------------------------------------------
# Count synapses
# ---------------------------------------------------------------------------

def _count_synapses(sim_config_eff, pre_gid, post_gid):
    bluecellulab.neuron.load_mechanisms(MECHANISMS_PATH)
    sim = bluecellulab.CircuitSimulation(sim_config_eff)
    sim.instantiate_gids(
        [(NODE_POP, post_gid)],
        add_synapses=True, add_minis=False, add_pulse_stimuli=False,
        intersect_pre_gids=[(NODE_POP, pre_gid)],
    )
    n = len(sim.cells[(NODE_POP, post_gid)].synapses)
    del sim
    return n


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compute EPSP basis for one pair using dhuruva_modified_edges.h5",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--pre-gid",     type=int, required=True)
    parser.add_argument("--post-gid",    type=int, required=True)
    parser.add_argument("--sim-config",  required=True,
                        help="Path to simulation_config.json (workdir must contain prespikes.h5)")
    parser.add_argument("--output-csv",  required=True)
    parser.add_argument("--num-trials",  type=int, default=5)
    parser.add_argument("--workers",     type=int, default=None)
    parser.add_argument("--circuit-config", default=None,
                        help="Override the 'network' field in sim_config (e.g. for new ion channels)")
    parser.add_argument("--waveforms", default=None, help="npz for each config's mean EPSP waveform (trial 0)")
    parser.add_argument("--configs", default="all", choices=["all", "ends"],
                        help="'ends' = only all-0 and all-1 (quick measure comparisons)")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    sim_config_eff = _strip_synapse_reports(args.sim_config, circuit_config=args.circuit_config)
    try:
        n_syn = _count_synapses(sim_config_eff, args.pre_gid, args.post_gid)
        if n_syn == 0:
            logger.error("No synapses found for %d->%d", args.pre_gid, args.post_gid)
            return
        logger.info("Pair %d->%d: %d synapses", args.pre_gid, args.post_gid, n_syn)

        configs = _generate_configs(n_syn)
        if args.configs == "ends":
            configs = [configs[0], configs[-1]]
        logger.info("%d rho configurations × %d trials = %d tasks",
                    len(configs), args.num_trials, len(configs) * args.num_trials)

        tasks = [
            (sim_config_eff, args.sim_config, args.pre_gid, args.post_gid, cfg, trial)
            for cfg in configs
            for trial in range(args.num_trials)
        ]

        n_workers = args.workers or min(multiprocessing.cpu_count(), len(tasks))
        with multiprocessing.Pool(n_workers) as pool:
            results = pool.map(_run_trial, tasks)

        df = pd.DataFrame([{k: v for k, v in r.items() if k != "wave"} for r in results])

        n_total = len(df)
        n_nan   = df["epsp"].isna().sum()
        if n_nan > 0:
            logger.warning("%d/%d trials returned NaN EPSP", n_nan, n_total)
        if n_nan == n_total:
            logger.error("All trials failed (all NaN) — not writing CSV. Check errors above.")
            raise SystemExit(1)

        # mean/std/count stay the raw-peak columns PairBasis reads; the Tazerart measures are extra
        g = df.groupby("config", sort=False)
        stats = g["epsp"].agg(["mean", "std", "count"])
        for m in ("pavg", "area"):
            stats[f"mean_{m}"] = g[f"epsp_{m}"].mean()
            stats[f"std_{m}"] = g[f"epsp_{m}"].std()
        stats = stats.reset_index()
        if args.waveforms:
            w = {r["config"]: r["wave"] for r in results if r["trial"] == 0 and r["wave"] is not None}
            np.savez(args.waveforms, t=WAVE_T, configs=np.array(list(w)), waves=np.array(list(w.values())))
        stats["pre_gid"]  = args.pre_gid
        stats["post_gid"] = args.post_gid

        os.makedirs(os.path.dirname(os.path.abspath(args.output_csv)), exist_ok=True)
        stats.to_csv(args.output_csv, index=False)
        logger.info("Saved %s", args.output_csv)

    finally:
        try:
            os.unlink(sim_config_eff)
        except Exception:
            pass


if __name__ == "__main__":
    main()
