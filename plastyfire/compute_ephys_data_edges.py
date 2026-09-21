"""
Edges-compatible EPSP test-pulse data generator.

Replaces compute_test_pulse_epsp_ratio.py for the edges fitting pipeline.

Runs test-pulse simulations using bluecellulab (which auto-loads Use_d/p,
gmax_d/p from dhuruva_modified_edges.h5) for a given rho configuration.
Plasticity is disabled during the test pulse (theta_d/theta_p set to 1e9).

Because synapse weights come from edges.h5 and plasticity is off, the output
ephys_data files are INDEPENDENT of the fitted parameters (gamma_d, gamma_p,
a-params). They only need to be generated once per (pre, post, rho_config)
triple and can be reused across all fitting iterations.

Output format (matches compute_test_pulse_epsp_ratio.py):
    ephys_data_{pre_gid}_{post_gid}_{rho_config_underscored}.pkl
    keys: t, v, rho_GB, pre_spikes, pre_gid, post_gid, rho_config, trial

Usage (CLI):
    python compute_ephys_data_edges.py \\
        /path/to/simulation_config.json \\
        <pre_gid> <post_gid> <rho_config e.g. 0,1,0,0> \\
        [--output-dir /path/to/output] [--trial 0]
"""

import argparse
import logging
import os
import pickle

import numpy as np

logger = logging.getLogger(__name__)

NODE_POP        = "S1nonbarrel_neurons"
MECHANISMS_PATH = "/project/ctb-emuller/dhuruva/DEES_cell_packages/"
DISABLE_THRESH  = 1e9   # set theta_d/theta_p to this to disable plasticity
SIM_DURATION_MS = 240000.0  # 4 min — same as compute_test_pulse_epsp_ratio.py


def get_epsp_value_edges(sim_config_path, pre_gid, post_gid, rho_config,
                         node_pop=NODE_POP, trial=0,
                         output_dir=None, no_cache=False):
    """
    Run a test-pulse simulation using bluecellulab + edges.h5 for the given
    rho configuration and save ephys_data pkl to disk.

    Synapse parameters (Use_d/p, gmax_d/p, rho0_GB) are auto-loaded from
    dhuruva_modified_edges.h5 by bluecellulab — no recipe/ParamsGenerator needed.

    :param sim_config_path: path to simulation_config.json (workdir must contain prespikes.h5)
    :param pre_gid: presynaptic GID
    :param post_gid: postsynaptic GID
    :param rho_config: comma-separated rho states, e.g. "0,1,0,0"
    :param node_pop: SONATA node population name
    :param trial: trial index (used as random seed and in filename lookup)
    :param output_dir: directory to write the pkl; defaults to 4 levels above workdir/ephys_data/
    :param no_cache: if True, overwrite existing file
    :returns: path to the saved pkl file
    """
    import bluecellulab
    from libsonata import SpikeReader

    safe_rho = rho_config.replace(",", "_")
    filename = f"ephys_data_{pre_gid}_{post_gid}_{safe_rho}.pkl"

    if output_dir:
        filepath = os.path.join(output_dir, filename)
    else:
        filepath = os.path.normpath(
            os.path.join(os.path.dirname(sim_config_path),
                         "..", "..", "..", "..", "ephys_data", filename)
        )

    if os.path.exists(filepath) and not no_cache:
        logger.info("Ephys file exists, skipping: %s", filepath)
        return filepath

    logger.info("Generating edges ephys file: %s", filepath)

    bluecellulab.neuron.load_mechanisms(MECHANISMS_PATH)
    bluecellulab.neuron.h.cao_CR_GluSynapse     = 2.0
    bluecellulab.neuron.h.minis_single_vesicle_GluSynapse = 0.0
    bluecellulab.neuron.h.init_depleted_GluSynapse        = 0.0

    np.random.seed(trial)
    sim = bluecellulab.CircuitSimulation(sim_config_path, base_seed=trial)

    workdir    = os.path.dirname(sim_config_path)
    pre_spikes = SpikeReader(os.path.join(workdir, "prespikes.h5"))[node_pop].get_dict()["timestamps"]

    sim.instantiate_gids(
        [(node_pop, post_gid)],
        add_synapses=True, add_minis=False, add_pulse_stimuli=True,
        intersect_pre_gids=[(node_pop, pre_gid)],
        pre_spike_trains={(node_pop, pre_gid): pre_spikes},
    )
    cell = sim.cells[(node_pop, post_gid)]

    for sec in cell.somatic + cell.axonal:
        sec.uninsert("SK_E2")

    # Set up rho_GB recordings
    rho_recorders = []
    for _syn_id, synapse in cell.synapses.items():
        rec = bluecellulab.neuron.h.Vector()
        rec.record(synapse.hsynapse._ref_rho_GB)
        rho_recorders.append(rec)

    rho_vals = [int(x) for x in rho_config.split(",")]
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
        # Disable plasticity during test pulse
        h.theta_d_GB = DISABLE_THRESH
        h.theta_p_GB = DISABLE_THRESH

    bluecellulab.neuron.h.cvode_active(1)
    sim.run(SIM_DURATION_MS, cvode=True)

    t = np.array(sim.get_time())
    v = np.array(sim.get_voltage_trace((node_pop, post_gid)))
    rho_gb = np.array([np.array(r) for r in rho_recorders])

    del sim

    data = {
        "t":          t,
        "v":          v,
        "rho_GB":     rho_gb,
        "pre_spikes": pre_spikes,
        "pre_gid":    pre_gid,
        "post_gid":   post_gid,
        "rho_config": rho_config,
        "trial":      trial,
    }

    os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)
    with open(filepath, "wb") as f:
        pickle.dump(data, f, protocol=-1)
    logger.info("Saved %s", filepath)
    return filepath


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
    )

    parser = argparse.ArgumentParser(
        description="Generate edges-compatible ephys_data pkl for one rho configuration"
    )
    parser.add_argument("sim_config_path",
                        help="Path to simulation_config.json (workdir must contain prespikes.h5)")
    parser.add_argument("pre_gid",  type=int, help="Presynaptic GID")
    parser.add_argument("post_gid", type=int, help="Postsynaptic GID")
    parser.add_argument("rho_config",
                        help="Comma-separated rho states, e.g. '0,1,0,0'")
    parser.add_argument("--node-pop", default=NODE_POP, help="SONATA node population name")
    parser.add_argument("--trial",    type=int, default=0, help="Trial index (used as RNG seed)")
    parser.add_argument("--output-dir", default=None,
                        help="Output directory (default: 4 levels above workdir/ephys_data/)")
    parser.add_argument("--no-cache", action="store_true",
                        help="Overwrite existing file")
    args = parser.parse_args()

    out = get_epsp_value_edges(
        args.sim_config_path,
        args.pre_gid,
        args.post_gid,
        args.rho_config,
        node_pop=args.node_pop,
        trial=args.trial,
        output_dir=args.output_dir,
        no_cache=args.no_cache,
    )
    print(f"Written: {out}")
