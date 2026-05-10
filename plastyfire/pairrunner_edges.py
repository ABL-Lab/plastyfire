"""
Runner for connected-pair STDP induction simulations using dhuruva_modified_edges.h5.
Thresholds (theta_d/theta_p) and synapse properties (rho0, Use_d/p, gmax_d/p) are read
directly from the edges HDF5 — no EPG calibration step.

Usage (run from the simulation workdir via sbatch --chdir=<workdir>):
    python pairrunner_edges.py --params chindemi --fastforward 280000
    python pairrunner_edges.py --params chindemi  # no fastforward
"""

import os
import sys
import time
import pickle
import argparse
import logging

import plastyfire.simulator_edges as sim_mod

# ── Named parameter presets ──────────────────────────────────────────────────

PARAM_PRESETS = {
    "chindemi": {
        "tau_effca_GB_GluSynapse": 278.3177658387,
        "gamma_d_GB_GluSynapse":   101.5387594661,
        "gamma_p_GB_GluSynapse":   216.1841700668,
    },
}


if __name__ == "__main__":
    workdir = os.getcwd()

    parser = argparse.ArgumentParser(
        description="Run induction sim; read synapse params from dhuruva_modified_edges.h5"
    )
    parser.add_argument(
        "--params", choices=list(PARAM_PRESETS.keys()),
        help="Named parameter preset for global GluSynapse HOC vars (e.g. chindemi)",
    )
    parser.add_argument(
        "--fastforward", type=float, default=None,
        help="Fast-forward start point in ms (default: no fast-forward)",
    )
    parser.add_argument(
        "--edges-h5", default=sim_mod.EDGES_H5_DEFAULT,
        help=f"Path to dhuruva_modified_edges.h5 (default: {sim_mod.EDGES_H5_DEFAULT})",
    )
    parser.add_argument(
        "--output-filename", default=None,
        help="Output pickle filename (default: simulation_edges.pkl)",
    )
    parser.add_argument(
        "--traces-output-dir", default=None,
        help=(
            "If given, simulation_traces.pkl is written to "
            "<traces-output-dir>/<pair_name>/<protocol>/simulation_traces.pkl "
            "instead of bluecellulab_results/ inside the workdir."
        ),
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
    )
    log = logging.getLogger(__name__)

    fit_params = dict(PARAM_PRESETS[args.params]) if args.params else {}

    log.info("workdir          : %s", workdir)
    log.info("params           : %s", args.params or "(none)")
    log.info("fastforward      : %s ms", args.fastforward or "disabled")
    log.info("edges_h5         : %s", args.edges_h5)
    log.info("traces_output_dir: %s", args.traces_output_dir or "(workdir/bluecellulab_results)")

    t_start = time.time()
    results = sim_mod.runconnectedpair_prefire_from_edges(
        workdir,
        fit_params=fit_params or None,
        edges_h5=args.edges_h5,
        fastforward=args.fastforward,
        traces_output_dir=args.traces_output_dir,
    )
    elapsed = time.strftime("%H:%M:%S", time.gmtime(time.time() - t_start))
    log.info("Simulation finished in %s", elapsed)

    filename = args.output_filename or "simulation_edges.pkl"
    out_path = os.path.join(workdir, filename)
    with open(out_path, "wb") as f:
        pickle.dump(results, f, protocol=-1)
    log.info("Results saved → %s", out_path)
