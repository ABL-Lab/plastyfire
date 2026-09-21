"""
CLI wrapper for fitting-loop use of the edges-based prefire simulation.

Accepts individual parameter values on the command line (like pairrunner.py does
for the recipe pipeline) so that evaluator_edges.py can call it via subprocess.run()
from inside a multiprocessing Pool worker without triggering nested-daemon errors.
The subprocess is a fresh OS process that is not a daemon, so it may freely spawn the
multiprocessing.Process used internally by simulator_edges.runconnectedpair_prefire_from_edges.

Usage (run from simulation workdir):
    python pairrunner_edges_fit.py \\
        --gamma_d_GB_GluSynapse=120.0 --gamma_p_GB_GluSynapse=210.0 \\
        --a00=2.0 --a01=1.5 --a10=2.0 --a11=2.0 --a20=3.0 --a21=2.0 --a30=4.0 --a31=2.0 \\
        --tau_effca_GB_GluSynapse=278.3177658387 \\
        --fastforward=280000 \\
        --param_hash=abc123def456 \\
        [--force]                    # keep runs that miss the induction spike count

Output: simulation_edges_{param_hash}.pkl written in workdir.
"""

import argparse
import hashlib
import logging
import os
import pickle
import time

import plastyfire.simulator_edges as sim_mod

# Parameters optimised by the fitting loop (tau_effca is fixed, not fitted)
FIT_PARAM_NAMES = [
    "gamma_d_GB_GluSynapse",
    "gamma_p_GB_GluSynapse",
    "a00", "a01",
    "a10", "a11",
    "a20", "a21",
    "a30", "a31",
]
FITTED_TAU = 278.3177658387  # fixed; injected via --tau_effca_GB_GluSynapse


if __name__ == "__main__":
    workdir = os.getcwd()

    parser = argparse.ArgumentParser(
        description="Run edges-based prefire induction sim with per-param arguments for fitting"
    )
    for name in FIT_PARAM_NAMES:
        parser.add_argument(f"--{name}", type=float, required=True,
                            help=f"GluSynapse model parameter {name}")
    # tau_effca is accepted but defaults to FITTED_TAU so callers may omit it
    parser.add_argument("--tau_effca_GB_GluSynapse", type=float, default=FITTED_TAU,
                        help="Ca* integration time constant (default: pre-optimised value)")
    parser.add_argument("--fastforward", type=float, default=None,
                        help="Fast-forward start point in ms")
    parser.add_argument("--edges-h5", default=sim_mod.EDGES_H5_DEFAULT,
                        help=f"Path to dhuruva_modified_edges.h5 (default: {sim_mod.EDGES_H5_DEFAULT})")
    parser.add_argument("--cpre-cpost-cache", default=None,
                        help="Path to precomputed c_pre/c_post cache pkl")
    parser.add_argument("--param_hash", type=str, default=None,
                        help="Hash tag used for output filename; auto-derived if omitted")
    parser.add_argument("--circuit-config", default=None,
                        help="Override the 'network' field in prefire_simulation_config.json "
                             "(e.g. for new ion channels)")
    parser.add_argument("--extracellular-calcium", type=float, default=None,
                        help="Override extracellular Ca (mM). Sets both "
                             "conditions.extracellular_calcium (release-probability "
                             "scaling) and GluSynapse.cao_CR (VDCC/NMDA driving force). "
                             "Default: leave the config value (2.0) alone.")
    parser.add_argument("--force", action="store_true",
                        help="Ignore the induction spike-count guardrail: write the "
                             "result pkl even when the post cell fired the wrong "
                             "number of APs during the pairing train. The pkl carries "
                             "guardrail_forced=True so the run stays identifiable.")
    parser.add_argument("--lean", action="store_true",
                        help="Record rho only: skip rho_timeseries.npy (~3.5 MB) and "
                             "simulation_traces.pkl (~275 MB) per workdir, and drop "
                             "t/v/spike trains from the output pkl (520 kB -> ~1 kB).")
    parser.add_argument("--out", default=None,
                        help="Explicit output pkl path. Default: "
                             "<workdir>/simulation_edges_<param_hash>.pkl")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
    )
    log = logging.getLogger(__name__)

    # Build fit_params dict passed to simulator_edges
    fit_params = {name: getattr(args, name) for name in FIT_PARAM_NAMES}
    fit_params["tau_effca_GB_GluSynapse"] = args.tau_effca_GB_GluSynapse

    # Determine output filename
    if args.param_hash:
        param_hash = args.param_hash
    else:
        param_hash = hashlib.md5(str(sorted(fit_params.items())).encode()).hexdigest()[:12]
    filename = f"simulation_edges_{param_hash}.pkl"

    log.info("workdir          : %s", workdir)
    log.info("param_hash       : %s", param_hash)
    log.info("fastforward      : %s ms", args.fastforward or "disabled")
    log.info("edges_h5         : %s", args.edges_h5)
    log.info("cpre_cpost_cache : %s", args.cpre_cpost_cache or "(none)")
    log.info("extracell. Ca    : %s mM", args.extracellular_calcium
             if args.extracellular_calcium is not None else "(config default)")
    log.info("lean             : %s", args.lean)
    log.info("force            : %s", args.force)
    if args.cpre_cpost_cache:
        import os as _os
        log.info("cpre_cpost_cache exists: %s", _os.path.exists(args.cpre_cpost_cache))
    log.info("fit_params : %s", fit_params)

    t_start = time.time()
    results = sim_mod.runconnectedpair_prefire_from_edges(
        workdir,
        fit_params=fit_params,
        edges_h5=args.edges_h5,
        fastforward=args.fastforward,
        cpre_cpost_cache=args.cpre_cpost_cache,
        circuit_config=args.circuit_config,
        bcl_subdir="bluecellulab_results_optimizer",
        extracellular_calcium=args.extracellular_calcium,
        lean=args.lean,
        force=args.force,
    )
    elapsed = time.strftime("%H:%M:%S", time.gmtime(time.time() - t_start))
    log.info("Simulation finished in %s", elapsed)

    out_path = args.out or os.path.join(workdir, filename)
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, "wb") as f:
        pickle.dump(results, f, protocol=-1)
    log.info("Results saved → %s", out_path)
