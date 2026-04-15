"""
Simple run script for `simulator/runconnectedpair()` that parses command line arguments, runs sim, and saves results
authors: Giuseppe Chindemi (12.2020) + minor modifications by András Ecker (06.2024)
"""

import os
import time
import argparse
import logging
import hashlib

import plastyfire.simulator as sim
from plastyfire.epg_dhuruva_custom_ratio import RHO_RATIO_ENV, parse_rho_ratio

FIT_PARAM_NAMES = ["enable_CICR_GluSynapse",
                   "gamma_d_GB_GluSynapse", "gamma_p_GB_GluSynapse",
                   "a00", "a01", "a10", "a11", "a20", "a21", "a30", "a31",
                   "delta_IP3_CICR_GluSynapse", "tau_IP3_CICR_GluSynapse",
                   "V_IP3R_CICR_GluSynapse", "V_RyR_CICR_GluSynapse",
                   "V_SERCA_CICR_GluSynapse", "K_SERCA_CICR_GluSynapse",
                   "V_leak_CICR_GluSynapse", "tau_extrusion_CICR_GluSynapse",
                   "tau_effca_GB_GluSynapse",
                   "tau_ref_CICR_GluSynapse", "K_h_ref_CICR_GluSynapse"]
FITTED_TAU = 278.318


if __name__ == "__main__":
    workdir = os.getcwd()
    # Parse command line
    parser = argparse.ArgumentParser()
    for param_name in FIT_PARAM_NAMES:
        parser.add_argument("--%s" % param_name, type=float, help="GluSynapse model parameter")
    parser.add_argument("--fastforward", type=float, help="Fastforward begin point (ms)")
    parser.add_argument("--epg-variant", "--variant",
                         choices=["epg", "epg_dhuruva", "epg_dhuruva_full", "epg_dhuruva_custom_ratio"],
                         default="epg_dhuruva",
                         help="Parameter generator module to use (default: epg_dhuruva)")
    parser.add_argument("--rho-ratio", "--rho_ratio", type=str, default=None,
                         help="Depressed:potentiated percentage split for epg_dhuruva_custom_ratio, e.g. 45:55")
    parser.add_argument("--param_hash", type=str, help="Hash of parameter values for unique output filename")
    parser.add_argument("--recipe-path", type=str, default=None, help="Path to custom recipe.csv file")
    parser.add_argument("--output-filename", type=str, default=None, help="Explicit output filename")
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="Verbose messages")
    parser.add_argument("--debug", default=False, action="store_true", help="Enable debug mode")
    args = parser.parse_args()
    # Create dictionary of fitted parameters, if needed
    fit_params = {param_name: getattr(args, param_name) for param_name in FIT_PARAM_NAMES if
                  getattr(args, param_name) is not None}
    if "tau_effca_GB_GluSynapse" not in fit_params:
        raise ValueError(f"Missing required parameter 'tau_effca_GB_GluSynapse'. The simulator must receive the optimized integration window from JAX.")
    if args.epg_variant == "epg_dhuruva_custom_ratio":
        parse_rho_ratio(args.rho_ratio)
    elif args.rho_ratio is not None:
        raise ValueError("--rho-ratio can only be used with --epg-variant=epg_dhuruva_custom_ratio")
    if args.rho_ratio is not None:
        os.environ[RHO_RATIO_ENV] = args.rho_ratio
    else:
        os.environ.pop(RHO_RATIO_ENV, None)
    # Configure logger
    logging.basicConfig(level=logging.DEBUG if args.debug else logging.INFO)
    logger = logging.getLogger(__name__)
    # Enable debugging mode
    if args.debug:
        sim.logger.setLevel(logging.DEBUG)
        sim.DEBUG = True
    if args.verbose:
        sim.logger.setLevel(logging.DEBUG)

    # Run simulation
    start_time = time.time()
    results = sim.runconnectedpair_induction(
        workdir,
        fit_params=fit_params,
        fastforward=args.fastforward,
        recipe_path=args.recipe_path,
        epg_variant=args.epg_variant,
    )
    logger.info("Simulation finished in: %s" % time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))
    # Store results
    import pickle
    if args.output_filename:
        filename = args.output_filename
    elif args.param_hash:
        filename = "simulation_%s.pkl" % args.param_hash
    else:
        filename = "simulation_results.pkl"
    
    filepath = os.path.join(workdir, filename)
    with open(filepath, "wb") as f:
        pickle.dump(results, f)
    logger.info("Data writing finished")

