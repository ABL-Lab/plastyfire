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

FIT_PARAM_NAMES = [  # "tau_effca_GB_GluSynapse",
                   "gamma_d_GB_GluSynapse", "gamma_p_GB_GluSynapse",
                   "a00", "a01", "a10", "a11", "a20", "a21", "a30", "a31"]
FITTED_TAU = 278.3177658387  # previously optimized time constant of Ca*


if __name__ == "__main__":
    workdir = os.getcwd()
    # Parse command line
    parser = argparse.ArgumentParser()
    for param_name in FIT_PARAM_NAMES:
        parser.add_argument("--%s" % param_name, type=float, help="GluSynapse model parameter")
    parser.add_argument("--fastforward", type=float, help="Fastforward begin point (ms)")
    parser.add_argument("--param_hash", type=str, help="Hash of parameter values for unique output filename")
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="Verbose messages")
    parser.add_argument("--debug", default=False, action="store_true", help="Enable debug mode")
    args = parser.parse_args()
    # Create dictionary of fitted parameters, if needed
    fit_params = {param_name: getattr(args, param_name) for param_name in FIT_PARAM_NAMES if
                  getattr(args, param_name) is not None}
    fit_params["tau_effca_GB_GluSynapse"] = FITTED_TAU
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
    results = sim.runconnectedpair(workdir, fit_params=fit_params, fastforward=args.fastforward)
    logger.info("Simulation finished in: %s" % time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))
    # Store results
    # Create unique filename based on parameter hash to avoid file locking conflicts
    if args.param_hash is not None:
        output_filename = f"simulation_{args.param_hash}.pkl"
    else:
        # Fallback: create hash from the actual parameter values
        param_values_list = [fit_params.get(name) for name in FIT_PARAM_NAMES if fit_params.get(name) is not None]
        param_hash = hashlib.md5(str(param_values_list).encode()).hexdigest()[:12]
        output_filename = f"simulation_{param_hash}.pkl"
    
    # Write results as pickle file in the format expected by Experiment class
    import pickle
    with open(os.path.join(workdir, output_filename), "wb") as f:
        pickle.dump(results, f)
    logger.info("Data writing finished")



