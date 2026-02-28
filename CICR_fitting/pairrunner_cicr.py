"""
Simple run script for `simulator/runconnectedpair()` that parses command line arguments, runs sim, and saves results
authors: Giuseppe Chindemi (12.2020) + minor modifications by András Ecker (06.2024)
"""

import os
import time
import sys
import argparse
import logging
import hashlib

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
os.environ["BLUECELLULAB_MOD_LIBRARY_PATH"] = "/home/dhuruva/projects/ctb-emuller/dhuruva/DEES_cell_packages/x86_64/libnrnmech.so"
# os.environ["BBP_MECH_DIR"] = "/home/dhuruva/projects/ctb-emuller/dhuruva/DEES_cell_packages"
import simulator_cicr as sim

FIT_PARAM_NAMES = [
    "gamma_d", "gamma_p",
    "a00", "a01", "a10", "a11", "a20", "a21", "a30", "a31",
    "tau_prime", "k_prime", "K_prime_limit",
    "K_trig", "n_trig", "Vmax_CICR", "tau_CICR", "tau_eff"
]
FITTED_TAU = 278.3177658387  # Not used directly if tau_eff is provided



if __name__ == "__main__":
    workdir = os.getcwd()
    # Parse command line
    parser = argparse.ArgumentParser()
    for param_name in FIT_PARAM_NAMES:
        parser.add_argument("--%s" % param_name, type=float, help="GluSynapse model parameter")
    parser.add_argument("--fastforward", type=float, help="Fastforward begin point (ms)")
    parser.add_argument("--param_hash", type=str, help="Hash of parameter values for unique output filename")
    parser.add_argument("--recipe-path", type=str, default=None, help="Path to custom recipe.csv file")
    parser.add_argument("--output-filename", type=str, default=None, help="Explicit output filename")
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="Verbose messages")
    parser.add_argument("--debug", default=False, action="store_true", help="Enable debug mode")
    args = parser.parse_args()
    # Create dictionary of fitted parameters, if needed
    fit_params = {param_name: getattr(args, param_name) for param_name in FIT_PARAM_NAMES if getattr(args, param_name) is not None}
    
    # Map back to the mod file terminology immediately before passing to simulator
    mapped_params = {}
    if "gamma_d" in fit_params: mapped_params["gamma_d_GB_GluSynapse"] = fit_params["gamma_d"]
    if "gamma_p" in fit_params: mapped_params["gamma_p_GB_GluSynapse"] = fit_params["gamma_p"]
    if "tau_eff" in fit_params: mapped_params["tau_effca_GB_GluSynapse"] = fit_params["tau_eff"]
    
    # Keep thresholds as aXX for simulator_cicr.py
    for a_param in ["a00", "a01", "a10", "a11", "a20", "a21", "a30", "a31"]:
        if a_param in fit_params: mapped_params[a_param] = fit_params[a_param]
        
    # Mapping new CICR params
    cicr_maps = {
        "tau_prime": "tau_prime_CICR_GluSynapse",
        "k_prime": "k_prime_CICR_GluSynapse",
        "K_prime_limit": "K_prime_limit_CICR_GluSynapse",
        "K_trig": "K_RyR_CICR_GluSynapse",
        "n_trig": "n_RyR_CICR_GluSynapse",
        "Vmax_CICR": "Vmax_CICR_GluSynapse",
        "tau_CICR": "tau_rel_CICR_GluSynapse"
    }
    for old_k, new_k in cicr_maps.items():
        if old_k in fit_params:
            mapped_params[new_k] = fit_params[old_k]
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
    results = sim.runconnectedpair(workdir, fit_params=mapped_params, fastforward=args.fastforward, recipe_path=args.recipe_path)
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



