#!/usr/bin/env python3
"""
Script to generate EPSP ratio comparison figure for 5 protocols.
Runs prefire simulations and ephys for mrk97_01 to mrk97_05 protocols
on L5TTPC_L5TTPC pairs and compares against biodata.

Author: Generated for evaluation
"""

import argparse
import hashlib
import logging
import multiprocessing as mp
import os
import pickle
import subprocess
import shutil
import sys
import time

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Add plastyfire to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "plastyfire"))

from plastyfire.config import OptConfig
from plastyfire.evaluator import compute_epsp_prefire

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Constants
CSVF_NAME = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/biodata/paired_recordings.csv"
CONFIGS_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/configs"

EVAL_PROTOCOLS = ["sjh06_02", "mrk97_01", "mrk97_02", "mrk97_07", "mrk97_08"]

# Best solution file path
BESTSOL_PKL = "/lustre06/project/6077694/dhuruva/plastyfire/bestsol.pkl"
CHECKPOINT_PKL = "/lustre06/project/6077694/dhuruva/plastyfire/plastyfire/checkpoint.pkl"

# Param names from modelfitter.py
FIT_PARAM_NAMES = [
    "gamma_d_GB_GluSynapse",
    "gamma_p_GB_GluSynapse",
    "a00",
    "a01",
    "a10",
    "a11",
    "a20",
    "a21",
    "a30",
    "a31",
]

MIN2MS = 60 * 1000.0

# Custom parameter sets for comparison
CUSTOM_PARAM_SETS = {
    "Chindemi_params": {
        "gamma_d_GB_GluSynapse": 101.5,
        "gamma_p_GB_GluSynapse": 216.2,
        "a00": 1.002,
        "a01": 1.954,
        "a10": 1.159,
        "a11": 2.483,
        "a20": 1.127,
        "a21": 2.456,
        "a30": 5.236,
        "a31": 1.782,
    },
    "Deep_learning_params": {
        "gamma_d_GB_GluSynapse": 163.46167854068054,
        "gamma_p_GB_GluSynapse": 221.93928040026307,
        "a00": 0.3296322326165755,
        "a01": 1.9167068577889812,
        "a10": 2.5252431110314606,
        "a11": 0.317159582514317,
        "a20": 7.2678533585713385,
        "a21": 4.264799660065239,
        "a30": 5.480241019952461,
        "a31": 0.3079526712830378,
    },
    "mrk97_07_mrk97_08": {
        "gamma_d_GB_GluSynapse": 199.887594498312,
        "gamma_p_GB_GluSynapse": 194.23843405523158,
        "a00": 1.0005007769557286,
        "a01": 2.706871658650942,
        "a10": 2.5132952169332246,
        "a11": 3.9826946063778337,
        "a20": 9.594554277367934,
        "a21": 3.0337323982995814,
        "a30": 1.6205328327559538,
        "a31": 2.692343865436747,
    },
}

def load_best_solution():
    """Load the best solution from bestsol.pkl or checkpoint.pkl"""
    if os.path.exists(BESTSOL_PKL):
        logger.info(f"Loading best solution from {BESTSOL_PKL}")
        with open(BESTSOL_PKL, "rb") as f:
            best_params = pickle.load(f)
        return best_params
    elif os.path.exists(CHECKPOINT_PKL):
        logger.warning(
            f"bestsol.pkl not found, loading from checkpoint: {CHECKPOINT_PKL}"
        )
        with open(CHECKPOINT_PKL, "rb") as f:
            checkpoint = pickle.load(f)
        if "halloffame" in checkpoint and len(checkpoint["halloffame"]) > 0:
            hof = checkpoint["halloffame"]
            best_ind = hof[0]
            best_params = dict(zip(FIT_PARAM_NAMES, best_ind))
            logger.info(f"Extracted best individual from checkpoint: {best_params}")
            return best_params
        else:
            raise ValueError("No hall of fame found in checkpoint")
    else:
        raise FileNotFoundError(f"Neither {BESTSOL_PKL} nor {CHECKPOINT_PKL} found")


def run_prefire_simulation(args):
    """Worker function to run a single prefire simulation"""
    sim_dict, fit_params, param_hash, results_dir, no_cache, recipe_path = args

    try:
        workdir = sim_dict["workdir"]
        pre_gid, post_gid = workdir.split("/")[-2].split("-")
        output_file = os.path.join(results_dir, f"simulation_{param_hash}_{sim_dict['protocol_id']}_{pre_gid}_{post_gid}.pkl")

        if os.path.exists(output_file):
            logger.info(f"Simulation already exists in results dir: {output_file}")
            return sim_dict["protocol_id"], output_file

        logger.info(f"Running prefire simulation for {sim_dict['protocol_id']} in {workdir}")

        param_args = [f"--{name}={value}" for name, value in fit_params.items()]
        param_args.append(f"--fastforward={sim_dict['fastforward']}")
        param_args.append(f"--param_hash={param_hash}")
        if recipe_path:
            param_args.append(f"--recipe-path={recipe_path}")
        
        workdir_output_file = os.path.join(workdir, f"simulation_{param_hash}_{sim_dict['protocol_id']}_{pre_gid}_{post_gid}.pkl")
        param_args.append(f"--output-filename={os.path.basename(workdir_output_file)}")

        script_path = os.path.join(os.path.dirname(__file__), "plastyfire", "pairrunner.py")
        cmd = [sys.executable, script_path] + param_args

        if no_cache and os.path.exists(workdir_output_file):
            try:
                os.remove(workdir_output_file)
            except OSError:
                pass

        result = subprocess.run(
            cmd,
            cwd=workdir,
            capture_output=True,
            text=True,
            timeout=1800,
        )

        if result.returncode != 0:
            logger.error(f"Simulation failed for {sim_dict['protocol_id']}: {result.stderr}")
            return None

        if not os.path.exists(workdir_output_file):
            logger.error(f"Output file not created for {sim_dict['protocol_id']}: {workdir_output_file}")
            return None

        shutil.move(workdir_output_file, output_file)
        return sim_dict["protocol_id"], output_file

    except Exception as e:
        logger.error(f"Simulation failed for {sim_dict['protocol_id']}: {e}")
        return None


def generate_ephys_data_worker(args):
    """Worker function to generate a single ephys file"""
    sim_config_path, pre_gid, post_gid, rho_str, ephys_file, no_cache, recipe_path, synapse_ids_str, fit_params = args

    if os.path.exists(ephys_file):
        return ephys_file

    try:
        script_path = os.path.join(os.path.dirname(__file__), "plastyfire", "compute_test_pulse_epsp_ratio.py")

        cmd = [
            sys.executable,
            script_path,
            sim_config_path,
            str(pre_gid),
            str(post_gid),
            rho_str,
            "--node_pop", "S1nonbarrel_neurons",
            "--trial", "0",
            "--output-dir", os.path.dirname(ephys_file),
        ]

        if recipe_path:
            cmd.extend(["--recipe-path", recipe_path])

        if synapse_ids_str:
            cmd.extend(["--synapse-ids", synapse_ids_str])

        if fit_params:
            for name, value in fit_params.items():
                cmd.append(f"--{name}={value}")

        if no_cache:
            cmd.append("--no-cache")

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

        if result.returncode != 0:
            logger.error(f"Ephys generation failed for {pre_gid}->{post_gid}: {result.stderr}")
            return None

        return ephys_file

    except Exception as e:
        logger.error(f"Error generating ephys for {pre_gid}->{post_gid}: {e}")
        return None


def generate_ephys_data(sim_dict, param_hash, fit_params, results_dir, no_cache, recipe_path):
    """Generate ephys data for initial and final rho states"""
    workdir = sim_dict["workdir"]
    pre_gid, post_gid = workdir.split("/")[-2].split("-")
    pkl_file = os.path.join(results_dir, f"simulation_{param_hash}_{sim_dict['protocol_id']}_{pre_gid}_{post_gid}.pkl")

    if not os.path.exists(pkl_file):
        return []

    with open(pkl_file, "rb") as f:
        raw_results = pickle.load(f)

    rho_data = raw_results["rho_GB"]
    if len(rho_data) > 100:
        rho_data = np.transpose(rho_data)

    synapse_ids_str = None
    if "synprop" in raw_results and "synapseID" in raw_results["synprop"]:
        synapse_ids = raw_results["synprop"]["synapseID"]
        synapse_ids_str = ",".join(map(str, synapse_ids))

    initial_rho = [0 if k[0] < 0.5 else 1 for k in rho_data]
    final_rho = [0 if k[-1] < 0.5 else 1 for k in rho_data]

    ephys_dir = results_dir
    os.makedirs(ephys_dir, exist_ok=True)
    sim_config_path = os.path.join(workdir, "simulation_config.json")
    jobs = []

    for rho_config in [initial_rho, final_rho]:
        rho_str = ",".join(map(str, rho_config))
        ephys_file = os.path.join(ephys_dir, f"ephys_data_{pre_gid}_{post_gid}_{rho_str.replace(',', '_')}.pkl")
        jobs.append((sim_config_path, pre_gid, post_gid, rho_str, ephys_file, no_cache, recipe_path, synapse_ids_str, fit_params))

    return jobs


def compute_epsp_ratios(sim_dict, param_hash, results_dir):
    """Compute EPSP ratios from ephys data"""
    workdir = sim_dict["workdir"]
    pre_gid, post_gid = workdir.split("/")[-2].split("-")
    pkl_file = os.path.join(results_dir, f"simulation_{param_hash}_{sim_dict['protocol_id']}_{pre_gid}_{post_gid}.pkl")
    
    with open(pkl_file, "rb") as f:
        raw_results = pickle.load(f)

    rho_data = raw_results["rho_GB"]
    if len(rho_data) > 100:
        rho_data = np.transpose(rho_data)

    initial_rho = [0 if k[0] < 0.5 else 1 for k in rho_data]
    final_rho = [0 if k[-1] < 0.5 else 1 for k in rho_data]

    ephys_dir = results_dir
    initial_rho_str = "_".join(map(str, initial_rho))
    final_rho_str = "_".join(map(str, final_rho))

    ephys_before_path = os.path.join(ephys_dir, f"ephys_data_{pre_gid}_{post_gid}_{initial_rho_str}.pkl")
    ephys_after_path = os.path.join(ephys_dir, f"ephys_data_{pre_gid}_{post_gid}_{final_rho_str}.pkl")

    if not os.path.exists(ephys_before_path) or not os.path.exists(ephys_after_path):
        return None

    try:
        epsp_before = compute_epsp_prefire(ephys_before_path)
        epsp_after = compute_epsp_prefire(ephys_after_path)
        epsp_ratio = epsp_after / epsp_before

        return {
            "protocol_id": sim_dict["protocol_id"],
            "epsp_ratio": epsp_ratio,
        }
    except Exception:
        return None


def main():
    parser = argparse.ArgumentParser(description="Generate EPSP ratio comparison figure")
    parser.add_argument("--sample-size", type=int, default=100, help="Number of pairs to evaluate")
    parser.add_argument("--workers", type=int, default=60, help="Number of parallel workers")
    parser.add_argument("--seed", type=int, default=1234, help="Random seed")
    parser.add_argument("--no-cache", action="store_true", help="Force re-computation")
    parser.add_argument("--recipe-path", type=str, default=None, help="Path to custom recipe.csv file")
    parser.add_argument(
        "--param-set",
        type=str,
        default="Chindemi_params",
        choices=["best"] + list(CUSTOM_PARAM_SETS.keys()),
        help=f"Parameter set to evaluate: 'best' or one of {list(CUSTOM_PARAM_SETS.keys())}",
    )
    args = parser.parse_args()

    if args.param_set == "best":
        try:
            params = load_best_solution()
            logger.info(f"Using best parameters: {params}")
        except Exception as e:
            logger.error(f"Failed to load best solution: {e}")
            return 1
    else:
        params = CUSTOM_PARAM_SETS[args.param_set]
        logger.info(f"Using custom parameter set '{args.param_set}': {params}")

    param_values = [params[name] for name in FIT_PARAM_NAMES]
    param_hash = hashlib.md5(str(param_values).encode()).hexdigest()[:12]
    
    results_base_dir = f"/lustre06/project/6077694/dhuruva/plastyfire/evaluation_results/figure_gen_{args.param_set}/sim_results"
    os.makedirs(results_base_dir, exist_ok=True)

    invitro_db = pd.read_csv(CSVF_NAME)
    invitro_db = invitro_db.loc[invitro_db["protocol_id"].isin(EVAL_PROTOCOLS)]

    all_sims = []
    for elem in invitro_db.itertuples():
        config = OptConfig(os.path.join(CONFIGS_DIR, f"{elem.pre_mtype}_{elem.post_mtype}.yaml"))
        fastforward = config.fastforward
        if fastforward is None:
            fastforward = config.C01_duration * MIN2MS + config.nreps * config.T
        fastforward = float(fastforward)

        index_file = os.path.join(os.path.split(os.path.split(config.out_dir)[0])[0], f"index_{elem.pre_mtype}_{elem.post_mtype}.csv")
        if not os.path.exists(index_file):
            continue

        sim_idx = pd.read_csv(index_file)
        sim_idx.set_index(["frequency", "dt"], inplace=True)
        sim_idx.sort_index(inplace=True)

        paths = sim_idx.loc[elem.frequency_train, elem.dt_train]["path"]
        if isinstance(paths, str):
            paths = pd.Series([paths])

        np.random.seed(args.seed)
        sample_size = min(args.sample_size, len(paths))
        paths = paths.sample(sample_size, random_state=args.seed)

        for path in paths:
            all_sims.append({
                "protocol_id": elem.protocol_id,
                "fastforward": fastforward,
                "workdir": "/".join(path.split("/")[:-1]),
            })

    logger.info(f"Starting {len(all_sims)} simulations...")

    # Step 1: Prefire
    sim_args = [(sim_dict, params, param_hash, results_base_dir, args.no_cache, args.recipe_path) for sim_dict in all_sims]
    with mp.Pool(processes=args.workers) as pool:
        sim_results = pool.map(run_prefire_simulation, sim_args)
    sim_results = [r for r in sim_results if r is not None]

    # Step 2: Ephys
    all_ephys_jobs = []
    for sim_dict in all_sims:
        jobs = generate_ephys_data(sim_dict, param_hash, params, results_base_dir, args.no_cache, args.recipe_path)
        all_ephys_jobs.extend(jobs)
    
    with mp.Pool(processes=args.workers) as pool:
        ephys_results = pool.map(generate_ephys_data_worker, all_ephys_jobs)

    # Step 3: EPSP Ratios
    epsp_args = [(sim_dict, param_hash, results_base_dir) for sim_dict in all_sims]
    with mp.Pool(processes=args.workers) as pool:
        epsp_results = pool.starmap(compute_epsp_ratios, epsp_args)
    epsp_results = [r for r in epsp_results if r is not None]

    results_df = pd.DataFrame(epsp_results)
    
    # Plotting
    plt.figure(figsize=(10, 6))
    
    frequencies = []
    sim_means = []
    sim_sems = []
    bio_means = []
    bio_sems = []

    for protocol in EVAL_PROTOCOLS:
        bio_data = invitro_db[invitro_db["protocol_id"] == protocol].iloc[0]
        sim_data = results_df[results_df["protocol_id"] == protocol]["epsp_ratio"]
        
        freq = bio_data["frequency_train"]
        frequencies.append(freq)
        
        bio_means.append(bio_data["mean_epsp_ratio"])
        bio_sems.append(bio_data["sem_epsp_ratio"])
        
        sim_means.append(sim_data.mean())
        sim_sems.append(sim_data.sem())

    plt.errorbar(frequencies, bio_means, yerr=bio_sems, fmt='o-', label='In Vitro', capsize=5)
    plt.errorbar(frequencies, sim_means, yerr=sim_sems, fmt='s-', label='In Silico', capsize=5)
    
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('EPSP Ratio')
    plt.title('EPSP Ratio vs Frequency Comparison')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    output_path = "epsp_ratio_comparison.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    logger.info(f"Figure saved to {output_path}")

    return 0

if __name__ == "__main__":
    sys.exit(main())
