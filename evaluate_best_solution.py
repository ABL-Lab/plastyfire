#!/usr/bin/env python3
"""
Evaluation script for the best solution from optimization.
Runs prefire simulations and ephys for mrk97_03, mrk97_07, mrk97_08 protocols
on 100 L5TTPC_L5TTPC pairs and compares against biodata.

Author: Generated for evaluation
"""

import argparse
import hashlib
import logging
import multiprocessing as mp
import os
import pickle
import subprocess
import sys
import time

import numpy as np
import pandas as pd

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
EVAL_PROTOCOLS = ["mrk97_03", "mrk97_07", "mrk97_08"]  # Evaluation protocols
# Note: tau_effca_GB_GluSynapse is hardcoded in pairrunner.py, not fitted

# Best solution file path
BESTSOL_PKL = "/lustre06/project/6077694/dhuruva/plastyfire/bestsol.pkl"
# Alternative: use checkpoint.pkl if bestsol.pkl doesn't exist
CHECKPOINT_PKL = (
    "/lustre06/project/6077694/dhuruva/plastyfire/plastyfire/checkpoint.pkl"
)

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
    "custom_v1": {
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
    # Add more custom parameter sets here if needed
    # "custom_v2": {...},
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
        # Extract best individual from hall of fame
        if "halloffame" in checkpoint and len(checkpoint["halloffame"]) > 0:
            hof = checkpoint["halloffame"]
            # Hall of fame is already sorted, best individual is at index 0
            best_ind = hof[0]
            # Convert to parameter dictionary
            best_params = dict(zip(FIT_PARAM_NAMES, best_ind))
            logger.info(f"Extracted best individual from checkpoint: {best_params}")
            logger.info(f"Fitness scores: {best_ind.fitness.values}")
            logger.info(
                f"Combined score (norm): {sum(x**2 for x in best_ind.fitness.values) ** 0.5}"
            )
            return best_params
        else:
            raise ValueError("No hall of fame found in checkpoint")
    else:
        raise FileNotFoundError(f"Neither {BESTSOL_PKL} nor {CHECKPOINT_PKL} found")


def run_prefire_simulation(args):
    """Worker function to run a single prefire simulation using subprocess to avoid daemon issues"""
    sim_dict, fit_params, param_hash = args

    try:
        workdir = sim_dict["workdir"]
        output_file = os.path.join(workdir, f"simulation_{param_hash}.pkl")

        # Check if already completed
        if os.path.exists(output_file):
            logger.info(
                f"Simulation already exists for {sim_dict['protocol_id']}: {output_file}"
            )
            return sim_dict["protocol_id"], output_file

        logger.info(
            f"Running prefire simulation for {sim_dict['protocol_id']} in {workdir}"
        )

        # Build parameter arguments (tau_effca_GB_GluSynapse is hardcoded in pairrunner.py)
        param_args = [f"--{name}={value}" for name, value in fit_params.items()]
        param_args.append(f"--fastforward={sim_dict['fastforward']}")
        param_args.append(f"--param_hash={param_hash}")

        # Run pairrunner.py in subprocess
        script_path = os.path.join(
            os.path.dirname(__file__), "plastyfire", "pairrunner.py"
        )
        cmd = [sys.executable, script_path] + param_args

        result = subprocess.run(
            cmd,
            cwd=workdir,
            capture_output=True,
            text=True,
            timeout=1800,  # 30 minute timeout per simulation
        )

        if result.returncode != 0:
            logger.error(
                f"Simulation failed for {sim_dict['protocol_id']}: {result.stderr}"
            )
            return None

        # Verify output file was created
        if not os.path.exists(output_file):
            logger.error(
                f"Output file not created for {sim_dict['protocol_id']}: {output_file}"
            )
            return None

        logger.info(
            f"Simulation completed for {sim_dict['protocol_id']}: {output_file}"
        )
        return sim_dict["protocol_id"], output_file

    except Exception as e:
        logger.error(f"Simulation failed for {sim_dict['protocol_id']}: {e}")
        import traceback

        traceback.print_exc()
        return None


def generate_ephys_data_worker(args):
    """Worker function to generate a single ephys file"""
    sim_config_path, pre_gid, post_gid, rho_str, ephys_file = args

    if os.path.exists(ephys_file):
        logger.info(f"Ephys file already exists: {ephys_file}")
        return ephys_file

    try:
        logger.info(f"Generating ephys data for {pre_gid}->{post_gid}: {rho_str}")

        script_path = os.path.join(
            os.path.dirname(__file__), "plastyfire", "compute_test_pulse_epsp_ratio.py"
        )

        # Run compute_test_pulse_epsp_ratio.py
        cmd = [
            sys.executable,
            script_path,
            sim_config_path,
            str(pre_gid),
            str(post_gid),
            rho_str,
            "--node_pop",
            "S1nonbarrel_neurons",
            "--trial",
            "0",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

        if result.returncode != 0:
            logger.error(
                f"Ephys generation failed for {pre_gid}->{post_gid}: {result.stderr}"
            )
            return None

        logger.info(f"Ephys file generated: {ephys_file}")
        return ephys_file

    except Exception as e:
        logger.error(f"Error generating ephys for {pre_gid}->{post_gid}: {e}")
        import traceback

        traceback.print_exc()
        return None


def generate_ephys_data(sim_dict, param_hash, fit_params):
    """Generate ephys data for initial and final rho states - returns job info for parallel processing"""
    workdir = sim_dict["workdir"]
    pkl_file = os.path.join(workdir, f"simulation_{param_hash}.pkl")

    if not os.path.exists(pkl_file):
        logger.error(f"Simulation pickle file not found: {pkl_file}")
        return []

    # Load simulation results to get rho states
    with open(pkl_file, "rb") as f:
        raw_results = pickle.load(f)

    # Get initial and final rho values
    rho_data = raw_results["rho_GB"]
    if len(rho_data) > 100:
        rho_data = np.transpose(rho_data)

    initial_rho = [0 if k[0] < 0.5 else 1 for k in rho_data]
    final_rho = [0 if k[-1] < 0.5 else 1 for k in rho_data]

    # Get pre_gid and post_gid from workdir path
    # Path structure: .../simulations/191917-186517/10Hz_-10ms
    # So GIDs are in the second-to-last component
    pre_gid, post_gid = workdir.split("/")[-2].split("-")

    # Generate ephys data directory
    ephys_dir = os.path.join(os.path.dirname(workdir), "..", "..", "..", "ephys_data")
    os.makedirs(ephys_dir, exist_ok=True)

    # Prepare jobs for both initial and final states
    sim_config_path = os.path.join(workdir, "simulation_config.json")
    jobs = []

    for rho_config, label in [(initial_rho, "initial"), (final_rho, "final")]:
        rho_str = ",".join(map(str, rho_config))
        ephys_file = os.path.join(
            ephys_dir,
            f"ephys_data_{pre_gid}_{post_gid}_{rho_str.replace(',', '_')}.pkl",
        )
        jobs.append((sim_config_path, pre_gid, post_gid, rho_str, ephys_file))

    return jobs


def compute_epsp_ratios(sim_dict, param_hash):
    """Compute EPSP ratios from ephys data"""
    workdir = sim_dict["workdir"]
    # Path structure: .../simulations/191917-186517/10Hz_-10ms
    # So GIDs are in the second-to-last component
    pre_gid, post_gid = workdir.split("/")[-2].split("-")

    # Load simulation results to get rho states
    pkl_file = os.path.join(workdir, f"simulation_{param_hash}.pkl")
    with open(pkl_file, "rb") as f:
        raw_results = pickle.load(f)

    rho_data = raw_results["rho_GB"]
    if len(rho_data) > 100:
        rho_data = np.transpose(rho_data)

    initial_rho = [0 if k[0] < 0.5 else 1 for k in rho_data]
    final_rho = [0 if k[-1] < 0.5 else 1 for k in rho_data]

    # Build ephys file paths
    ephys_dir = os.path.join(os.path.dirname(workdir), "..", "..", "..", "ephys_data")
    initial_rho_str = "_".join(map(str, initial_rho))
    final_rho_str = "_".join(map(str, final_rho))

    ephys_before_path = os.path.join(
        ephys_dir, f"ephys_data_{pre_gid}_{post_gid}_{initial_rho_str}.pkl"
    )
    ephys_after_path = os.path.join(
        ephys_dir, f"ephys_data_{pre_gid}_{post_gid}_{final_rho_str}.pkl"
    )

    if not os.path.exists(ephys_before_path) or not os.path.exists(ephys_after_path):
        logger.error(
            f"Ephys files not found: {ephys_before_path} or {ephys_after_path}"
        )
        return None

    try:
        # Compute EPSP values
        epsp_before = compute_epsp_prefire(ephys_before_path)
        epsp_after = compute_epsp_prefire(ephys_after_path)
        epsp_ratio = epsp_after / epsp_before

        logger.info(f"EPSP ratio for {sim_dict['protocol_id']}: {epsp_ratio:.4f}")

        return {
            "protocol_id": sim_dict["protocol_id"],
            "epsp_ratio": epsp_ratio,
            "epsp_before": epsp_before,
            "epsp_after": epsp_after,
            "ephys_before_path": ephys_before_path,
            "ephys_after_path": ephys_after_path,
            "pre_gid": pre_gid,
            "post_gid": post_gid,
        }
    except RuntimeError as e:
        logger.warning(
            f"Failed to compute EPSP ratio for {sim_dict['protocol_id']} ({pre_gid}-{post_gid}): {e}"
        )
        return None
    except Exception as e:
        logger.error(
            f"Unexpected error computing EPSP ratio for {sim_dict['protocol_id']} ({pre_gid}-{post_gid}): {e}"
        )
        import traceback

        traceback.print_exc()
        return None


def main():
    parser = argparse.ArgumentParser(description="Evaluate best solution")
    parser.add_argument(
        "--sample-size", type=int, default=100, help="Number of pairs to evaluate"
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=60,
        help="Number of parallel workers (default: 60)",
    )
    parser.add_argument("--seed", type=int, default=1234, help="Random seed")
    parser.add_argument(
        "--sim-workers",
        type=int,
        default=None,
        help="Workers for simulations (default: workers)",
    )
    parser.add_argument(
        "--ephys-workers",
        type=int,
        default=None,
        help="Workers for ephys generation (default: workers)",
    )
    parser.add_argument(
        "--partial",
        action="store_true",
        help="Skip ephys generation and use only existing ephys files",
    )
    parser.add_argument(
        "--param-set",
        type=str,
        default="best",
        choices=["best"] + list(CUSTOM_PARAM_SETS.keys()),
        help=f"Parameter set to evaluate: 'best' for optimized solution, or one of {list(CUSTOM_PARAM_SETS.keys())}",
    )
    parser.add_argument(
        "--compare-all",
        action="store_true",
        help="Evaluate all parameter sets (best + all custom sets) and compare results",
    )
    args = parser.parse_args()

    # Set worker counts for each step
    sim_workers = args.sim_workers if args.sim_workers else args.workers
    ephys_workers = args.ephys_workers if args.ephys_workers else args.workers

    logger.info(
        f"Configuration: {args.workers} default workers, {sim_workers} sim workers, {ephys_workers} ephys workers"
    )

    # Determine which parameter sets to evaluate
    if args.compare_all:
        # Evaluate all parameter sets
        param_sets_to_eval = {}
        try:
            param_sets_to_eval["best"] = load_best_solution()
        except Exception as e:
            logger.warning(f"Failed to load best solution: {e}")

        param_sets_to_eval.update(CUSTOM_PARAM_SETS)
        logger.info(f"Comparing all parameter sets: {list(param_sets_to_eval.keys())}")
    else:
        # Evaluate single parameter set
        if args.param_set == "best":
            try:
                current_params = load_best_solution()
                logger.info(f"Using best parameters: {current_params}")
            except Exception as e:
                logger.error(f"Failed to load best solution: {e}")
                return 1
        else:
            current_params = CUSTOM_PARAM_SETS[args.param_set]
            logger.info(f"Using custom parameter set '{args.param_set}': {current_params}")

        param_sets_to_eval = {args.param_set: current_params}

    # Load biodata
    invitro_db = pd.read_csv(CSVF_NAME)
    invitro_db = invitro_db.loc[invitro_db["protocol_id"].isin(EVAL_PROTOCOLS)]

    logger.info(f"Evaluating protocols: {EVAL_PROTOCOLS}")
    logger.info(
        f"In vitro data:\n{invitro_db[['protocol_id', 'mean_epsp_ratio', 'sem_epsp_ratio']]}"
    )

    # Find all simulations for evaluation protocols (do this once for all parameter sets)
    all_sims = []
    for elem in invitro_db.itertuples():
        # Load simulation config
        config = OptConfig(
            os.path.join(CONFIGS_DIR, f"{elem.pre_mtype}_{elem.post_mtype}.yaml")
        )

        fastforward = config.fastforward
        if fastforward is None:
            fastforward = config.C01_duration * MIN2MS + config.nreps * config.T
        fastforward = float(fastforward)

        # Load simulation index
        index_file = os.path.join(
            os.path.split(os.path.split(config.out_dir)[0])[0],
            f"index_{elem.pre_mtype}_{elem.post_mtype}.csv",
        )

        if not os.path.exists(index_file):
            logger.error(f"Index file not found: {index_file}")
            continue

        sim_idx = pd.read_csv(index_file)
        sim_idx.set_index(["frequency", "dt"], inplace=True)
        sim_idx.sort_index(inplace=True)

        paths = sim_idx.loc[elem.frequency_train, elem.dt_train]["path"]
        if isinstance(paths, str):
            paths = pd.Series([paths])

        # Sample paths
        np.random.seed(args.seed)
        sample_size = min(args.sample_size, len(paths))
        paths = paths.sample(sample_size, random_state=args.seed)

        logger.info(f"Protocol {elem.protocol_id}: Found {len(paths)} simulation paths")

        for path in paths:
            all_sims.append(
                {
                    "protocol_id": elem.protocol_id,
                    "fastforward": fastforward,
                    "workdir": "/".join(path.split("/")[:-1]),
                }
            )

    logger.info(f"Total simulations to run: {len(all_sims)}")

    # Store all results for comparison
    all_param_set_results = {}
    all_comparisons = {}

    # Evaluate each parameter set
    for param_set_name, params in param_sets_to_eval.items():
        logger.info("\n" + "=" * 80)
        logger.info(f"EVALUATING PARAMETER SET: {param_set_name}")
        logger.info("=" * 80)
        logger.info(f"Parameters: {params}")

        # Create hash for these params
        param_values = [params[name] for name in FIT_PARAM_NAMES]
        param_hash = hashlib.md5(str(param_values).encode()).hexdigest()[:12]
        logger.info(f"Parameter hash: {param_hash}")

        # Step 1: Run prefire simulations in parallel using subprocess
        logger.info("=" * 60)
        logger.info(f"STEP 1: Running prefire simulations ({sim_workers} workers)")
        logger.info("=" * 60)

        sim_args = [(sim_dict, params, param_hash) for sim_dict in all_sims]

        logger.info(
            f"Starting {len(sim_args)} prefire simulations with {sim_workers} workers..."
        )
        step1_start = time.time()
        with mp.Pool(processes=sim_workers) as pool:
            sim_results = pool.map(run_prefire_simulation, sim_args)
        step1_time = time.time() - step1_start

        # Filter out failed simulations
        sim_results = [r for r in sim_results if r is not None]
        logger.info(
            f"Completed {len(sim_results)} / {len(all_sims)} simulations in {step1_time / 60:.2f} minutes"
        )

        if len(sim_results) == 0:
            logger.error(f"All simulations failed for {param_set_name}! Skipping.")
            continue

        # Step 2: Generate ephys data in parallel (unless --partial mode)
        if args.partial:
            logger.info("=" * 60)
            logger.info("STEP 2: Skipping ephys generation (--partial mode)")
            logger.info("=" * 60)
            logger.info("Will use only existing ephys files for EPSP ratio computation")
            step2_time = 0
        else:
            logger.info("=" * 60)
            logger.info(f"STEP 2: Generating ephys data ({ephys_workers} workers)")
            logger.info("=" * 60)

            # Collect all ephys jobs from all simulations
            all_ephys_jobs = []
            for sim_dict in all_sims:
                jobs = generate_ephys_data(sim_dict, param_hash, params)
                all_ephys_jobs.extend(jobs)

            logger.info(
                f"Starting {len(all_ephys_jobs)} ephys generation jobs with {ephys_workers} workers..."
            )
            step2_start = time.time()
            with mp.Pool(processes=ephys_workers) as pool:
                ephys_results = pool.map(generate_ephys_data_worker, all_ephys_jobs)
            step2_time = time.time() - step2_start

            # Filter out failed ephys generations
            ephys_results = [r for r in ephys_results if r is not None]
            logger.info(
                f"Completed {len(ephys_results)} / {len(all_ephys_jobs)} ephys files in {step2_time / 60:.2f} minutes"
            )

        # Step 3: Compute EPSP ratios in parallel
        logger.info("=" * 60)
        logger.info(f"STEP 3: Computing EPSP ratios ({args.workers} workers)")
        logger.info("=" * 60)

        epsp_args = [(sim_dict, param_hash) for sim_dict in all_sims]

        logger.info(f"Computing EPSP ratios for {len(epsp_args)} simulations...")
        step3_start = time.time()
        with mp.Pool(processes=args.workers) as pool:
            epsp_results = pool.starmap(compute_epsp_ratios, epsp_args)
        step3_time = time.time() - step3_start

        # Filter out failed computations
        epsp_results = [r for r in epsp_results if r is not None]
        logger.info(
            f"Computed {len(epsp_results)} / {len(all_sims)} EPSP ratios in {step3_time / 60:.2f} minutes"
        )

        # Create results dataframe
        results_df = pd.DataFrame(epsp_results)

        # Compute statistics per protocol
        logger.info("=" * 60)
        logger.info(f"STEP 4: Results Summary for {param_set_name}")
        logger.info("=" * 60)

        insilico_stats = results_df.groupby("protocol_id")["epsp_ratio"].agg(
            ["mean", "sem", "std", "count"]
        )
        insilico_stats.columns = [
            "mean_epsp_ratio_sim",
            "sem_epsp_ratio_sim",
            "std_epsp_ratio_sim",
            "count",
        ]

        # Merge with in vitro data
        comparison = pd.merge(
            invitro_db[
                [
                    "protocol_id",
                    "mean_epsp_ratio",
                    "sem_epsp_ratio",
                    "frequency_train",
                    "dt_train",
                ]
            ],
            insilico_stats,
            on="protocol_id",
        )

        # Compute error
        comparison["error"] = np.abs(
            (comparison["mean_epsp_ratio"] - comparison["mean_epsp_ratio_sim"])
            / comparison["sem_epsp_ratio"]
        )
        comparison["param_set"] = param_set_name

        logger.info(f"\nComparison of In Vitro vs In Silico Results ({param_set_name}):")
        logger.info("\n" + comparison.to_string())

        # Store results for this parameter set
        all_param_set_results[param_set_name] = results_df
        all_comparisons[param_set_name] = comparison

        # Save results for this parameter set
        output_dir = f"/lustre06/project/6077694/dhuruva/plastyfire/evaluation_results/{param_set_name}"
        os.makedirs(output_dir, exist_ok=True)

        results_df.to_csv(os.path.join(output_dir, "detailed_results.csv"), index=False)
        comparison.to_csv(os.path.join(output_dir, "comparison_summary.csv"), index=False)

        # Save params
        with open(os.path.join(output_dir, "params.pkl"), "wb") as f:
            pickle.dump(params, f)

        logger.info(f"\nResults for {param_set_name} saved to {output_dir}")

        # Print timing summary for this parameter set
        total_time = step1_time + step2_time + step3_time
        logger.info(f"\nTiming for {param_set_name}:")
        logger.info(f"  Step 1 (Prefire Simulations): {step1_time / 60:.2f} minutes")
        logger.info(f"  Step 2 (Ephys Generation):    {step2_time / 60:.2f} minutes")
        logger.info(f"  Step 3 (EPSP Computation):    {step3_time / 60:.2f} minutes")
        logger.info(f"  Total Time:                   {total_time / 60:.2f} minutes")

    # Final comparison across all parameter sets
    if len(all_comparisons) > 1:
        logger.info("\n" + "=" * 80)
        logger.info("FINAL COMPARISON ACROSS ALL PARAMETER SETS")
        logger.info("=" * 80)

        combined_comparison = pd.concat(all_comparisons.values(), ignore_index=True)

        # Create pivot table for easy comparison
        pivot_comparison = combined_comparison.pivot_table(
            index=["protocol_id", "frequency_train", "dt_train"],
            columns="param_set",
            values=["mean_epsp_ratio_sim", "error"]
        )

        logger.info("\nMean EPSP Ratios by Parameter Set:")
        logger.info("\n" + pivot_comparison["mean_epsp_ratio_sim"].to_string())
        logger.info("\nErrors (normalized by SEM) by Parameter Set:")
        logger.info("\n" + pivot_comparison["error"].to_string())

        # Save combined comparison
        combined_output_dir = "/lustre06/project/6077694/dhuruva/plastyfire/evaluation_results"
        combined_comparison.to_csv(
            os.path.join(combined_output_dir, "combined_comparison.csv"), index=False
        )
        logger.info(f"\nCombined comparison saved to {combined_output_dir}/combined_comparison.csv")

    # Plot comparisons
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # Plot for each parameter set individually
        for param_set_name, results_df in all_param_set_results.items():
            fig, axes = plt.subplots(1, 3, figsize=(15, 5))

            for idx, protocol_id in enumerate(EVAL_PROTOCOLS):
                ax = axes[idx]

                # Get data for this protocol
                protocol_results = results_df[results_df["protocol_id"] == protocol_id][
                    "epsp_ratio"
                ]
                invitro_value = invitro_db[invitro_db["protocol_id"] == protocol_id][
                    "mean_epsp_ratio"
                ].iloc[0]

                # Plot histogram
                ax.hist(
                    protocol_results, bins=20, alpha=0.7, color="blue", edgecolor="black"
                )

                # Add vertical lines
                ax.axvline(
                    invitro_value,
                    color="red",
                    linestyle="--",
                    linewidth=2,
                    label=f"In vitro: {invitro_value:.3f}",
                )
                ax.axvline(
                    protocol_results.mean(),
                    color="green",
                    linestyle="-",
                    linewidth=2,
                    label=f"In silico: {protocol_results.mean():.3f}",
                )

                # Labels
                freq = invitro_db[invitro_db["protocol_id"] == protocol_id][
                    "frequency_train"
                ].iloc[0]
                dt = invitro_db[invitro_db["protocol_id"] == protocol_id]["dt_train"].iloc[
                    0
                ]
                ax.set_title(f"{protocol_id}\nf={freq}Hz, dt={dt}ms")
                ax.set_xlabel("EPSP Ratio")
                ax.set_ylabel("Count")
                ax.legend()
                ax.grid(alpha=0.3)

            plt.suptitle(f"Parameter Set: {param_set_name}", fontsize=14, fontweight="bold")
            plt.tight_layout()

            param_output_dir = f"/lustre06/project/6077694/dhuruva/plastyfire/evaluation_results/{param_set_name}"
            plt.savefig(
                os.path.join(param_output_dir, "evaluation_comparison.png"),
                dpi=150,
                bbox_inches="tight",
            )
            logger.info(f"Plot saved to {param_output_dir}/evaluation_comparison.png")
            plt.close()

        # If comparing multiple parameter sets, create a comparison plot
        if len(all_param_set_results) > 1:
            fig, axes = plt.subplots(1, 3, figsize=(18, 6))

            colors = plt.cm.Set1(np.linspace(0, 1, len(all_param_set_results)))

            for idx, protocol_id in enumerate(EVAL_PROTOCOLS):
                ax = axes[idx]

                # Get in vitro value
                invitro_value = invitro_db[invitro_db["protocol_id"] == protocol_id][
                    "mean_epsp_ratio"
                ].iloc[0]

                # Plot vertical line for in vitro
                ax.axvline(
                    invitro_value,
                    color="black",
                    linestyle="--",
                    linewidth=2.5,
                    label=f"In vitro: {invitro_value:.3f}",
                    zorder=10
                )

                # Plot mean values for each parameter set
                for color_idx, (param_set_name, results_df) in enumerate(all_param_set_results.items()):
                    protocol_results = results_df[results_df["protocol_id"] == protocol_id]["epsp_ratio"]
                    mean_val = protocol_results.mean()

                    ax.axvline(
                        mean_val,
                        color=colors[color_idx],
                        linestyle="-",
                        linewidth=2,
                        label=f"{param_set_name}: {mean_val:.3f}",
                        alpha=0.8
                    )

                # Labels
                freq = invitro_db[invitro_db["protocol_id"] == protocol_id]["frequency_train"].iloc[0]
                dt = invitro_db[invitro_db["protocol_id"] == protocol_id]["dt_train"].iloc[0]
                ax.set_title(f"{protocol_id}\nf={freq}Hz, dt={dt}ms")
                ax.set_xlabel("EPSP Ratio")
                ax.set_ylabel("Density")
                ax.legend(fontsize=8)
                ax.grid(alpha=0.3)

            plt.suptitle("Comparison Across All Parameter Sets", fontsize=14, fontweight="bold")
            plt.tight_layout()

            combined_output_dir = "/lustre06/project/6077694/dhuruva/plastyfire/evaluation_results"
            plt.savefig(
                os.path.join(combined_output_dir, "combined_comparison.png"),
                dpi=150,
                bbox_inches="tight",
            )
            logger.info(f"Combined plot saved to {combined_output_dir}/combined_comparison.png")
            plt.close()

    except Exception as e:
        logger.warning(f"Failed to create plot: {e}")
        import traceback
        traceback.print_exc()

    logger.info("\n" + "=" * 80)
    logger.info("EVALUATION COMPLETED SUCCESSFULLY!")
    logger.info("=" * 80)
    return 0


if __name__ == "__main__":
    sys.exit(main())
