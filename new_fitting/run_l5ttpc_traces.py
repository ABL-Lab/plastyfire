#!/usr/bin/env python3
"""
Run L5TTPC_L5TTPC plasticity simulations with Chindemi parameters,
recording synaptic traces (effcai_GB, cai_CR, vsyn, rho_GB) during induction.

Usage:
    python run_l5ttpc_traces.py --workers 30
    python run_l5ttpc_traces.py --max-pairs 2 --workers 2  # Test run

Author: Generated for trace recording
"""

import argparse
import logging
import multiprocessing as mp
import os
import pickle
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Base directories
SIMULATIONS_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC/simulations"
OUTPUT_BASE_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/trace_results"

# Chindemi parameters (from evaluate_best_solution.py)
CHINDEMI_PARAMS = {
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
}

# Traces that will be saved (already recorded by simulator.py SYNREC)
TRACE_KEYS = ["rho_GB", "effcai_GB", "cai_CR", "vsyn", "Use_GB", "gmax_AMPA", "ica_NMDA", "ica_VDCC"]


def run_single_simulation(args):
    """Worker function to run a single simulation"""
    pair_dir, freq_dt_dir, output_dir = args
    
    workdir = os.path.join(pair_dir, freq_dt_dir)
    pair_name = os.path.basename(pair_dir)
    output_file = os.path.join(output_dir, pair_name, freq_dt_dir, "simulation_traces.pkl")
    
    # Skip if already completed
    if os.path.exists(output_file):
        logger.info(f"Skipping {pair_name}/{freq_dt_dir} - already exists")
        return {"status": "skipped", "pair": pair_name, "freq_dt": freq_dt_dir}
    
    # Create output directory
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    logger.info(f"Running {pair_name}/{freq_dt_dir}")
    
    try:
        # Build pairrunner command with Chindemi parameters
        script_path = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/plastyfire/pairrunner.py"
        
        param_args = [f"--{name}={value}" for name, value in CHINDEMI_PARAMS.items()]
        
        # Temporary output file in workdir
        temp_output = os.path.join(workdir, "simulation_traces_temp.pkl")
        param_args.append(f"--output-filename=simulation_traces_temp.pkl")
        
        cmd = [sys.executable, script_path] + param_args
        
        result = subprocess.run(
            cmd,
            cwd=workdir,
            capture_output=True,
            text=True,
            timeout=7200  # 2 hour timeout
        )
        
        if result.returncode != 0:
            logger.error(f"Simulation failed for {pair_name}/{freq_dt_dir}: {result.stderr[:500]}")
            return {"status": "failed", "pair": pair_name, "freq_dt": freq_dt_dir, "error": result.stderr[:500]}
        
        # Load results and extract only the traces we need
        if os.path.exists(temp_output):
            with open(temp_output, "rb") as f:
                full_results = pickle.load(f)
            
            # Extract traces and essential data
            trace_results = {
                "t": full_results.get("t"),
                "v": full_results.get("v"),
                "prespikes": full_results.get("prespikes"),
                "postspikes": full_results.get("postspikes"),
                "synprop": full_results.get("synprop"),
                "params": CHINDEMI_PARAMS,
            }
            
            # Add all trace keys
            for key in TRACE_KEYS:
                if key in full_results:
                    trace_results[key] = full_results[key]
            
            # Save to output directory
            with open(output_file, "wb") as f:
                pickle.dump(trace_results, f)
            
            # Clean up temp file
            os.remove(temp_output)
            
            logger.info(f"Completed {pair_name}/{freq_dt_dir} - saved to {output_file}")
            return {"status": "success", "pair": pair_name, "freq_dt": freq_dt_dir, "output": output_file}
        else:
            logger.error(f"Output file not created for {pair_name}/{freq_dt_dir}")
            return {"status": "failed", "pair": pair_name, "freq_dt": freq_dt_dir, "error": "No output file"}
    
    except subprocess.TimeoutExpired:
        logger.error(f"Timeout for {pair_name}/{freq_dt_dir}")
        return {"status": "timeout", "pair": pair_name, "freq_dt": freq_dt_dir}
    except Exception as e:
        logger.error(f"Error for {pair_name}/{freq_dt_dir}: {e}")
        return {"status": "error", "pair": pair_name, "freq_dt": freq_dt_dir, "error": str(e)}


def main():
    parser = argparse.ArgumentParser(description="Run L5TTPC simulations with Chindemi params and trace recording")
    parser.add_argument("--workers", type=int, default=30, help="Number of parallel workers")
    parser.add_argument("--max-pairs", type=int, default=None, help="Max number of pairs to process (for testing)")
    parser.add_argument("--freq", type=str, default=None, help="Filter by frequency, e.g., '10Hz'")
    parser.add_argument("--output-dir", type=str, default=None, help="Custom output directory")
    args = parser.parse_args()
    
    # Output directory
    output_dir = args.output_dir or os.path.join(OUTPUT_BASE_DIR, "Chindemi_params")
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Output directory: {output_dir}")
    
    # Find all simulation directories
    simulations_path = Path(SIMULATIONS_DIR)
    pair_dirs = sorted([d for d in simulations_path.iterdir() if d.is_dir() and d.name != "single_cells"])
    
    if args.max_pairs:
        pair_dirs = pair_dirs[:args.max_pairs]
    
    logger.info(f"Found {len(pair_dirs)} pair directories")
    
    # Build job list
    jobs = []
    for pair_dir in pair_dirs:
        freq_dt_dirs = [d.name for d in pair_dir.iterdir() if d.is_dir()]
        
        if args.freq:
            freq_dt_dirs = [d for d in freq_dt_dirs if args.freq in d]
        
        for freq_dt in freq_dt_dirs:
            jobs.append((str(pair_dir), freq_dt, output_dir))
    
    logger.info(f"Total jobs: {len(jobs)}")
    
    # Run in parallel
    start_time = time.time()
    
    with mp.Pool(processes=args.workers) as pool:
        results = pool.map(run_single_simulation, jobs)
    
    elapsed = time.time() - start_time
    
    # Summary
    success = sum(1 for r in results if r["status"] == "success")
    skipped = sum(1 for r in results if r["status"] == "skipped")
    failed = sum(1 for r in results if r["status"] in ["failed", "error", "timeout"])
    
    logger.info("=" * 60)
    logger.info(f"SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Total jobs:    {len(jobs)}")
    logger.info(f"Successful:    {success}")
    logger.info(f"Skipped:       {skipped}")
    logger.info(f"Failed:        {failed}")
    logger.info(f"Time elapsed:  {elapsed/60:.2f} minutes")
    logger.info(f"Output dir:    {output_dir}")
    
    # Save summary
    summary_file = os.path.join(output_dir, "run_summary.pkl")
    with open(summary_file, "wb") as f:
        pickle.dump({
            "results": results,
            "params": CHINDEMI_PARAMS,
            "elapsed_time": elapsed,
        }, f)
    logger.info(f"Summary saved to {summary_file}")


if __name__ == "__main__":
    main()
