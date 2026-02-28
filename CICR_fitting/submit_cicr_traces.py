#!/usr/bin/env python3
"""
Submit L5TTPC_L5TTPC plasticity simulations with Chindemi parameters to SLURM,
recording synaptic traces (effcai_GB, cai_CR, vsyn, rho_GB) during induction.

Usage:
    # Submit all 100 pairs as SLURM jobs (overnight run)
    python submit_l5ttpc_traces.py
    
    # Test with 2 pairs first
    python submit_l5ttpc_traces.py --max-pairs 2
    
    # Use multiprocessing instead of SLURM
    python submit_l5ttpc_traces.py --execution-mode cpu --workers 30

Author: Generated for trace recording
"""

import argparse
import logging
import os
import pickle
import subprocess
import sys
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

import json

# Base directories
SIMULATIONS_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC/simulations"
OUTPUT_BASE_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/trace_results_CICR"
PLASTYFIRE_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire"

# Primed Junctional parameters
PARAM_FILE = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/jax_fitting/best_params_primed_junctional.json"
with open(PARAM_FILE, "r") as f:
    CHINDEMI_PARAMS = json.load(f)

def _submitit_worker(batch, output_dir, cpus_per_task):
    """Worker function executed by submitit on the SLURM node."""
    from multiprocessing import Pool
    import time
    
    # Add output_dir to each job tuple
    jobs_with_output = [(pair_dir, freq_dt, output_dir) for pair_dir, freq_dt in batch]
    
    logger.info(f"Starting submitit batch of {len(batch)} simulations with {cpus_per_task} parallel workers...")
    
    import submit_cicr_traces
    with Pool(processes=cpus_per_task) as pool:
        results = pool.map(submit_cicr_traces._cpu_worker, jobs_with_output)
        
    success = sum(1 for r in results if r["status"] == "success")
    skipped = sum(1 for r in results if r["status"] == "skipped")
    failed = sum(1 for r in results if r["status"] in ["failed", "error"])
    
    logger.info(f"Batch completed. Success: {success}, Skipped: {skipped}, Failed: {failed}")
    return results


def submit_slurm_jobs(jobs, output_dir, timeout_hours=12, mem_gb=8, cpus_per_task=4,
                       slurm_account="ctb-emuller", batch_size=10):
    """Submit jobs to SLURM using submitit."""
    try:
        import submitit
    except ImportError:
        logger.error("submitit module not found. Please install it: pip install submitit")
        sys.exit(1)
        
    log_folder = Path(output_dir) / "logs"
    log_folder.mkdir(parents=True, exist_ok=True)
    
    # Create batches
    batches = []
    for i in range(0, len(jobs), batch_size):
        batch = jobs[i:i + batch_size]
        batches.append(batch)
    
    logger.info(f"Submitting {len(batches)} SLURM jobs (each processing up to {batch_size} traces) for {len(jobs)} total simulations...")
    
    executor = submitit.AutoExecutor(folder=str(log_folder))
    executor.update_parameters(
        timeout_min=timeout_hours * 60,
        cpus_per_task=cpus_per_task,
        mem_gb=mem_gb,
        slurm_account=slurm_account,
        slurm_setup=[f"source {PLASTYFIRE_DIR}/setupenv.sh"]
    )
    
    with executor.batch():
        submitted_jobs = [executor.submit(_submitit_worker, batch, output_dir, cpus_per_task) for batch in batches]
    
    logger.info(f"\n{'='*60}")
    logger.info(f"Submitted {len(submitted_jobs)} SLURM jobs via submitit")
    logger.info(f"Logs: {log_folder}/")
    logger.info(f"Results: {output_dir}/")
    logger.info(f"\nMonitor with: squeue -u $USER (or submitit tools)")
    
    return [j.job_id for j in submitted_jobs]


def _cpu_worker(args):
    """Worker function for CPU multiprocessing mode (must be at module level for pickling)."""
    pair_dir, freq_dt, output_dir = args
    workdir = os.path.join(pair_dir, freq_dt)
    pair_name = os.path.basename(pair_dir)
    output_file = os.path.join(output_dir, pair_name, freq_dt, "simulation_traces.pkl")
    
    if os.path.exists(output_file):
        return {"status": "skipped", "pair": pair_name}
    
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    param_args = [f"--{name}={value}" for name, value in CHINDEMI_PARAMS.items()]
    param_args.append("--output-filename=simulation_traces_temp.pkl")
    
    env = os.environ.copy()
    current_paths = ":".join(sys.path)
    env["PYTHONPATH"] = f"{current_paths}:{env.get('PYTHONPATH', '')}:{PLASTYFIRE_DIR}"
    # env["BLUECELLULAB_MOD_LIBRARY_PATH"] = "/home/dhuruva/projects/ctb-emuller/dhuruva/DEES_cell_packages/x86_64/libnrnmech.so"
    # env["BBP_MECH_DIR"] = "/home/dhuruva/projects/ctb-emuller/dhuruva/DEES_cell_packages"

    cmd = [sys.executable, f"{PLASTYFIRE_DIR}/CICR_fitting/pairrunner_cicr.py"] + param_args
    
    try:
        result = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True, timeout=7200, env=env)
        
        temp_file = os.path.join(workdir, "simulation_traces_temp.pkl")
        if os.path.exists(temp_file):
            import shutil
            shutil.move(temp_file, output_file)
            return {"status": "success", "pair": pair_name, "freq_dt": freq_dt}
        else:
            # Log the error for debugging
            return {
                "status": "failed", 
                "pair": pair_name, 
                "freq_dt": freq_dt,
                "returncode": result.returncode,
                "stderr": result.stderr[-2500:] if result.stderr else "",
                "stdout": result.stdout[-2500:] if result.stdout else "",
                "cmd": " ".join(cmd[:3]) + "..."  # First 3 parts of command
            }
    except Exception as e:
        return {"status": "error", "pair": pair_name, "freq_dt": freq_dt, "error": str(e)}


def run_cpu_mode(jobs, output_dir, workers=30):
    """Run simulations locally using multiprocessing."""
    from multiprocessing import Pool
    import time
    
    logger.info(f"Running {len(jobs)} simulations with {workers} workers...")
    
    # Add output_dir to each job tuple
    jobs_with_output = [(pair_dir, freq_dt, output_dir) for pair_dir, freq_dt in jobs]
    
    start = time.time()
    with Pool(workers) as pool:
        results = pool.map(_cpu_worker, jobs_with_output)
    elapsed = time.time() - start
    
    success = sum(1 for r in results if r["status"] == "success")
    skipped = sum(1 for r in results if r["status"] == "skipped")
    failed = sum(1 for r in results if r["status"] in ["failed", "error"])
    
    logger.info(f"\n{'='*60}")
    logger.info(f"Completed in {elapsed/60:.1f} minutes")
    logger.info(f"Success: {success}, Skipped: {skipped}, Failed: {failed}")

    if failed > 0:
        logger.error("-" * 40)
        for r in results:
            if r["status"] == "failed":
                logger.error(f"Failed {r.get('pair')}/{r.get('freq_dt')}: \n{r.get('stderr', '')[:500]}")
            elif r["status"] == "error":
                logger.error(f"Error {r.get('pair')}/{r.get('freq_dt')}: {r.get('error')}")



def main():
    parser = argparse.ArgumentParser(
        description="Submit L5TTPC trace simulations with Chindemi params",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Submit all to SLURM (overnight run)
  python submit_l5ttpc_traces.py
  
  # Test with 2 pairs
  python submit_l5ttpc_traces.py --max-pairs 2
  
  # Run locally with multiprocessing
  python submit_l5ttpc_traces.py --execution-mode cpu --workers 30
        """
    )
    
    parser.add_argument("--execution-mode", choices=["slurm", "cpu"], default="slurm",
                        help="Execution mode (default: slurm)")
    parser.add_argument("--max-pairs", type=int, default=None,
                        help="Max pairs to process (for testing)")
    parser.add_argument("--freq", type=str, default=None,
                        help="Filter by frequency, e.g., '10Hz'")
    parser.add_argument("--output-dir", type=str, default=None,
                        help="Custom output directory")
    
    # CPU mode options
    parser.add_argument("--workers", type=int, default=30,
                        help="Workers for CPU mode (default: 30)")
    
    # SLURM options
    parser.add_argument("--timeout-hours", type=int, default=1,
                        help="SLURM timeout in hours (default: 1)")
    parser.add_argument("--mem-gb", type=int, default=16,
                        help="Memory per job in GB (default: 16)")
    parser.add_argument("--cpus-per-task", type=int, default=0,
                        help="CPUs per SLURM job (0 to match batch-size)")
    parser.add_argument("--batch-size", type=int, default=18,
                        help="Simulations per SLURM job (default: 18)")
    parser.add_argument("--slurm-account", default="ctb-emuller",
                        help="SLURM account (default: ctb-emuller)")
    
    args = parser.parse_args()
    
    # Output directory
    output_dir = args.output_dir or os.path.join(OUTPUT_BASE_DIR, "PrimedJunctional_params")
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Output directory: {output_dir}")
    
    # Find all simulation directories
    simulations_path = Path(SIMULATIONS_DIR)
    pair_dirs = sorted([d for d in simulations_path.iterdir() 
                       if d.is_dir() and d.name != "single_cells"])
    
    if args.max_pairs:
        pair_dirs = pair_dirs[:args.max_pairs]
    
    logger.info(f"Found {len(pair_dirs)} pair directories")
    
    # Filter protocols
    allowed_freqs = ["2Hz_", "5Hz_", "10Hz_", "20Hz_", "30Hz_", "40Hz_"]
    allowed_dts = ["_-10ms", "_5ms", "_10ms"]
    
    # Build job list
    jobs = []
    for pair_dir in pair_dirs:
        freq_dt_dirs = [d.name for d in pair_dir.iterdir() if d.is_dir()]
        
        if args.freq:
            freq_dt_dirs = [d for d in freq_dt_dirs if args.freq in d]
        
        for freq_dt in freq_dt_dirs:
            # Filter exactly 2, 5, 10, 20, 30, 40Hz with -10, 5, 10ms
            matches_freq = any(freq_dt.startswith(f) for f in allowed_freqs)
            matches_dt = any(dt in freq_dt for dt in allowed_dts)
            
            if matches_freq and matches_dt:
                jobs.append((str(pair_dir), freq_dt))
    
    logger.info(f"Total simulations: {len(jobs)}")
    
    # Save job metadata
    metadata = {
        "params": CHINDEMI_PARAMS,
        "jobs": jobs,
        "output_dir": output_dir,
    }
    with open(os.path.join(output_dir, "job_metadata.pkl"), "wb") as f:
        pickle.dump(metadata, f)
    
    # Execute
    final_cpus = args.cpus_per_task if args.cpus_per_task > 0 else args.batch_size
    
    if args.execution_mode == "slurm":
        submit_slurm_jobs(
            jobs, output_dir,
            timeout_hours=args.timeout_hours,
            mem_gb=args.mem_gb,
            cpus_per_task=final_cpus,
            batch_size=args.batch_size,
            slurm_account=args.slurm_account
        )
    else:
        run_cpu_mode(jobs, output_dir, workers=args.workers)


if __name__ == "__main__":
    main()
