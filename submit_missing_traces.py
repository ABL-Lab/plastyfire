#!/usr/bin/env python3
"""
Submit trace collection for the 26 missing L5TTPC pairs.
These pairs have basis files but no traces in trace_results/Chindemi_params/

Usage:
    python submit_missing_traces.py --execution-mode slurm
    python submit_missing_traces.py --execution-mode cpu --workers 30
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

# Base directories
SIMULATIONS_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC/simulations"
OUTPUT_BASE_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/trace_results"
PLASTYFIRE_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire"

# Missing pairs (have basis but no traces)
MISSING_PAIRS = [
    "202571-195398", "203038-207594", "203171-197644", "203300-206195", "203710-206592",
    "203919-187736", "203964-196350", "204120-188861", "204194-205064", "204883-200512",
    "205300-191467", "205327-192314", "205559-199162", "205643-189193", "206314-189443",
    "206930-184305", "207423-186224", "207475-192618", "207836-184831", "208093-203811",
    "208904-200676", "208966-196199", "208976-188476", "209495-189342", "209709-184968",
    "210209-180231"
]

# Chindemi parameters (from submit_l5ttpc_traces.py)
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


def _cpu_worker(args):
    """Worker function for CPU multiprocessing mode."""
    pair_dir, freq_dt, output_dir = args
    workdir = os.path.join(pair_dir, freq_dt)
    pair_name = os.path.basename(pair_dir)
    output_file = os.path.join(output_dir, pair_name, freq_dt, "simulation_traces.pkl")
    
    if os.path.exists(output_file):
        return {"status": "skipped", "pair": pair_name, "freq_dt": freq_dt}
    
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    param_args = [f"--{name}={value}" for name, value in CHINDEMI_PARAMS.items()]
    param_args.append("--output-filename=simulation_traces_temp.pkl")
    
    cmd = [sys.executable, f"{PLASTYFIRE_DIR}/plastyfire/pairrunner.py"] + param_args
    
    try:
        result = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True, timeout=7200)
        
        temp_file = os.path.join(workdir, "simulation_traces_temp.pkl")
        if os.path.exists(temp_file):
            import shutil
            shutil.move(temp_file, output_file)
            return {"status": "success", "pair": pair_name, "freq_dt": freq_dt}
        else:
            return {
                "status": "failed", 
                "pair": pair_name, 
                "freq_dt": freq_dt,
                "returncode": result.returncode,
                "stderr": result.stderr[:500] if result.stderr else "",
            }
    except Exception as e:
        return {"status": "error", "pair": pair_name, "freq_dt": freq_dt, "error": str(e)}


def run_cpu_mode(jobs, output_dir, workers=30):
    """Run simulations locally using multiprocessing."""
    from multiprocessing import Pool
    import time
    
    logger.info(f"Running {len(jobs)} simulations with {workers} workers...")
    
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
    
    # Log failed jobs
    for r in results:
        if r["status"] in ["failed", "error"]:
            logger.error(f"Failed: {r}")


def submit_slurm_jobs(jobs, output_dir, timeout_hours=12, mem_gb=8, cpus_per_task=4,
                       slurm_account="ctb-emuller", batch_size=10):
    """Submit jobs to SLURM."""
    log_folder = Path(output_dir) / "logs_missing"
    log_folder.mkdir(parents=True, exist_ok=True)
    
    batches = [jobs[i:i + batch_size] for i in range(0, len(jobs), batch_size)]
    
    logger.info(f"Submitting {len(batches)} SLURM jobs for {len(jobs)} simulations...")
    
    time_str = f"{timeout_hours:02d}:00:00"
    job_ids = []
    
    for batch_idx, batch in enumerate(batches):
        job_name = f"missing_trace_{batch_idx:04d}"
        workdirs = [os.path.join(pair_dir, freq_dt) for pair_dir, freq_dt in batch]
        workdirs_str = " ".join(f'"{w}"' for w in workdirs)
        param_args = " ".join(f"--{name}={value}" for name, value in CHINDEMI_PARAMS.items())
        
        sbatch_script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --account={slurm_account}
#SBATCH --time={time_str}
#SBATCH --mem={mem_gb}G
#SBATCH --cpus-per-task={cpus_per_task}
#SBATCH --output={log_folder}/{job_name}_%j.out
#SBATCH --error={log_folder}/{job_name}_%j.err

source {PLASTYFIRE_DIR}/setupenv.sh

for workdir in {workdirs_str}; do
    echo "Processing: $workdir"
    
    pair_name=$(basename $(dirname $workdir))
    freq_dt=$(basename $workdir)
    output_file="{output_dir}/$pair_name/$freq_dt/simulation_traces.pkl"
    
    if [ -f "$output_file" ]; then
        echo "Skipping $pair_name/$freq_dt - already exists"
        continue
    fi
    
    mkdir -p "{output_dir}/$pair_name/$freq_dt"
    
    cd $workdir
    {sys.executable} {PLASTYFIRE_DIR}/plastyfire/pairrunner.py {param_args} --output-filename=simulation_traces_temp.pkl
    
    if [ -f "$workdir/simulation_traces_temp.pkl" ]; then
        mv "$workdir/simulation_traces_temp.pkl" "$output_file"
        echo "Completed: $pair_name/$freq_dt"
    else
        echo "Failed: $pair_name/$freq_dt"
    fi
done

echo "Batch {batch_idx} completed"
"""
        
        result = subprocess.run(["sbatch"], input=sbatch_script, capture_output=True, text=True)
        
        if result.returncode == 0:
            job_id = result.stdout.strip().split()[-1]
            job_ids.append(job_id)
            logger.info(f"Submitted batch {batch_idx}: job {job_id} ({len(batch)} sims)")
        else:
            logger.error(f"Failed to submit batch {batch_idx}: {result.stderr}")
    
    logger.info(f"\n{'='*60}")
    logger.info(f"Submitted {len(job_ids)} SLURM jobs")
    logger.info(f"Logs: {log_folder}/")
    logger.info(f"\nMonitor with: squeue -u $USER")
    
    return job_ids


def main():
    parser = argparse.ArgumentParser(description="Submit trace collection for missing L5TTPC pairs")
    parser.add_argument("--execution-mode", choices=["slurm", "cpu"], default="slurm")
    parser.add_argument("--workers", type=int, default=30, help="Workers for CPU mode")
    parser.add_argument("--timeout-hours", type=int, default=12)
    parser.add_argument("--mem-gb", type=int, default=8)
    parser.add_argument("--batch-size", type=int, default=10)
    args = parser.parse_args()
    
    output_dir = os.path.join(OUTPUT_BASE_DIR, "Chindemi_params")
    os.makedirs(output_dir, exist_ok=True)
    
    simulations_path = Path(SIMULATIONS_DIR)
    
    # Build job list for missing pairs only
    jobs = []
    for pair_str in MISSING_PAIRS:
        pair_dir = simulations_path / pair_str
        if not pair_dir.exists():
            logger.warning(f"Pair directory not found: {pair_dir}")
            continue
        
        freq_dt_dirs = [d.name for d in pair_dir.iterdir() if d.is_dir()]
        for freq_dt in freq_dt_dirs:
            jobs.append((str(pair_dir), freq_dt))
    
    logger.info(f"Found {len(jobs)} simulations for {len(MISSING_PAIRS)} missing pairs")
    
    if args.execution_mode == "slurm":
        submit_slurm_jobs(jobs, output_dir, 
                          timeout_hours=args.timeout_hours,
                          mem_gb=args.mem_gb,
                          batch_size=args.batch_size)
    else:
        run_cpu_mode(jobs, output_dir, workers=args.workers)


if __name__ == "__main__":
    main()
