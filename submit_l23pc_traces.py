#!/usr/bin/env python3
"""
Submit L23PC_L5TTPC plasticity simulations with Chindemi parameters to SLURM,
recording synaptic traces (effcai_GB, cai_CR, vsyn, rho_GB) during induction.

Usage:
    # Submit all pairs as SLURM jobs (overnight run)
    python submit_l23pc_traces.py

    # Test with 2 pairs first
    python submit_l23pc_traces.py --max-pairs 2

    # Use multiprocessing instead of SLURM
    python submit_l23pc_traces.py --execution-mode cpu --workers 30

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

# Base directories
SIMULATIONS_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/refitting_results/fitting/n100/seed19091997/L23PC_L5TTPC/simulations"
OUTPUT_BASE_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/trace_results"
PLASTYFIRE_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire"

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


def submit_slurm_jobs(jobs, output_dir, timeout_hours=12, mem_gb=8, cpus_per_task=4,
                       slurm_account="ctb-emuller", batch_size=10):
    """Submit jobs to SLURM as batch jobs.

    Args:
        jobs: List of (pair_dir, freq_dt_dir) tuples
        output_dir: Output directory for results
        timeout_hours: Job timeout in hours (default: 12 for overnight)
        mem_gb: Memory per job in GB
        cpus_per_task: CPUs per job
        slurm_account: SLURM account
        batch_size: Number of simulations per SLURM job
    """
    log_folder = Path(output_dir) / "logs"
    log_folder.mkdir(parents=True, exist_ok=True)

    # Create batches
    batches = []
    for i in range(0, len(jobs), batch_size):
        batch = jobs[i:i + batch_size]
        batches.append(batch)

    logger.info(f"Submitting {len(batches)} SLURM jobs for {len(jobs)} simulations...")

    # Time format
    hours = timeout_hours
    time_str = f"{hours:02d}:00:00"

    job_ids = []

    for batch_idx, batch in enumerate(batches):
        job_name = f"l23_trace_{batch_idx:04d}"

        # Build the list of workdirs for this batch
        workdirs = [os.path.join(pair_dir, freq_dt) for pair_dir, freq_dt in batch]
        workdirs_str = " ".join(f'"{w}"' for w in workdirs)

        # Build parameter args string
        param_args = " ".join(f"--{name}={value}" for name, value in CHINDEMI_PARAMS.items())

        # Create sbatch script
        sbatch_script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --account={slurm_account}
#SBATCH --time={time_str}
#SBATCH --mem={mem_gb}G
#SBATCH --cpus-per-task={cpus_per_task}
#SBATCH --output={log_folder}/{job_name}_%j.out
#SBATCH --error={log_folder}/{job_name}_%j.err

source {PLASTYFIRE_DIR}/setupenv.sh

# Run each simulation in this batch
for workdir in {workdirs_str}; do
    echo "Processing: $workdir"

    pair_name=$(basename $(dirname $workdir))
    freq_dt=$(basename $workdir)
    output_file="{output_dir}/$pair_name/$freq_dt/simulation_traces.pkl"

    # Skip if already exists
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

        # Submit via stdin
        result = subprocess.run(
            ["sbatch"],
            input=sbatch_script,
            capture_output=True,
            text=True
        )

        if result.returncode == 0:
            job_id = result.stdout.strip().split()[-1]
            job_ids.append(job_id)
            logger.info(f"Submitted batch {batch_idx}: job {job_id} ({len(batch)} sims)")
        else:
            logger.error(f"Failed to submit batch {batch_idx}: {result.stderr}")

    logger.info(f"\n{'='*60}")
    logger.info(f"Submitted {len(job_ids)} SLURM jobs")
    logger.info(f"Logs: {log_folder}/")
    logger.info(f"Results: {output_dir}/")
    logger.info(f"\nMonitor with: squeue -u $USER")

    return job_ids


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

    cmd = [sys.executable, f"{PLASTYFIRE_DIR}/plastyfire/pairrunner.py"] + param_args

    try:
        result = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True, timeout=7200)

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
                "stderr": result.stderr[:500] if result.stderr else "",
                "stdout": result.stdout[:500] if result.stdout else "",
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



def main():
    parser = argparse.ArgumentParser(
        description="Submit L23PC_L5TTPC trace simulations with Chindemi params",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Submit all to SLURM (overnight run)
  python submit_l23pc_traces.py

  # Test with 2 pairs
  python submit_l23pc_traces.py --max-pairs 2

  # Run locally with multiprocessing
  python submit_l23pc_traces.py --execution-mode cpu --workers 30
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
    parser.add_argument("--timeout-hours", type=int, default=12,
                        help="SLURM timeout in hours (default: 12)")
    parser.add_argument("--mem-gb", type=int, default=8,
                        help="Memory per job in GB (default: 8)")
    parser.add_argument("--cpus-per-task", type=int, default=4,
                        help="CPUs per SLURM job (default: 4)")
    parser.add_argument("--batch-size", type=int, default=10,
                        help="Simulations per SLURM job (default: 10)")
    parser.add_argument("--slurm-account", default="ctb-emuller",
                        help="SLURM account (default: ctb-emuller)")

    args = parser.parse_args()

    # Output directory - use L23PC_Chindemi_params subfolder
    output_dir = args.output_dir or os.path.join(OUTPUT_BASE_DIR, "L23PC_Chindemi_params")
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Output directory: {output_dir}")

    # Find all simulation directories
    simulations_path = Path(SIMULATIONS_DIR)
    pair_dirs = sorted([d for d in simulations_path.iterdir()
                       if d.is_dir() and d.name != "single_cells"])

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
            jobs.append((str(pair_dir), freq_dt))

    logger.info(f"Total simulations: {len(jobs)}")

    # Save job metadata
    metadata = {
        "params": CHINDEMI_PARAMS,
        "jobs": jobs,
        "output_dir": output_dir,
        "pathway": "L23PC_L5TTPC",
    }
    with open(os.path.join(output_dir, "job_metadata.pkl"), "wb") as f:
        pickle.dump(metadata, f)

    # Execute
    if args.execution_mode == "slurm":
        submit_slurm_jobs(
            jobs, output_dir,
            timeout_hours=args.timeout_hours,
            mem_gb=args.mem_gb,
            cpus_per_task=args.cpus_per_task,
            batch_size=args.batch_size,
            slurm_account=args.slurm_account
        )
    else:
        run_cpu_mode(jobs, output_dir, workers=args.workers)


if __name__ == "__main__":
    main()
