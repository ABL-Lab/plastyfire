#!/usr/bin/env python3
"""
Submit STDP trace simulations to SLURM or run locally.

Usage:
    # Submit all to SLURM (overnight run)
    python submit_stdp_traces.py
    
    # Test with 2 pairs
    python submit_stdp_traces.py --max-pairs 2
    
    # Run locally
    python submit_stdp_traces.py --execution-mode cpu --workers 30
"""

import argparse
import logging
import os
import pickle
import subprocess
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

SIMULATIONS_BASE_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/refitting_results/fitting/n100/seed19091997"
OUTPUT_BASE_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/trace_results"
PLASTYFIRE_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire"

DHURUVA_PARAMS = {
    "tau_effca_GB_GluSynapse": 111.4966,
    "gamma_d_GB_GluSynapse": 88.70172982129458,
    "gamma_p_GB_GluSynapse": 268.32068636660046,
    "a00": 0.4777981934076907,
    "a01": 0.6934394338463974,
    "a10": 1.0586973073436239,
    "a11": 3.465063729118567,
    "a20": 0.7322653809967383,
    "a21": 0.14746008318800907,
    "a30": 1.88136100637871,
    "a31": 4.447189652961518,
}

DHURUVA_PARAMS_V2 = {
    "tau_effca_GB_GluSynapse": 73.21360951908201,
    "gamma_d_GB_GluSynapse": 109.4093089925484,
    "gamma_p_GB_GluSynapse": 596.8444303473575,
    "a00": 0.2933119882453923,
    "a01": 1.5969307382598152,
    "a10": 0.8408881959013523,
    "a11": 4.8170123961535785,
    "a20": 0.6607260497664416,
    "a21": 0.010296203791785974,
    "a30": 3.161785481914619,
    "a31": 3.6158549749054814,
    "tau_prime_CICR_GluSynapse": 1109.1151367901289,
    "k_prime_CICR_GluSynapse": 17516.914653528733,
    "K_prime_limit_CICR_GluSynapse": 0.00046795676802132444,
    "K_RyR_CICR_GluSynapse": 0.0020996724134375957,
    "n_RyR_CICR_GluSynapse": 4.457570893183441,
    "Vmax_CICR_GluSynapse": 0.0005044554223042098,
    "tau_rel_CICR_GluSynapse": 223.5880115766854,
}

pair_types = [
    "L5TTPC_L5TTPC_STDP",
    "L23PC_L5TTPC_STDP",
    "L23PC_L23PC_STDP",
    "L4PC_L23PC_STDP"
]

def submit_slurm_jobs(jobs, output_dir, timeout_hours=4, mem_gb=8, cpus_per_task=10,
                       slurm_account="ctb-emuller", batch_size=10, parallel_sims=4):
    log_folder = Path(output_dir) / "logs"
    log_folder.mkdir(parents=True, exist_ok=True)
    
    batches = [jobs[i:i + batch_size] for i in range(0, len(jobs), batch_size)]
    logger.info(f"Submitting {len(batches)} SLURM jobs for {len(jobs)} simulations...")
    
    time_str = f"{timeout_hours:02d}:00:00"
    job_ids = []
    
    for batch_idx, batch in enumerate(batches):
        job_name = f"stdp_trace_{batch_idx:04d}"
        
        workdirs = [os.path.join(t_dir, pair_dir, freq_dt) for t_dir, pair_dir, freq_dt in batch]
        workdirs_str = " ".join(f'"{w}"' for w in workdirs)
        
        param_args = " ".join(f"--{name}={value}" for name, value in DHURUVA_PARAMS.items())
        
        sbatch_script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --account={slurm_account}
#SBATCH --time={time_str}
#SBATCH --mem={mem_gb}G
#SBATCH --cpus-per-task={cpus_per_task}
#SBATCH --output={log_folder}/{job_name}_%j.out
#SBATCH --error={log_folder}/{job_name}_%j.err

source {PLASTYFIRE_DIR}/setupenv.sh

run_sim() {{
    local workdir=$1
    echo "Processing: $workdir"
    
    local freq_dt=$(basename "$workdir")
    local pair_name=$(basename $(dirname "$workdir"))
    local type_name=$(basename $(dirname $(dirname $(dirname "$workdir"))))
    
    local output_file="{output_dir}/$type_name/$pair_name/$freq_dt/simulation_traces.pkl"
    
    if [ -f "$output_file" ]; then
        echo "Skipping $type_name/$pair_name/$freq_dt - already exists"
        return 0
    fi
    
    mkdir -p "{output_dir}/$type_name/$pair_name/$freq_dt"
    
    cd "$workdir"
    {sys.executable} {PLASTYFIRE_DIR}/plastyfire/pairrunner.py {param_args} --output-filename=simulation_traces_temp.pkl
    
    if [ -f "$workdir/simulation_traces_temp.pkl" ]; then
        mv "$workdir/simulation_traces_temp.pkl" "$output_file"
        echo "Completed: $type_name/$pair_name/$freq_dt"
    else
        echo "Failed: $type_name/$pair_name/$freq_dt"
    fi
}}

for workdir in {workdirs_str}; do
    run_sim "$workdir" &
    while [ $(jobs -r -p | wc -l) -ge {parallel_sims} ]; do
        wait -n
    done
done
wait
echo "Batch {batch_idx} completed"
"""
        result = subprocess.run(["sbatch"], input=sbatch_script, capture_output=True, text=True)
        if result.returncode == 0:
            job_id = result.stdout.strip().split()[-1]
            job_ids.append(job_id)
            logger.info(f"Submitted batch {batch_idx}: job {job_id} ({len(batch)} sims)")
        else:
            logger.error(f"Failed to submit batch {batch_idx}: {result.stderr}")
    return job_ids

def _cpu_worker(args):
    type_dir, pair_dir, freq_dt, output_dir = args
    workdir = os.path.join(type_dir, pair_dir, freq_dt)
    pair_name = os.path.basename(pair_dir)
    type_name = os.path.basename(os.path.dirname(type_dir))
    
    output_file = os.path.join(output_dir, type_name, pair_name, freq_dt, "simulation_traces.pkl")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    param_args = [f"--{name}={value}" for name, value in DHURUVA_PARAMS.items()]
    param_args.append("--output-filename=simulation_traces_temp.pkl")
    cmd = [sys.executable, f"{PLASTYFIRE_DIR}/plastyfire/pairrunner.py"] + param_args
    
    try:
        result = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True, timeout=7200)
        temp_file = os.path.join(workdir, "simulation_traces_temp.pkl")
        if os.path.exists(temp_file):
            import shutil
            shutil.move(temp_file, output_file)
            return {"status": "success", "pair": pair_name, "freq_dt": freq_dt}
        return {"status": "failed", "pair": pair_name, "freq_dt": freq_dt}
    except Exception as e:
        return {"status": "error", "pair": pair_name, "freq_dt": freq_dt, "error": str(e)}

def run_cpu_mode(jobs, output_dir, workers=30):
    from multiprocessing import Pool
    import time
    logger.info(f"Running {len(jobs)} simulations with {workers} workers...")
    jobs_with_output = [(t, p, f, output_dir) for t, p, f in jobs]
    
    start = time.time()
    with Pool(workers) as pool:
        results = pool.map(_cpu_worker, jobs_with_output)
    elapsed = time.time() - start
    
    success = sum(1 for r in results if r["status"] == "success")
    failed = sum(1 for r in results if r["status"] in ["failed", "error"])
    logger.info(f"Completed in {elapsed/60:.1f} minutes. Success: {success}, Failed: {failed}")

def main():
    parser = argparse.ArgumentParser(description="Submit STDP trace simulations")
    parser.add_argument("--params", choices=["v1", "v2"], default="v1")
    parser.add_argument("--execution-mode", choices=["slurm", "cpu"], default="slurm")
    parser.add_argument("--max-pairs", type=int, default=None)
    parser.add_argument("--workers", type=int, default=30)
    args = parser.parse_args()
    
    global DHURUVA_PARAMS
    output_subdir = "DHURUVA_PARAMS_STDP"
    if args.params == "v2":
        DHURUVA_PARAMS = DHURUVA_PARAMS_V2
        output_subdir = "DHURUVA_PARAMS_STDP_V2"
        
    output_dir = os.path.join(OUTPUT_BASE_DIR, output_subdir)
    os.makedirs(output_dir, exist_ok=True)
    
    jobs = []
    
    for t_type in pair_types:
        sims_path = Path(SIMULATIONS_BASE_DIR) / t_type / "simulations"
        if not sims_path.exists():
            continue
            
        pair_dirs = sorted([d for d in sims_path.iterdir() if d.is_dir() and d.name != "single_cells"])
        if args.max_pairs:
            pair_dirs = pair_dirs[:args.max_pairs]
            
        for pair_dir in pair_dirs:
            for freq_dt_dir in [d.name for d in pair_dir.iterdir() if d.is_dir()]:
                jobs.append((str(sims_path), pair_dir.name, freq_dt_dir))
                
    logger.info(f"Total simulations: {len(jobs)}")
    
    if args.execution_mode == "slurm":
        submit_slurm_jobs(jobs, output_dir)
    else:
        run_cpu_mode(jobs, output_dir, workers=args.workers)

if __name__ == "__main__":
    main()
