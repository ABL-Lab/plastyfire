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

# Base directories
SIMULATIONS_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC/simulations"
SIMULATIONS_DIR_STDP = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations"
OUTPUT_BASE_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/trace_results"
PLASTYFIRE_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire"

# Chindemi parameters (from evaluate_best_solution.py)
# CHINDEMI_PARAMS = {
#     "gamma_d_GB_GluSynapse": 101.5,
#     "gamma_p_GB_GluSynapse": 216.2,
#     "a00": 1.002,
#     "a01": 1.954,
#     "a10": 1.159,
#     "a11": 2.483,
#     "a20": 1.127,
#     "a21": 2.456,
#     "a30": 5.236,
#     "a31": 1.782,
# }

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

DHURUVA_PARAMS_V3 = {
    "gamma_d_GB_GluSynapse": 50.4690,
    "gamma_p_GB_GluSynapse": 161.2494,
    "a00": 0.5169,
    "a01": 4.6079,
    "a10": 1.8423,
    "a11": 4.9230,
    "a20": 0.1282,
    "a21": 0.4070,
    "a30": 0.3265,
    "a31": 4.9275,
    "tau_prime_CICR_GluSynapse": 436.0976,
    "k_prime_CICR_GluSynapse": 44310.3178,
    "K_prime_limit_CICR_GluSynapse": 0.0003,
    "K_RyR_CICR_GluSynapse": 0.0013,
    "n_RyR_CICR_GluSynapse": 5.0744,
    "Vmax_CICR_GluSynapse": 0.0022,
    "tau_rel_CICR_GluSynapse": 313.2375,
    "tau_effca_GB_GluSynapse": 100.7488,
}

DHURUVA_PARAMS_V4 = {
    "gamma_d_GB_GluSynapse": 171.5880,
    "gamma_p_GB_GluSynapse": 598.7721,
    "a00": 2.3404,
    "a01": 2.2347,
    "a10": 2.7113,
    "a11": 4.0093,
    "a20": 2.6263,
    "a21": 2.6543,
    "a30": 4.4278,
    "a31": 3.9129,
    "tau_prime_CICR_GluSynapse": 1989.2330,
    "k_prime_CICR_GluSynapse": 41683.2555,
    "K_prime_limit_CICR_GluSynapse": 0.0010,
    "K_RyR_CICR_GluSynapse": 0.0049,
    "n_RyR_CICR_GluSynapse": 4.5044,
    "Vmax_CICR_GluSynapse": 0.0019,
    "tau_rel_CICR_GluSynapse": 450.8465,
    "tau_effca_GB_GluSynapse": 344.6220,
}

DHURUVA_PARAMS_V5 = {
    "gamma_d_GB_GluSynapse": 69.35582333947616,
    "gamma_p_GB_GluSynapse": 172.8162926754338,
    "a00": 1.0338281805504534,
    "a01": 2.0824027912480663,
    "a10": 1.50981772882111,
    "a11": 4.937595300732932,
    "a20": 1.024923299559302,
    "a21": 1.0040979194963504,
    "a30": 4.869465155545409,
    "a31": 4.37081939321435,
    "tau_prime_CICR_GluSynapse": 870.7718136117252,
    "k_prime_CICR_GluSynapse": 17653.30973320602,
    "K_prime_limit_CICR_GluSynapse": 0.0001738667305935812,
    "K_RyR_CICR_GluSynapse": 0.002333447308131663,
    "n_RyR_CICR_GluSynapse": 4.637471212566512,
    "Vmax_CICR_GluSynapse": 0.0009454229634058555,
    "tau_rel_CICR_GluSynapse": 71.00509272517166,
    "tau_effca_GB_GluSynapse": 279.72967330670116,
}

DHURUVA_PARAMS_V6 = {
    "gamma_d_GB_GluSynapse": 88.37748728897184,
    "gamma_p_GB_GluSynapse": 75.21696897014064,
    "a00": 1.0016365565804994,
    "a01": 1.4995212257656416,
    "a10": 1.1327719903735114,
    "a11": 3.077770462888927,
    "a20": 1.370297324550421,
    "a21": 2.4040852180572143,
    "a30": 1.9262572961622118,
    "a31": 3.508641458270441,
    # CICR params — confirmed HOC global names (all end in _GluSynapse in compiled .so)
    "k_fill_CICR_GluSynapse": 1.6465937576950695,
    "tau_leak_CICR_GluSynapse": 96157.34999520841,
    "K_clip_CICR_GluSynapse": 0.0008326648885177493,
    "K_trig_CICR_GluSynapse": 0.0019013479472819725,
    "n_trig_CICR_GluSynapse": 5.737730388740266,
    "Vmax_CICR_GluSynapse": 1.6609943109456529,
    "tau_CICR_rel_GluSynapse": 807.4956732356754,
    "g_cicr_GB_GluSynapse": 47234.51951479494,
    "tau_effca_GB_GluSynapse": 484.9361541639738,
}

DHURUVA_PARAMS_V7 = {
    "enable_CICR_GluSynapse": 1,
    "gamma_d_GB_GluSynapse": 69.61242611028226,
    "gamma_p_GB_GluSynapse": 166.96578019613747,
    "a00": 1.0003490121729528,
    "a01": 1.0877438002708641,
    "a10": 1.3917207001380245,
    "a11": 4.129505808304705,
    "a20": 2.3079274077684175,
    "a21": 1.9619384872398133,
    "a30": 4.456560780131431,
    "a31": 4.80485148600644,
    "delta_IP3_CICR_GluSynapse": 2.3564846028724626,
    "tau_IP3_CICR_GluSynapse": 4198.8339980611045,
    "phi_serca_CICR_GluSynapse": 0.00010926586342370828,
    "k_leak_CICR_GluSynapse": 0.00038983484712417705,
    "V_RyR_CICR_GluSynapse": 0.0008831061403601626,
    "K_T_CICR_GluSynapse": 0.0016997377463884833,
    "n_T_CICR_GluSynapse": 9.677747356489755,
    "V_IP3R_CICR_GluSynapse": 0.009883023779226525,
    "g_RyR_CICR_GluSynapse": 14060.428640146583,
    "tau_effca_GB_GluSynapse": 306.6936809019515,
}

def submit_slurm_jobs(jobs, output_dir, timeout_hours=12, mem_gb=8, cpus_per_task=10,
                       slurm_account="ctb-emuller", batch_size=10, parallel_sims=4):
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
    time_str = f"00:30:00"
    
    job_ids = []
    
    for batch_idx, batch in enumerate(batches):
        job_name = f"trace_{batch_idx:04d}"
        
        # Build the list of workdirs for this batch
        workdirs = [os.path.join(pair_dir, freq_dt) for pair_dir, freq_dt in batch]
        workdirs_str = " ".join(f'"{w}"' for w in workdirs)
        
        # Build parameter args string
        param_args = " ".join(f"--{name}={value}" for name, value in DHURUVA_PARAMS.items())
        
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

run_sim() {{
    local workdir=$1
    echo "Processing: $workdir"
    
    local pair_name=$(basename $(dirname "$workdir"))
    local freq_dt=$(basename "$workdir")
    local output_file="{output_dir}/$pair_name/$freq_dt/simulation_traces.pkl"
    
    # Skip if already exists
    if [ -f "$output_file" ]; then
        echo "Skipping $pair_name/$freq_dt - already exists"
        return 0
    fi
    
    mkdir -p "{output_dir}/$pair_name/$freq_dt"
    
    cd "$workdir"
    {sys.executable} {PLASTYFIRE_DIR}/plastyfire/pairrunner.py {param_args} --output-filename=simulation_traces_temp.pkl
    
    if [ -f "$workdir/simulation_traces_temp.pkl" ]; then
        mv "$workdir/simulation_traces_temp.pkl" "$output_file"
        echo "Completed: $pair_name/$freq_dt"
    else
        echo "Failed: $pair_name/$freq_dt"
    fi
}}

# Run each simulation in this batch
for workdir in {workdirs_str}; do
    run_sim "$workdir" &
    
    # Limit number of parallel jobs
    while [ $(jobs -r -p | wc -l) -ge {parallel_sims} ]; do
        wait -n
    done
done

wait
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
    
    if False:
        return {"status": "skipped", "pair": pair_name}
    
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
    
    parser.add_argument("--params", choices=["v1", "v2", "v3", "v4", "v5", "v6", "v7"], default="v1",
                        help="Parameter set to use: v1=DHURUVA_PARAMS, ..., v7=DHURUVA_PARAMS_V7 (default: v1)")
    parser.add_argument("--execution-mode", choices=["slurm", "cpu"], default="slurm",
                        help="Execution mode (default: slurm)")
    parser.add_argument("--max-pairs", type=int, default=None,
                        help="Max pairs to process (for testing)")
    parser.add_argument("--freq", type=str, default=None,
                        help="Filter by frequency, e.g., '10Hz'")
    parser.add_argument("--delays", type=str, nargs="+", default=None,
                        help="Filter by delays, e.g., '0' '-30' '10'")
    parser.add_argument("--output-dir", type=str, default=None,
                        help="Custom output directory")
    parser.add_argument("--simulations-dir", type=str, default=None,
                        help="Custom simulations directory (default depends on --params)")
    
    # CPU mode options
    parser.add_argument("--workers", type=int, default=30,
                        help="Workers for CPU mode (default: 30)")
    
    # SLURM options
    parser.add_argument("--timeout-hours", type=int, default=2,
                        help="SLURM timeout in hours (default: 2)")
    parser.add_argument("--mem-gb", type=int, default=4,
                        help="Memory per job in GB (default: 8)")
    parser.add_argument("--cpus-per-task", type=int, default=10,
                        help="CPUs per SLURM job (default: 10)")
    parser.add_argument("--batch-size", type=int, default=10,
                        help="Simulations per SLURM job (default: 10)")
    parser.add_argument("--parallel-sims", type=int, default=4,
                        help="Parallel simulations to pack inside one SLURM job (default: 4)")
    parser.add_argument("--slurm-account", default="ctb-emuller",
                        help="SLURM account (default: ctb-emuller)")
    
    args = parser.parse_args()
    
    # Select parameter set
    global DHURUVA_PARAMS
    default_simulations_dir = SIMULATIONS_DIR
    if args.params == "v7":
        DHURUVA_PARAMS = DHURUVA_PARAMS_V7
        default_subdir = "DHURUVA_PARAMS_V7"
        default_simulations_dir = SIMULATIONS_DIR_STDP
    elif args.params == "v6":
        DHURUVA_PARAMS = DHURUVA_PARAMS_V6
        default_subdir = "DHURUVA_PARAMS_V6"
        default_simulations_dir = SIMULATIONS_DIR_STDP
    elif args.params == "v5":
        DHURUVA_PARAMS = DHURUVA_PARAMS_V5
        default_subdir = "DHURUVA_PARAMS_V5"
    elif args.params == "v4":
        DHURUVA_PARAMS = DHURUVA_PARAMS_V4
        default_subdir = "DHURUVA_PARAMS_V4"
    elif args.params == "v3":
        DHURUVA_PARAMS = DHURUVA_PARAMS_V3
        default_subdir = "DHURUVA_PARAMS_V3"
    elif args.params == "v2":
        DHURUVA_PARAMS = DHURUVA_PARAMS_V2
        default_subdir = "DHURUVA_PARAMS_V2"
    else:
        default_subdir = "DHURUVA_PARAMS"

    # Output directory
    output_dir = args.output_dir or os.path.join(OUTPUT_BASE_DIR, default_subdir)
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Output directory: {output_dir}")
    
    # Find all simulation directories
    simulations_path = Path(args.simulations_dir if args.simulations_dir is not None else default_simulations_dir)
    pair_dirs = sorted([d for d in simulations_path.iterdir() 
                       if d.is_dir() and d.name != "single_cells"])
    
    if args.max_pairs:
        pair_dirs = pair_dirs[:args.max_pairs]
    
    logger.info(f"Found {len(pair_dirs)} pair directories in {simulations_path}")
    
    # Build job list
    jobs = []
    for pair_dir in pair_dirs:
        freq_dt_dirs = [d.name for d in pair_dir.iterdir() if d.is_dir()]
        
        if args.freq:
            freq_dt_dirs = [d for d in freq_dt_dirs if args.freq in d]
            
        if args.delays:
            # Match {delay}ms in the directory name, typically {freq}Hz_{delay}ms
            target_delays = [f"{d}ms" for d in args.delays]
            freq_dt_dirs = [d for d in freq_dt_dirs if any(d.endswith(f"_{delay}") for delay in target_delays)]
        
        for freq_dt in freq_dt_dirs:
            jobs.append((str(pair_dir), freq_dt))
    
    logger.info(f"Total simulations: {len(jobs)}")
    
    # Save job metadata
    metadata = {
        "params": DHURUVA_PARAMS,
        "jobs": jobs,
        "output_dir": output_dir,
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
            slurm_account=args.slurm_account,
            parallel_sims=args.parallel_sims
        )
    else:
        run_cpu_mode(jobs, output_dir, workers=args.workers)


if __name__ == "__main__":
    main()
