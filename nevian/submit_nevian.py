#!/usr/bin/env python3
"""
Submit Nevian & Sakmann 2006 protocol simulations to SLURM.
===========================================================
Submits pairrunner.py jobs for all simulation workdirs created by
setup_nevian_simulations.py.

Usage:
    # First create the simulation folders (if not done yet)
    python setup_nevian_simulations.py --protocols control

    # Submit all control protocols for all pairs
    python submit_nevian.py

    # Specific protocols / pairs
    python submit_nevian.py --protocols LTP_3ap_50hz_dt+10ms,LTD_3ap_50hz_dt-10ms
    python submit_nevian.py --max-pairs 3 --protocols control

    # Run locally (no SLURM)
    python submit_nevian.py --execution-mode cpu --workers 8
"""

import argparse
import logging
import os
import pickle
import subprocess
import sys
from pathlib import Path

# Add plastyfire to path
PLASTYFIRE_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PLASTYFIRE_ROOT))

from nevian_protocols import ALL_PROTOCOLS, CONTROL_PROTOCOLS, PROTOCOLS_BY_ID

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# --- Directories ---
PLASTYFIRE_DIR  = str(PLASTYFIRE_ROOT)
NEVIAN_DIR      = str(Path(__file__).parent)
DEFAULT_SIMS_DIR = str(Path(__file__).parent / "simulations")
DEFAULT_OUT_DIR  = str(Path(__file__).parent / "results")
PAIRRUNNER      = str(PLASTYFIRE_ROOT / "plastyfire" / "pairrunner.py")

# =============================================================================
# Parameter sets  (copy the one you are currently fitting / validating against)
# =============================================================================

# V7 parameters (latest CICR model)
PARAMS_V7 = {
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

# No-CICR baseline (original Chindemi-style)
PARAMS_NO_CICR = {
    "tau_effca_GB_GluSynapse": 278.318,
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

PARAM_SETS = {
    "v7": PARAMS_V7,
    "no_cicr": PARAMS_NO_CICR,
}

DEFAULT_PARAMS = "v7"

# =============================================================================

def build_param_args_str(params: dict) -> str:
    """Build --key=value string for pairrunner.py"""
    return " ".join(f"--{k}={v}" for k, v in params.items())


def build_param_args_list(params: dict) -> list:
    """Build [--key=value, ...] list for subprocess"""
    return [f"--{k}={v}" for k, v in params.items()]


def discover_workdirs(sims_dir: str,
                      protocol_filter: list = None,
                      max_pairs: int = None) -> list:
    """
    Discover all (pair_dir, protocol_id) workdir paths under sims_dir.
    Returns list of (pair_dir_path, protocol_id) tuples.
    """
    sims_path = Path(sims_dir)
    if not sims_path.is_dir():
        logger.error(f"Simulations directory not found: {sims_dir}")
        return []

    pair_dirs = sorted([
        d for d in sims_path.iterdir()
        if d.is_dir() and "-" in d.name and d.name != "single_cells"
    ])
    if max_pairs:
        pair_dirs = pair_dirs[:max_pairs]

    jobs = []
    for pair_dir in pair_dirs:
        proto_dirs = sorted([d for d in pair_dir.iterdir() if d.is_dir()])
        for proto_dir in proto_dirs:
            pid = proto_dir.name
            if protocol_filter and pid not in protocol_filter:
                continue
            # Sanity check: workdir must have simulation_config.json
            if not (proto_dir / "simulation_config.json").exists():
                logger.warning(f"Missing simulation_config.json in {proto_dir}, skipping.")
                continue
            jobs.append((str(pair_dir), pid))

    return jobs


def submit_slurm_jobs(jobs: list,
                      output_dir: str,
                      params: dict,
                      param_label: str,
                      timeout_hours: int = 2,
                      mem_gb: int = 8,
                      cpus_per_task: int = 10,
                      slurm_account: str = "ctb-emuller",
                      batch_size: int = 5,
                      parallel_sims: int = 2) -> list:
    """Submit SLURM batch jobs for the given list of (pair_dir, protocol_id) jobs."""
    log_folder = Path(output_dir) / "logs"
    log_folder.mkdir(parents=True, exist_ok=True)

    param_args_str = build_param_args_str(params)
    time_str = f"{timeout_hours:02d}:00:00"

    # Build batches
    batches = [jobs[i:i+batch_size] for i in range(0, len(jobs), batch_size)]
    logger.info(f"Submitting {len(batches)} SLURM batch jobs for {len(jobs)} simulations...")

    job_ids = []
    for batch_idx, batch in enumerate(batches):
        job_name = f"nevian_{param_label}_{batch_idx:04d}"

        workdirs_lines = "\n".join(
            f'  run_sim "{pair_dir}/{protocol_id}"' for pair_dir, protocol_id in batch
        )

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
    local pair_name=$(basename $(dirname "$workdir"))
    local protocol_id=$(basename "$workdir")
    local output_file="{output_dir}/$pair_name/$protocol_id/simulation_traces.pkl"

    if [ -f "$output_file" ]; then
        echo "[SKIP] $pair_name/$protocol_id already done"
        return 0
    fi

    mkdir -p "{output_dir}/$pair_name/$protocol_id"
    cd "$workdir"
    {sys.executable} {PAIRRUNNER} {param_args_str} --output-filename=simulation_traces.pkl

    if [ -f "$workdir/simulation_traces.pkl" ]; then
        mv "$workdir/simulation_traces.pkl" "$output_file"
        echo "[OK]   $pair_name/$protocol_id"
    else
        echo "[FAIL] $pair_name/$protocol_id"
    fi
}}

# Run batch (up to {parallel_sims} in parallel)
{workdirs_lines.replace('run_sim', 'run_sim') }

wait
echo "Batch {batch_idx} done"
"""
        # Replace the lines to add & and wait logic
        sim_lines = []
        for pair_dir, protocol_id in batch:
            sim_lines.append(f'  run_sim "{pair_dir}/{protocol_id}" &')
            sim_lines.append(f'  while [ $(jobs -r -p | wc -l) -ge {parallel_sims} ]; do wait -n; done')

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
    local pair_name=$(basename $(dirname "$workdir"))
    local protocol_id=$(basename "$workdir")
    local output_file="{output_dir}/$pair_name/$protocol_id/simulation_traces.pkl"

    if [ -f "$output_file" ]; then
        echo "[SKIP] $pair_name/$protocol_id"
        return 0
    fi

    mkdir -p "{output_dir}/$pair_name/$protocol_id"
    cd "$workdir"
    {sys.executable} {PAIRRUNNER} {param_args_str} --output-filename=simulation_traces.pkl

    if [ -f "$workdir/simulation_traces.pkl" ]; then
        mv "$workdir/simulation_traces.pkl" "$output_file"
        echo "[OK]   $pair_name/$protocol_id"
    else
        echo "[FAIL] $pair_name/$protocol_id"
    fi
}}

{chr(10).join(sim_lines)}

wait
echo "Batch {batch_idx} completed ({len(batch)} sims)"
"""
        result = subprocess.run(["sbatch"], input=sbatch_script, capture_output=True, text=True)
        if result.returncode == 0:
            jid = result.stdout.strip().split()[-1]
            job_ids.append(jid)
            logger.info(f"  Submitted batch {batch_idx:4d}: job {jid} ({len(batch)} sims)")
        else:
            logger.error(f"  Failed batch {batch_idx}: {result.stderr.strip()}")

    logger.info(f"\nSubmitted {len(job_ids)}/{len(batches)} jobs")
    logger.info(f"Logs   : {log_folder}/")
    logger.info(f"Results: {output_dir}/")
    logger.info(f"Monitor: squeue -u $USER")
    return job_ids


def _cpu_worker(args):
    """Multiprocessing worker for local CPU execution."""
    pair_dir, protocol_id, output_dir, params = args
    workdir = os.path.join(pair_dir, protocol_id)
    pair_name = os.path.basename(pair_dir)
    output_file = os.path.join(output_dir, pair_name, protocol_id, "simulation_traces.pkl")

    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    param_args = build_param_args_list(params)
    param_args.append("--output-filename=simulation_traces_temp.pkl")
    cmd = [sys.executable, PAIRRUNNER] + param_args

    try:
        result = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True, timeout=7200)
        temp = os.path.join(workdir, "simulation_traces_temp.pkl")
        if os.path.exists(temp):
            import shutil
            shutil.move(temp, output_file)
            return {"status": "success", "pair": pair_name, "protocol": protocol_id}
        return {"status": "failed", "pair": pair_name, "protocol": protocol_id,
                "returncode": result.returncode,
                "stderr": result.stderr[:500]}
    except Exception as e:
        return {"status": "error", "pair": pair_name, "protocol": protocol_id, "error": str(e)}


def run_cpu_mode(jobs: list, output_dir: str, params: dict, workers: int = 8):
    """Run simulations locally using multiprocessing."""
    from multiprocessing import Pool
    import time

    jobs_with_args = [(pair_dir, protocol_id, output_dir, params)
                      for pair_dir, protocol_id in jobs]
    logger.info(f"Running {len(jobs)} simulations with {workers} workers...")
    start = time.time()
    with Pool(workers) as pool:
        results = pool.map(_cpu_worker, jobs_with_args)
    elapsed = time.time() - start

    success = sum(1 for r in results if r["status"] == "success")
    failed  = sum(1 for r in results if r["status"] in ["failed", "error"])
    logger.info(f"Done in {elapsed/60:.1f} min  |  Success: {success}  Failed: {failed}")

    for r in results:
        if r["status"] != "success":
            logger.warning(f"  [{r['status'].upper()}] {r['pair']}/{r.get('protocol', '?')}: "
                           f"{r.get('stderr', r.get('error', ''))}")


def main():
    parser = argparse.ArgumentParser(
        description="Submit Nevian & Sakmann 2006 STDP simulations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--sims-dir", default=DEFAULT_SIMS_DIR,
                        help=f"Root of simulation workdirs (default: {DEFAULT_SIMS_DIR})")
    parser.add_argument("--output-dir", default=None,
                        help="Output directory for results (default: nevian/results/<param_label>/)")
    parser.add_argument("--params", choices=list(PARAM_SETS.keys()), default=DEFAULT_PARAMS,
                        help=f"Parameter set to use (default: {DEFAULT_PARAMS})")
    parser.add_argument("--protocols", default="control",
                        help="'all', 'control', or comma-separated protocol IDs")
    parser.add_argument("--max-pairs", type=int, default=None,
                        help="Limit to N pairs (for testing)")
    parser.add_argument("--execution-mode", choices=["slurm", "cpu"], default="slurm")
    # SLURM options
    parser.add_argument("--timeout-hours", type=int, default=2)
    parser.add_argument("--mem-gb", type=int, default=8)
    parser.add_argument("--cpus-per-task", type=int, default=10)
    parser.add_argument("--batch-size", type=int, default=5,
                        help="Simulations per SLURM job (default: 5)")
    parser.add_argument("--parallel-sims", type=int, default=2,
                        help="Parallel sims inside one SLURM job (default: 2)")
    parser.add_argument("--slurm-account", default="ctb-emuller")
    # CPU options
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    # Protocol filter
    if args.protocols == "all":
        protocol_filter = None  # no filter = all
    elif args.protocols == "control":
        protocol_filter = [p.protocol_id for p in CONTROL_PROTOCOLS]
    else:
        protocol_filter = [p.strip() for p in args.protocols.split(",")]

    # Params
    params = PARAM_SETS[args.params]
    param_label = args.params

    # Output dir
    output_dir = args.output_dir or os.path.join(NEVIAN_DIR, "results", param_label)
    os.makedirs(output_dir, exist_ok=True)

    # Discover workdirs
    jobs = discover_workdirs(args.sims_dir, protocol_filter, args.max_pairs)
    if not jobs:
        logger.error("No simulation workdirs found. Run setup_nevian_simulations.py first.")
        sys.exit(1)

    logger.info(f"Found {len(jobs)} simulations to run")
    logger.info(f"Params : {param_label}")
    logger.info(f"Output : {output_dir}")

    # Save metadata
    meta = {"params": params, "param_label": param_label, "jobs": jobs, "output_dir": output_dir}
    with open(os.path.join(output_dir, "job_metadata.pkl"), "wb") as f:
        pickle.dump(meta, f)

    if args.execution_mode == "slurm":
        submit_slurm_jobs(
            jobs, output_dir, params, param_label,
            timeout_hours=args.timeout_hours,
            mem_gb=args.mem_gb,
            cpus_per_task=args.cpus_per_task,
            slurm_account=args.slurm_account,
            batch_size=args.batch_size,
            parallel_sims=args.parallel_sims,
        )
    else:
        run_cpu_mode(jobs, output_dir, params, workers=args.workers)


if __name__ == "__main__":
    main()
