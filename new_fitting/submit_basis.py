#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import submitit
from pathlib import Path

def get_pairs_and_paths(csv_path):
    df = pd.read_csv(csv_path)
    # We need unique pairs, but also their paths.
    # Assuming one path per pair? Or picking the first one?
    # The CSV might have multiple entries for different frequencies/dts.
    # We just need A valid simulation directory for the pair.
    # So we groupby pair and take the first 'path'.
    grouped = df.groupby(['pregid', 'postgid'])['path'].first().reset_index()
    return grouped.values.tolist()

def main():
    parser = argparse.ArgumentParser(description="Submit Basis Simulations")
    parser.add_argument("--csv", required=True, help="Path to index CSV")
    # parser.add_argument("--sim-config", required=True, help="Path to simulation_config.json") # Removed
    parser.add_argument("--output-dir", default="basis_results", help="Dir to save results")
    parser.add_argument("--partition", default=None, help="SLURM partition")
    parser.add_argument("--account", default="ctb-emuller", help="SLURM account")
    parser.add_argument("--timeout", default=60, type=int, help="Timeout in minutes")
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Executor setup
    log_folder = output_dir / "logs"
    executor = submitit.AutoExecutor(folder=log_folder)
    
    submitit_params = {
        "timeout_min": args.timeout,
        "slurm_account": args.account,
        "nodes": 1,
        "cpus_per_task": 20, 
        "mem_gb": 64,
        "slurm_array_parallelism": 20
    }
    
    if args.partition:
        submitit_params["slurm_partition"] = args.partition
        
    executor.update_parameters(**submitit_params)
    
    # Get pairs and their simulation paths
    pairs_data = get_pairs_and_paths(args.csv)
    print(f"Found {len(pairs_data)} pairs to submit.")
    
    worker_script = Path("run_basis_pair.py").absolute()
    
    jobs = []
    with executor.batch():
        for pre, post, sim_path in pairs_data:
            pre = int(pre)
            post = int(post)
            
            # Resolve config path
            # sim_path points to simulation.batch file
            sim_dir = Path(sim_path).parent
            sim_config = sim_dir / "simulation_config.json"
            
            # Wrapper function to run the command
            def run_cmd(cmd_list):
                import subprocess
                subprocess.check_call(cmd_list)

            # Construct the command
            out_csv = str(output_dir / f"basis_{pre}_{post}.csv")
            cmd = [
                "python", str(worker_script),
                "--pre-gid", str(pre),
                "--post-gid", str(post),
                "--sim-config", str(sim_config),
                "--output-csv", out_csv,
                "--workers", "18"
            ]
            
            job = executor.submit(run_cmd, cmd)
            jobs.append(job)
    
    for i, job in enumerate(jobs):
        print(f"Submitted job {job.job_id} for pair {i}")
            
    print(f"Submitted {len(jobs)} jobs. Logs in {log_folder}")

if __name__ == "__main__":
    main()
