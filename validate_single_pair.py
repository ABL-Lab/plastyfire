#!/usr/bin/env python3
import os
import subprocess
import argparse
import pandas as pd
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Validate Basis Run for Single Pair")
    parser.add_argument("--csv", required=True, help="Path to index CSV")
    # parser.add_argument("--sim-config", required=True, help="Path to simulation_config.json") # Removed per user request
    parser.add_argument("--index", type=int, default=0, help="Index of pair in CSV to test")
    parser.add_argument("--output-dir", default="validation_results", help="Dir to save results")
    args = parser.parse_args()
    
    # Read CSV to get a pair
    df = pd.read_csv(args.csv)
    # We need path as well now
    try:
        pair_data = df.iloc[args.index]
        pre = int(pair_data['pregid'])
        post = int(pair_data['postgid'])
        sim_path = pair_data['path']
        print(f"Testing pair {args.index}: {pre} -> {post}")
        print(f"Simulation path from CSV: {sim_path}")
    except IndexError:
        print(f"Index {args.index} out of range. Found {len(df)} pairs.")
        return
        
    # Resolve simulation_config.json from the CSV path
    # CSV path points to simulation.batch in the specific sim dir
    sim_dir = Path(sim_path).parent
    sim_config_candidate = sim_dir / "simulation_config.json"
    
    if not sim_config_candidate.exists():
        print(f"Warning: Config not found at {sim_config_candidate}")
        # Try to find it in the directory or parents?
        # Revert to the heuristic used in run_basis_pair if needed, 
        # But for now, let's assume the CSV path is reliable for the directory structure
        # User showed path: .../10Hz_-10ms/simulation_config.json
        # CSV path: .../40Hz_5ms/simulation.batch
        # So sim_dir / simulation_config.json SHOULD exist.
        pass
    else:
        print(f"Resolved config: {sim_config_candidate}")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    out_csv = output_dir / f"basis_{pre}_{post}.csv"
    
    worker_script = Path("run_basis_pair.py").absolute()
    if not worker_script.exists():
        print(f"Error: {worker_script} not found.")
        return

    cmd = [
        "python", str(worker_script),
        "--pre-gid", str(pre),
        "--post-gid", str(post),
        "--sim-config", str(sim_config_candidate),
        "--output-csv", str(out_csv),
        # Use few trials/workers for validation
        "--num-trials", "2",
        "--workers", "2"
    ]
    
    print("Running command:")
    print(" ".join(cmd))
    
    try:
        subprocess.check_call(cmd)
        print("\nSUCCESS: Validation run completed.")
        if out_csv.exists():
            print(f"Output CSV created at: {out_csv}")
            res_df = pd.read_csv(out_csv)
            print("\nPreview of results:")
            print(res_df)
        else:
            print("WARNING: Output CSV not found (maybe no synapses?). check logs.")
            
    except subprocess.CalledProcessError as e:
        print(f"\nFAILURE: Command exited with error code {e.returncode}")
        
if __name__ == "__main__":
    main()
