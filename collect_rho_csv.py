#!/usr/bin/env python3
"""
Collect per-synapse initial/final rho values from simulation_traces.pkl files
for a given parameter set and write them to a CSV.

Usage:
    python collect_rho_csv.py --params v17
    python collect_rho_csv.py --params v17 --output rho_v17.csv
"""

import os
import sys
import argparse
import pickle
import glob
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from submit_l5ttpc_traces import (
    CHINDEMI_PARAMS, DHURUVA_PARAMS, DHURUVA_PARAMS_V2, DHURUVA_PARAMS_V3, DHURUVA_PARAMS_V4,
    DHURUVA_PARAMS_V5, DHURUVA_PARAMS_V6, DHURUVA_PARAMS_V7, DHURUVA_PARAMS_V8,
    DHURUVA_PARAMS_V9, DHURUVA_PARAMS_V10, DHURUVA_PARAMS_V11, DHURUVA_PARAMS_V12,
    DHURUVA_PARAMS_V13, DHURUVA_PARAMS_V14, DHURUVA_PARAMS_V15, DHURUVA_PARAMS_V16,
    DHURUVA_PARAMS_V17, DHURUVA_PARAMS_V18, DHURUVA_PARAMS_V19
)

TRACE_RESULTS_BASE = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/trace_results"
TRACE_RESULTS_DIRS = {
    "v1":  os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS"),
    "v2":  os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V2"),
    "v3":  os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V3"),
    "v4":  os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V4"),
    "v5":  os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_STDP_V5"),
    "v6":  os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V6"),
    "v7":  os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V7"),
    "v8":  os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V8"),
    "v9":  os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V9"),
    "v10": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V10"),
    "v11": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V11"),
    "v12": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V12"),
    "v13": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V13"),
    "v14": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V14"),
    "v15": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V15"),
    "v16": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V16"),
    "v17": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V17"),
    "v18": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V18"),
    "v19": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V19"),
    "chindemi": os.path.join(TRACE_RESULTS_BASE, "CHINDEMI_PARAMS"),
}

PARAMS_MAP = {
    "v1": DHURUVA_PARAMS, "v2": DHURUVA_PARAMS_V2, "v3": DHURUVA_PARAMS_V3,
    "v4": DHURUVA_PARAMS_V4, "v5": DHURUVA_PARAMS_V5, "v6": DHURUVA_PARAMS_V6,
    "v7": DHURUVA_PARAMS_V7, "v8": DHURUVA_PARAMS_V8, "v9": DHURUVA_PARAMS_V9,
    "v10": DHURUVA_PARAMS_V10, "v11": DHURUVA_PARAMS_V11, "v12": DHURUVA_PARAMS_V12,
    "v13": DHURUVA_PARAMS_V13, "v14": DHURUVA_PARAMS_V14, "v15": DHURUVA_PARAMS_V15,
    "v16": DHURUVA_PARAMS_V16, "v17": DHURUVA_PARAMS_V17, "v18": DHURUVA_PARAMS_V18,
    "v19": DHURUVA_PARAMS_V19,
    "chindemi": CHINDEMI_PARAMS,
}


def collect_rho(trace_results_dir):
    """Walk the results tree and collect per-synapse rho rows.

    Returns a list of dicts with keys:
        pregid, postgid, freq, delay, syn_id,
        initial_rho_raw, final_rho_raw, initial_rho, final_rho
    """
    rows = []

    for pair_dir in sorted(glob.glob(os.path.join(trace_results_dir, "*"))):
        if not os.path.isdir(pair_dir) or os.path.basename(pair_dir) == "logs":
            continue

        pair_name = os.path.basename(pair_dir)
        try:
            pre_gid, post_gid = map(int, pair_name.split("-"))
        except ValueError:
            continue

        for freq_dt_dir in sorted(glob.glob(os.path.join(pair_dir, "*Hz_*"))):
            freq_dt = os.path.basename(freq_dt_dir)   # e.g. "10Hz_10ms"
            try:
                freq_str, dt_str = freq_dt.split("_")
                freq  = float(freq_str.replace("Hz", ""))
                delay = float(dt_str.replace("ms", ""))
            except ValueError:
                continue

            pkl_path = os.path.join(freq_dt_dir, "simulation_traces.pkl")
            if not os.path.exists(pkl_path):
                continue

            try:
                with open(pkl_path, "rb") as f:
                    sim_data = pickle.load(f)
            except Exception as e:
                print(f"  Skipping {pair_name}/{freq_dt} — corrupted pickle ({e})")
                continue

            if not sim_data or "rho_GB" not in sim_data:
                print(f"  Skipping {pair_name}/{freq_dt} — missing rho_GB")
                continue

            rho_gb_raw = sim_data["rho_GB"]

            # Normalise to dict {syn_id: 1d array}
            if isinstance(rho_gb_raw, dict):
                rho_dict = {k: np.asarray(v, dtype=float) for k, v in rho_gb_raw.items()}
            else:
                arr = np.asarray(rho_gb_raw, dtype=float)
                if arr.ndim == 2:
                    if arr.shape[0] > arr.shape[1]:
                        arr = arr.T          # (synapses, timepoints)
                else:
                    arr = arr.reshape(1, -1)
                rho_dict = {i: arr[i] for i in range(arr.shape[0])}

            for syn_id, trace in rho_dict.items():
                if len(trace) == 0:
                    continue
                irho_raw = float(trace[0])
                frho_raw = float(trace[-1])
                rows.append({
                    "pregid":          pre_gid,
                    "postgid":         post_gid,
                    "freq":            freq,
                    "delay":           delay,
                    "syn_id":          syn_id,
                    "initial_rho_raw": irho_raw,
                    "final_rho_raw":   frho_raw,
                    "initial_rho":     1 if irho_raw >= 0.5 else 0,
                    "final_rho":       1 if frho_raw >= 0.5 else 0,
                })

    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Collect per-synapse rho values into a CSV")
    parser.add_argument("--params", choices=list(TRACE_RESULTS_DIRS.keys()), default="v17",
                        help="Parameter set to collect (default: v17)")
    parser.add_argument("--results-dir", default=None,
                        help="Override default results directory")
    parser.add_argument("--output", default=None,
                        help="Output CSV filename (default: rho_<params>.csv)")
    args = parser.parse_args()

    trace_results_dir = args.results_dir or TRACE_RESULTS_DIRS[args.params]
    output_csv        = args.output or f"rho_{args.params}.csv"

    print(f"Collecting rho values from: {trace_results_dir}")
    rows = collect_rho(trace_results_dir)

    if not rows:
        print("No data found.")
        sys.exit(1)

    df = pd.DataFrame(rows, columns=[
        "pregid", "postgid", "freq", "delay", "syn_id",
        "initial_rho_raw", "final_rho_raw", "initial_rho", "final_rho",
    ])
    df.to_csv(output_csv, index=False)
    print(f"Wrote {len(df)} rows ({df[['pregid','postgid','freq','delay']].drop_duplicates().shape[0]} simulations) → {output_csv}")
    print(df.groupby("delay")[["initial_rho", "final_rho"]].mean().rename(
        columns={"initial_rho": "mean_initial_rho", "final_rho": "mean_final_rho"}
    ).to_string())
