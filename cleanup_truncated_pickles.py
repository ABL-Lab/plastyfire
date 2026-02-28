#!/usr/bin/env python3
"""
Find and delete truncated/corrupted simulation_traces.pkl files.
Run with --dry-run first to see what would be deleted.

Usage:
    python cleanup_truncated_pickles.py --version v7 --dry-run
    python cleanup_truncated_pickles.py --version v7
"""
import os
import glob
import pickle
import argparse

TRACE_RESULTS_BASE = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/trace_results"
TRACE_RESULTS_DIRS = {
    "v1": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS"),
    "v2": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V2"),
    "v3": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V3"),
    "v4": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V4"),
    "v5": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_STDP_V5"),
    "v6": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V6"),
    "v7": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V7"),
}

parser = argparse.ArgumentParser()
parser.add_argument("--version", choices=list(TRACE_RESULTS_DIRS.keys()), default="v7")
parser.add_argument("--dry-run", action="store_true", help="Print what would be deleted without deleting")
parser.add_argument("--results-dir", type=str, default=None)
args = parser.parse_args()

results_dir = args.results_dir or TRACE_RESULTS_DIRS[args.version]
pattern = os.path.join(results_dir, "*", "*", "simulation_traces.pkl")
all_pkls = sorted(glob.glob(pattern))

print(f"Scanning {len(all_pkls)} pickle files in {results_dir}...")

ok, bad = 0, 0
for pkl_path in all_pkls:
    try:
        with open(pkl_path, "rb") as f:
            pickle.load(f)
        ok += 1
    except Exception as e:
        bad += 1
        rel = os.path.relpath(pkl_path, results_dir)
        if args.dry_run:
            print(f"  [DRY-RUN] Would delete: {rel}  ({type(e).__name__}: {e})")
        else:
            print(f"  Deleting: {rel}  ({type(e).__name__}: {e})")
            os.remove(pkl_path)
            # Also remove the parent dir if empty so SLURM re-runs it cleanly
            parent = os.path.dirname(pkl_path)
            if not os.listdir(parent):
                os.rmdir(parent)

print(f"\nDone. OK: {ok}, Corrupted {'(would delete)' if args.dry_run else '(deleted)'}: {bad}")
