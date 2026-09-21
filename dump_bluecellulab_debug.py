#!/usr/bin/env python3
"""
Dump bluecellulab simulation results to a CSV for detailed debugging.

For every bluecellulab_results/rho.h5 found under the refitting_results tree,
extracts per-pair per-protocol:
  pre_gid, post_gid, freq, delay,
  initial_rho_array, final_rho_array,
  mean_epsp_before, mean_epsp_after, epsp_ratio

Arrays are stored as JSON strings so the CSV remains single-line per row.

Usage:
    python dump_bluecellulab_debug.py
    python dump_bluecellulab_debug.py --freq 10 --output debug_bluecellulab_10Hz.csv
    python dump_bluecellulab_debug.py --freq 10 --workers 8 --ratio-method delta_method
    python dump_bluecellulab_debug.py --results-dir /path/to/refitting_results \
                                      --basis-dir /path/to/basis_results_edges_mini
"""

import argparse
import concurrent.futures
import glob
import json
import logging
import multiprocessing
import os
import sys

import h5py
import numpy as np
import pandas as pd

PLASTYFIRE_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(PLASTYFIRE_ROOT, "new_fitting"))
from effcai_to_epsp import fetch_epsp_ratio

logging.basicConfig(
    level=logging.WARNING,          # suppress per-pair INFO noise in worker procs
    format="%(asctime)s %(levelname)s %(message)s",
)
logger = logging.getLogger(__name__)

DEFAULT_RESULTS_DIR = os.path.join(PLASTYFIRE_ROOT, "refitting_results")
DEFAULT_BASIS_DIR   = os.path.join(PLASTYFIRE_ROOT, "basis_results_edges_mini")
NODE_POP            = "S1nonbarrel_neurons"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _find_basis_path(pre_gid, post_gid, basis_dir):
    path = os.path.join(basis_dir, f"basis_{pre_gid}_{post_gid}.csv")
    return path if os.path.exists(path) else None


def _read_rho_h5(rho_h5_path):
    """
    Return (initial_rho_list, final_rho_list) as raw float arrays (not binarised),
    plus the binarised lists used for EPSP extrapolation.
    rho.h5 layout:  report/<node_pop>/data  shape (n_timepoints, n_synapses)
                    first row = initial state, last row = final state.
    """
    with h5py.File(rho_h5_path, "r") as f:
        data = f[f"report/{NODE_POP}/data"][()]   # (n_timepoints, n_synapses)

    raw_initial = data[0].tolist()
    raw_final   = data[-1].tolist()
    # Binarise for EPSP extrapolation (same rule as plot_compare_ndamus_bluecellulab.py)
    bin_initial = [1 if v >= 0.5 else 0 for v in raw_initial]
    bin_final   = [1 if v >= 0.5 else 0 for v in raw_final]
    return raw_initial, raw_final, bin_initial, bin_final


# ---------------------------------------------------------------------------
# Per-simulation worker (top-level for pickling)
# ---------------------------------------------------------------------------

def _process_one(args):
    """
    Process a single (sim_dir, pre_gid, post_gid, freq, delay_ms, basis_path)
    tuple.  Returns a dict row or None on failure.
    """
    sim_dir, pre_gid, post_gid, freq, delay_ms, basis_path, ratio_method = args

    rho_h5 = os.path.join(sim_dir, "bluecellulab_results", "rho.h5")
    if not os.path.exists(rho_h5):
        return None

    try:
        raw_initial, raw_final, bin_initial, bin_final = _read_rho_h5(rho_h5)
    except Exception as e:
        logger.warning("rho.h5 read error %d-%d dt=%g: %s", pre_gid, post_gid, delay_ms, e)
        return None

    try:
        res = fetch_epsp_ratio(
            pre_gid, post_gid,
            freq, delay_ms,
            bin_initial, bin_final,
            basis_dir=os.path.dirname(basis_path),
            ratio_method=ratio_method,
        )
    except Exception as e:
        logger.warning("EPSP ratio error %d-%d dt=%g: %s", pre_gid, post_gid, delay_ms, e)
        res = {
            "epsp_before_mean": float("nan"),
            "epsp_after_mean":  float("nan"),
            "ratio_mean":       float("nan"),
        }

    n_syn    = len(raw_initial)
    n_pot_i  = int(sum(bin_initial))
    n_pot_f  = int(sum(bin_final))
    n_0to1   = int(sum(1 for i, f in zip(bin_initial, bin_final) if i == 0 and f == 1))
    n_1to0   = int(sum(1 for i, f in zip(bin_initial, bin_final) if i == 1 and f == 0))

    return {
        "pre_gid":          pre_gid,
        "post_gid":         post_gid,
        "freq_hz":          freq,
        "delay_ms":         delay_ms,
        "n_synapses":       n_syn,
        "n_potentiated_initial": n_pot_i,
        "n_potentiated_final":   n_pot_f,
        "n_0to1":           n_0to1,
        "n_1to0":           n_1to0,
        # raw float rho values (not binarised) — useful for spotting borderline cases
        "initial_rho_array": json.dumps([round(v, 6) for v in raw_initial]),
        "final_rho_array":   json.dumps([round(v, 6) for v in raw_final]),
        # binarised arrays used in EPSP extrapolation
        "initial_rho_bin":   json.dumps(bin_initial),
        "final_rho_bin":     json.dumps(bin_final),
        # EPSP results
        "mean_epsp_before":  res["epsp_before_mean"],
        "mean_epsp_after":   res["epsp_after_mean"],
        "epsp_ratio":        res["ratio_mean"],
    }


# ---------------------------------------------------------------------------
# Collect jobs
# ---------------------------------------------------------------------------

def collect_jobs(results_dir, basis_dir, freq, ratio_method):
    jobs = []
    pattern = os.path.join(
        results_dir, "fitting", "*", "seed*", "*_STDP",
        "simulations", "*-*", f"{int(freq)}Hz_*",
    )
    sim_dirs = sorted(glob.glob(pattern))
    missing_basis = 0

    for sim_dir in sim_dirs:
        # Skip if no bluecellulab result yet
        if not os.path.exists(os.path.join(sim_dir, "bluecellulab_results", "rho.h5")):
            continue

        pair_name = os.path.basename(os.path.dirname(sim_dir))
        freq_dt   = os.path.basename(sim_dir)           # e.g. "10Hz_-10ms"

        try:
            pre_gid, post_gid = [int(x) for x in pair_name.split("-")]
        except ValueError:
            continue

        dt_str = freq_dt.split("_", 1)[1].replace("ms", "")
        try:
            delay_ms = float(dt_str)
        except ValueError:
            continue

        basis_path = _find_basis_path(pre_gid, post_gid, basis_dir)
        if basis_path is None:
            missing_basis += 1
            logger.warning("No basis CSV for %s — EPSP will be NaN", pair_name)
            # Still include the row with NaN EPSP so rho arrays are visible
            basis_path = "__missing__"

        jobs.append((sim_dir, pre_gid, post_gid, float(freq), delay_ms,
                     basis_path, ratio_method))

    if missing_basis:
        print(f"  Warning: {missing_basis} pairs have no basis CSV — EPSP columns will be NaN")

    return jobs


# ---------------------------------------------------------------------------
# Top-level dispatch (must be module-level to be picklable by ProcessPoolExecutor)
# ---------------------------------------------------------------------------

def _dispatch(job_args):
    """Route to _process_one, or return NaN-EPSP row when basis CSV is missing."""
    if job_args[5] == "__missing__":
        sim_dir, pre_gid, post_gid, freq, delay_ms, _, _ratio_method = job_args
        rho_h5 = os.path.join(sim_dir, "bluecellulab_results", "rho.h5")
        try:
            raw_initial, raw_final, bin_initial, bin_final = _read_rho_h5(rho_h5)
        except Exception:
            return None
        return {
            "pre_gid":               pre_gid,
            "post_gid":              post_gid,
            "freq_hz":               freq,
            "delay_ms":              delay_ms,
            "n_synapses":            len(raw_initial),
            "n_potentiated_initial": sum(1 if v >= 0.5 else 0 for v in raw_initial),
            "n_potentiated_final":   sum(1 if v >= 0.5 else 0 for v in raw_final),
            "n_0to1":                sum(1 for i, f in zip(raw_initial, raw_final)
                                         if i < 0.5 and f >= 0.5),
            "n_1to0":                sum(1 for i, f in zip(raw_initial, raw_final)
                                         if i >= 0.5 and f < 0.5),
            "initial_rho_array":     json.dumps([round(v, 6) for v in raw_initial]),
            "final_rho_array":       json.dumps([round(v, 6) for v in raw_final]),
            "initial_rho_bin":       json.dumps([1 if v >= 0.5 else 0 for v in raw_initial]),
            "final_rho_bin":         json.dumps([1 if v >= 0.5 else 0 for v in raw_final]),
            "mean_epsp_before":      float("nan"),
            "mean_epsp_after":       float("nan"),
            "epsp_ratio":            float("nan"),
        }
    return _process_one(job_args)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Dump bluecellulab rho + EPSP debug info to CSV",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--results-dir", default=DEFAULT_RESULTS_DIR,
                        help=f"Root refitting results dir (default: {DEFAULT_RESULTS_DIR})")
    parser.add_argument("--basis-dir",   default=DEFAULT_BASIS_DIR,
                        help=f"Basis results dir (default: {DEFAULT_BASIS_DIR})")
    parser.add_argument("--freq",        type=float, default=10.,
                        help="Induction frequency in Hz (default: 10)")
    parser.add_argument("--workers",     type=int, default=None,
                        help="Parallel workers (default: cpu_count)")
    parser.add_argument("--ratio-method", default="delta_method",
                        choices=["simple", "delta_method"],
                        help="EPSP ratio method (default: delta_method)")
    parser.add_argument("--output",      default=None,
                        help="Output CSV path (default: bluecellulab_debug_<freq>Hz.csv)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable DEBUG logging")
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    output = args.output or f"bluecellulab_debug_{int(args.freq)}Hz.csv"

    print(f"Results dir  : {args.results_dir}")
    print(f"Basis dir    : {args.basis_dir}")
    print(f"Freq         : {int(args.freq)} Hz")
    print(f"Ratio method : {args.ratio_method}")
    print(f"Output       : {output}")

    jobs = collect_jobs(args.results_dir, args.basis_dir, args.freq, args.ratio_method)
    if not jobs:
        print("No bluecellulab_results/rho.h5 files found. Nothing to do.")
        sys.exit(0)

    print(f"Found {len(jobs)} simulation results to process...")

    n_workers = args.workers or multiprocessing.cpu_count()
    rows = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as pool:
        for i, result in enumerate(pool.map(_dispatch, jobs), 1):
            if result is not None:
                rows.append(result)
            if i % 50 == 0 or i == len(jobs):
                print(f"  Processed {i}/{len(jobs)} ...", flush=True)

    if not rows:
        print("No results collected.")
        sys.exit(1)

    df = pd.DataFrame(rows, columns=[
        "pre_gid", "post_gid", "freq_hz", "delay_ms",
        "n_synapses",
        "n_potentiated_initial", "n_potentiated_final",
        "n_0to1", "n_1to0",
        "initial_rho_array", "final_rho_array",
        "initial_rho_bin",   "final_rho_bin",
        "mean_epsp_before", "mean_epsp_after", "epsp_ratio",
    ])
    df = df.sort_values(["delay_ms", "pre_gid", "post_gid"]).reset_index(drop=True)

    df.to_csv(output, index=False)
    print(f"\nSaved {len(df)} rows → {output}")

    # Quick summary
    valid = df.dropna(subset=["epsp_ratio"])
    print(f"\nSummary (freq={int(args.freq)} Hz, {len(valid)} pairs with valid EPSP):")
    summary = (valid.groupby("delay_ms")["epsp_ratio"]
               .agg(["mean", "sem", "count"])
               .rename(columns={"mean": "epsp_ratio_mean", "sem": "epsp_ratio_sem",
                                 "count": "n_pairs"}))
    print(summary.to_string())


if __name__ == "__main__":
    main()
