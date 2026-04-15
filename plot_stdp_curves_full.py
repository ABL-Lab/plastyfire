#!/usr/bin/env python3
"""
Plot STDP curves from FULL simulation traces.

Unlike plot_stdp_curves.py, this script does NOT use basis-function extrapolation.
Instead it reads the voltage trace from each simulation_traces.pkl and computes
the EPSP ratio directly via plastyfire.ephysutils.Experiment.compute_epsp_ratio(),
exactly as done in Plot_Chindemi_Fig_6.ipynb.

The simulation is assumed to have the standard L5TTPC_L5TTPC structure:
    - C01 connectivity test : 4 min of 0.25 Hz test pulses  (first 60 pre-spikes)
    - Induction protocol    : 60 pre-post pairs at 10 Hz
    - C02 connectivity test : 4 min of 0.25 Hz test pulses  (last 60 pre-spikes)

EPSP ratio = mean(C02 EPSPs[-n:]) / mean(C01 EPSPs[-n:])  with n=60

Usage:
    python plot_stdp_curves_full.py --params chindemi
    python plot_stdp_curves_full.py --params chindemi --results-dir /path/to/full_trace_results/CHINDEMI_PARAMS
    python plot_stdp_curves_full.py --params v19 --workers 8
"""

import os
import sys
import argparse
import pickle
import glob
import time
import multiprocessing
import concurrent.futures
from concurrent.futures.process import BrokenProcessPool

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from plastyfire.ephysutils import Experiment
from submit_l5ttpc_traces import (
    CHINDEMI_PARAMS, DHURUVA_PARAMS, DHURUVA_PARAMS_V2, DHURUVA_PARAMS_V3, DHURUVA_PARAMS_V4,
    DHURUVA_PARAMS_V5, DHURUVA_PARAMS_V6, DHURUVA_PARAMS_V7, DHURUVA_PARAMS_V8,
    DHURUVA_PARAMS_V9, DHURUVA_PARAMS_V10, DHURUVA_PARAMS_V11, DHURUVA_PARAMS_V12,
    DHURUVA_PARAMS_V13, DHURUVA_PARAMS_V14, DHURUVA_PARAMS_V15, DHURUVA_PARAMS_V16,
    DHURUVA_PARAMS_V17, DHURUVA_PARAMS_V18, DHURUVA_PARAMS_V19
)

# Default output dirs — mirrors plot_stdp_curves.py but points at full_trace_results
FULL_TRACE_RESULTS_BASE = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/full_trace_results"
TRACE_RESULTS_BASE      = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/trace_results"
 
TRACE_RESULTS_DIRS = {
    "v1":       os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS"),
    "v2":       os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V2"),
    "v3":       os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V3"),
    "v4":       os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V4"),
    "v5":       os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_STDP_V5"),
    "v6":       os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V6"),
    "v7":       os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V7"),
    "v8":       os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V8"),
    "v9":       os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V9"),
    "v10":      os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V10"),
    "v11":      os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V11"),
    "v12":      os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V12"),
    "v13":      os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V13"),
    "v14":      os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V14"),
    "v15":      os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V15"),
    "v16":      os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V16"),
    "v17":      os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V17"),
    "v18":      os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V18"),
    "v19":      os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V19"),
    "chindemi": os.path.join(FULL_TRACE_RESULTS_BASE, "CHINDEMI_PARAMS"),
}

# L5TTPC_L5TTPC protocol parameters (from configs/L5TTPC_L5TTPC.yaml)
C01_DURATION_MIN = 4.   # minutes
C02_DURATION_MIN = 4.   # minutes
TEST_PULSE_PERIOD_S = 4. # seconds between test pulses  → 0.25 Hz
N_EPSP = 60              # test pulses in each C0x window (4 min / 4 s = 60)


def _process_single_sim(args):
    """Worker: load one pkl and compute EPSP ratio directly from voltage trace.

    Parameters
    ----------
    args : tuple
        (pair_name, freq_dt, pkl_path, dt)

    Returns
    -------
    dict or None
        {"pair", "dt", "epsp_ratio", "n_0to1", "n_0to0", "n_1to1", "n_1to0"} on success, None on failure.
    """
    pair_name, freq_dt, pkl_path, dt = args

    try:
        with open(pkl_path, "rb") as f:
            sim_data = pickle.load(f)
    except Exception as e:
        print(f"Skipping {pair_name} {freq_dt} — corrupted pickle ({type(e).__name__}: {e}). Deleting {pkl_path}")
        try:
            os.remove(pkl_path)
        except OSError:
            pass
        return None

    if not sim_data or "v" not in sim_data or "t" not in sim_data:
        print(f"Skipping {pair_name} {freq_dt} — missing 't'/'v' in pickle.")
        return None

    if "prespikes" not in sim_data:
        print(f"Skipping {pair_name} {freq_dt} — missing 'prespikes' in pickle.")
        return None

    # Extract rho_GB initial/final states
    rho_gb_raw = sim_data.get("rho_GB")
    n_0to1 = n_0to0 = n_1to1 = n_1to0 = 0
    if rho_gb_raw is not None:
        if isinstance(rho_gb_raw, dict):
            rho_trace = np.array(list(rho_gb_raw.values()))  # (n_syn, n_t)
        else:
            rho_trace = np.asarray(rho_gb_raw)
            if rho_trace.ndim == 2 and rho_trace.shape[0] > rho_trace.shape[1]:
                rho_trace = rho_trace.T
            elif rho_trace.ndim == 1:
                rho_trace = rho_trace.reshape(1, -1)
        initial = np.array([1 if r[0]  >= 0.5 else 0 for r in rho_trace])
        final   = np.array([1 if r[-1] >= 0.5 else 0 for r in rho_trace])
        n_0to1 = int(np.sum((initial == 0) & (final == 1)))
        n_0to0 = int(np.sum((initial == 0) & (final == 0)))
        n_1to1 = int(np.sum((initial == 1) & (final == 1)))
        n_1to0 = int(np.sum((initial == 1) & (final == 0)))

    try:
        exp = Experiment(
            data=sim_data,
            c01duration=C01_DURATION_MIN,
            c02duration=C02_DURATION_MIN,
            period=TEST_PULSE_PERIOD_S,
        )
        epsp_before, epsp_after, epsp_ratio, _, _ = exp.compute_epsp_ratio(n=N_EPSP, full=True)
        return {"pair": pair_name, "dt": dt,
                "epsp_before": epsp_before, "epsp_after": epsp_after, "epsp_ratio": epsp_ratio,
                "n_0to1": n_0to1, "n_0to0": n_0to0, "n_1to1": n_1to1, "n_1to0": n_1to0}
    except Exception as e:
        print(f"Error computing EPSP ratio for {pair_name} {freq_dt}: {e}")
        return None


def _run_jobs_serial(jobs):
    data = []
    total_jobs = len(jobs)
    progress_every = max(1, total_jobs // 20)
    start_time = time.time()
    for idx, job in enumerate(jobs, start=1):
        result = _process_single_sim(job)
        if result is not None:
            data.append(result)
        if idx == 1 or idx % progress_every == 0 or idx == total_jobs:
            elapsed = time.time() - start_time
            print(
                f"Completed {idx}/{total_jobs} simulations "
                f"({100.0 * idx / total_jobs:5.1f}%) in {elapsed:7.1f}s "
                f"[valid results: {len(data)}]",
                flush=True,
            )
    return data


def _run_jobs_parallel(jobs, n_workers):
    data = []
    total_jobs = len(jobs)
    progress_every = max(1, total_jobs // 20)
    start_time = time.time()
    mp_context = multiprocessing.get_context("spawn")
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers, mp_context=mp_context) as pool:
        futures = [pool.submit(_process_single_sim, job) for job in jobs]
        for idx, future in enumerate(concurrent.futures.as_completed(futures), start=1):
            result = future.result()
            if result is not None:
                data.append(result)
            if idx == 1 or idx % progress_every == 0 or idx == total_jobs:
                elapsed = time.time() - start_time
                print(
                    f"Completed {idx}/{total_jobs} simulations "
                    f"({100.0 * idx / total_jobs:5.1f}%) in {elapsed:7.1f}s "
                    f"[valid results: {len(data)}]",
                    flush=True,
                )
    return data


def process_results(trace_results_dir, n_workers=None):
    """Scan directory tree, dispatch all pkl files to worker pool.

    Parameters
    ----------
    trace_results_dir : str
    n_workers : int or None
    """
    jobs = []

    for pair_dir in sorted(glob.glob(os.path.join(trace_results_dir, "*"))):
        if not os.path.isdir(pair_dir) or os.path.basename(pair_dir) in ("logs", "figures"):
            continue

        pair_name = os.path.basename(pair_dir)
        if "-" not in pair_name:
            continue

        for freq_dt_dir in sorted(glob.glob(os.path.join(pair_dir, "10Hz_*"))):
            freq_dt = os.path.basename(freq_dt_dir)
            dt_str  = freq_dt.replace("10Hz_", "").replace("ms", "")
            try:
                dt = float(dt_str)
            except ValueError:
                continue

            pkl_path = os.path.join(freq_dt_dir, "simulation_traces.pkl")
            if not os.path.exists(pkl_path):
                continue

            jobs.append((pair_name, freq_dt, pkl_path, dt))

    if not jobs:
        return pd.DataFrame()

    requested_workers = n_workers or multiprocessing.cpu_count()
    n_workers = max(1, min(int(requested_workers), len(jobs)))
    print(
        f"Processing {len(jobs)} simulations from {trace_results_dir} "
        f"across up to {n_workers} CPU cores...",
        flush=True,
    )

    current_workers = n_workers
    while True:
        if current_workers == 1:
            if current_workers != n_workers:
                print("Retrying serially after repeated worker-pool failures...", flush=True)
            data = _run_jobs_serial(jobs)
            break

        try:
            print(
                f"Running with {current_workers} worker processes (spawn context)...",
                flush=True,
            )
            data = _run_jobs_parallel(jobs, current_workers)
            break
        except BrokenProcessPool as e:
            next_workers = max(1, current_workers // 2)
            print(
                f"Worker pool failed abruptly with {current_workers} workers "
                f"({type(e).__name__}: {e}). Retrying with {next_workers} workers..."
                ,
                flush=True,
            )
            current_workers = next_workers

    return pd.DataFrame(data)


def plot_stdp_curve(df, version, output_filename):
    if df.empty:
        print("No data collected!")
        return

    summary = df.groupby("dt")["epsp_ratio"].agg(["mean", "sem"]).reset_index()
    summary = summary.sort_values("dt")

    # In vitro reference data — Markram 1997 (10 Hz)
    invitro_dt   = [-10,   5,      10]
    invitro_mean = [0.7922, 1.2038, 1.2013]
    invitro_sem  = [0.0259, 0.0644, 0.0626]

    plt.figure(figsize=(6, 5))

    plt.errorbar(summary["dt"], summary["mean"], yerr=summary["sem"],
                 fmt="o-", color="#66b3e6", label=f"in silico ({version})", capsize=0)

    plt.errorbar(invitro_dt, invitro_mean, yerr=invitro_sem,
                 fmt="o-", color="#ffaa55", label="in vitro (Markram 1997)", capsize=0)

    plt.axhline(1.0, color="k", linestyle="--", alpha=0.5)
    plt.axvline(0.0, color="k", linestyle="--", alpha=0.5)

    plt.xlabel(r"$\Delta t$ (ms)", fontsize=12)
    plt.ylabel("EPSP ratio", fontsize=12)
    plt.title(f"Frequency = 10 Hz ({version}) — full simulation", fontsize=14)
    plt.legend(frameon=False, loc="upper left")

    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_filename, dpi=300)
    print(f"Saved {output_filename}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot STDP curves from full simulation traces (no basis extrapolation)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python plot_stdp_curves_full.py --params chindemi
  python plot_stdp_curves_full.py --params chindemi \\
      --results-dir /path/to/full_trace_results/CHINDEMI_PARAMS
  python plot_stdp_curves_full.py --params v19 --workers 8
        """,
    )
    parser.add_argument(
        "--params",
        choices=["v1","v2","v3","v4","v5","v6","v7","v8","v9","v10",
                 "v11","v12","v13","v14","v15","v16","v17","v18","v19","chindemi"],
        default="chindemi",
        help="Parameter set (default: chindemi)",
    )
    parser.add_argument("--results-dir", type=str, default=None,
                        help="Override default results directory")
    parser.add_argument("--workers", type=int, default=None,
                        help="Worker processes (default: all CPU cores)")
    parser.add_argument("--output", type=str, default=None,
                        help="Output PNG filename (default: stdp_curve_full_10Hz_<params>.png)")
    args = parser.parse_args()

    version           = args.params
    trace_results_dir = args.results_dir or TRACE_RESULTS_DIRS[version]
    output_filename   = args.output or f"stdp_curve_full_10Hz_{version}.png"

    print("\n" + "=" * 55)
    print(f"STDP Curve (full simulation) — {version}")
    print(f"Results dir : {trace_results_dir}")
    print(f"Output      : {output_filename}")
    print("=" * 55)

    df = process_results(trace_results_dir, n_workers=args.workers)

    if not df.empty:
        n_pairs   = df["pair"].nunique()
        total_sim = len(df)
        print(f"Pairs processed  : {n_pairs}")
        print(f"Total simulations: {total_sim}")
        print("\nSimulation counts by Δt:")
        print(df.groupby("dt").size().to_string())

        summary = df.groupby("dt")["epsp_ratio"].agg(["mean", "sem", "count"]).reset_index()
        summary = summary.sort_values("dt")
        print("\nStatistics by Δt (mean of per-pair ratios):")
        print(summary.to_string(index=False))

        invitro = {-50.0: 1.0, -10.0: 0.7922, 5.0: 1.2038, 10.0: 1.2013}
        overlap = summary[summary["dt"].isin(invitro.keys())]
        if not overlap.empty:
            errors = [(row["mean"] - invitro[row["dt"]]) ** 2 for _, row in overlap.iterrows()]
            rmse = np.sqrt(np.mean(errors))
            print(f"\nRMSE vs in-vitro (Δt = -50, -10, 5, 10 ms): {rmse:.4f}")

        # Rho state transition summary
        rho_cols = ["n_0to1", "n_0to0", "n_1to1", "n_1to0"]
        if all(c in df.columns for c in rho_cols):
            rho_summary = df.groupby("dt")[rho_cols].sum().reset_index()
            rho_summary["n_total"] = rho_summary[rho_cols].sum(axis=1)
            for c in rho_cols:
                rho_summary[c + "_pct"] = 100.0 * rho_summary[c] / rho_summary["n_total"]

            print("\n" + "=" * 65)
            print("Rho state transitions  (per Δt, summed over all pairs & synapses)")
            print("=" * 65)
            print(f"{'Δt':>8}  {'0→1 (LTP)':>12}  {'0→0 (dep)':>12}  {'1→1 (pot)':>12}  {'1→0 (LTD)':>12}  {'total':>7}")
            print("-" * 65)
            for _, row in rho_summary.sort_values("dt").iterrows():
                print(
                    f"{row['dt']:>8.0f}  "
                    f"{int(row['n_0to1']):>6} ({row['n_0to1_pct']:5.1f}%)  "
                    f"{int(row['n_0to0']):>6} ({row['n_0to0_pct']:5.1f}%)  "
                    f"{int(row['n_1to1']):>6} ({row['n_1to1_pct']:5.1f}%)  "
                    f"{int(row['n_1to0']):>6} ({row['n_1to0_pct']:5.1f}%)  "
                    f"{int(row['n_total']):>7}"
                )
            print("=" * 65)

        plot_stdp_curve(df, version, output_filename)
        print("\nDone.")
    else:
        print(f"No results found in {trace_results_dir}")

    print("=" * 55)
