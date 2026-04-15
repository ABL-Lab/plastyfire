#!/usr/bin/env python3
"""
Extract EPSP ratios from completed neurodamus SONATA simulations and plot STDP curves.

Reads voltage traces directly from out/out.h5 (SONATA HDF5) and prespikes from
prespikes.h5.  Does NOT use the basis-function extrapolation.

Directory structure expected:
    <results_dir>/fitting/<n_label>/seed<seed>/<pre>_<post>_STDP/simulations/
        <pre_gid>-<post_gid>/<freq>Hz_<dt>ms/
            prespikes.h5
            out/out.h5

Usage:
    python plot_refitting_stdp.py
    python plot_refitting_stdp.py --results-dir /path/to/refitting_results
    python plot_refitting_stdp.py --freq 10 --workers 8 --output my_stdp.png
"""

import os
import sys
import glob
import argparse
import multiprocessing
import concurrent.futures

import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plastyfire.ephysutils import Experiment

DEFAULT_RESULTS_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "refitting_results"
)

# Protocol parameters (from configs/L5TTPC_L5TTPC_STDP.yaml)
C01_DURATION_MIN = 4.
C02_DURATION_MIN = 4.
TEST_PULSE_PERIOD_S = 4.
N_EPSP = 60

# In-vitro reference — Markram et al. 1997 (10 Hz)
INVITRO_DT   = [-10,    5,      10]
INVITRO_MEAN = [0.7922, 1.2038, 1.2013]
INVITRO_SEM  = [0.0259, 0.0644, 0.0626]

NODE_POP = "S1nonbarrel_neurons"


def _read_soma_voltage(soma_h5_path):
    """Return (t, v) arrays from a SONATA soma report."""
    with h5py.File(soma_h5_path, "r") as h:
        t_info = h[f"report/{NODE_POP}/mapping/time"][()]  # [start, end, dt]
        t_start, t_end, dt = float(t_info[0]), float(t_info[1]), float(t_info[2])
        n = round((t_end - t_start) / dt)
        t = np.linspace(t_start, t_start + (n - 1) * dt, n)
        v = h[f"report/{NODE_POP}/data"][()].ravel().astype(np.float64)
    min_len = min(len(t), len(v))
    return t[:min_len], v[:min_len]


def _read_prespikes(prespikes_h5_path):
    """Return sorted spike timestamps (ms) from a SONATA spikes file."""
    with h5py.File(prespikes_h5_path, "r") as h:
        ts = h[f"spikes/{NODE_POP}/timestamps"][()]
    return np.sort(ts.astype(np.float64))


def _process_single_sim(args):
    sim_dir, pre_gid, post_gid, dt_ms = args
    pair_name = f"{pre_gid}-{post_gid}"

    soma_h5      = os.path.join(sim_dir, "out", "soma.h5")
    prespikes_h5 = os.path.join(sim_dir, "prespikes.h5")

    if not os.path.exists(soma_h5):
        return None
    if not os.path.exists(prespikes_h5):
        return None

    try:
        t, v      = _read_soma_voltage(soma_h5)
        prespikes = _read_prespikes(prespikes_h5)
    except Exception as e:
        print(f"Read error {pair_name} dt={dt_ms}: {e}")
        return None

    try:
        exp = Experiment(
            data={"t": t, "v": v, "prespikes": prespikes},
            c01duration=C01_DURATION_MIN,
            c02duration=C02_DURATION_MIN,
            period=TEST_PULSE_PERIOD_S,
        )
        epsp_before, epsp_after, epsp_ratio, epsp_before_std, epsp_after_std = \
            exp.compute_epsp_ratio(n=N_EPSP, full=True)
    except Exception as e:
        print(f"EPSP error {pair_name} dt={dt_ms}: {e}")
        return None

    return {
        "pair":            pair_name,
        "pre_gid":         pre_gid,
        "post_gid":        post_gid,
        "dt":              dt_ms,
        "epsp_before":     epsp_before,
        "epsp_after":      epsp_after,
        "epsp_ratio":      epsp_ratio,
        "epsp_before_std": epsp_before_std,
        "epsp_after_std":  epsp_after_std,
    }


def collect_jobs(results_dir, freq=10):
    """Scan results_dir for completed simulations (soma.h5 present)."""
    jobs = []
    freq_pattern = f"{int(freq)}Hz_*"

    sim_dirs = glob.glob(
        os.path.join(results_dir, "fitting", "*", "seed*", "*_STDP", "simulations",
                     "*-*", freq_pattern)
    )

    for sim_dir in sorted(sim_dirs):
        if not os.path.exists(os.path.join(sim_dir, "out", "soma.h5")):
            continue

        pair_name = os.path.basename(os.path.dirname(sim_dir))
        freq_dt   = os.path.basename(sim_dir)

        try:
            pre_gid, post_gid = [int(x) for x in pair_name.split("-")]
        except ValueError:
            continue

        dt_str = freq_dt.split("_", 1)[1].replace("ms", "")
        try:
            dt_ms = float(dt_str)
        except ValueError:
            continue

        jobs.append((sim_dir, pre_gid, post_gid, dt_ms))

    return jobs


def process_results(results_dir, freq=10, n_workers=None):
    """Collect all jobs and process them in parallel."""
    jobs = collect_jobs(results_dir, freq=freq)
    if not jobs:
        print(f"No completed sims found in {results_dir} (freq={freq} Hz)")
        return pd.DataFrame()

    n_workers = n_workers or multiprocessing.cpu_count()
    print(f"Processing {len(jobs)} simulations with {n_workers} workers...")

    data = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as pool:
        for result in pool.map(_process_single_sim, jobs):
            if result is not None:
                data.append(result)

    return pd.DataFrame(data)


def plot_stdp_curve(df, output_filename, freq=10):
    """Plot STDP curve (in-silico vs Markram 1997) and save."""
    summary = df.groupby("dt")["epsp_ratio"].agg(["mean", "sem"]).reset_index()
    summary = summary.sort_values("dt")

    fig, ax = plt.subplots(figsize=(6, 5))

    ax.errorbar(summary["dt"], summary["mean"], yerr=summary["sem"],
                fmt="o-", color="#66b3e6", label="in silico (refitting)", capsize=0)

    ax.errorbar(INVITRO_DT, INVITRO_MEAN, yerr=INVITRO_SEM,
                fmt="o-", color="#ffaa55", label="in vitro (Markram 1997)", capsize=0)

    ax.axhline(1.0, color="k", linestyle="--", alpha=0.5)
    ax.axvline(0.0, color="k", linestyle="--", alpha=0.5)

    ax.set_xlabel(r"$\Delta t$ (ms)", fontsize=12)
    ax.set_ylabel("EPSP ratio", fontsize=12)
    ax.set_title(f"Frequency = {freq} Hz — refitting (full simulation)", fontsize=13)
    ax.legend(frameon=False, loc="upper left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_filename, dpi=300)
    print(f"Saved {output_filename}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract EPSP ratios from neurodamus SONATA output and plot STDP curves",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python plot_refitting_stdp.py
  python plot_refitting_stdp.py --results-dir /path/to/refitting_results
  python plot_refitting_stdp.py --freq 10 --workers 8 --output refitting_10Hz.png
        """,
    )
    parser.add_argument("--results-dir", default=DEFAULT_RESULTS_DIR,
                        help=f"Root refitting results directory (default: {DEFAULT_RESULTS_DIR})")
    parser.add_argument("--freq", type=float, default=10.,
                        help="Induction frequency in Hz (default: 10)")
    parser.add_argument("--workers", type=int, default=None,
                        help="Worker processes (default: all CPU cores)")
    parser.add_argument("--output", default=None,
                        help="Output PNG file (default: stdp_curve_refitting_<freq>Hz.png)")
    parser.add_argument("--csv", default=None,
                        help="Also save results DataFrame as CSV")
    args = parser.parse_args()

    output = args.output or f"stdp_curve_refitting_{int(args.freq)}Hz.png"

    print("\n" + "=" * 60)
    print(f"STDP Curve — refitting results ({int(args.freq)} Hz, full sims)")
    print(f"Results dir : {args.results_dir}")
    print(f"Output      : {output}")
    print("=" * 60)

    df = process_results(args.results_dir, freq=args.freq, n_workers=args.workers)

    if df.empty:
        print("No data collected — exiting.")
        return

    n_pairs = df["pair"].nunique()
    print(f"\nPairs processed  : {n_pairs}")
    print(f"Total simulations: {len(df)}")
    print("\nSim counts by Δt:")
    print(df.groupby("dt").size().to_string())

    summary = df.groupby("dt")["epsp_ratio"].agg(["mean", "sem", "count"]).reset_index()
    print("\nEPSP ratio statistics by Δt:")
    print(summary.sort_values("dt").to_string(index=False))

    invitro_map = dict(zip(INVITRO_DT, INVITRO_MEAN))
    overlap = summary[summary["dt"].isin(invitro_map.keys())]
    if not overlap.empty:
        errors = [(row["mean"] - invitro_map[row["dt"]]) ** 2 for _, row in overlap.iterrows()]
        rmse = np.sqrt(np.mean(errors))
        print(f"\nRMSE vs in-vitro (Δt = {INVITRO_DT}): {rmse:.4f}")

    if args.csv:
        df.to_csv(args.csv, index=False)
        print(f"Saved CSV: {args.csv}")

    plot_stdp_curve(df, output, freq=args.freq)
    print("\nDone.")
    print("=" * 60)


if __name__ == "__main__":
    main()
