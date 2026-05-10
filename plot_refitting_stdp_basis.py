#!/usr/bin/env python3
"""
Extract EPSP ratios from completed neurodamus SONATA simulations and plot STDP curves.

Reads rho traces directly from out/rho.h5 (SONATA HDF5), binarizes initial and
final rho per synapse, then extrapolates EPSP ratios via the basis linear
superposition (same approach as plot_stdp_curves.py).

Directory structure expected:
    <results_dir>/fitting/<n_label>/seed<seed>/<pre>_<post>_STDP/simulations/
        <pre_gid>-<post_gid>/<freq>Hz_<dt>ms/
            out/rho.h5

Usage:
    python plot_refitting_stdp_basis.py
    python plot_refitting_stdp_basis.py --results-dir /path/to/refitting_results
    python plot_refitting_stdp_basis.py --freq 10 --workers 8 --output my_stdp.png
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

PLASTYFIRE_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(PLASTYFIRE_ROOT, "new_fitting"))
from effcai_to_epsp import fetch_epsp_ratio

DEFAULT_RESULTS_DIR = os.path.join(PLASTYFIRE_ROOT, "refitting_results")
DEFAULT_BASIS_DIRS = [
    os.path.join(PLASTYFIRE_ROOT, "basis_results"),
    os.path.join(PLASTYFIRE_ROOT, "basis_results_old"),
]

NODE_POP = "S1nonbarrel_neurons"

# In-vitro reference — Markram et al. 1997 (10 Hz)
INVITRO_DT   = [-10,    5,      10]
INVITRO_MEAN = [0.7922, 1.2038, 1.2013]
INVITRO_SEM  = [0.0259, 0.0644, 0.0626]


def _resolve_basis_dirs(basis_dir=None):
    if basis_dir is None:
        candidates = DEFAULT_BASIS_DIRS
    elif isinstance(basis_dir, (list, tuple)):
        candidates = list(basis_dir)
    else:
        candidates = [basis_dir]
    existing = [p for p in candidates if os.path.isdir(p)]
    if not existing:
        raise FileNotFoundError(f"No basis directory found. Checked: {', '.join(candidates)}")
    return existing


def _find_basis_path(pre_gid, post_gid, basis_dirs):
    filename = f"basis_{pre_gid}_{post_gid}.csv"
    for basis_dir in basis_dirs:
        path = os.path.join(basis_dir, filename)
        if os.path.exists(path):
            return path
    return None


def _read_rho_h5(rho_h5_path):
    """Return (initial_rho_binary, final_rho_binary, n_0to1, n_0to0, n_1to1, n_1to0)."""
    with h5py.File(rho_h5_path, "r") as f:
        data = f[f"report/{NODE_POP}/data"][()]  # (n_timepoints, n_synapses)
        # ND element_ids are local synapse IDs in non-sequential order; sort ascending
        # to match basis CSV positions (BCL iterates in ascending local synapse ID order).
        eids     = f[f"report/{NODE_POP}/mapping/element_ids"][()]
        sort_idx = np.argsort(eids)
        data     = data[:, sort_idx]

    initial = [1 if v >= 0.5 else 0 for v in data[0]]
    final   = [1 if v >= 0.5 else 0 for v in data[-1]]

    n_0to1 = sum(1 for i, f in zip(initial, final) if i == 0 and f == 1)
    n_0to0 = sum(1 for i, f in zip(initial, final) if i == 0 and f == 0)
    n_1to1 = sum(1 for i, f in zip(initial, final) if i == 1 and f == 1)
    n_1to0 = sum(1 for i, f in zip(initial, final) if i == 1 and f == 0)

    return initial, final, n_0to1, n_0to0, n_1to1, n_1to0


def _process_single_sim(args):
    sim_dir, pre_gid, post_gid, dt_ms, basis_path = args
    pair_name = f"{pre_gid}-{post_gid}"

    rho_h5 = os.path.join(sim_dir, "out", "rho.h5")
    if not os.path.exists(rho_h5):
        return None

    try:
        initial_rho, final_rho, n_0to1, n_0to0, n_1to1, n_1to0 = _read_rho_h5(rho_h5)
    except Exception as e:
        print(f"Read error {pair_name} dt={dt_ms}: {e}")
        return None

    try:
        res = fetch_epsp_ratio(
            pre_gid, post_gid,
            10.0, dt_ms,
            initial_rho, final_rho,
            basis_dir=os.path.dirname(basis_path),
            ratio_method="delta_method",
        )
    except Exception as e:
        print(f"EPSP error {pair_name} dt={dt_ms}: {e}")
        return None

    return {
        "pair":      pair_name,
        "pre_gid":   pre_gid,
        "post_gid":  post_gid,
        "dt":        dt_ms,
        "epsp_ratio": res["ratio_mean"],
        "n_0to1":    n_0to1,
        "n_0to0":    n_0to0,
        "n_1to1":    n_1to1,
        "n_1to0":    n_1to0,
    }


def collect_jobs(results_dir, basis_dirs, freq=10):
    jobs = []
    freq_pattern = f"{int(freq)}Hz_*"

    sim_dirs = glob.glob(
        os.path.join(results_dir, "fitting", "*", "seed*", "*_STDP",
                     "simulations", "*-*", freq_pattern)
    )

    for sim_dir in sorted(sim_dirs):
        if not os.path.exists(os.path.join(sim_dir, "out", "rho.h5")):
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

        basis_path = _find_basis_path(pre_gid, post_gid, basis_dirs)
        if basis_path is None:
            print(f"Skipping {pair_name}, no basis file.")
            continue

        jobs.append((sim_dir, pre_gid, post_gid, dt_ms, basis_path))

    return jobs


def process_results(results_dir, freq=10, n_workers=None, basis_dir=None):
    basis_dirs = _resolve_basis_dirs(basis_dir)
    jobs = collect_jobs(results_dir, basis_dirs, freq=freq)
    if not jobs:
        print(f"No completed sims with rho.h5 found in {results_dir} (freq={freq} Hz)")
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
    summary = df.groupby("dt")["epsp_ratio"].agg(["mean", "sem"]).reset_index()
    summary = summary.sort_values("dt")

    fig, ax = plt.subplots(figsize=(6, 5))

    ax.errorbar(summary["dt"], summary["mean"], yerr=summary["sem"],
                fmt="o-", color="#66b3e6", label="in silico (refitting, basis)", capsize=0)

    ax.errorbar(INVITRO_DT, INVITRO_MEAN, yerr=INVITRO_SEM,
                fmt="o-", color="#ffaa55", label="in vitro (Markram 1997)", capsize=0)

    ax.axhline(1.0, color="k", linestyle="--", alpha=0.5)
    ax.axvline(0.0, color="k", linestyle="--", alpha=0.5)

    ax.set_xlabel(r"$\Delta t$ (ms)", fontsize=12)
    ax.set_ylabel("EPSP ratio", fontsize=12)
    ax.set_title(f"Frequency = {freq} Hz — refitting (basis extrapolation)", fontsize=13)
    ax.legend(frameon=False, loc="upper left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_filename, dpi=300)
    print(f"Saved {output_filename}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract EPSP ratios via basis extrapolation from refitting rho.h5 and plot STDP curves",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python plot_refitting_stdp_basis.py
  python plot_refitting_stdp_basis.py --results-dir /path/to/refitting_results
  python plot_refitting_stdp_basis.py --freq 10 --workers 8 --output refitting_basis_10Hz.png
        """,
    )
    parser.add_argument("--results-dir", default=DEFAULT_RESULTS_DIR,
                        help=f"Root refitting results directory (default: {DEFAULT_RESULTS_DIR})")
    parser.add_argument("--basis-dir", default=None,
                        help="Basis directory (default: auto-detect basis_results, then basis_results_old)")
    parser.add_argument("--freq", type=float, default=10.,
                        help="Induction frequency in Hz (default: 10)")
    parser.add_argument("--workers", type=int, default=None,
                        help="Worker processes (default: all CPU cores)")
    parser.add_argument("--output", default=None,
                        help="Output PNG file (default: stdp_curve_refitting_basis_<freq>Hz.png)")
    parser.add_argument("--csv", default=None,
                        help="Also save results DataFrame as CSV")
    args = parser.parse_args()

    output = args.output or f"stdp_curve_refitting_basis_{int(args.freq)}Hz.png"
    basis_dirs = _resolve_basis_dirs(args.basis_dir)

    print("\n" + "=" * 60)
    print(f"STDP Curve — refitting results, basis extrapolation ({int(args.freq)} Hz)")
    print(f"Results dir : {args.results_dir}")
    print(f"Basis dirs  : {', '.join(basis_dirs)}")
    print(f"Output      : {output}")
    print("=" * 60)

    df = process_results(args.results_dir, freq=args.freq,
                         n_workers=args.workers, basis_dir=args.basis_dir)

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

    if args.csv:
        df.to_csv(args.csv, index=False)
        print(f"Saved CSV: {args.csv}")

    plot_stdp_curve(df, output, freq=args.freq)
    print("\nDone.")
    print("=" * 60)


if __name__ == "__main__":
    main()
