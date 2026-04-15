"""
Plot EPSP ratio vs frequency at fixed Δt = +5 ms and Δt = −10 ms,
mirroring the two-panel figure from Sjöström et al. 2001.

Usage:
    python plot_freq_curves.py --params v19
    python plot_freq_curves.py --params v19 --output freq_curves_v19.png
"""

import os
import sys
import argparse
import pickle
import glob
import multiprocessing
import concurrent.futures
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.append("/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/new_fitting")
from effcai_to_epsp import fetch_epsp_ratio
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
BASIS_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/basis_results"

# ---------------------------------------------------------------------------
# In vitro reference data — Sjöström et al. 2001 (approximate digitised values)
# ---------------------------------------------------------------------------
# Panel d: Δt = +5 ms
INVITRO_5MS = {
    "freq": [0.1,  5.0,  10.0, 20.0, 40.0, 50.0],
    "mean": [0.97, 1.02, 1.14, 1.42, 1.50, 1.52],
    "sem":  [0.05, 0.04, 0.05, 0.07, 0.06, 0.06],
}
# Panel e: Δt = −10 ms
INVITRO_NEG10MS = {
    "freq": [10.0],
    "mean": [0.79],
    "sem":  [0.07],
}
# Panel f: Δt = +10 ms
INVITRO_10MS = {
    "freq": [10.0],
    "mean": [1.20],
    "sem":  [0.06],
}


def _process_single_sim(args):
    """Worker: load one pickle at a fixed dt and compute the EPSP ratio."""
    pre_gid, post_gid, pair_name, freq_dt, pkl_path, freq, dt, params_dict = args

    try:
        with open(pkl_path, "rb") as f:
            sim_data = pickle.load(f)
    except Exception as e:
        print(f"Skipping {pair_name} {freq_dt} — corrupted pickle ({type(e).__name__}: {e}). Deleting.")
        try:
            os.remove(pkl_path)
        except OSError:
            pass
        return None

    if not sim_data or "rho_GB" not in sim_data:
        print(f"Skipping {pair_name} {freq_dt} — missing rho_GB.")
        return None

    rho_gb_raw = sim_data["rho_GB"]
    if isinstance(rho_gb_raw, dict):
        rho_trace = np.array(list(rho_gb_raw.values()))
    else:
        rho_trace = np.asarray(rho_gb_raw)
        if rho_trace.ndim == 2:
            if rho_trace.shape[0] > rho_trace.shape[1]:
                rho_trace = rho_trace.T
        else:
            rho_trace = rho_trace.reshape(1, -1)

    initial_rho = [1 if r[0]  >= 0.5 else 0 for r in rho_trace]
    final_rho   = [1 if r[-1] >= 0.5 else 0 for r in rho_trace]

    try:
        res = fetch_epsp_ratio(pre_gid, post_gid, freq, dt, initial_rho, final_rho,
                               basis_dir=BASIS_DIR, params=params_dict)
        return {"pair": pair_name, "freq": freq, "dt": dt, "epsp_ratio": res["ratio_mean"]}
    except Exception as e:
        print(f"Error processing {pair_name} {freq_dt}: {e}")
        return None


def process_results(trace_results_dir, fixed_dts, params_dict=None, n_workers=None):
    """Collect EPSP ratios for all frequencies at each fixed dt in *fixed_dts*.

    Parameters
    ----------
    trace_results_dir : str
    fixed_dts : list[float]  e.g. [5.0, -10.0]
    params_dict : dict or None
    n_workers : int or None

    Returns
    -------
    pd.DataFrame  with columns [pair, freq, dt, epsp_ratio]
    """
    # Build a set of target dt strings for fast dir-name matching
    dt_strings = {dt: f"{int(dt)}ms" if dt == int(dt) else f"{dt}ms" for dt in fixed_dts}
    # e.g. {5.0: "5ms", -10.0: "-10ms"}

    jobs = []

    for pair_dir in glob.glob(os.path.join(trace_results_dir, "*")):
        if not os.path.isdir(pair_dir) or os.path.basename(pair_dir) == "logs":
            continue
        pair_name = os.path.basename(pair_dir)
        try:
            pre_gid, post_gid = map(int, pair_name.split("-"))
        except ValueError:
            continue

        basis_path = os.path.join(BASIS_DIR, f"basis_{pre_gid}_{post_gid}.csv")
        if not os.path.exists(basis_path):
            print(f"Skipping {pair_name}, no basis file.")
            continue

        for freq_dt_dir in glob.glob(os.path.join(pair_dir, "*Hz_*")):
            freq_dt = os.path.basename(freq_dt_dir)   # e.g. "10Hz_5ms"
            try:
                freq_str, dt_str = freq_dt.split("_", 1)
                freq = float(freq_str.replace("Hz", ""))
                dt   = float(dt_str.replace("ms", ""))
            except ValueError:
                continue

            if dt not in fixed_dts:
                continue

            pkl_path = os.path.join(freq_dt_dir, "simulation_traces.pkl")
            if not os.path.exists(pkl_path):
                continue

            jobs.append((pre_gid, post_gid, pair_name, freq_dt, pkl_path, freq, dt, params_dict))

    if not jobs:
        return pd.DataFrame()

    n_workers = n_workers or multiprocessing.cpu_count()
    print(f"Processing {len(jobs)} simulations across {n_workers} cores...")

    data = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as pool:
        for result in pool.map(_process_single_sim, jobs):
            if result is not None:
                data.append(result)

    return pd.DataFrame(data)


def plot_freq_curves(df, version, output_filename):
    """Reproduce the three-panel frequency-dependence figure."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    panels = [
        (5.0,   INVITRO_5MS,      r"$\Delta t = 5$ ms",    "d"),
        (-10.0, INVITRO_NEG10MS,  r"$\Delta t = -10$ ms",  "e"),
        (10.0,  INVITRO_10MS,     r"$\Delta t = +10$ ms",  "f"),
    ]

    for ax, (dt_val, invitro, title, label) in zip(axes, panels):
        sub = df[df["dt"] == dt_val]

        # In silico
        if not sub.empty:
            summary = sub.groupby("freq")["epsp_ratio"].agg(["mean", "sem"]).reset_index()
            summary = summary.sort_values("freq")
            ax.errorbar(summary["freq"], summary["mean"], yerr=summary["sem"],
                        fmt="o-", color="#6baed6", label=f"in silico ({version})",
                        capsize=0, lw=1.5, ms=5)

        # In vitro
        ax.errorbar(invitro["freq"], invitro["mean"], yerr=invitro["sem"],
                    fmt="o-", color="#fd8d3c", label="in vitro",
                    capsize=0, lw=1.5, ms=5)

        ax.axhline(1.0, color="k", ls="--", lw=0.8, alpha=0.6)
        ax.set_xlabel("Frequency (Hz)", fontsize=11)
        ax.set_ylabel("EPSP ratio", fontsize=11)
        ax.set_title(title, fontsize=12)
        ax.text(-0.08, 1.05, label, transform=ax.transAxes,
                fontsize=13, fontweight="bold", va="top")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(frameon=False, fontsize=9, loc="upper left")

    plt.tight_layout()
    plt.savefig(output_filename, dpi=300)
    print(f"Saved {output_filename}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot EPSP ratio vs frequency at fixed Δt")
    parser.add_argument("--params",
                        choices=list(TRACE_RESULTS_DIRS.keys()), default="v19",
                        help="Parameter set to analyse (default: v19)")
    parser.add_argument("--results-dir", default=None,
                        help="Override default results directory")
    parser.add_argument("--output", default=None,
                        help="Output PNG filename")
    parser.add_argument("--workers", type=int, default=None,
                        help="Parallel worker processes (default: all cores)")
    args = parser.parse_args()

    version           = args.params
    trace_results_dir = args.results_dir or TRACE_RESULTS_DIRS[version]
    output_filename   = args.output or f"freq_curves_{version}.png"

    params_mapping = {
        "v1": DHURUVA_PARAMS, "v2": DHURUVA_PARAMS_V2, "v3": DHURUVA_PARAMS_V3,
        "v4": DHURUVA_PARAMS_V4, "v5": DHURUVA_PARAMS_V5, "v6": DHURUVA_PARAMS_V6,
        "v7": DHURUVA_PARAMS_V7, "v8": DHURUVA_PARAMS_V8, "v9": DHURUVA_PARAMS_V9,
        "v10": DHURUVA_PARAMS_V10, "v11": DHURUVA_PARAMS_V11, "v12": DHURUVA_PARAMS_V12,
        "v13": DHURUVA_PARAMS_V13, "v14": DHURUVA_PARAMS_V14, "v15": DHURUVA_PARAMS_V15,
        "v16": DHURUVA_PARAMS_V16, "v17": DHURUVA_PARAMS_V17, "v18": DHURUVA_PARAMS_V18,
        "v19": DHURUVA_PARAMS_V19, "chindemi": CHINDEMI_PARAMS,
    }
    selected_params = params_mapping.get(version, DHURUVA_PARAMS_V19)

    print("\n" + "=" * 55)
    print(f"Frequency-Curve Analysis  ({version})")
    print(f"Results dir : {trace_results_dir}")
    print("=" * 55)

    df = process_results(trace_results_dir,
                         fixed_dts=[5.0, -10.0, 10.0],
                         params_dict=selected_params,
                         n_workers=args.workers)

    if not df.empty:
        print(f"\nTotal rows : {len(df)}")
        print("\nMean EPSP ratio by (dt, freq):")
        summary = df.groupby(["dt", "freq"])["epsp_ratio"].agg(["mean", "sem", "count"])
        print(summary.to_string())
        plot_freq_curves(df, version, output_filename)
    else:
        print("No results found — make sure simulations at the target frequencies exist.")

    print("=" * 55)
