#!/usr/bin/env python3
"""
Plot STDP curve for completed fitted-params simulations.

Reads rho.h5 from bluecellulab_results_fitted/<pair>/<protocol>/bluecellulab_results/
and computes EPSP ratios via linear superposition from basis_results_edges_mini/.
Optionally overlays the chindemi baseline from bluecellulab_results/.

Usage:
    python plot_stdp_fitted.py
    python plot_stdp_fitted.py --no-baseline --output fitted_only.png
    python plot_stdp_fitted.py --fitted-dir /path/to/bluecellulab_results_fitted
"""

import argparse
import os
import sys

import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PLASTYFIRE_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(PLASTYFIRE_ROOT, "new_fitting"))
from effcai_to_epsp import fetch_epsp_ratio

DEFAULT_FITTED_DIR    = "/project/ctb-emuller/dhuruva/plastyfire/bluecellulab_results_fitted"
DEFAULT_FITTING2_DIR  = "/project/ctb-emuller/dhuruva/plastyfire/bluecellulab_results_fitting2"
DEFAULT_BASELINE_DIR  = "/project/ctb-emuller/dhuruva/plastyfire/bluecellulab_results"
DEFAULT_BASIS_DIR    = os.path.join(PLASTYFIRE_ROOT, "basis_results_edges_mini")
NODE_POP             = "S1nonbarrel_neurons"

# In-vitro reference — Markram et al. 1997 / Chindemi et al. (10 Hz)
INVITRO_DT   = [-10,    5,      10]
INVITRO_MEAN = [0.7922, 1.2038, 1.2013]
INVITRO_SEM  = [0.0259, 0.0644, 0.0626]

# Fitted params (JAX best-fit, loss = 0.000289)
FITTED_PARAMS = {
    "gamma_d": 60.6964,  "gamma_p": 178.2431,
    "a00": 0.7967, "a01": 1.1547, "a10": 1.1382, "a11": 1.9027,
    "a20": 2.2637, "a21": 7.0244, "a30": 3.5947, "a31": 6.2696,
    "tau_eff": 200.4818,
}


# ── HDF5 reader (same SONATA format as _write_rho_h5 in simulator_edges.py) ──

def _read_rho_h5(rho_h5_path):
    with h5py.File(rho_h5_path, "r") as f:
        pop = f["report"][NODE_POP]
        data = pop["data"][:]          # shape (2, n_synapses): [initial, final]
        initial_rho = data[0].tolist()
        final_rho   = data[1].tolist()
    n_0to1 = sum(1 for i, f in zip(initial_rho, final_rho) if i < 0.5 and f >= 0.5)
    n_0to0 = sum(1 for i, f in zip(initial_rho, final_rho) if i < 0.5 and f < 0.5)
    n_1to1 = sum(1 for i, f in zip(initial_rho, final_rho) if i >= 0.5 and f >= 0.5)
    n_1to0 = sum(1 for i, f in zip(initial_rho, final_rho) if i >= 0.5 and f < 0.5)
    return initial_rho, final_rho, n_0to1, n_0to0, n_1to1, n_1to0


def _parse_protocol(protocol):
    """Parse '10Hz_-10ms' → (freq=10, dt=-10)."""
    freq_part, dt_part = protocol.split("_")
    freq = int(freq_part.replace("Hz", ""))
    dt   = int(dt_part.replace("ms", ""))
    return freq, dt


def _compute_ratio(rho_h5_path, pre_gid, post_gid, freq, dt, basis_dir):
    initial_rho, final_rho, *_ = _read_rho_h5(rho_h5_path)
    initial_bin = [1 if r >= 0.5 else 0 for r in initial_rho]
    final_bin   = [1 if r >= 0.5 else 0 for r in final_rho]
    try:
        result = fetch_epsp_ratio(
            pre_gid, post_gid, freq, dt,
            initial_rho_config=initial_bin,
            final_rho_config=final_bin,
            basis_dir=basis_dir,
            n_trials=10,
            seed=42,
        )
        return result["ratio_mean"], result.get("ratio_sem", result.get("ratio_std", np.nan))
    except Exception as e:
        print(f"  WARN: fetch_epsp_ratio failed for {pre_gid}-{post_gid} {freq}Hz_{dt}ms: {e}")
        return np.nan, np.nan


def collect_results(results_dir, basis_dir, label=""):
    """
    Scan results_dir for completed rho.h5 files and compute EPSP ratios.
    Expected structure: <results_dir>/<pair>/<protocol>/bluecellulab_results/rho.h5
    Returns a DataFrame with columns: pair, pre_gid, post_gid, freq, dt, ratio_mean, ratio_sem
    """
    rows = []
    rho_files = []
    for pair_dir in sorted(os.listdir(results_dir)):
        pair_path = os.path.join(results_dir, pair_dir)
        if not os.path.isdir(pair_path):
            continue
        try:
            pre_gid, post_gid = [int(x) for x in pair_dir.split("-")]
        except ValueError:
            continue
        for proto in sorted(os.listdir(pair_path)):
            rho_path = os.path.join(pair_path, proto, "bluecellulab_results", "rho.h5")
            if os.path.isfile(rho_path):
                rho_files.append((rho_path, pair_dir, pre_gid, post_gid, proto))

    print(f"[{label}] Found {len(rho_files)} completed rho.h5 files")

    for rho_path, pair_dir, pre_gid, post_gid, proto in rho_files:
        try:
            freq, dt = _parse_protocol(proto)
        except Exception:
            continue
        basis_csv = os.path.join(basis_dir, f"basis_{pre_gid}_{post_gid}.csv")
        if not os.path.isfile(basis_csv):
            continue
        ratio_mean, ratio_sem = _compute_ratio(rho_path, pre_gid, post_gid, freq, dt, basis_dir)
        rows.append(dict(pair=pair_dir, pre_gid=pre_gid, post_gid=post_gid,
                         freq=freq, dt=dt, ratio_mean=ratio_mean, ratio_sem=ratio_sem))

    return pd.DataFrame(rows)


def collect_results_workdirs(refitting_results_dir, basis_dir, label=""):
    """
    Collect from the deep refitting_results workdir tree.
    Expected: <refitting_results_dir>/fitting/*/seed*/*_STDP/simulations/<pair>/<proto>/
                bluecellulab_results/rho.h5
    """
    import glob as _glob
    pattern = os.path.join(
        refitting_results_dir, "fitting", "*", "seed*", "*_STDP",
        "simulations", "*-*", "*Hz_*ms", "bluecellulab_results", "rho.h5",
    )
    rho_files = sorted(_glob.glob(pattern))
    print(f"[{label}] Found {len(rho_files)} completed rho.h5 files in workdirs")

    rows = []
    for rho_path in rho_files:
        parts = rho_path.split(os.sep)
        proto    = parts[-3]   # e.g. 10Hz_10ms
        pair_dir = parts[-4]   # e.g. 180164-197248
        try:
            pre_gid, post_gid = [int(x) for x in pair_dir.split("-")]
            freq, dt = _parse_protocol(proto)
        except Exception:
            continue
        basis_csv = os.path.join(basis_dir, f"basis_{pre_gid}_{post_gid}.csv")
        if not os.path.isfile(basis_csv):
            continue
        ratio_mean, ratio_sem = _compute_ratio(rho_path, pre_gid, post_gid, freq, dt, basis_dir)
        rows.append(dict(pair=pair_dir, pre_gid=pre_gid, post_gid=post_gid,
                         freq=freq, dt=dt, ratio_mean=ratio_mean, ratio_sem=ratio_sem))

    return pd.DataFrame(rows)


def _aggregate(df):
    """Group by dt, compute mean ± SEM of ratio_mean across pairs."""
    out = (
        df.dropna(subset=["ratio_mean"])
          .groupby("dt")["ratio_mean"]
          .agg(mean="mean", sem=lambda x: x.std(ddof=1) / np.sqrt(len(x)), count="count")
          .reset_index()
          .sort_values("dt")
    )
    return out


def plot(fitted_df, fitting2_df, baseline_df, output_path):
    fig, ax = plt.subplots(figsize=(7, 5))

    # In-vitro reference
    ax.errorbar(INVITRO_DT, INVITRO_MEAN, yerr=INVITRO_SEM,
                fmt="ks", markersize=7, capsize=4, label="In vitro (Markram 1997)", zorder=5)

    # Baseline chindemi (if provided)
    if baseline_df is not None and len(baseline_df):
        agg = _aggregate(baseline_df)
        ax.errorbar(agg["dt"], agg["mean"], yerr=agg["sem"],
                    fmt="o--", color="#66b3e6", capsize=3,
                    label=f"Chindemi params (n={len(baseline_df['pair'].unique())} pairs)")

    # Fitted params (round 1)
    if fitted_df is not None and len(fitted_df):
        agg = _aggregate(fitted_df)
        ax.errorbar(agg["dt"], agg["mean"], yerr=agg["sem"],
                    fmt="o-", color="#e07b39", capsize=3, linewidth=2,
                    label=f"Fitted (γ_d=60.7, γ_p=178.2, τ=200.5 ms) n={len(fitted_df['pair'].unique())}")

    # Fitting2 params (round 2)
    if fitting2_df is not None and len(fitting2_df):
        agg = _aggregate(fitting2_df)
        ax.errorbar(agg["dt"], agg["mean"], yerr=agg["sem"],
                    fmt="s-", color="#6abf69", capsize=3, linewidth=2,
                    label=f"Fitting2 (γ_d=58.5, γ_p=246.3, τ=401.8 ms) n={len(fitting2_df['pair'].unique())}")

    ax.axhline(1.0, color="gray", linestyle=":", linewidth=0.8)
    ax.set_xlabel("Pre→Post spike delay (ms)", fontsize=12)
    ax.set_ylabel("EPSP ratio (after / before)", fontsize=12)
    ax.set_title("STDP curve — Fitted vs Chindemi params (10 Hz)", fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    print(f"Saved → {output_path}")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fitted-dir", default=DEFAULT_FITTED_DIR,
                        help=f"bluecellulab_results_fitted directory (default: {DEFAULT_FITTED_DIR})")
    parser.add_argument("--fitting2-dir", default=DEFAULT_FITTING2_DIR,
                        help=f"bluecellulab_results_fitting2 directory (default: {DEFAULT_FITTING2_DIR})")
    parser.add_argument("--baseline-dir", default=DEFAULT_BASELINE_DIR,
                        help=f"bluecellulab_results (chindemi) directory (default: {DEFAULT_BASELINE_DIR})")
    parser.add_argument("--baseline-workdirs",
                        default=os.path.join(PLASTYFIRE_ROOT, "refitting_results"),
                        help="refitting_results root for chindemi baseline (deep workdir scan)")
    parser.add_argument("--basis-dir", default=DEFAULT_BASIS_DIR,
                        help=f"basis CSV directory (default: {DEFAULT_BASIS_DIR})")
    parser.add_argument("--no-baseline", action="store_true",
                        help="Skip chindemi overlay")
    parser.add_argument("--output", default="stdp_fitted_vs_chindemi.png",
                        help="Output PNG path (default: stdp_fitted_vs_chindemi.png)")
    args = parser.parse_args()

    print(f"Fitted dir  : {args.fitted_dir}")
    print(f"Basis dir   : {args.basis_dir}")

    fitted_df = collect_results(args.fitted_dir, args.basis_dir, label="fitted")

    fitting2_df = None
    if os.path.isdir(args.fitting2_dir):
        fitting2_df = collect_results(args.fitting2_dir, args.basis_dir, label="fitting2")
    else:
        print(f"fitting2 dir not found yet: {args.fitting2_dir}")

    baseline_df = None
    if not args.no_baseline:
        print(f"Baseline workdirs: {args.baseline_workdirs}")
        baseline_df = collect_results_workdirs(args.baseline_workdirs, args.basis_dir, label="chindemi")

    print("\nFitted results summary by dt:")
    print(_aggregate(fitted_df).to_string(index=False) if fitted_df is not None and len(fitted_df) else "  (none)")

    if fitting2_df is not None and len(fitting2_df):
        print("\nFitting2 results summary by dt:")
        print(_aggregate(fitting2_df).to_string(index=False))

    if baseline_df is not None and len(baseline_df):
        print("\nChindemi results summary by dt:")
        print(_aggregate(baseline_df).to_string(index=False))

    plot(fitted_df, fitting2_df, baseline_df, args.output)


if __name__ == "__main__":
    main()
