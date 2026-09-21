"""
Plot STDP curve for gen-30 best params: original HOC vs new ion channels vs in vitro.

Reads simulation_edges_{hash}.pkl from every pair/protocol workdir,
computes EPSP ratio via basis method, groups by dt, plots mean ± SEM.

Usage:
    python plot_stdp_ion_channels.py [--out stdp_ion_channels.png]
"""

import argparse
import os
import pickle
from concurrent.futures import ProcessPoolExecutor, as_completed

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plastyfire.evaluator_edges import compute_epsp_from_basis, _load_basis

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(
    SCRIPT_DIR,
    "refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations",
)

# Two basis dirs — one per HOC variant
BASIS_DIR_ORIG = os.path.join(SCRIPT_DIR, "basis_results_edges_mini")
BASIS_DIR_IC   = os.path.join(SCRIPT_DIR, "basis_results_edges_ion_channels")

HASH_ORIG = "d157569160e0"       # gen-30 best, original HOC
HASH_IC   = "d157569160e0_ic"    # same params, new ion channels

PROTO_TO_DT = {
    "10Hz_5ms":   5,
    "10Hz_10ms":  10,
    "10Hz_30ms":  30,
    "10Hz_50ms":  50,
    "10Hz_-10ms": -10,
    "10Hz_-30ms": -30,
    "10Hz_-50ms": -50,
}

INVITRO_DT   = [-10,    5,     10]
INVITRO_MEAN = [0.7922, 1.2038, 1.2013]
INVITRO_SEM  = [0.0259, 0.0644, 0.0626]


def _worker(args):
    pair, proto_dir, dt, param_hash, basis_dir = args
    pre_gid, post_gid = pair.split("-")
    workdir = os.path.join(RESULTS_DIR, pair, proto_dir)
    pkl     = os.path.join(workdir, f"simulation_edges_{param_hash}.pkl")

    if not os.path.exists(pkl):
        return None

    try:
        with open(pkl, "rb") as f:
            raw = pickle.load(f)
        initial_rho = list(raw["initial_rho"])
        final_rho   = list(raw["final_rho"])
        basis_df    = _load_basis(basis_dir, pre_gid, post_gid)
        ep_b_m, ep_b_s = compute_epsp_from_basis(basis_df, initial_rho)
        ep_a_m, _      = compute_epsp_from_basis(basis_df, final_rho)
        if ep_b_m == 0:
            return None
        cv2   = min((ep_b_s / ep_b_m) ** 2, 0.25)
        ratio = (ep_a_m / ep_b_m) * (1.0 + cv2)
        return {"pair": pair, "dt": dt, "epsp_ratio": ratio, "variant": param_hash}
    except Exception as e:
        print(f"  Error {pair} {proto_dir} [{param_hash}]: {e}")
        return None


def collect(param_hash, basis_dir):
    jobs = []
    for pair_dir in sorted(os.scandir(RESULTS_DIR), key=lambda e: e.name):
        if not pair_dir.is_dir():
            continue
        pair = pair_dir.name
        for proto_dir, dt in PROTO_TO_DT.items():
            workdir = os.path.join(RESULTS_DIR, pair, proto_dir)
            if os.path.isfile(os.path.join(workdir, "prefire_simulation_config.json")):
                jobs.append((pair, proto_dir, dt, param_hash, basis_dir))

    rows, missing = [], 0
    with ProcessPoolExecutor() as ex:
        for r in as_completed(ex.submit(_worker, j) for j in jobs):
            result = r.result()
            if result is None:
                missing += 1
            else:
                rows.append(result)
    return pd.DataFrame(rows), missing


def summarise(df):
    return (
        df.groupby("dt")["epsp_ratio"]
        .agg(["mean", "sem", "count"])
        .reset_index()
        .sort_values("dt")
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="stdp_ion_channels.png")
    args = parser.parse_args()

    print("Collecting original HOC results …")
    df_orig, miss_orig = collect(HASH_ORIG, BASIS_DIR_ORIG)
    print(f"  {len(df_orig)} ratios  ({miss_orig} missing)")

    print("Collecting new ion channel results …")
    df_ic, miss_ic = collect(HASH_IC, BASIS_DIR_IC)
    print(f"  {len(df_ic)} ratios  ({miss_ic} missing)")

    if df_orig.empty and df_ic.empty:
        print("No data at all — run submit scripts first.")
        return

    fig, ax = plt.subplots(figsize=(7, 5))

    if not df_orig.empty:
        s = summarise(df_orig)
        print("\nOriginal HOC by dt:")
        for _, r in s.iterrows():
            print(f"  dt={int(r['dt']):+4d} ms: {r['mean']:.4f} ± {r['sem']:.4f}  (n={int(r['count'])})")
        ax.errorbar(s["dt"], s["mean"], yerr=s["sem"],
                    fmt="o-", color="#4c9be8", lw=2, ms=6, capsize=4,
                    label="Gen-30 best — original HOC")

    if not df_ic.empty:
        s = summarise(df_ic)
        print("\nNew ion channels by dt:")
        for _, r in s.iterrows():
            print(f"  dt={int(r['dt']):+4d} ms: {r['mean']:.4f} ± {r['sem']:.4f}  (n={int(r['count'])})")
        ax.errorbar(s["dt"], s["mean"], yerr=s["sem"],
                    fmt="s-", color="#2ca02c", lw=2, ms=6, capsize=4,
                    label="Gen-30 best — new ion channels")

    ax.errorbar(INVITRO_DT, INVITRO_MEAN, yerr=INVITRO_SEM,
                fmt="^--", color="#ff7f2a", lw=2, ms=6, capsize=4,
                label="In vitro (Markram 1997, 10 Hz)")

    ax.axhline(1.0, color="gray", ls="--", lw=1, alpha=0.5)
    ax.axvline(0.0, color="gray", ls="--", lw=1, alpha=0.5)
    ax.set_xlabel(r"$\Delta t$ (ms)", fontsize=13)
    ax.set_ylabel("Mean EPSP ratio", fontsize=13)
    ax.set_title(
        "STDP curve — original HOC vs new ion channels\n"
        r"(gen-30 best params, $\gamma_d$=86.2, $\gamma_p$=185.1, 10 Hz)",
        fontsize=12,
    )
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, alpha=0.2, axis="y")

    plt.tight_layout()
    out_path = os.path.join(SCRIPT_DIR, args.out)
    plt.savefig(out_path, dpi=300)
    print(f"\nSaved → {out_path}")


if __name__ == "__main__":
    main()
