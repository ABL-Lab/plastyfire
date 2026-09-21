"""
Plot STDP curve for optimizer best params vs in vitro data (Markram 1997, 10 Hz).

Reads simulation_edges_23cd68d13dc3.pkl from every pair/protocol workdir,
computes EPSP ratio via basis method, groups by dt, and plots mean ± SEM.

In vitro reference (10 Hz only):
    dt = -10 ms : 0.7922 ± 0.0259  (mrk97_08)
    dt =  +5 ms : 1.2038 ± 0.0644  (mrk97_03)
    dt = +10 ms : 1.2013 ± 0.0626  (mrk97_07)

Usage:
    python plot_stdp_optimizer_best.py [--out stdp_optimizer_best.png]
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

# ---------------------------------------------------------------------------
# Paths & constants
# ---------------------------------------------------------------------------
SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(
    SCRIPT_DIR,
    "refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations",
)
BASIS_DIR   = os.path.join(SCRIPT_DIR, "basis_results_edges_mini")
PARAM_HASH  = "d157569160e0"  # gen 30 best

# Protocol directory name → dt in ms
PROTO_TO_DT = {
    "10Hz_5ms":   5,
    "10Hz_10ms":  10,
    "10Hz_30ms":  30,
    "10Hz_50ms":  50,
    "10Hz_-10ms": -10,
    "10Hz_-30ms": -30,
    "10Hz_-50ms": -50,
}

# In vitro reference (Markram 1997, 10 Hz, mrk97 protocols)
INVITRO_DT   = [-10,    5,     10]
INVITRO_MEAN = [0.7922, 1.2038, 1.2013]
INVITRO_SEM  = [0.0259, 0.0644, 0.0626]


# ---------------------------------------------------------------------------
# Worker (top-level for multiprocessing)
# ---------------------------------------------------------------------------

def _worker(args):
    """(pair, proto_dir, dt) → dict or None."""
    pair, proto_dir, dt = args
    pre_gid, post_gid = pair.split("-")
    workdir = os.path.join(RESULTS_DIR, pair, proto_dir)
    pkl     = os.path.join(workdir, f"simulation_edges_{PARAM_HASH}.pkl")

    if not os.path.exists(pkl):
        return None

    try:
        with open(pkl, "rb") as f:
            raw = pickle.load(f)
        initial_rho = list(raw["initial_rho"])
        final_rho   = list(raw["final_rho"])
        basis_df    = _load_basis(BASIS_DIR, pre_gid, post_gid)
        ep_b_m, ep_b_s = compute_epsp_from_basis(basis_df, initial_rho)
        ep_a_m, _      = compute_epsp_from_basis(basis_df, final_rho)
        if ep_b_m == 0:
            return None
        cv2   = min((ep_b_s / ep_b_m) ** 2, 0.25)
        ratio = (ep_a_m / ep_b_m) * (1.0 + cv2)
        return {"pair": pair, "dt": dt, "epsp_ratio": ratio}
    except Exception as e:
        print(f"  Error {pair} {proto_dir}: {e}")
        return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="stdp_optimizer_best.png")
    args = parser.parse_args()

    # Build job list
    jobs = []
    for pair_dir in sorted(os.scandir(RESULTS_DIR), key=lambda e: e.name):
        if not pair_dir.is_dir():
            continue
        pair = pair_dir.name
        for proto_dir, dt in PROTO_TO_DT.items():
            workdir = os.path.join(RESULTS_DIR, pair, proto_dir)
            if os.path.isfile(os.path.join(workdir, "prefire_simulation_config.json")):
                jobs.append((pair, proto_dir, dt))

    print(f"Processing {len(jobs)} pair/protocol combinations …")
    rows, missing = [], 0
    with ProcessPoolExecutor() as ex:
        for r in as_completed(ex.submit(_worker, j) for j in jobs):
            result = r.result()
            if result is None:
                missing += 1
            else:
                rows.append(result)

    df = pd.DataFrame(rows)
    print(f"Collected {len(df)} EPSP ratios  ({missing} missing pkls)")

    if df.empty:
        print("No data — run submit_optimizer_best_700.sh and wait for jobs to finish.")
        return

    summary = (
        df.groupby("dt")["epsp_ratio"]
        .agg(["mean", "sem", "count"])
        .reset_index()
        .sort_values("dt")
    )
    print("\nSummary by dt:")
    for _, row in summary.iterrows():
        print(f"  dt={int(row['dt']):+4d} ms: mean={row['mean']:.4f} ± {row['sem']:.4f}  (n={int(row['count'])})")

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(7, 5))

    ax.errorbar(
        summary["dt"], summary["mean"], yerr=summary["sem"],
        fmt="o-", color="#4c9be8", lw=2, ms=6, capsize=4,
        label="Optimizer best (BCL prefire, basis EPSP)",
    )
    ax.errorbar(
        INVITRO_DT, INVITRO_MEAN, yerr=INVITRO_SEM,
        fmt="s--", color="#ff7f2a", lw=2, ms=6, capsize=4,
        label="In vitro (Markram 1997, 10 Hz)",
    )

    ax.axhline(1.0, color="gray", ls="--", lw=1, alpha=0.5)
    ax.axvline(0.0, color="gray", ls="--", lw=1, alpha=0.5)

    ax.set_xlabel(r"$\Delta t$ (ms)", fontsize=13)
    ax.set_ylabel("Mean EPSP ratio", fontsize=13)
    ax.set_title(
        "STDP curve — Optimizer best vs In vitro\n"
        r"(10 Hz, $\gamma_d$=86.2, $\gamma_p$=185.1, optimised a-params, gen 30)",
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
