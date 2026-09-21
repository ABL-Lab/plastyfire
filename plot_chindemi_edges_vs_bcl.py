#!/usr/bin/env python3
"""
Validate the optimizer against bluecellulab-edges for Chindemi params.

Panel 1 (STDP curve): in-vitro vs bluecellulab-edges vs optimizer
  - Bluecellulab-edges : bluecellulab_results/rho.h5, binary rho>=0.5, delta_method
                         (same as plot_compare_ndamus_bluecellulab.py — validated reference)
  - Optimizer          : simulation_edges_{hash}.pkl, binary rho>=0.5, delta_method
                         (same method as BCL — matches optimizer fitness signal)

Panel 2 (rho subplot): continuous rho for one pair / one synapse
  - Bluecellulab-edges : simulation_edges.pkl  initial_rho → final_rho step
  - Optimizer          : simulation_edges_{hash}.pkl  initial_rho → final_rho step
"""

import hashlib
import os
import pickle
import sys
import glob
import multiprocessing
import concurrent.futures

import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
PLASTYFIRE_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(PLASTYFIRE_ROOT, "new_fitting"))
sys.path.insert(0, PLASTYFIRE_ROOT)

from effcai_to_epsp import fetch_epsp_ratio
from plastyfire.evaluator_edges import compute_epsp_from_basis, _load_basis

RESULTS_DIR = os.path.join(PLASTYFIRE_ROOT, "refitting_results")
BASIS_DIR   = os.path.join(PLASTYFIRE_ROOT, "basis_results_edges_mini")
OUTPUT      = os.path.join(PLASTYFIRE_ROOT, "stdp_compare_10Hz.png")

NODE_POP = "S1nonbarrel_neurons"

CHINDEMI_PARAMS = [101.5, 216.2, 1.002, 1.954, 1.159, 2.483, 1.127, 2.456, 5.236, 1.782]
CHINDEMI_HASH   = hashlib.md5(str(CHINDEMI_PARAMS).encode()).hexdigest()[:12]

# Example pair/dt for the rho subplot
EXAMPLE_PAIR = "181002-188173"
EXAMPLE_DT   = "10Hz_10ms"

# In-vitro reference (Markram et al. 1997, 10 Hz)
INVITRO_DT   = [-10,    5,      10]
INVITRO_MEAN = [0.7922, 1.2038, 1.2013]
INVITRO_SEM  = [0.0259, 0.0644, 0.0626]


# ---------------------------------------------------------------------------
# BCL EPSP: bluecellulab_results/rho.h5, binary rho, delta_method
# (identical to plot_compare_ndamus_bluecellulab.py)
# ---------------------------------------------------------------------------

def _read_rho_h5_binary(rho_h5_path):
    with h5py.File(rho_h5_path, "r") as f:
        data = f[f"report/{NODE_POP}/data"][()]
    return ([1 if v >= 0.5 else 0 for v in data[0]],
            [1 if v >= 0.5 else 0 for v in data[-1]])


def _process_bcl(args):
    sim_dir, pre_gid, post_gid, dt_ms, basis_path = args
    rho_h5 = os.path.join(sim_dir, "bluecellulab_results", "rho.h5")
    if not os.path.exists(rho_h5):
        return None
    try:
        initial_rho, final_rho = _read_rho_h5_binary(rho_h5)
        res = fetch_epsp_ratio(pre_gid, post_gid, 10.0, dt_ms,
                               initial_rho, final_rho,
                               basis_dir=os.path.dirname(basis_path),
                               ratio_method="delta_method")
        return {"pre_gid": pre_gid, "post_gid": post_gid,
                "dt": dt_ms, "epsp_ratio": res["ratio_mean"]}
    except Exception as e:
        print(f"  BCL error {pre_gid}-{post_gid} dt={dt_ms}: {e}")
        return None


def collect_bcl(results_dir=RESULTS_DIR, basis_dir=BASIS_DIR, freq=10, n_workers=None):
    sim_dirs = glob.glob(os.path.join(results_dir, "fitting", "*", "seed*",
                                      "*_STDP", "simulations", "*-*", f"{int(freq)}Hz_*"))
    jobs = []
    for sd in sorted(sim_dirs):
        if not os.path.exists(os.path.join(sd, "bluecellulab_results", "rho.h5")):
            continue
        pair   = os.path.basename(os.path.dirname(sd))
        dt_str = os.path.basename(sd).split("_", 1)[1].replace("ms", "")
        try:
            pre_gid, post_gid = [int(x) for x in pair.split("-")]
            dt_ms = float(dt_str)
        except ValueError:
            continue
        bp = os.path.join(basis_dir, f"basis_{pre_gid}_{post_gid}.csv")
        if os.path.exists(bp):
            jobs.append((sd, pre_gid, post_gid, dt_ms, bp))

    n_workers = n_workers or multiprocessing.cpu_count()
    print(f"  [BCL] {len(jobs)} sims, {n_workers} workers...")
    data = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as pool:
        for r in pool.map(_process_bcl, jobs):
            if r is not None:
                data.append(r)
    return pd.DataFrame(data)


# ---------------------------------------------------------------------------
# Optimizer EPSP: simulation_edges_{hash}.pkl, binary rho>=0.5, delta_method
# (matches BCL and the fixed optimizer fitness signal)
# ---------------------------------------------------------------------------

def _process_optimizer(args):
    sim_dir, pre_gid, post_gid, dt_ms, basis_path, param_hash = args
    pkl = os.path.join(sim_dir, f"simulation_edges_{param_hash}.pkl")
    if not os.path.exists(pkl):
        return None
    try:
        with open(pkl, "rb") as f:
            raw = pickle.load(f)
        basis_df    = _load_basis(basis_path.replace(f"basis_{pre_gid}_{post_gid}.csv", ""),
                                  str(pre_gid), str(post_gid))
        epsp_before_mean, epsp_before_std = compute_epsp_from_basis(basis_df, list(raw["initial_rho"]))
        epsp_after_mean,  _               = compute_epsp_from_basis(basis_df, list(raw["final_rho"]))
        if epsp_before_mean == 0:
            return None
        cv2 = min((epsp_before_std / epsp_before_mean) ** 2, 0.25) if epsp_before_mean > 0 else 0.0
        epsp_ratio = (epsp_after_mean / epsp_before_mean) * (1.0 + cv2)
        return {"pre_gid": pre_gid, "post_gid": post_gid,
                "dt": dt_ms, "epsp_ratio": epsp_ratio}
    except Exception as e:
        print(f"  Optimizer error {pre_gid}-{post_gid} dt={dt_ms}: {e}")
        return None


def collect_optimizer(param_hash=CHINDEMI_HASH, results_dir=RESULTS_DIR,
                      basis_dir=BASIS_DIR, freq=10, n_workers=None):
    sim_dirs = glob.glob(os.path.join(results_dir, "fitting", "*", "seed*",
                                      "*_STDP", "simulations", "*-*", f"{int(freq)}Hz_*"))
    jobs = []
    for sd in sorted(sim_dirs):
        if not os.path.exists(os.path.join(sd, f"simulation_edges_{param_hash}.pkl")):
            continue
        pair   = os.path.basename(os.path.dirname(sd))
        dt_str = os.path.basename(sd).split("_", 1)[1].replace("ms", "")
        try:
            pre_gid, post_gid = [int(x) for x in pair.split("-")]
            dt_ms = float(dt_str)
        except ValueError:
            continue
        bp = os.path.join(basis_dir, f"basis_{pre_gid}_{post_gid}.csv")
        if os.path.exists(bp):
            jobs.append((sd, pre_gid, post_gid, dt_ms, bp, param_hash))

    n_workers = n_workers or multiprocessing.cpu_count()
    print(f"  [Optimizer] {len(jobs)} sims, {n_workers} workers...")
    data = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as pool:
        for r in pool.map(_process_optimizer, jobs):
            if r is not None:
                data.append(r)
    return pd.DataFrame(data)


# ---------------------------------------------------------------------------
# Rho subplot: step function from initial→final at last prespike, one synapse
# ---------------------------------------------------------------------------

def load_rho_steps(sim_dir, param_hash, synapse_idx=0):
    """
    Load rho step data for one synapse from both BCL and optimizer pkl files.
    Returns dict with keys 'bcl' and 'optimizer', each containing:
      t_start, t_end, t_prespike, initial_rho, final_rho
    """
    out = {}

    # BCL: simulation_edges.pkl (validated run)
    bcl_pkl = os.path.join(sim_dir, "simulation_edges.pkl")
    if os.path.exists(bcl_pkl):
        with open(bcl_pkl, "rb") as f:
            d = pickle.load(f)
        t_arr = np.asarray(d["t"])
        prespikes = np.asarray(d["prespikes"])
        out["bcl"] = {
            "t_start":    float(t_arr[0]),
            "t_end":      float(t_arr[-1]),
            "t_prespike": float(prespikes.max()) if len(prespikes) else float(t_arr[-1]),
            "initial_rho": float(np.asarray(d["initial_rho"])[synapse_idx]),
            "final_rho":   float(np.asarray(d["final_rho"])[synapse_idx]),
        }

    # Optimizer: simulation_edges_{hash}.pkl
    opt_pkl = os.path.join(sim_dir, f"simulation_edges_{param_hash}.pkl")
    if os.path.exists(opt_pkl):
        with open(opt_pkl, "rb") as f:
            d = pickle.load(f)
        t_arr = np.asarray(d["t"])
        prespikes = np.asarray(d["prespikes"])
        out["optimizer"] = {
            "t_start":    float(t_arr[0]),
            "t_end":      float(t_arr[-1]),
            "t_prespike": float(prespikes.max()) if len(prespikes) else float(t_arr[-1]),
            "initial_rho": float(np.asarray(d["initial_rho"])[synapse_idx]),
            "final_rho":   float(np.asarray(d["final_rho"])[synapse_idx]),
        }

    return out


def _rho_step(ax, rho_data, color, label):
    """Draw rho as a step: initial → (holds) → final_rho at last prespike."""
    t0, t1 = rho_data["t_start"], rho_data["t_end"]
    tp      = rho_data["t_prespike"]
    r0, r1  = rho_data["initial_rho"], rho_data["final_rho"]
    t = np.array([t0, tp, tp, t1])
    r = np.array([r0, r0, r1, r1])
    ax.plot(t, r, color=color, lw=1.8, label=label)
    ax.scatter([t0, t1], [r0, r1], color=color, s=40, zorder=5)


# ---------------------------------------------------------------------------
# Main plot
# ---------------------------------------------------------------------------

def plot(df_bcl, df_opt, rho_steps, output=OUTPUT):
    def summ(df):
        return (df.groupby("dt")["epsp_ratio"]
                  .agg(["mean", "sem"])
                  .reset_index().sort_values("dt"))

    s_bcl = summ(df_bcl)
    s_opt = summ(df_opt)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(
        f"Chindemi params — Optimizer validation  (hash: {CHINDEMI_HASH})",
        fontsize=12
    )

    # ── Panel 1: STDP curve ──────────────────────────────────────────────────
    ax = axes[0]
    ax.errorbar(s_bcl["dt"], s_bcl["mean"], yerr=s_bcl["sem"],
                fmt="o-", color="#66c97f", capsize=4,
                label="Bluecellulab-edges (binary rho, delta_method)", zorder=3)
    ax.errorbar(s_opt["dt"], s_opt["mean"], yerr=s_opt["sem"],
                fmt="s--", color="#e06c75", capsize=4, markersize=7,
                label="Optimizer (continuous rho, basis ratio)", zorder=4)
    ax.errorbar(INVITRO_DT, INVITRO_MEAN, yerr=INVITRO_SEM,
                fmt="o-", color="#e5c07b", capsize=4,
                label="In vitro (Markram 1997)", zorder=2)
    ax.axhline(1.0, color="k", lw=0.8, ls="--", alpha=0.4)
    ax.axvline(0.0, color="k", lw=0.8, ls="--", alpha=0.4)
    ax.set_xlabel(r"$\Delta t$ (ms)", fontsize=11)
    ax.set_ylabel("EPSP ratio", fontsize=11)
    ax.set_title("STDP curve  (10 Hz)", fontsize=11)
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # ── Panel 2: rho step for one synapse ───────────────────────────────────
    ax = axes[1]
    if "bcl" in rho_steps:
        _rho_step(ax, rho_steps["bcl"],       "#66c97f", "Bluecellulab-edges")
    if "optimizer" in rho_steps:
        _rho_step(ax, rho_steps["optimizer"], "#e06c75", "Optimizer")
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("t (ms)", fontsize=11)
    ax.set_ylabel("rho (synapse 0)", fontsize=11)
    ax.set_title(f"Rho — {EXAMPLE_PAIR}  {EXAMPLE_DT}  synapse 0", fontsize=11)
    ax.legend(frameon=False, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(output, dpi=200, bbox_inches="tight")
    print(f"\nSaved: {output}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print(f"Chindemi param hash: {CHINDEMI_HASH}")

    print("\nCollecting BCL results (bluecellulab_results/rho.h5, binary rho, delta_method)...")
    df_bcl = collect_bcl()
    print(f"  → {len(df_bcl)} rows, {df_bcl['dt'].nunique()} dt values, "
          f"{df_bcl[['pre_gid','post_gid']].drop_duplicates().shape[0]} pairs")

    print("\nCollecting optimizer results (simulation_edges_{hash}.pkl, continuous rho, basis)...")
    df_opt = collect_optimizer()
    print(f"  → {len(df_opt)} rows, {df_opt['dt'].nunique()} dt values, "
          f"{df_opt[['pre_gid','post_gid']].drop_duplicates().shape[0]} pairs")

    def print_summary(df, label):
        s = (df.groupby("dt")["epsp_ratio"]
               .agg(["mean", "sem", "count"])
               .reset_index().sort_values("dt"))
        print(f"\n{label}:")
        print(s.to_string(index=False))

    print_summary(df_bcl, "Bluecellulab-edges (all dt)")
    print_summary(df_opt, "Optimizer (Chindemi, ±10ms)")

    # Load rho steps for example pair
    example_sim_dir = os.path.join(
        RESULTS_DIR, "fitting", "n100", "seed19091997",
        "L5TTPC_L5TTPC_STDP", "simulations", EXAMPLE_PAIR, EXAMPLE_DT
    )
    rho_steps = load_rho_steps(example_sim_dir, CHINDEMI_HASH, synapse_idx=0)
    print(f"\nRho steps for {EXAMPLE_PAIR} {EXAMPLE_DT} synapse 0:")
    for src, v in rho_steps.items():
        print(f"  {src}: initial={v['initial_rho']:.4f}  final={v['final_rho']:.4f}  "
              f"last_prespike={v['t_prespike']:.1f} ms")

    plot(df_bcl, df_opt, rho_steps)


if __name__ == "__main__":
    main()
