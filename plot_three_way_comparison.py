#!/usr/bin/env python3
"""
Three-way STDP comparison for Chindemi params:
  1. BCL full       -- bluecellulab_results/rho.h5  (full 520s sim, basis EPSP)
  2. BCL prefire    -- bluecellulab_results_prefire/rho.h5  (42s prefire, basis EPSP)
  3. Optimizer      -- simulation_edges_{chindemi_hash}.pkl  (42s prefire, basis EPSP)

All three use the same basis method (compute_epsp_from_basis + delta-method CV² correction)
so any difference is due to the simulation, not the EPSP computation.

Usage:
    python plot_three_way_comparison.py [--workers N] [--output FILE]
"""

import os
import sys
import glob
import pickle
import hashlib
import argparse
import concurrent.futures

import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PLASTYFIRE_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PLASTYFIRE_ROOT)

from plastyfire.evaluator_edges import compute_epsp_from_basis, _load_basis

RESULTS_DIR = os.path.join(PLASTYFIRE_ROOT, "refitting_results")
BASIS_DIR   = os.path.join(PLASTYFIRE_ROOT, "basis_results_edges_mini")
NODE_POP    = "S1nonbarrel_neurons"

# Chindemi params list (same order as in evaluator_edges EDGES_FIT_PARAM_NAMES,
# but chindemi has no a-params — hash computed from full param dict for optimizer pkl)
CHINDEMI_PARAMS = [101.5, 216.2, 1.002, 1.954, 1.159, 2.483, 1.127, 2.456, 5.236, 1.782]
CHINDEMI_HASH   = hashlib.md5(str(CHINDEMI_PARAMS).encode()).hexdigest()[:12]

# Hash for pairrunner_edges_fit.py run (WITH a-params, theta computed from c_pre/c_post).
# Auto-computed by pairrunner_edges_fit.py from sorted(fit_params.items()) incl. tau.
CHINDEMI_FIT_PARAMS = {
    "gamma_d_GB_GluSynapse": 101.5, "gamma_p_GB_GluSynapse": 216.2,
    "a00": 1.002, "a01": 1.954, "a10": 1.159, "a11": 2.483,
    "a20": 1.127, "a21": 2.456, "a30": 5.236, "a31": 1.782,
    "tau_effca_GB_GluSynapse": 278.3177658387,
}
CHINDEMI_FIT_HASH = hashlib.md5(str(sorted(CHINDEMI_FIT_PARAMS.items())).encode()).hexdigest()[:12]

# In-vitro reference (Markram et al. 1997, 10 Hz)
INVITRO_DT   = [-10,    5,    10]
INVITRO_MEAN = [0.7922, 1.2038, 1.2013]
INVITRO_SEM  = [0.0259, 0.0644, 0.0626]

TARGET_DT = {10.0: (1.2013, 0.0626), -10.0: (0.7922, 0.0259)}


# ---------------------------------------------------------------------------
# Shared EPSP computation from rho.h5
# ---------------------------------------------------------------------------

def _epsp_from_rho_h5(rho_h5, pre_gid, post_gid, sort_by_element_id=False):
    """Read rho.h5, compute EPSP ratio via basis method. Returns ratio or None."""
    with h5py.File(rho_h5, "r") as f:
        data = f[f"report/{NODE_POP}/data"][()]
        if sort_by_element_id:
            eids = f[f"report/{NODE_POP}/mapping/element_ids"][()]
            sort_idx = np.argsort(eids)
            data = data[:, sort_idx]
    initial_rho = list(data[0])
    final_rho   = list(data[-1])
    try:
        basis_df = _load_basis(BASIS_DIR, str(pre_gid), str(post_gid))
    except Exception:
        return None
    ep_b_m, ep_b_s = compute_epsp_from_basis(basis_df, initial_rho)
    ep_a_m, _      = compute_epsp_from_basis(basis_df, final_rho)
    if ep_b_m == 0:
        return None
    cv2 = min((ep_b_s / ep_b_m) ** 2, 0.25)
    return (ep_a_m / ep_b_m) * (1.0 + cv2)


def _epsp_from_pkl(pkl_path, pre_gid, post_gid):
    """Read simulation_edges_{hash}.pkl, compute EPSP ratio via basis method."""
    with open(pkl_path, "rb") as f:
        raw = pickle.load(f)
    initial_rho = list(raw["initial_rho"])
    final_rho   = list(raw["final_rho"])
    try:
        basis_df = _load_basis(BASIS_DIR, str(pre_gid), str(post_gid))
    except Exception:
        return None
    ep_b_m, ep_b_s = compute_epsp_from_basis(basis_df, initial_rho)
    ep_a_m, _      = compute_epsp_from_basis(basis_df, final_rho)
    if ep_b_m == 0:
        return None
    cv2 = min((ep_b_s / ep_b_m) ** 2, 0.25)
    return (ep_a_m / ep_b_m) * (1.0 + cv2)


# ---------------------------------------------------------------------------
# Workers (top-level for ProcessPoolExecutor pickling)
# ---------------------------------------------------------------------------

def _worker_bcl_full(args):
    sd, pre_gid, post_gid, dt_ms = args
    # out/rho.h5 = Neurodamus (pydamus) chindemi run, same 42s induction protocol
    # element_ids are local synapse IDs written in non-sequential order → must sort
    rho_h5 = os.path.join(sd, "out", "rho.h5")
    if not os.path.exists(rho_h5):
        return None
    try:
        ratio = _epsp_from_rho_h5(rho_h5, pre_gid, post_gid, sort_by_element_id=True)
        if ratio is None:
            return None
        return {"pre_gid": pre_gid, "post_gid": post_gid, "dt": dt_ms, "epsp_ratio": ratio}
    except Exception as e:
        print(f"  BCL full error {pre_gid}-{post_gid} dt={dt_ms}: {e}")
        return None


def _worker_bcl_prefire(args):
    sd, pre_gid, post_gid, dt_ms = args
    rho_h5 = os.path.join(sd, "bluecellulab_results_edges_chindemi_params_prefire", "rho.h5")
    if not os.path.exists(rho_h5):
        return None
    try:
        ratio = _epsp_from_rho_h5(rho_h5, pre_gid, post_gid)
        if ratio is None:
            return None
        return {"pre_gid": pre_gid, "post_gid": post_gid, "dt": dt_ms, "epsp_ratio": ratio}
    except Exception as e:
        print(f"  BCL prefire error {pre_gid}-{post_gid} dt={dt_ms}: {e}")
        return None


def _worker_optimizer(args):
    sd, pre_gid, post_gid, dt_ms, param_hash = args
    pkl = os.path.join(sd, f"simulation_edges_{param_hash}.pkl")
    if not os.path.exists(pkl):
        return None
    try:
        ratio = _epsp_from_pkl(pkl, pre_gid, post_gid)
        if ratio is None:
            return None
        return {"pre_gid": pre_gid, "post_gid": post_gid, "dt": dt_ms, "epsp_ratio": ratio}
    except Exception as e:
        print(f"  Optimizer error {pre_gid}-{post_gid} dt={dt_ms}: {e}")
        return None


def _worker_optimizer_fit(args):
    """Fourth leg: pairrunner_edges_fit.py with a-params (theta from c_pre/c_post × a-params)."""
    sd, pre_gid, post_gid, dt_ms, param_hash = args
    pkl = os.path.join(sd, f"simulation_edges_{param_hash}.pkl")
    if not os.path.exists(pkl):
        return None
    try:
        ratio = _epsp_from_pkl(pkl, pre_gid, post_gid)
        if ratio is None:
            return None
        return {"pre_gid": pre_gid, "post_gid": post_gid, "dt": dt_ms, "epsp_ratio": ratio}
    except Exception as e:
        print(f"  Optimizer-fit error {pre_gid}-{post_gid} dt={dt_ms}: {e}")
        return None


# ---------------------------------------------------------------------------
# Collection
# ---------------------------------------------------------------------------

def _find_sim_dirs(results_dir=RESULTS_DIR, freq=10):
    pattern = os.path.join(results_dir, "fitting", "*", "seed*",
                           "*_STDP", "simulations", "*-*", f"{int(freq)}Hz_*")
    rows = []
    for sd in sorted(glob.glob(pattern)):
        pair   = os.path.basename(os.path.dirname(sd))
        dt_str = os.path.basename(sd).split("_", 1)[1].replace("ms", "")
        try:
            pre_gid, post_gid = [int(x) for x in pair.split("-")]
            dt_ms = float(dt_str)
        except ValueError:
            continue
        basis_csv = os.path.join(BASIS_DIR, f"basis_{pre_gid}_{post_gid}.csv")
        if not os.path.exists(basis_csv):
            continue
        rows.append((sd, pre_gid, post_gid, dt_ms))
    return rows


def collect(source, n_workers=4):
    """
    source: "bcl_full", "bcl_prefire", or "optimizer"
    Returns DataFrame with columns [pre_gid, post_gid, dt, epsp_ratio].
    """
    sim_dirs = _find_sim_dirs()
    if source == "optimizer":
        jobs = [(sd, pre, post, dt, CHINDEMI_HASH) for sd, pre, post, dt in sim_dirs]
        worker = _worker_optimizer
    elif source == "optimizer_fit":
        jobs = [(sd, pre, post, dt, CHINDEMI_FIT_HASH) for sd, pre, post, dt in sim_dirs]
        worker = _worker_optimizer_fit
    elif source == "bcl_prefire":
        jobs = sim_dirs
        worker = _worker_bcl_prefire
    else:  # bcl_full
        jobs = sim_dirs
        worker = _worker_bcl_full

    print(f"  [{source}] {len(jobs)} jobs, {n_workers} workers ...")
    rows = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as pool:
        for r in pool.map(worker, jobs):
            if r is not None:
                rows.append(r)
    df = pd.DataFrame(rows)
    print(f"  [{source}] → {len(df)} rows, "
          f"{df['dt'].nunique() if not df.empty else 0} dt values, "
          f"{df[['pre_gid','post_gid']].drop_duplicates().shape[0] if not df.empty else 0} pairs")
    return df


# ---------------------------------------------------------------------------
# Summary helpers
# ---------------------------------------------------------------------------

def summ(df):
    if df.empty:
        return pd.DataFrame(columns=["dt", "mean", "sem", "count"])
    return (df.groupby("dt")["epsp_ratio"]
              .agg(mean="mean",
                   sem=lambda x: x.std() / np.sqrt(len(x)),
                   count="count")
              .reset_index().sort_values("dt"))


def print_summary(label, df):
    if df.empty or "dt" not in df.columns:
        print(f"\n{label}: no data")
        return
    s = summ(df[df["dt"].isin(TARGET_DT.keys())])
    print(f"\n{label}:")
    for _, row in s.iterrows():
        tgt_mean, tgt_sem = TARGET_DT[row["dt"]]
        print(f"  dt={row['dt']:+.0f}ms  mean={row['mean']:.4f}±{row['sem']:.4f} "
              f"(n={int(row['count'])})  target={tgt_mean:.4f}  "
              f"diff={100*(row['mean']-tgt_mean)/tgt_mean:+.1f}%")


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot(df_bcl, df_prefire, df_opt, df_fit, output):
    all_dts = []
    for df in [df_bcl, df_prefire, df_opt, df_fit]:
        if not df.empty and "dt" in df.columns:
            all_dts += df["dt"].tolist()
    dts_all = sorted(set(all_dts)) if all_dts else []

    s_bcl     = summ(df_bcl)
    s_prefire = summ(df_prefire)
    s_opt     = summ(df_opt)
    s_fit     = summ(df_fit)

    # Training protocols only for bar chart
    dts_train = sorted(TARGET_DT.keys())

    fig, axes = plt.subplots(1, 2, figsize=(15, 5))
    fig.suptitle(
        "Chindemi params — 4-way comparison  (10 Hz)\n"
        "Neurodamus (100p)  |  BCL prefire (100p)  |  pairrunner no-a (100p)  |  pairrunner fit a-params (30p, ±10ms)\n"
        "All: same 42s induction, chindemi γ, basis EPSP",
        fontsize=9
    )

    # ── Panel 1: STDP curve across all dt ─────────────────────────────────
    ax = axes[0]
    ax.errorbar(INVITRO_DT, INVITRO_MEAN, yerr=INVITRO_SEM,
                fmt="o-", color="#e5c07b", capsize=4, lw=2,
                label="In vitro (Markram 1997)", zorder=2)
    if not s_bcl.empty:
        ax.errorbar(s_bcl["dt"], s_bcl["mean"], yerr=s_bcl["sem"],
                    fmt="o-", color="#66c97f", capsize=4, lw=2,
                    label=f"Neurodamus (n≤{int(s_bcl['count'].max())})", zorder=4)
    if not s_prefire.empty:
        ax.errorbar(s_prefire["dt"], s_prefire["mean"], yerr=s_prefire["sem"],
                    fmt="s--", color="#61afef", capsize=4, lw=2, markersize=7,
                    label=f"BCL prefire (n≤{int(s_prefire['count'].max())})", zorder=5)
    if not s_opt.empty:
        ax.errorbar(s_opt["dt"], s_opt["mean"], yerr=s_opt["sem"],
                    fmt="^:", color="#e06c75", capsize=4, lw=2, markersize=7,
                    label=f"pairrunner no-a (n≤{int(s_opt['count'].max())})", zorder=3)
    if not s_fit.empty:
        ax.errorbar(s_fit["dt"], s_fit["mean"], yerr=s_fit["sem"],
                    fmt="D-.", color="#c678dd", capsize=4, lw=2, markersize=7,
                    label=f"pairrunner fit a-params (n≤{int(s_fit['count'].max())})", zorder=6)
    ax.axhline(1.0, color="k", lw=0.8, ls="--", alpha=0.4)
    ax.axvline(0.0, color="k", lw=0.8, ls="--", alpha=0.4)
    ax.set_xlabel(r"$\Delta t$ (ms)", fontsize=11)
    ax.set_ylabel("EPSP ratio", fontsize=11)
    ax.set_title("STDP curve (10 Hz, all dt)", fontsize=11)
    ax.legend(frameon=False, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # ── Panel 2: bar chart for training protocols ──────────────────────────
    ax2 = axes[1]
    dt_labels = [f"dt={int(dt):+d} ms" for dt in dts_train]
    x     = np.arange(len(dts_train))
    width = 0.16  # 5 bars: offsets -2, -1, 0, +1, +2 × width

    def _bar_vals(s, dts):
        means = [s[s["dt"]==dt]["mean"].values[0] if dt in s["dt"].values else np.nan for dt in dts]
        sems  = [s[s["dt"]==dt]["sem"].values[0]  if dt in s["dt"].values else np.nan for dt in dts]
        return means, sems

    # In vitro
    tgt_means = [TARGET_DT[dt][0] for dt in dts_train]
    tgt_sems  = [TARGET_DT[dt][1] for dt in dts_train]
    ax2.bar(x - 2*width, tgt_means, width, yerr=tgt_sems,
            capsize=4, color="#e5c07b", alpha=0.9, label="In vitro target")

    if not s_bcl.empty:
        m, s = _bar_vals(s_bcl, dts_train)
        ax2.bar(x - 1*width, m, width, yerr=s,
                capsize=4, color="#66c97f", alpha=0.9, label="Neurodamus")

    if not s_prefire.empty:
        m, s = _bar_vals(s_prefire, dts_train)
        ax2.bar(x + 0*width, m, width, yerr=s,
                capsize=4, color="#61afef", alpha=0.9, label="BCL prefire (100p)")

    if not s_opt.empty:
        m, s = _bar_vals(s_opt, dts_train)
        ax2.bar(x + 1*width, m, width, yerr=s,
                capsize=4, color="#e06c75", alpha=0.9, label="pairrunner no-a (100p)")

    if not s_fit.empty:
        m, s = _bar_vals(s_fit, dts_train)
        ax2.bar(x + 2*width, m, width, yerr=s,
                capsize=4, color="#c678dd", alpha=0.9, label="pairrunner fit a-params (30p)")

    ax2.axhline(1.0, color="k", lw=0.8, ls="--", alpha=0.4)
    ax2.set_xticks(x)
    ax2.set_xticklabels(dt_labels)
    ax2.set_ylabel("Mean EPSP ratio", fontsize=11)
    ax2.set_title("Training protocols: 3-way comparison", fontsize=11)
    ax2.legend(frameon=False, fontsize=8)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(output, dpi=200, bbox_inches="tight")
    print(f"\nSaved: {output}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--output",  default="stdp_three_way_chindemi.png")
    args = parser.parse_args()

    print(f"Chindemi hash (no a-params): {CHINDEMI_HASH}")
    print(f"Chindemi fit hash (a-params): {CHINDEMI_FIT_HASH}")

    print("\nCollecting Neurodamus (out/rho.h5) ...")
    df_bcl = collect("bcl_full", n_workers=args.workers)

    print("\nCollecting BCL prefire (bluecellulab_results_edges_chindemi_params_prefire/rho.h5) ...")
    df_prefire = collect("bcl_prefire", n_workers=args.workers)

    print("\nCollecting pairrunner no-a (simulation_edges_{chindemi_hash}.pkl) ...")
    df_opt = collect("optimizer", n_workers=args.workers)

    print("\nCollecting pairrunner fit a-params (simulation_edges_{fit_hash}.pkl, 30 pairs) ...")
    df_fit = collect("optimizer_fit", n_workers=args.workers)

    print_summary("Neurodamus (chindemi, 42s)", df_bcl)
    print_summary("BCL prefire chindemi (42s, basis)", df_prefire)
    print_summary("pairrunner no-a (100p, basis)", df_opt)
    print_summary("pairrunner fit a-params (30p, basis)", df_fit)

    if df_bcl.empty and df_prefire.empty:
        print("ERROR: no data at all — check paths")
        sys.exit(1)

    plot(df_bcl, df_prefire, df_opt, df_fit, args.output)


if __name__ == "__main__":
    main()
