#!/usr/bin/env python3
"""
Plot STDP curve for GA best-fit params validated against full BCL simulations.

Reads simulation_edges.pkl (written by pairrunner_edges.py --params ga_best) which
contains continuous final_rho, then computes EPSP ratios via compute_epsp_from_basis
(binary threshold at 0.5 + delta-method bias correction) — exactly matching the
optimizer's fitness evaluation.

Usage:
    python plot_ga_best_bcl.py [--workers N] [--output FILE]
"""

import os
import sys
import glob
import pickle
import argparse
import concurrent.futures
import multiprocessing

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PLASTYFIRE_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(PLASTYFIRE_ROOT, "new_fitting"))
sys.path.insert(0, PLASTYFIRE_ROOT)

from plastyfire.evaluator_edges import compute_epsp_from_basis, _load_basis

RESULTS_DIR    = os.path.join(PLASTYFIRE_ROOT, "refitting_results")
BASIS_DIR      = os.path.join(PLASTYFIRE_ROOT, "basis_results_edges_mini")
CHECKPOINT_PKL = os.path.join(PLASTYFIRE_ROOT, "checkpoint_edges.pkl")
CACHE_DIR      = os.path.join(PLASTYFIRE_ROOT, ".cache")
# simulation_edges.pkl was written (overwritten) by the ga_best BCL run
SIM_PKL_NAME = "simulation_edges.pkl"

# Training protocols (mrk97_07 = +10ms, mrk97_08 = -10ms)
TARGET_DT = {10.0: (1.2013, 0.0626), -10.0: (0.7922, 0.0259)}

# In-vitro reference (Markram 1997, 10 Hz) — all measured dt values
INVITRO_DT   = [-10,    5,    10]
INVITRO_MEAN = [0.7922, 1.2038, 1.2013]
INVITRO_SEM  = [0.0259, 0.0644, 0.0626]

# GA best params (hall of fame, gen 28, fitness=0.3942)
GA_BEST_PARAMS = {
    "tau_effca_GB_GluSynapse": 278.3177658387,
    "gamma_d_GB_GluSynapse":    94.698202,
    "gamma_p_GB_GluSynapse":   211.404335,
    "a00": 1.000406, "a01": 1.954000,
    "a10": 1.160326, "a11": 2.590750,
    "a20": 1.164831, "a21": 2.489699,
    "a30": 3.670398, "a31": 1.566718,
}


def _process_one(args):
    sim_dir, pre_gid, post_gid, dt_ms = args
    pkl_path  = os.path.join(sim_dir, SIM_PKL_NAME)
    basis_csv = os.path.join(BASIS_DIR, f"basis_{pre_gid}_{post_gid}.csv")
    if not os.path.exists(pkl_path) or not os.path.exists(basis_csv):
        return None
    try:
        with open(pkl_path, "rb") as f:
            raw = pickle.load(f)
        # Use continuous rho — same as evaluator_edges.compute_epsp_ratio_basis_batch
        initial_rho = list(raw["initial_rho"])
        final_rho   = list(raw["final_rho"])
        basis_df = _load_basis(BASIS_DIR, str(pre_gid), str(post_gid))
        epsp_before_mean, epsp_before_std = compute_epsp_from_basis(basis_df, initial_rho)
        epsp_after_mean,  _               = compute_epsp_from_basis(basis_df, final_rho)
        if epsp_before_mean == 0:
            return None
        cv2 = min((epsp_before_std / epsp_before_mean) ** 2, 0.25)
        epsp_ratio = (epsp_after_mean / epsp_before_mean) * (1.0 + cv2)
        # Binary transitions for diagnostic info
        ini_bin = [1 if r >= 0.5 else 0 for r in initial_rho]
        fin_bin = [1 if r >= 0.5 else 0 for r in final_rho]
        return {"pre_gid": pre_gid, "post_gid": post_gid,
                "dt": dt_ms, "epsp_ratio": epsp_ratio,
                "n_syn": len(initial_rho),
                "n_0to1": sum(1 for i, f in zip(ini_bin, fin_bin) if i == 0 and f == 1),
                "n_1to0": sum(1 for i, f in zip(ini_bin, fin_bin) if i == 1 and f == 0)}
    except Exception as e:
        print(f"  ERROR {pre_gid}-{post_gid} dt={dt_ms}: {e}")
        return None


def get_optimizer_pkl_hash():
    """Return (pkl_hash_12, resdb) from the checkpoint hall-of-fame best individual."""
    import hashlib as _hl
    with open(CHECKPOINT_PKL, "rb") as f:
        cp = pickle.load(f)
    hof  = cp["halloffame"][0]
    h12  = _hl.md5(str(list(hof)).encode()).hexdigest()[:12]
    hfull = _hl.md5(str(list(hof)).encode()).hexdigest()
    cache_pkl = os.path.join(CACHE_DIR, f"{hfull}.pkl")
    resdb = None
    if os.path.exists(cache_pkl):
        with open(cache_pkl, "rb") as f:
            resdb = pickle.load(f).get("resdb")
    return h12, resdb


def collect_bcl_ga_best(results_dir=RESULTS_DIR, freq=10, n_workers=4,
                        opt_pairs_only=False, opt_pkl_hash=None):
    """Collect BCL ga_best EPSP ratios from simulation_edges.pkl.

    opt_pairs_only: if True, restrict to pairs where simulation_edges_{opt_pkl_hash}.pkl exists
                    (i.e. the same 30 pairs the optimizer evaluated).
    """
    pattern = os.path.join(results_dir, "fitting", "*", "seed*",
                           "*_STDP", "simulations", "*-*", f"{int(freq)}Hz_*")
    sim_dirs = sorted(glob.glob(pattern))
    jobs = []
    for sd in sim_dirs:
        pair   = os.path.basename(os.path.dirname(sd))
        dt_str = os.path.basename(sd).split("_", 1)[1].replace("ms", "")
        try:
            pre_gid, post_gid = [int(x) for x in pair.split("-")]
            dt_ms = float(dt_str)
        except ValueError:
            continue
        if not os.path.exists(os.path.join(sd, SIM_PKL_NAME)):
            continue
        if not os.path.exists(os.path.join(BASIS_DIR, f"basis_{pre_gid}_{post_gid}.csv")):
            continue
        if opt_pairs_only and opt_pkl_hash:
            if not os.path.exists(os.path.join(sd, f"simulation_edges_{opt_pkl_hash}.pkl")):
                continue
        jobs.append((sd, pre_gid, post_gid, dt_ms))

    label = f"{'optimizer 30' if opt_pairs_only else 'all 100'} pairs"
    print(f"  Collecting BCL ga_best ({label}): {len(jobs)} sim dirs, {n_workers} workers...")
    rows = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as pool:
        for r in pool.map(_process_one, jobs):
            if r is not None:
                rows.append(r)
    return pd.DataFrame(rows)


def summ_df(df):
    return (df.groupby("dt")["epsp_ratio"]
              .agg(mean="mean", sem=lambda x: x.std() / np.sqrt(len(x)), count="count")
              .reset_index().sort_values("dt"))


def print_summary(label, df):
    s = summ_df(df[df["dt"].isin(TARGET_DT.keys())])
    print(f"\n{label}:")
    for _, row in s.iterrows():
        tgt_mean, tgt_sem = TARGET_DT[row["dt"]]
        print(f"  dt={row['dt']:+.0f}ms  mean={row['mean']:.4f}±{row['sem']:.4f} "
              f"(n={int(row['count'])})  target={tgt_mean:.4f}  "
              f"diff={100*(row['mean']-tgt_mean)/tgt_mean:+.1f}%")


def plot(df_all, df_opt30, resdb_opt, output):
    """
    df_all    : BCL ga_best results for all 100 pairs (all dt)
    df_opt30  : BCL ga_best results for the same 30 pairs the optimizer evaluated
    resdb_opt : optimizer's own cached EPSP ratios (30 pairs, basis method)
    """
    # Training protocols only
    dts = sorted(TARGET_DT.keys())

    s_all   = summ_df(df_all[df_all["dt"].isin(dts)])
    s_opt30 = summ_df(df_opt30[df_opt30["dt"].isin(dts)])
    s_opt   = summ_df(resdb_opt.rename(columns={"protocol_id": "dt"}).assign(
                  dt=resdb_opt["protocol_id"].map({"mrk97_07": 10.0, "mrk97_08": -10.0})
              )) if resdb_opt is not None else None

    # Full STDP curve (all dt)
    s_all_full = summ_df(df_all)

    # ── Figure ───────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))
    fig.suptitle(
        "GA best-fit params: BCL simulation vs optimizer (training protocols)\n"
        r"$\gamma_d$=94.70, $\gamma_p$=211.40, $\tau$=278.32 (fixed)",
        fontsize=11
    )

    # Panel 1: STDP curve — all 100 BCL pairs across all dt
    ax = axes[0]
    ax.errorbar(INVITRO_DT, INVITRO_MEAN, yerr=INVITRO_SEM,
                fmt="o-", color="#e5c07b", capsize=4, lw=2,
                label="In vitro (Markram 1997)", zorder=2)
    ax.errorbar(s_all_full["dt"], s_all_full["mean"], yerr=s_all_full["sem"],
                fmt="s-", color="#61afef", capsize=4, lw=2, markersize=7,
                label=f"BCL sim — ga_best (n=100 pairs)", zorder=3)
    for dt_val in dts:
        ax.axvline(dt_val, color="gray", lw=0.8, ls=":", alpha=0.5)
    ax.axhline(1.0, color="k", lw=0.8, ls="--", alpha=0.4)
    ax.axvline(0.0, color="k", lw=0.8, ls="--", alpha=0.4)
    ax.set_xlabel(r"$\Delta t$ (ms)", fontsize=11)
    ax.set_ylabel("EPSP ratio", fontsize=11)
    ax.set_title("STDP curve  (10 Hz, all 100 pairs)", fontsize=11)
    ax.legend(frameon=False, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Panel 2: bar chart for training protocols — 4 bars: invitro / optimizer / BCL-30 / BCL-100
    ax2 = axes[1]
    dt_labels = [f"dt={int(dt):+d} ms" for dt in dts]
    x     = np.arange(len(dts))
    width = 0.20

    tgt_means = [TARGET_DT[dt][0] for dt in dts]
    tgt_sems  = [TARGET_DT[dt][1] for dt in dts]

    ax2.bar(x - 1.5*width, tgt_means, width, yerr=tgt_sems,
            capsize=4, color="#e5c07b", alpha=0.9, label="In vitro target")

    if s_opt is not None:
        opt_means = [s_opt[s_opt["dt"]==dt]["mean"].values[0] if dt in s_opt["dt"].values else np.nan for dt in dts]
        opt_sems  = [s_opt[s_opt["dt"]==dt]["sem"].values[0]  if dt in s_opt["dt"].values else np.nan for dt in dts]
        ax2.bar(x - 0.5*width, opt_means, width, yerr=opt_sems,
                capsize=4, color="#98c379", alpha=0.9, label="Optimizer basis (n=30)")

    bcl30_means = [s_opt30[s_opt30["dt"]==dt]["mean"].values[0] if dt in s_opt30["dt"].values else np.nan for dt in dts]
    bcl30_sems  = [s_opt30[s_opt30["dt"]==dt]["sem"].values[0]  if dt in s_opt30["dt"].values else np.nan for dt in dts]
    ax2.bar(x + 0.5*width, bcl30_means, width, yerr=bcl30_sems,
            capsize=4, color="#61afef", alpha=0.9, label="BCL full sim (n=30, same pairs)")

    bcl100_means = [s_all[s_all["dt"]==dt]["mean"].values[0] if dt in s_all["dt"].values else np.nan for dt in dts]
    bcl100_sems  = [s_all[s_all["dt"]==dt]["sem"].values[0]  if dt in s_all["dt"].values else np.nan for dt in dts]
    ax2.bar(x + 1.5*width, bcl100_means, width, yerr=bcl100_sems,
            capsize=4, color="#c678dd", alpha=0.9, label="BCL full sim (n=100, all pairs)")

    ax2.axhline(1.0, color="k", lw=0.8, ls="--", alpha=0.4)
    ax2.set_xticks(x)
    ax2.set_xticklabels(dt_labels)
    ax2.set_ylabel("Mean EPSP ratio", fontsize=11)
    ax2.set_title("Training protocols: optimizer vs BCL validation", fontsize=11)
    ax2.legend(frameon=False, fontsize=8)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(output, dpi=200, bbox_inches="tight")
    print(f"\nSaved: {output}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--output",  default="stdp_ga_best_bcl.png")
    args = parser.parse_args()

    # Get optimizer HOF hash and cached resdb
    opt_hash, resdb_opt = get_optimizer_pkl_hash()
    print(f"Optimizer HOF pkl_hash: {opt_hash}")
    if resdb_opt is not None:
        print(f"Optimizer resdb: {len(resdb_opt)} rows")

    print("\nCollecting BCL ga_best — all 100 pairs (all dt)...")
    df_all = collect_bcl_ga_best(n_workers=args.workers,
                                 opt_pairs_only=False, opt_pkl_hash=opt_hash)
    print(f"  → {len(df_all)} rows, {df_all['dt'].nunique()} dt values, "
          f"{df_all[['pre_gid','post_gid']].drop_duplicates().shape[0]} pairs")

    print("\nCollecting BCL ga_best — same 30 pairs as optimizer...")
    df_opt30 = collect_bcl_ga_best(n_workers=args.workers,
                                   opt_pairs_only=True, opt_pkl_hash=opt_hash)
    print(f"  → {len(df_opt30)} rows, "
          f"{df_opt30[['pre_gid','post_gid']].drop_duplicates().shape[0]} pairs")

    if df_all.empty:
        print("ERROR: no data — check SIM_PKL_NAME and results path")
        sys.exit(1)

    print_summary("Optimizer basis (n=30)", resdb_opt.rename(
        columns={"protocol_id":"dt"}).assign(
        dt=resdb_opt["protocol_id"].map({"mrk97_07":10.0,"mrk97_08":-10.0}))
    ) if resdb_opt is not None else None
    print_summary("BCL full sim — same 30 pairs", df_opt30)
    print_summary("BCL full sim — all 100 pairs", df_all)

    plot(df_all, df_opt30, resdb_opt, args.output)


if __name__ == "__main__":
    main()
