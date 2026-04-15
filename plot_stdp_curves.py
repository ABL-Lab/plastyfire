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

PLASTYFIRE_ROOT = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire"

# Add new_fitting to path to import effcai_to_epsp
sys.path.append(os.path.join(PLASTYFIRE_ROOT, "new_fitting"))
from effcai_to_epsp import fetch_epsp_ratio
from submit_l5ttpc_traces import (
    CHINDEMI_PARAMS, DHURUVA_PARAMS, DHURUVA_PARAMS_V2, DHURUVA_PARAMS_V3, DHURUVA_PARAMS_V4,
    DHURUVA_PARAMS_V5, DHURUVA_PARAMS_V6, DHURUVA_PARAMS_V7, DHURUVA_PARAMS_V8,
    DHURUVA_PARAMS_V9, DHURUVA_PARAMS_V10, DHURUVA_PARAMS_V11, DHURUVA_PARAMS_V12,
    DHURUVA_PARAMS_V13, DHURUVA_PARAMS_V14, DHURUVA_PARAMS_V15, DHURUVA_PARAMS_V16,
    DHURUVA_PARAMS_V17, DHURUVA_PARAMS_V18, DHURUVA_PARAMS_V19
)

TRACE_RESULTS_BASE = os.path.join(PLASTYFIRE_ROOT, "full_trace_results")
TRACE_RESULTS_DIRS = {
    "v1": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS"),
    "v2": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V2"),
    "v3": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V3"),
    "v4": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V4"),
    "v5": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_STDP_V5"),
    "v6": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V6"),
    "v7": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V7"),
    "v8": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V8"),
    "v9": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V9"),
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
DEFAULT_BASIS_DIRS = [
    os.path.join(PLASTYFIRE_ROOT, "basis_results"),
    os.path.join(PLASTYFIRE_ROOT, "basis_results_old"),
]


def _resolve_basis_dirs(basis_dir=None):
    """Return existing basis directories, preferring the current basis_results path."""
    if basis_dir is None:
        candidates = DEFAULT_BASIS_DIRS
    elif isinstance(basis_dir, (list, tuple)):
        candidates = list(basis_dir)
    else:
        candidates = [basis_dir]

    existing = [path for path in candidates if os.path.isdir(path)]
    if not existing:
        checked = ", ".join(candidates)
        raise FileNotFoundError(f"No basis directory found. Checked: {checked}")
    return existing


def _find_basis_path(pre_gid, post_gid, basis_dirs):
    filename = f"basis_{pre_gid}_{post_gid}.csv"
    for basis_dir in basis_dirs:
        basis_path = os.path.join(basis_dir, filename)
        if os.path.exists(basis_path):
            return basis_path
    return None


def _process_single_sim(args):
    """Worker function: load one full-trace pickle and compute the basis EPSP ratio.

    This function is called in a separate process by ProcessPoolExecutor.
    It must be defined at module level so it can be pickled.

    Parameters
    ----------
    args : tuple
        (pre_gid, post_gid, pair_name, freq_dt, pkl_path, dt, params_dict, basis_path)

    Returns
    -------
    dict or None
        {"pair": ..., "dt": ..., "epsp_ratio": ...} on success, None on failure.
    """
    pre_gid, post_gid, pair_name, freq_dt, pkl_path, dt, params_dict, basis_path = args

    try:
        with open(pkl_path, "rb") as f:
            sim_data = pickle.load(f)
    except Exception as e:
        print(f"Skipping {pair_name} {freq_dt} - corrupted pickle ({type(e).__name__}: {e}). Deleting {pkl_path}")
        try:
            os.remove(pkl_path)
        except OSError as rm_err:
            print(f"  Warning: could not delete {pkl_path}: {rm_err}")
        return None

    if not sim_data or "rho_GB" not in sim_data:
        print(f"Skipping {pair_name} {freq_dt} - missing data in pickle.")
        return None

    rho_gb_raw = sim_data["rho_GB"]
    # _runconnectedpair_prefire_process stores SYNREC as dict {syn_id: array};
    # old _runconnectedpair_process stored it as a 2D numpy array.
    if isinstance(rho_gb_raw, dict):
        rho_trace = np.array(list(rho_gb_raw.values()))
    else:
        rho_trace = np.asarray(rho_gb_raw)
        if rho_trace.ndim == 2:
            if rho_trace.shape[0] > rho_trace.shape[1]:
                rho_trace = rho_trace.T
        else:
            rho_trace = rho_trace.reshape(1, -1)

    initial_rho = [1 if r[0] >= 0.5 else 0 for r in rho_trace]
    final_rho = [1 if r[-1] >= 0.5 else 0 for r in rho_trace]

    initial = np.array(initial_rho)
    final = np.array(final_rho)
    n_0to1 = int(np.sum((initial == 0) & (final == 1)))
    n_0to0 = int(np.sum((initial == 0) & (final == 0)))
    n_1to1 = int(np.sum((initial == 1) & (final == 1)))
    n_1to0 = int(np.sum((initial == 1) & (final == 0)))

    try:
        basis_dir = os.path.dirname(basis_path)
        res = fetch_epsp_ratio(
            pre_gid,
            post_gid,
            10.0,
            dt,
            initial_rho,
            final_rho,
            basis_dir=basis_dir,
            params=params_dict,
            ratio_method="delta_method",
        )
        return {
            "pair": pair_name,
            "dt": dt,
            "epsp_ratio": res["ratio_mean"],
            "n_0to1": n_0to1,
            "n_0to0": n_0to0,
            "n_1to1": n_1to1,
            "n_1to0": n_1to0,
        }
    except Exception as e:
        print(f"Error processing {pair_name} {freq_dt}: {e}")
        return None


def process_results(trace_results_dir, params_dict=None, n_workers=None, basis_dir=None):
    """Process all full-trace simulation results in parallel across CPU cores."""
    jobs = []
    basis_dirs = _resolve_basis_dirs(basis_dir)

    for pair_dir in glob.glob(os.path.join(trace_results_dir, "*")):
        if not os.path.isdir(pair_dir) or os.path.basename(pair_dir) == "logs":
            continue

        pair_name = os.path.basename(pair_dir)
        try:
            pre_gid, post_gid = map(int, pair_name.split("-"))
        except ValueError:
            continue

        basis_path = _find_basis_path(pre_gid, post_gid, basis_dirs)
        if basis_path is None:
            print(f"Skipping {pair_name}, no basis file.")
            continue

        for freq_dt_dir in glob.glob(os.path.join(pair_dir, "10Hz_*")):
            freq_dt = os.path.basename(freq_dt_dir)
            dt_str = freq_dt.replace("10Hz_", "").replace("ms", "")
            try:
                dt = float(dt_str)
            except ValueError:
                continue

            pkl_path = os.path.join(freq_dt_dir, "simulation_traces.pkl")
            if not os.path.exists(pkl_path):
                continue

            jobs.append((pre_gid, post_gid, pair_name, freq_dt, pkl_path, dt, params_dict, basis_path))

    if not jobs:
        return pd.DataFrame()

    n_workers = n_workers or multiprocessing.cpu_count()
    print(f"Processing {len(jobs)} full-trace simulations across {n_workers} CPU cores...")

    data = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as pool:
        for result in pool.map(_process_single_sim, jobs):
            if result is not None:
                data.append(result)

    return pd.DataFrame(data)


def plot_stdp_curve(df, version, output_filename):
    if df.empty:
        print("No data collected!")
        return

    summary = df.groupby("dt")["epsp_ratio"].agg(["mean", "sem"]).reset_index()
    summary = summary.sort_values("dt")

    invitro_dt = [-10, 5, 10]
    invitro_mean = [0.7922, 1.2038, 1.2013]
    invitro_sem = [0.0259, 0.0644, 0.0626]

    plt.figure(figsize=(6, 5))

    plt.errorbar(summary["dt"], summary["mean"], yerr=summary["sem"],
                 fmt="o-", color="#66b3e6", label=f"in silico ({version})", capsize=0)

    plt.errorbar(invitro_dt, invitro_mean, yerr=invitro_sem,
                 fmt="o-", color="#ffaa55", label="in vitro", capsize=0)

    plt.axhline(1.0, color="k", linestyle="--", alpha=0.5)
    plt.axvline(0.0, color="k", linestyle="--", alpha=0.5)

    plt.xlabel(r"$\Delta t$ (ms)", fontsize=12)
    plt.ylabel("EPSP ratio", fontsize=12)
    plt.title(f"Frequency = 10 Hz ({version})", fontsize=14)
    plt.legend(frameon=False, loc="upper left")

    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_filename, dpi=300)
    print(f"Saved {output_filename}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot STDP curves from full-trace results using basis EPSP extrapolation")
    parser.add_argument("--params", choices=["v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8", "v9", "v10", "v11", "v12", "v13", "v14", "v15", "v16", "v17", "v18", "v19", "chindemi"], default="v1",
                        help="Parameter set results to analyze (default: v1)")
    parser.add_argument("--results-dir", type=str, default=None,
                        help="Override default full-trace results directory")
    parser.add_argument("--basis-dir", type=str, default=None,
                        help="Basis directory to use for EPSP extrapolation (default: auto-detect basis_results, then basis_results_old)")
    parser.add_argument("--workers", type=int, default=None,
                        help="Number of parallel worker processes (default: all CPU cores)")
    args = parser.parse_args()

    version = args.params
    trace_results_dir = args.results_dir or TRACE_RESULTS_DIRS[version]
    basis_dirs = _resolve_basis_dirs(args.basis_dir)

    print("\n" + "=" * 50)
    print(f"STDP Curve Analysis Summary ({version}, basis from full traces)")
    print(f"Full-trace results dir: {trace_results_dir}")
    print(f"Basis dirs: {', '.join(basis_dirs)}")
    print("=" * 50)

    params_mapping = {
        "v1": DHURUVA_PARAMS, "v2": DHURUVA_PARAMS_V2, "v3": DHURUVA_PARAMS_V3,
        "v4": DHURUVA_PARAMS_V4, "v5": DHURUVA_PARAMS_V5, "v6": DHURUVA_PARAMS_V6,
        "v7": DHURUVA_PARAMS_V7, "v8": DHURUVA_PARAMS_V8, "v9": DHURUVA_PARAMS_V9,
        "v10": DHURUVA_PARAMS_V10, "v11": DHURUVA_PARAMS_V11, "v12": DHURUVA_PARAMS_V12,
        "v13": DHURUVA_PARAMS_V13, "v14": DHURUVA_PARAMS_V14, "v15": DHURUVA_PARAMS_V15,
        "v16": DHURUVA_PARAMS_V16, "v17": DHURUVA_PARAMS_V17, "v18": DHURUVA_PARAMS_V18,
        "v19": DHURUVA_PARAMS_V19,
        "chindemi": CHINDEMI_PARAMS,
    }
    selected_params = params_mapping.get(version, DHURUVA_PARAMS)

    df = process_results(trace_results_dir, params_dict=selected_params, n_workers=args.workers,
                         basis_dir=args.basis_dir)

    if not df.empty:
        n_pairs = df["pair"].nunique()
        total_sims = len(df)
        print(f"Total pairs processed: {n_pairs}")
        print(f"Total simulations:     {total_sims}")
        print("\nSimulation Counts by Delta T:")
        print(df.groupby("dt").size().to_string())

        summary = df.groupby("dt")["epsp_ratio"].agg(["mean", "sem", "count"]).reset_index()
        summary = summary.sort_values("dt")

        print("\nStatistics by Delta T:")
        print(summary.to_string(index=False))

        invitro = {-50.0: 1.0, -10.0: 0.7922, 5.0: 1.2038, 10.0: 1.2013}

        overlap = summary[summary["dt"].isin(invitro.keys())]
        if not overlap.empty:
            errors = []
            for _, row in overlap.iterrows():
                target = invitro[row["dt"]]
                error = row["mean"] - target
                errors.append(error ** 2)

            rmse = np.sqrt(np.mean(errors))
            print(f"\nRMSE against In-Vitro (10Hz, -50, -10, 5, 10ms): {rmse:.4f}")

        output_filename = f"stdp_curve_10Hz_{version}.png"
        plot_stdp_curve(df, version, output_filename)

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

        print("\nAnalysis complete.")
    else:
        print("No results found in " + trace_results_dir)
    print("=" * 50)
