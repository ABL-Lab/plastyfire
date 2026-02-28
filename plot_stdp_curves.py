import os
import sys
import argparse
import pickle
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Add new_fitting to path to import effcai_to_epsp
sys.path.append("/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/new_fitting")
from effcai_to_epsp import fetch_epsp_ratio

TRACE_RESULTS_BASE = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/trace_results"
TRACE_RESULTS_DIRS = {
    "v1": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS"),
    "v2": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V2"),
    "v3": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V3"),
    "v4": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V4"),
    "v5": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_STDP_V5"),
    "v6": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V6"),
    "v7": os.path.join(TRACE_RESULTS_BASE, "DHURUVA_PARAMS_V7"),
}
BASIS_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/basis_results"

def process_results(trace_results_dir):
    data = []

    # Iterate over pairs
    for pair_dir in glob.glob(os.path.join(trace_results_dir, "*")):
        if not os.path.isdir(pair_dir) or os.path.basename(pair_dir) == "logs":
            continue

        pair_name = os.path.basename(pair_dir)
        try:
            pre_gid, post_gid = map(int, pair_name.split("-"))
        except ValueError:
            continue

        # Ensure basis exists
        basis_path = os.path.join(BASIS_DIR, f"basis_{pre_gid}_{post_gid}.csv")
        if not os.path.exists(basis_path):
            print(f"Skipping {pair_name}, no basis file.")
            continue

        # Process 10Hz directories
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

            try:
                with open(pkl_path, "rb") as f:
                    sim_data = pickle.load(f)
            except Exception as e:
                print(f"Skipping {pair_name} {freq_dt} - corrupted pickle ({type(e).__name__}: {e})")
                continue

            if not sim_data or "rho_GB" not in sim_data:
                print(f"Skipping {pair_name} {freq_dt} - missing data in pickle.")
                continue

            rho_trace = np.asarray(sim_data["rho_GB"])
            if rho_trace.ndim == 2:
                # Transpose if shape is (timepoints, synapses) instead of (synapses, timepoints)
                if rho_trace.shape[0] > 100 and rho_trace.shape[1] < 100:
                    rho_trace = rho_trace.T
            else:
                rho_trace = rho_trace.reshape(1, -1)

            initial_rho = [1 if r[0] >= 0.5 else 0 for r in rho_trace]
            final_rho = [1 if r[-1] >= 0.5 else 0 for r in rho_trace]

            try:
                res = fetch_epsp_ratio(pre_gid, post_gid, 10.0, dt, initial_rho, final_rho, basis_dir=BASIS_DIR)
                epsp_ratio = res["ratio_mean"]
                data.append({
                    "pair": pair_name,
                    "dt": dt,
                    "epsp_ratio": epsp_ratio
                })
            except Exception as e:
                print(f"Error processing {pair_name} {freq_dt}: {e}")

    return pd.DataFrame(data)

def plot_stdp_curve(df, version, output_filename):
    if df.empty:
        print("No data collected!")
        return

    # Group by dt
    summary = df.groupby("dt")["epsp_ratio"].agg(["mean", "sem"]).reset_index()
    summary = summary.sort_values("dt")

    # In vitro data for 10Hz (Markram 1997)
    invitro_dt = [-10, 5, 10]
    invitro_mean = [0.7922, 1.2038, 1.2013]
    invitro_sem = [0.0259, 0.0644, 0.0626]

    plt.figure(figsize=(6, 5))

    # Plot in silico
    plt.errorbar(summary["dt"], summary["mean"], yerr=summary["sem"],
                 fmt="o-", color="#66b3e6", label=f"in silico ({version})", capsize=0)

    # Plot in vitro
    plt.errorbar(invitro_dt, invitro_mean, yerr=invitro_sem,
                 fmt="o-", color="#ffaa55", label="in vitro", capsize=0)

    plt.axhline(1.0, color="k", linestyle="--", alpha=0.5)
    plt.axvline(0.0, color="k", linestyle="--", alpha=0.5)

    plt.xlabel(r"$\Delta t$ (ms)", fontsize=12)
    plt.ylabel("EPSP ratio", fontsize=12)
    plt.title(f"Frequency = 10 Hz ({version})", fontsize=14)
    plt.legend(frameon=False, loc="upper left")

    # Spines
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_filename, dpi=300)
    print(f"Saved {output_filename}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot STDP curves from simulation results")
    parser.add_argument("--params", choices=["v1", "v2", "v3", "v4", "v5", "v6", "v7"], default="v1",
                        help="Parameter set results to analyze (default: v1)")
    parser.add_argument("--results-dir", type=str, default=None,
                        help="Override default results directory")
    args = parser.parse_args()

    version = args.params
    trace_results_dir = args.results_dir or TRACE_RESULTS_DIRS[version]

    print("\n" + "="*50)
    print(f"STDP Curve Analysis Summary ({version})")
    print(f"Results dir: {trace_results_dir}")
    print("="*50)

    # Process all directories and calculate means instead of replacing `df` completely
    df = process_results(trace_results_dir)

    if not df.empty:
        n_pairs = df["pair"].nunique()
        total_sims = len(df)
        print(f"Total pairs processed: {n_pairs}")
        print(f"Total simulations:     {total_sims}")
        print("\nSimulation Counts by Delta T:")
        print(df.groupby("dt").size().to_string())

        # Calculate Stats
        summary = df.groupby("dt")["epsp_ratio"].agg(["mean", "sem", "count"]).reset_index()
        summary = summary.sort_values("dt")

        print("\nStatistics by Delta T:")
        print(summary.to_string(index=False))

        # Calculate RMSE against in-vitro data (provided in plot_stdp_curve)
        invitro = { -50.0: 1.0, -10.0: 0.7922, 5.0: 1.2038, 10.0: 1.2013 }

        overlap = summary[summary["dt"].isin(invitro.keys())]
        if not overlap.empty:
            errors = []
            for _, row in overlap.iterrows():
                target = invitro[row["dt"]]
                error = row["mean"] - target
                errors.append(error**2)

            rmse = np.sqrt(np.mean(errors))
            print(f"\nRMSE against In-Vitro (10Hz, -50, -10, 5, 10ms): {rmse:.4f}")

        output_filename = f"stdp_curve_10Hz_{version}.png"
        plot_stdp_curve(df, version, output_filename)
        print("\nAnalysis complete.")
    else:
        print("No results found in " + trace_results_dir)
    print("="*50)
