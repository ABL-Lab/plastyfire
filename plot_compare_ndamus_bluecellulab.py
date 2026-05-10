#!/usr/bin/env python3
"""
Side-by-side STDP curve comparison: neurodamus vs bluecellulab.

Both simulators write rho.h5 in the same refitting_results tree but under
different subdirectories:
  neurodamus   : <sim_dir>/out/rho.h5
  bluecellulab : <sim_dir>/bluecellulab_results/rho.h5

Usage:
    python plot_compare_ndamus_bluecellulab.py
    python plot_compare_ndamus_bluecellulab.py --freq 10 --workers 8
    python plot_compare_ndamus_bluecellulab.py --output comparison_10Hz.png
    python plot_compare_ndamus_bluecellulab.py --basis-dir basis_results_edges_mini
    python plot_compare_ndamus_bluecellulab.py --trace-figs --figs-dir neurodamus_bluecellulab_figures
"""

import os
import sys
import glob
import pickle
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
DEFAULT_BASIS_DIR   = os.path.join(PLASTYFIRE_ROOT, "basis_results_edges_mini")

NODE_POP = "S1nonbarrel_neurons"

INVITRO_DT   = [-10,    5,      10]
INVITRO_MEAN = [0.7922, 1.2038, 1.2013]
INVITRO_SEM  = [0.0259, 0.0644, 0.0626]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _find_basis_path(pre_gid, post_gid, basis_dir):
    path = os.path.join(basis_dir, f"basis_{pre_gid}_{post_gid}.csv")
    return path if os.path.exists(path) else None


def _read_rho_h5(rho_h5_path, sort_by_element_id=False):
    with h5py.File(rho_h5_path, "r") as f:
        data = f[f"report/{NODE_POP}/data"][()]   # (n_timepoints, n_synapses)
        if sort_by_element_id:
            eids     = f[f"report/{NODE_POP}/mapping/element_ids"][()]
            sort_idx = np.argsort(eids)
            data     = data[:, sort_idx]
    initial = [1 if v >= 0.5 else 0 for v in data[0]]
    final   = [1 if v >= 0.5 else 0 for v in data[-1]]
    n_0to1 = sum(1 for i, f_ in zip(initial, final) if i == 0 and f_ == 1)
    n_0to0 = sum(1 for i, f_ in zip(initial, final) if i == 0 and f_ == 0)
    n_1to1 = sum(1 for i, f_ in zip(initial, final) if i == 1 and f_ == 1)
    n_1to0 = sum(1 for i, f_ in zip(initial, final) if i == 1 and f_ == 0)
    return initial, final, n_0to1, n_0to0, n_1to1, n_1to0


def _process_single_sim(args):
    sim_dir, pre_gid, post_gid, dt_ms, basis_path, rho_subdir = args
    rho_h5 = os.path.join(sim_dir, rho_subdir, "rho.h5")
    if not os.path.exists(rho_h5):
        return None
    # ND element_ids are local synapse IDs (non-sequential); sort ascending to match
    # the basis CSV positions (BCL iterates in ascending local synapse ID order).
    # BCL element_ids are global circuit IDs — data is already in BCL iteration order.
    sort_elements = (rho_subdir == "out")
    try:
        initial_rho, final_rho, n_0to1, n_0to0, n_1to1, n_1to0 = _read_rho_h5(rho_h5, sort_by_element_id=sort_elements)
    except Exception as e:
        print(f"Read error {pre_gid}-{post_gid} dt={dt_ms}: {e}")
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
        print(f"EPSP error {pre_gid}-{post_gid} dt={dt_ms}: {e}")
        return None
    return {
        "pair":       f"{pre_gid}-{post_gid}",
        "pre_gid":    pre_gid,
        "post_gid":   post_gid,
        "dt":         dt_ms,
        "epsp_ratio": res["ratio_mean"],
        "n_0to1":     n_0to1,
        "n_0to0":     n_0to0,
        "n_1to1":     n_1to1,
        "n_1to0":     n_1to0,
    }


def collect_jobs(results_dir, basis_dir, rho_subdir, freq=10):
    jobs = []
    sim_dirs = glob.glob(
        os.path.join(results_dir, "fitting", "*", "seed*", "*_STDP",
                     "simulations", "*-*", f"{int(freq)}Hz_*")
    )
    for sim_dir in sorted(sim_dirs):
        if not os.path.exists(os.path.join(sim_dir, rho_subdir, "rho.h5")):
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
        basis_path = _find_basis_path(pre_gid, post_gid, basis_dir)
        if basis_path is None:
            print(f"No basis for {pair_name}, skipping.")
            continue
        jobs.append((sim_dir, pre_gid, post_gid, dt_ms, basis_path, rho_subdir))
    return jobs


def process_results(results_dir, basis_dir, rho_subdir, freq=10, n_workers=None):
    jobs = collect_jobs(results_dir, basis_dir, rho_subdir, freq=freq)
    if not jobs:
        print(f"No rho.h5 found under '{rho_subdir}/' (freq={freq} Hz)")
        return pd.DataFrame()
    n_workers = n_workers or multiprocessing.cpu_count()
    print(f"  [{rho_subdir}] {len(jobs)} sims, {n_workers} workers...")
    data = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as pool:
        for result in pool.map(_process_single_sim, jobs):
            if result is not None:
                data.append(result)
    return pd.DataFrame(data)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _plot_one(ax, df, title, color):
    if df.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(title)
        return

    summary = df.groupby("dt")["epsp_ratio"].agg(["mean", "sem"]).reset_index().sort_values("dt")

    ax.errorbar(summary["dt"], summary["mean"], yerr=summary["sem"],
                fmt="o-", color=color, capsize=3, label="in silico")
    ax.errorbar(INVITRO_DT, INVITRO_MEAN, yerr=INVITRO_SEM,
                fmt="o-", color="#ffaa55", capsize=3, label="in vitro (Markram 1997)")

    ax.axhline(1.0, color="k", linestyle="--", alpha=0.4, lw=0.8)
    ax.axvline(0.0, color="k", linestyle="--", alpha=0.4, lw=0.8)

    n_pairs = df["pair"].nunique()
    ax.set_title(f"{title}\n(n={n_pairs} pairs)", fontsize=11)
    ax.set_xlabel(r"$\Delta t$ (ms)", fontsize=10)
    ax.set_ylabel("EPSP ratio", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    invitro_map = dict(zip(INVITRO_DT, INVITRO_MEAN))
    overlap = summary[summary["dt"].isin(invitro_map)]
    if not overlap.empty:
        rmse = np.sqrt(np.mean([(r["mean"] - invitro_map[r["dt"]]) ** 2
                                 for _, r in overlap.iterrows()]))
        ax.text(0.97, 0.05, f"RMSE={rmse:.3f}", ha="right", va="bottom",
                transform=ax.transAxes, fontsize=8, color="gray")


def plot_comparison(df_ndamus, df_bcl, output, freq=10):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    fig.suptitle(f"STDP curve comparison — {int(freq)} Hz", fontsize=13, y=1.01)

    _plot_one(axes[0], df_ndamus, "Neurodamus",    "#4c9be8")
    _plot_one(axes[1], df_bcl,    "Bluecellulab",  "#66c97f")

    # shared y-axis label only on left
    axes[1].set_ylabel("")

    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches="tight")
    print(f"Saved {output}")


# ---------------------------------------------------------------------------
# Per-pair trace figures
# ---------------------------------------------------------------------------

def _load_ndamus_traces(workdir):
    """Load neurodamus voltage and rho time series from out/.

    Synapses are reordered by ascending element_id so the order matches the
    bluecellulab global_ids sort order.
    """
    soma_path = os.path.join(workdir, "out", "soma.h5")
    rho_path  = os.path.join(workdir, "out",  "rho.h5")
    if not (os.path.exists(soma_path) and os.path.exists(rho_path)):
        return None
    with h5py.File(soma_path, "r") as f:
        soma_data = f[f"report/{NODE_POP}/data"][()]       # (n_t, 1)
        t_meta    = f[f"report/{NODE_POP}/mapping/time"][()] # [t_start, t_stop, dt]
    with h5py.File(rho_path, "r") as f:
        rho_data    = f[f"report/{NODE_POP}/data"][()]         # (n_t, n_syn)
        element_ids = f[f"report/{NODE_POP}/mapping/element_ids"][()]
    # Sort columns by element_id so synapses are in a canonical order
    sort_idx = np.argsort(element_ids)
    rho_data = rho_data[:, sort_idx]
    t = np.linspace(t_meta[0], t_meta[1], len(soma_data), endpoint=False)
    return {"t": t, "v": soma_data[:, 0], "rho": rho_data, "tstop": t_meta[1]}


def _load_bcl_traces(workdir):
    """Load bluecellulab voltage and rho from simulation_edges.pkl.

    Synapses are reordered by ascending global_id so the order matches the
    neurodamus element_ids sort order.
    """
    pkl_path = os.path.join(workdir, "simulation_edges.pkl")
    if not os.path.exists(pkl_path):
        return None
    with open(pkl_path, "rb") as f:
        r = pickle.load(f)
    global_ids  = np.asarray(r["global_ids"])
    initial_rho = np.asarray(r["initial_rho"])
    final_rho   = np.asarray(r["final_rho"])
    # Sort by global_id for canonical synapse ordering
    sort_idx    = np.argsort(global_ids)
    return {
        "t":           np.asarray(r["t"]),
        "v":           np.asarray(r["v"]),
        "prespikes":   np.asarray(r["prespikes"]),
        "postspikes":  np.asarray(r["postspikes"]),
        "initial_rho": initial_rho[sort_idx],
        "final_rho":   final_rho[sort_idx],
    }


def plot_pair_traces(workdir, figs_dir):
    """
    Generate one stacked figure per workdir with N+1 rows:
      Row 0       : soma voltage overlay (blue=neurodamus, red=bluecellulab), zoomed 990–3000 ms
      Rows 1..N   : rho(t) per synapse overlay (blue=neurodamus time series,
                    red=bluecellulab step from initial→final at last prespike)
    Saved as neurodamus_bluecellulab_figures/{pair}/{freq_dt}.png
    """
    nd  = _load_ndamus_traces(workdir)
    bcl = _load_bcl_traces(workdir)
    if nd is None and bcl is None:
        return None

    pair    = os.path.basename(os.path.dirname(workdir))
    freq_dt = os.path.basename(workdir)

    n_syn = nd["rho"].shape[1] if nd is not None else len(bcl["initial_rho"])
    n_rows = 1 + n_syn

    fig, axes = plt.subplots(n_rows, 1, figsize=(10, 3 * n_rows), sharex=False)
    if n_rows == 1:
        axes = [axes]
    fig.suptitle(f"Pair {pair}  —  {freq_dt}", fontsize=12, y=1.01)

    V_ZOOM = (990., 3000.)   # ms

    # ── Row 0: soma voltage zoomed 990–3000 ms ────────────────────────────────
    ax = axes[0]
    if nd is not None:
        mask = (nd["t"] >= V_ZOOM[0]) & (nd["t"] <= V_ZOOM[1])
        ax.plot(nd["t"][mask], nd["v"][mask], lw=0.8, color="#4c9be8",
                label="Neurodamus", zorder=2)
    if bcl is not None:
        mask = (bcl["t"] >= V_ZOOM[0]) & (bcl["t"] <= V_ZOOM[1])
        ax.plot(bcl["t"][mask], bcl["v"][mask], lw=0.8, color="#e05252",
                alpha=0.8, label="Bluecellulab", zorder=1)
    ax.set_xlim(V_ZOOM)
    ax.set_xlabel("t (ms)")
    ax.set_ylabel("V (mV)")
    ax.set_title("Soma voltage (990–3000 ms)")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # ── Rows 1..N: rho per synapse ────────────────────────────────────────────
    # Build bluecellulab rho step function: rho=initial until last prespike, then final
    if bcl is not None:
        t_last_pre = float(bcl["prespikes"].max()) if len(bcl["prespikes"]) else 42000.
    else:
        t_last_pre = 42000.

    for i in range(n_syn):
        ax = axes[1 + i]

        if nd is not None:
            ax.plot(nd["t"], nd["rho"][:, i], lw=0.8, color="#4c9be8",
                    label="Neurodamus", zorder=2)

        if bcl is not None:
            rho0 = float(bcl["initial_rho"][i])
            rho1 = float(bcl["final_rho"][i])
            t_bcl = bcl["t"]
            t_end = float(t_bcl[-1])
            # Step: initial → final at last prespike
            t_step = np.array([0., t_last_pre, t_last_pre, t_end])
            r_step = np.array([rho0, rho0, rho1, rho1])
            ax.plot(t_step, r_step, lw=1.2, color="#e05252",
                    linestyle="--", label="Bluecellulab (initial→final)", zorder=1)
            ax.scatter([0., t_end], [rho0, rho1], color="#e05252", s=30, zorder=3)

        ax.set_ylim(-0.05, 1.05)
        ax.set_xlabel("t (ms)")
        ax.set_ylabel("rho")
        ax.set_title(f"Synapse {i}  rho(t)")
        ax.legend(frameon=False, fontsize=8, loc="upper right")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout()

    pair_dir = os.path.join(figs_dir, pair)
    os.makedirs(pair_dir, exist_ok=True)
    out_path = os.path.join(pair_dir, f"{freq_dt}.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _trace_worker(args):
    workdir, figs_dir = args
    try:
        return plot_pair_traces(workdir, figs_dir)
    except Exception as e:
        print(f"Error {workdir}: {e}")
        return None


def generate_trace_figures(results_dir, figs_dir, freq=10, n_workers=None):
    """Generate per-workdir trace figures for all pairs with both simulators done."""
    freq_str = f"{int(freq)}Hz"
    pattern  = os.path.join(results_dir, "fitting", "*", "seed*", "*_STDP",
                             "simulations", "*-*", f"{freq_str}_*")
    workdirs = sorted(glob.glob(pattern))

    # Keep only dirs that have at least one of neurodamus or bluecellulab results
    valid = [d for d in workdirs
             if os.path.exists(os.path.join(d, "out", "rho.h5")) or
                os.path.exists(os.path.join(d, "simulation_edges.pkl"))]
    print(f"Found {len(valid)} workdirs with at least one simulator result.")

    os.makedirs(figs_dir, exist_ok=True)
    n_workers = n_workers or min(multiprocessing.cpu_count(), 8)
    print(f"Generating trace figures with {n_workers} workers → {figs_dir}")

    args = [(d, figs_dir) for d in valid]
    saved = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as pool:
        for result in pool.map(_trace_worker, args):
            if result is not None:
                saved.append(result)

    print(f"Saved {len(saved)} figures.")
    return saved


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Side-by-side STDP comparison: neurodamus vs bluecellulab",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--results-dir", default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--basis-dir",   default=DEFAULT_BASIS_DIR)
    parser.add_argument("--freq",        type=float, default=10.)
    parser.add_argument("--workers",     type=int, default=None)
    parser.add_argument("--output",      default=None)
    parser.add_argument("--csv-ndamus",  default=None, help="Save neurodamus DataFrame as CSV")
    parser.add_argument("--csv-bcl",     default=None, help="Save bluecellulab DataFrame as CSV")
    parser.add_argument("--trace-figs",  action="store_true",
                        help="Generate per-pair voltage+rho trace figures in --figs-dir")
    parser.add_argument("--figs-dir",    default=None,
                        help="Output directory for per-pair trace figures (default: neurodamus_bluecellulab_figures/ next to --results-dir)")
    parser.add_argument("--no-stdp-curve", action="store_true",
                        help="Skip the STDP summary curve (useful with --trace-figs only)")
    args = parser.parse_args()

    output = args.output or f"stdp_compare_{int(args.freq)}Hz.png"

    figs_dir = args.figs_dir or os.path.join(
        os.path.dirname(os.path.abspath(args.results_dir)),
        "neurodamus_bluecellulab_figures",
    )

    print(f"Results dir : {args.results_dir}")
    print(f"Basis dir   : {args.basis_dir}")
    print(f"Freq        : {int(args.freq)} Hz")
    print(f"Figs dir    : {figs_dir}")

    if args.trace_figs:
        print(f"\n--- Generating per-pair trace figures → {figs_dir} ---")
        generate_trace_figures(args.results_dir, figs_dir,
                               freq=args.freq, n_workers=args.workers)

    # Always plot one example pair trace alongside the STDP curve
    example_pair = "180164-197248"
    example_dt   = f"{int(args.freq)}Hz_5ms"
    example_workdir = os.path.join(
        args.results_dir, "fitting", "n100", "seed19091997",
        "L5TTPC_L5TTPC_STDP", "simulations", example_pair, example_dt,
    )
    if os.path.exists(example_workdir):
        print(f"\n--- Example pair trace: {example_pair} {example_dt} ---")
        saved = plot_pair_traces(example_workdir, figs_dir)
        if saved:
            print(f"Saved trace figure: {saved}")

    if not args.no_stdp_curve:
        print("\n--- Neurodamus (out/) ---")
        df_ndamus = process_results(args.results_dir, args.basis_dir,
                                    "out", args.freq, args.workers)

        print("\n--- Bluecellulab (bluecellulab_results/) ---")
        df_bcl = process_results(args.results_dir, args.basis_dir,
                                 "bluecellulab_results", args.freq, args.workers)

        if args.csv_ndamus and not df_ndamus.empty:
            df_ndamus.to_csv(args.csv_ndamus, index=False)
            print(f"Saved {args.csv_ndamus}")
        if args.csv_bcl and not df_bcl.empty:
            df_bcl.to_csv(args.csv_bcl, index=False)
            print(f"Saved {args.csv_bcl}")

        plot_comparison(df_ndamus, df_bcl, output, freq=args.freq)

    print("Done.")


if __name__ == "__main__":
    main()
