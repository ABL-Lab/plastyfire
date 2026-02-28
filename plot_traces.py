#!/usr/bin/env python3
"""
Plot synaptic traces (effcai_GB, cai_CR, vsyn, rho_GB) from L5TTPC trace simulations.

Usage:
    # Plot a single simulation
    python plot_traces.py --pkl trace_results/Chindemi_params/180164-197248/2Hz_-10ms/simulation_traces.pkl
    
    # Plot all simulations for a pair
    python plot_traces.py --pair 180164-197248
    
    # Plot first N synapses only
    python plot_traces.py --pkl <path> --max-synapses 5

Author: Generated for trace visualization
"""

import argparse
import glob
import os
import pickle
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Output directory
RESULTS_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/trace_results/DHURUVA_PARAMS_V7"
FIGURES_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/figures/plastyfire/traces"


def plot_single_simulation(pkl_path, output_dir=None, max_synapses=10, show_legend=True):
    """Plot traces from a single simulation pickle file.
    
    Args:
        pkl_path: Path to simulation_traces.pkl
        output_dir: Output directory for figures
        max_synapses: Max number of synapses to plot (for clarity)
        show_legend: Whether to show legend
    """
    # Load data
    with open(pkl_path, "rb") as f:
        data = pickle.load(f)
    
    # Extract path info for naming
    path_parts = Path(pkl_path).parts
    pair_name = path_parts[-3]  # e.g., 180164-197248
    freq_dt = path_parts[-2]    # e.g., 2Hz_-10ms
    
    # Output directory
    if output_dir is None:
        output_dir = FIGURES_DIR
    os.makedirs(output_dir, exist_ok=True)
    
    # Get time vector
    t = data.get("t")

    # Helper to ensure arrays are (n_synapses, n_timepoints)
    def ensure_synapse_first(arr, t_len):
        """Transpose array if needed so shape is (n_synapses, n_timepoints)."""
        if arr is None:
            return None
        arr = np.asarray(arr)
        if arr.ndim == 1:
            return arr.reshape(1, -1)
        # If first dimension matches time length, transpose
        if arr.shape[0] == t_len and arr.shape[1] != t_len:
            return arr.T
        return arr

    # Get rho and determine time length
    rho = data.get("rho_GB")
    if rho is not None:
        rho = np.asarray(rho)
        # Determine time length from t or from rho shape
        if t is None:
            t_len = max(rho.shape)
            t = np.arange(t_len)
        else:
            t_len = len(t)
        rho = ensure_synapse_first(rho, t_len)
    elif t is not None:
        t_len = len(t)
    else:
        print(f"No time or rho data found in {pkl_path}")
        return

    # Convert time to seconds for readability
    t_sec = t / 1000.0

    # Number of synapses
    n_synapses = min(max_synapses, rho.shape[0] if rho is not None else 0)
    if n_synapses == 0:
        print(f"No synaptic data found in {pkl_path}")
        return
    
    # Create figure with subplots
    fig, axes = plt.subplots(5, 1, figsize=(14, 12), sharex=True)
    fig.suptitle(f"{pair_name} / {freq_dt}\nChindemi Parameters", fontsize=14, fontweight="bold")
    
    # Color map for synapses
    colors = plt.cm.tab10(np.linspace(0, 1, n_synapses))
    
    # 1. Soma voltage
    ax = axes[0]
    v = data.get("v")
    if v is not None:
        ax.plot(t_sec, v, "k-", linewidth=0.8)
    ax.set_ylabel("V$_{soma}$ (mV)")
    ax.set_title("Somatic Membrane Potential")
    ax.grid(alpha=0.3)
    
    # Mark pre and post spikes
    if "prespikes" in data:
        for spike in data["prespikes"][:50]:  # Limit markers
            ax.axvline(spike/1000, color="blue", alpha=0.1, linewidth=0.5)
    if "postspikes" in data:
        for spike in data["postspikes"]:
            ax.axvline(spike/1000, color="red", alpha=0.3, linewidth=1)
    
    # 2. rho_GB (plasticity state)
    ax = axes[1]
    if rho is not None:
        for i in range(n_synapses):
            ax.plot(t_sec, rho[i], color=colors[i], linewidth=0.8, alpha=0.7, label=f"Syn {i}")
    ax.set_ylabel("ρ$_{GB}$")
    ax.set_ylim(-0.1, 1.1)
    ax.set_title("Plasticity State (0=depressed, 1=potentiated)")
    ax.grid(alpha=0.3)
    if show_legend and n_synapses <= 10:
        ax.legend(loc="upper right", fontsize=7, ncol=2)
    
    # 3. effcai_GB (effective calcium)
    ax = axes[2]
    effcai = data.get("effcai_GB")
    if effcai is not None:
        effcai = ensure_synapse_first(np.asarray(effcai), t_len)
        for i in range(min(n_synapses, effcai.shape[0])):
            ax.plot(t_sec, effcai[i], color=colors[i], linewidth=0.8, alpha=0.7)
    ax.set_ylabel("effCa$_i$ (μM)")
    ax.set_title("Effective Intracellular Calcium")
    ax.grid(alpha=0.3)

    # 4. cai_CR (calcium concentration)
    ax = axes[3]
    cai = data.get("cai_CR")
    if cai is not None:
        cai = ensure_synapse_first(np.asarray(cai), t_len)
        for i in range(min(n_synapses, cai.shape[0])):
            ax.plot(t_sec, cai[i], color=colors[i], linewidth=0.8, alpha=0.7)
    ax.set_ylabel("Ca$_i$ (μM)")
    ax.set_title("Calcium Concentration")
    ax.grid(alpha=0.3)

    # 5. vsyn (synaptic voltage)
    ax = axes[4]
    vsyn = data.get("vsyn")
    if vsyn is not None:
        vsyn = ensure_synapse_first(np.asarray(vsyn), t_len)
        for i in range(min(n_synapses, vsyn.shape[0])):
            ax.plot(t_sec, vsyn[i], color=colors[i], linewidth=0.8, alpha=0.7)
    ax.set_ylabel("V$_{syn}$ (mV)")
    ax.set_xlabel("Time (s)")
    ax.set_title("Synaptic Voltage")
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    fig_name = f"{pair_name}_{freq_dt}_traces.png"
    fig_path = os.path.join(output_dir, fig_name)
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    
    print(f"Saved: {fig_path}")
    return fig_path


def plot_rho_summary(pkl_path, output_dir=None):
    """Plot summary of rho evolution for all synapses."""
    with open(pkl_path, "rb") as f:
        data = pickle.load(f)

    path_parts = Path(pkl_path).parts
    pair_name = path_parts[-3]
    freq_dt = path_parts[-2]

    if output_dir is None:
        output_dir = FIGURES_DIR
    os.makedirs(output_dir, exist_ok=True)

    rho = data.get("rho_GB")
    if rho is None:
        return

    rho = np.asarray(rho)
    t = data.get("t")

    # Determine time length and ensure rho is (n_synapses, n_timepoints)
    if t is not None:
        t_len = len(t)
    else:
        t_len = max(rho.shape)
        t = np.arange(t_len)

    # Transpose if needed
    if rho.ndim == 2 and rho.shape[0] == t_len and rho.shape[1] != t_len:
        rho = rho.T

    t_sec = t / 1000.0
    
    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f"{pair_name} / {freq_dt} - Plasticity Summary", fontsize=12)
    
    # Left: All rho traces
    ax = axes[0]
    for i, trace in enumerate(rho):
        ax.plot(t_sec, trace, alpha=0.3, linewidth=0.5)
    ax.plot(t_sec, np.mean(rho, axis=0), "k-", linewidth=2, label="Mean")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("ρ$_{GB}$")
    ax.set_ylim(-0.1, 1.1)
    ax.set_title(f"All {len(rho)} Synapses")
    ax.legend()
    ax.grid(alpha=0.3)
    
    # Right: Initial vs Final distribution
    ax = axes[1]
    initial_rho = [r[0] for r in rho]
    final_rho = [r[-1] for r in rho]
    
    ax.hist(initial_rho, bins=20, alpha=0.5, label=f"Initial (mean={np.mean(initial_rho):.2f})", color="blue")
    ax.hist(final_rho, bins=20, alpha=0.5, label=f"Final (mean={np.mean(final_rho):.2f})", color="red")
    ax.set_xlabel("ρ$_{GB}$")
    ax.set_ylabel("Count")
    ax.set_title("Plasticity State Distribution")
    ax.legend()
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    
    fig_name = f"{pair_name}_{freq_dt}_rho_summary.png"
    fig_path = os.path.join(output_dir, fig_name)
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    
    print(f"Saved: {fig_path}")
    return fig_path


def plot_all_for_pair(pair_name, results_dir=None, output_dir=None, max_synapses=10):
    """Plot all simulations for a given pair."""
    if results_dir is None:
        results_dir = RESULTS_DIR
    
    pair_dir = os.path.join(results_dir, pair_name)
    if not os.path.isdir(pair_dir):
        print(f"Pair directory not found: {pair_dir}")
        return
    
    pkl_files = glob.glob(os.path.join(pair_dir, "*/simulation_traces.pkl"))
    print(f"Found {len(pkl_files)} simulations for pair {pair_name}")
    
    for pkl_path in sorted(pkl_files):
        try:
            plot_single_simulation(pkl_path, output_dir, max_synapses)
            plot_rho_summary(pkl_path, output_dir)
        except Exception as e:
            print(f"Error plotting {pkl_path}: {e}")


def main():
    parser = argparse.ArgumentParser(description="Plot synaptic traces from L5TTPC simulations")
    
    parser.add_argument("--pkl", type=str, help="Path to a single simulation_traces.pkl file")
    parser.add_argument("--pair", type=str, help="Pair name (e.g., 180164-197248) to plot all its simulations")
    parser.add_argument("--all", action="store_true", help="Plot all available simulations")
    parser.add_argument("--max-synapses", type=int, default=10, help="Max synapses to plot (default: 10)")
    parser.add_argument("--output-dir", type=str, default=None, help="Output directory for figures")
    parser.add_argument("--results-dir", type=str, default=RESULTS_DIR, help="Results directory")
    
    args = parser.parse_args()
    
    if args.pkl:
        plot_single_simulation(args.pkl, args.output_dir, args.max_synapses)
        plot_rho_summary(args.pkl, args.output_dir)
    elif args.pair:
        plot_all_for_pair(args.pair, args.results_dir, args.output_dir, args.max_synapses)
    elif args.all:
        # Find all pairs
        pairs = [d for d in os.listdir(args.results_dir) 
                 if os.path.isdir(os.path.join(args.results_dir, d))]
        print(f"Plotting {len(pairs)} pairs...")
        for pair in sorted(pairs):
            if pair == "logs":
                continue
            plot_all_for_pair(pair, args.results_dir, args.output_dir, args.max_synapses)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
