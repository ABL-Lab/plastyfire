#!/usr/bin/env python3
"""
Plot a histogram of pairwise EPSP values for all pairs in basis_results/,
using the initial rho0_GB values from the modified edges.h5 and linear
superposition over the per-pair basis CSVs.

Usage:
    python plot_basis_epsp_histogram.py [--output epsp_histogram.png]
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

# ── Paths ────────────────────────────────────────────────────────────────────
REPO_DIR    = os.path.dirname(os.path.abspath(__file__))
BASIS_DIR   = os.path.join(REPO_DIR, "basis_results")
EDGES_H5    = "/lustre06/project/6077694/dhuruva/plastyfire/data/dhuruva_modified_edges.h5"
EDGE_POP    = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"
POP_BASE    = f"edges/{EDGE_POP}"
POP_PATH    = f"edges/{EDGE_POP}/0"
EXC_TYPE_MIN = 100

sys.path.insert(0, os.path.join(REPO_DIR, "new_fitting"))
from effcai_to_epsp import load_basis, extrapolate_epsp


# ── Edge lookup ──────────────────────────────────────────────────────────────

def get_rho0_for_pair(h5file, pre_gid, post_gid):
    """Return rho0_GB values for EXC synapses from pre_gid → post_gid."""
    n2r = h5file[f"{POP_BASE}/indices/target_to_source/node_id_to_ranges"]
    r2e = h5file[f"{POP_BASE}/indices/target_to_source/range_to_edge_id"]

    r0, r1 = int(n2r[post_gid][0]), int(n2r[post_gid][1])
    if r0 == r1:
        return None

    edge_ids = []
    for e_start, e_end in r2e[r0:r1]:
        edge_ids.extend(range(int(e_start), int(e_end)))
    edge_ids = np.array(edge_ids, dtype=np.int64)

    src   = h5file[f"{POP_BASE}/source_node_id"][edge_ids]
    stype = h5file[f"{POP_PATH}/syn_type_id"][edge_ids]

    mask = (src == pre_gid) & (stype >= EXC_TYPE_MIN)
    if not mask.any():
        return None

    syn_ids = edge_ids[mask]
    rho0 = h5file[f"{POP_PATH}/rho0_GB"][syn_ids]
    return rho0.astype(int).tolist()


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Plot pairwise EPSP histogram from basis + edges.h5 rho0")
    parser.add_argument("--output", default="epsp_initial_histogram.png")
    parser.add_argument("--basis-dir", default=BASIS_DIR)
    parser.add_argument("--edges-h5", default=EDGES_H5)
    args = parser.parse_args()

    # Parse pairs from basis filenames
    pairs = []
    for fname in sorted(os.listdir(args.basis_dir)):
        if not fname.startswith("basis_") or not fname.endswith(".csv"):
            continue
        stem = fname[len("basis_"):-len(".csv")]
        parts = stem.split("_")
        if len(parts) == 2:
            try:
                pairs.append((int(parts[0]), int(parts[1])))
            except ValueError:
                continue

    print(f"Found {len(pairs)} pairs in {args.basis_dir}")

    epsp_values = []
    skipped = []

    with h5py.File(args.edges_h5, "r") as hf:
        for pre_gid, post_gid in pairs:
            rho0 = get_rho0_for_pair(hf, pre_gid, post_gid)
            if rho0 is None:
                skipped.append((pre_gid, post_gid, "no edges found"))
                continue

            try:
                baseline_mean, baseline_std, singleton_means, singleton_stds, n_syn = \
                    load_basis(pre_gid, post_gid, basis_dir=args.basis_dir)
            except (FileNotFoundError, ValueError) as exc:
                skipped.append((pre_gid, post_gid, str(exc)))
                continue

            if len(rho0) != n_syn:
                skipped.append((pre_gid, post_gid,
                                 f"rho0 length {len(rho0)} != basis n_syn {n_syn}"))
                continue

            epsp_mean, epsp_std, _ = extrapolate_epsp(
                rho0, baseline_mean, baseline_std,
                singleton_means, singleton_stds, n_trials=0
            )
            epsp_values.append(epsp_mean)

    if skipped:
        print(f"Skipped {len(skipped)} pairs:")
        for pre, post, reason in skipped:
            print(f"  {pre}->{post}: {reason}")

    epsp_values = np.array(epsp_values)
    n = len(epsp_values)
    mean_val = epsp_values.mean()
    std_val  = epsp_values.std()

    print(f"\nEPSP summary: n={n}, mean={mean_val:.4f} mV, std={std_val:.4f} mV")

    # ── Plot ─────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 5))

    n_bins = min(20, max(5, n // 4))
    ax.hist(epsp_values, bins=n_bins, color="steelblue", edgecolor="white", linewidth=0.5)

    annotation = (
        f"n = {n} pairs\n"
        f"mean = {mean_val:.3f} mV\n"
        f"std  = {std_val:.3f} mV"
    )
    ax.text(
        0.97, 0.95, annotation,
        transform=ax.transAxes,
        ha="right", va="top",
        fontsize=11,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="gray", alpha=0.9),
    )

    ax.set_xlabel("Initial EPSP amplitude (mV)", fontsize=12)
    ax.set_ylabel("Number of pairs", fontsize=12)
    ax.set_title("Pairwise initial EPSP distribution", fontsize=12)
    ax.axvline(mean_val, color="crimson", linestyle="--", linewidth=1.5, label=f"mean = {mean_val:.3f} mV")
    ax.legend(fontsize=10)

    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
