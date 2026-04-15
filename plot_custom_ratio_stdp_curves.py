#!/usr/bin/env python3
"""
Overlay full-simulation STDP curves for all CHINDEMI custom-ratio result folders.

This reuses the direct EPSP-ratio extraction from plot_stdp_curves_full.py and
plots every CHINDEMI_PARAMS_custom_ratio_<depressed>_<potentiated> directory on
the same axes with a different color per ratio.

Example:
    python plot_custom_ratio_stdp_curves.py \
        --results-base-dir /project/ctb-emuller/dhuruva/plastyfire/full_trace_results \
        --workers 40 \
        --output stdp_curve_custom_ratio_stack.png
"""

import argparse
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

from plot_stdp_curves_full import FULL_TRACE_RESULTS_BASE, process_results


CUSTOM_RATIO_DIR_RE = re.compile(r"^CHINDEMI_PARAMS_custom_ratio_(\d+)_(\d+)$")


def discover_custom_ratio_dirs(results_base_dir):
    """Return sorted [(label, dir_path, depressed_pct, potentiated_pct), ...]."""
    ratio_dirs = []
    for entry in sorted(os.listdir(results_base_dir)):
        entry_path = os.path.join(results_base_dir, entry)
        match = CUSTOM_RATIO_DIR_RE.match(entry)
        if not match or not os.path.isdir(entry_path):
            continue
        depressed_pct = int(match.group(1))
        potentiated_pct = int(match.group(2))
        label = f"{depressed_pct}:{potentiated_pct}"
        ratio_dirs.append((label, entry_path, depressed_pct, potentiated_pct))

    ratio_dirs.sort(key=lambda item: (item[2], item[3], item[0]))
    return ratio_dirs


def summarize_curve(df):
    """Return mean/sem/count by delta-t for one curve."""
    summary = df.groupby("dt")["epsp_ratio"].agg(["mean", "sem", "count"]).reset_index()
    return summary.sort_values("dt")


def plot_custom_ratio_curves(curves, output_filename, title):
    if not curves:
        raise ValueError("No custom-ratio curves available to plot.")

    plt.figure(figsize=(7.5, 5.5))
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(curves)))

    for color, curve in zip(colors, curves):
        plt.errorbar(
            curve["summary"]["dt"],
            curve["summary"]["mean"],
            yerr=curve["summary"]["sem"],
            fmt="o-",
            color=color,
            linewidth=1.8,
            markersize=5,
            capsize=0,
            label=curve["label"],
        )

    invitro_dt = [-10, 5, 10]
    invitro_mean = [0.7922, 1.2038, 1.2013]
    invitro_sem = [0.0259, 0.0644, 0.0626]
    plt.errorbar(
        invitro_dt,
        invitro_mean,
        yerr=invitro_sem,
        fmt="o--",
        color="black",
        linewidth=1.6,
        markersize=5,
        capsize=0,
        label="in vitro",
        alpha=0.85,
    )

    plt.axhline(1.0, color="k", linestyle="--", alpha=0.4)
    plt.axvline(0.0, color="k", linestyle="--", alpha=0.4)
    plt.xlabel(r"$\Delta t$ (ms)", fontsize=12)
    plt.ylabel("EPSP ratio", fontsize=12)
    plt.title(title, fontsize=14)
    plt.legend(frameon=False, loc="best", title="Depressed:Potentiated")

    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_filename, dpi=300)
    print(f"Saved {output_filename}")


def main():
    parser = argparse.ArgumentParser(
        description="Overlay full-simulation STDP curves for all CHINDEMI custom-ratio folders."
    )
    parser.add_argument(
        "--results-base-dir",
        default=FULL_TRACE_RESULTS_BASE,
        help="Base directory containing CHINDEMI_PARAMS_custom_ratio_*_* folders.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Worker processes to use per custom-ratio folder.",
    )
    parser.add_argument(
        "--output",
        default="stdp_curve_custom_ratio_stack.png",
        help="Output PNG filename.",
    )
    parser.add_argument(
        "--title",
        default="Frequency = 10 Hz (Chindemi custom ratios) - full simulation",
        help="Plot title.",
    )
    args = parser.parse_args()

    ratio_dirs = discover_custom_ratio_dirs(args.results_base_dir)
    if not ratio_dirs:
        raise FileNotFoundError(
            f"No CHINDEMI custom-ratio result folders found in {args.results_base_dir}"
        )

    print("\n" + "=" * 72, flush=True)
    print("Stacked STDP curves for CHINDEMI custom ratios", flush=True)
    print(f"Base dir : {args.results_base_dir}", flush=True)
    print(f"Output   : {args.output}", flush=True)
    print(f"Ratios   : {', '.join(label for label, *_ in ratio_dirs)}", flush=True)
    print("=" * 72, flush=True)

    curves = []
    total_ratios = len(ratio_dirs)
    for idx, (label, ratio_dir, depressed_pct, potentiated_pct) in enumerate(ratio_dirs, start=1):
        print(
            f"\n[{idx}/{total_ratios}] Processing custom ratio {label} from {ratio_dir}",
            flush=True,
        )
        df = process_results(ratio_dir, n_workers=args.workers)
        if df.empty:
            print(f"Skipping {label}: no usable simulations found.", flush=True)
            continue

        summary = summarize_curve(df)
        print(
            f"Finished {label}: {df['pair'].nunique()} pairs, {len(df)} valid simulations",
            flush=True,
        )
        print(summary.to_string(index=False), flush=True)
        curves.append(
            {
                "label": label,
                "ratio_dir": ratio_dir,
                "depressed_pct": depressed_pct,
                "potentiated_pct": potentiated_pct,
                "summary": summary,
            }
        )

    if not curves:
        raise RuntimeError("All custom-ratio result folders were empty or unusable.")

    plot_custom_ratio_curves(curves, args.output, args.title)


if __name__ == "__main__":
    main()
