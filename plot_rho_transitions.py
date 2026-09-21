"""
Per-synapse rho transition census, split by apical / basal.

For every pair x protocol result pkl, classify each synapse by what its rho did:
    DEP->POT : initial rho < 0.5, final >= 0.5   (depressed -> potentiated)
    POT->DEP : initial rho >= 0.5, final < 0.5   (potentiated -> depressed)
    stayed   : no crossing of the 0.5 boundary

The 0.5 boundary is the same binarisation compute_epsp_from_basis applies when
turning rho into an EPSP (rho >= 0.5 -> 1 else 0), so these counts are exactly
the transitions the STDP curve can actually see. Because that binarisation
discards sub-threshold movement, the plot also reports mean |final - initial|
per group, which does not threshold.

Section type comes from edges.h5 afferent_section_type (3 = apical), matching
simulator_edges._apply_theta_from_a_params.

Usage:
    python plot_rho_transitions.py --param-hash b8c7ff3ecf0a \
        --label "DE fit2" --out rho_transitions_de_fit2.png

`--param-hash` selects which simulation_edges_<hash>.pkl set to read; it defaults
to 0ce64fa83b85 (the original ic run) and will silently find nothing if the hash
does not match any pkl on disk, so check the "Missing/unreadable pkls" line.
"""

import argparse
import glob
import os
import pickle
from collections import defaultdict

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SIMS_DIR = os.path.join(
    SCRIPT_DIR,
    "refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations")
EDGES_H5 = os.path.join(SCRIPT_DIR, "data/dhuruva_modified_edges.h5")
EDGE_POP = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"
PARAM_HASH = "0ce64fa83b85"   # overridden by --param-hash

# dt order matching the STDP curve
PROTOCOLS = [("10Hz_-50ms", -50), ("10Hz_-30ms", -30), ("10Hz_-10ms", -10),
             ("10Hz_5ms", 5), ("10Hz_10ms", 10), ("10Hz_30ms", 30),
             ("10Hz_50ms", 50)]

RHO_THRESH = 0.5


def load_section_types():
    with h5py.File(EDGES_H5, "r") as f:
        return f[f"edges/{EDGE_POP}/0/afferent_section_type"][:]


def main():
    global PARAM_HASH
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="rho_transitions.png")
    ap.add_argument("--param-hash", default=PARAM_HASH,
                    help=f"simulation_edges_<hash>.pkl to read (default: {PARAM_HASH})")
    ap.add_argument("--label", default="IC gen-4 best",
                    help="text appended to the figure title")
    args = ap.parse_args()
    PARAM_HASH = args.param_hash
    print(f"param_hash : {PARAM_HASH}")

    print("Loading section types from edges.h5 …")
    sec_types = load_section_types()

    # counts[dt][loc][category] and magnitude accumulator
    counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    deltas = defaultdict(lambda: defaultdict(list))
    missing = 0

    for pair in sorted(os.listdir(SIMS_DIR)):
        pair_dir = os.path.join(SIMS_DIR, pair)
        if not os.path.isdir(pair_dir):
            continue
        for proto, dt in PROTOCOLS:
            pkl = os.path.join(pair_dir, proto, f"simulation_edges_{PARAM_HASH}.pkl")
            if not os.path.exists(pkl):
                missing += 1
                continue
            try:
                with open(pkl, "rb") as f:
                    d = pickle.load(f)
            except Exception:
                missing += 1
                continue
            gids = d["global_ids"]
            r0 = np.asarray(d["initial_rho"], dtype=float)
            r1 = np.asarray(d["final_rho"], dtype=float)
            for gid, a, b in zip(gids, r0, r1):
                loc = "apical" if sec_types[gid] == 3 else "basal"
                was_pot = a >= RHO_THRESH
                is_pot = b >= RHO_THRESH
                if not was_pot and is_pot:
                    cat = "DEP->POT"
                elif was_pot and not is_pot:
                    cat = "POT->DEP"
                else:
                    cat = "stayed"
                counts[dt][loc][cat] += 1
                deltas[dt][loc].append(b - a)

    print(f"Missing/unreadable pkls: {missing}")

    dts = [dt for _, dt in PROTOCOLS]
    cats = ["DEP->POT", "POT->DEP", "stayed"]
    colors = {"DEP->POT": "#2e7d32", "POT->DEP": "#c62828", "stayed": "#9e9e9e"}

    fig, axes = plt.subplots(2, 2, figsize=(15, 9))

    for col, loc in enumerate(["apical", "basal"]):
        # ── top row: percentage stacked bars ──
        ax = axes[0][col]
        totals = [sum(counts[dt][loc].values()) for dt in dts]
        bottom = np.zeros(len(dts))
        for cat in cats:
            vals = np.array([
                100.0 * counts[dt][loc][cat] / t if t else 0.0
                for dt, t in zip(dts, totals)])
            ax.bar([str(d) for d in dts], vals, bottom=bottom,
                   label=cat, color=colors[cat], edgecolor="white", linewidth=0.6)
            for x, (v, b) in enumerate(zip(vals, bottom)):
                if v >= 3.0:
                    ax.text(x, b + v / 2, f"{v:.0f}%", ha="center", va="center",
                            fontsize=8, color="white", fontweight="bold")
            bottom += vals
        n_syn = totals[0] if totals else 0
        ax.set_title(f"{loc}  (n≈{n_syn} synapses per Δt)", fontsize=12)
        ax.set_ylabel("% of synapses")
        ax.set_ylim(0, 100)
        ax.axvline(2.5, color="k", ls="--", lw=0.8, alpha=0.5)
        if col == 0:
            ax.legend(loc="lower left", fontsize=9)

        # ── bottom row: mean rho change (no thresholding) ──
        ax2 = axes[1][col]
        means = [np.mean(deltas[dt][loc]) if deltas[dt][loc] else 0.0 for dt in dts]
        sems = [np.std(deltas[dt][loc]) / np.sqrt(len(deltas[dt][loc]))
                if deltas[dt][loc] else 0.0 for dt in dts]
        bar_colors = ["#2e7d32" if m > 0 else "#c62828" for m in means]
        ax2.bar([str(d) for d in dts], means, yerr=sems, capsize=3,
                color=bar_colors, alpha=0.85)
        ax2.axhline(0, color="k", lw=0.8)
        ax2.axvline(2.5, color="k", ls="--", lw=0.8, alpha=0.5)
        ax2.set_xlabel("Δt (ms)")
        ax2.set_ylabel("mean Δρ (final − initial)")
        ax2.set_title(f"{loc} — mean rho change (unthresholded)", fontsize=11)

    fig.suptitle(f"Per-synapse rho transitions across the 0.5 boundary — {args.label}",
                 fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out_path = os.path.join(SCRIPT_DIR, args.out)
    plt.savefig(out_path, dpi=200)
    print(f"Saved → {out_path}")

    # ── text summary ──
    print("\nTransition census (counts):")
    hdr = f"{'dt':>5} {'loc':>7} {'n':>6} {'DEP->POT':>9} {'POT->DEP':>9} {'stayed':>8} {'meanDrho':>9}"
    print(hdr)
    print("-" * len(hdr))
    for dt in dts:
        for loc in ["apical", "basal"]:
            c = counts[dt][loc]
            t = sum(c.values())
            md = np.mean(deltas[dt][loc]) if deltas[dt][loc] else 0.0
            print(f"{dt:>5} {loc:>7} {t:>6} {c['DEP->POT']:>9} "
                  f"{c['POT->DEP']:>9} {c['stayed']:>8} {md:>+9.4f}")


if __name__ == "__main__":
    main()
