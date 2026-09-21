"""
Plot fitness evolution + mean EPSP ratio vs in vitro target from checkpoint_edges.pkl.

Usage:
    python plot_edges_fitness.py [--checkpoint checkpoint_edges.pkl] [--out fitness.png]
"""
import argparse
import hashlib
import os
import pickle
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CSVF_NAME = os.path.join(os.path.dirname(__file__),
                         "biodata/paired_recordings.csv")
PROTOCOL_IDX = ["mrk97_07", "mrk97_08"]
CACHE_DIR = ".cache"


def load_best_resdb(best_individual):
    """Load the resdb (simulated EPSP ratios) for the best individual from .cache."""
    cachekey = hashlib.md5(str(list(best_individual)).encode()).hexdigest()
    pkl_path = os.path.join(CACHE_DIR, f"{cachekey}.pkl")
    if not os.path.exists(pkl_path):
        return None
    with open(pkl_path, "rb") as f:
        data = pickle.load(f)
    return data.get("resdb")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", default="checkpoint_edges.pkl")
    parser.add_argument("--out", default="fitness_evolution.png")
    args = parser.parse_args()

    try:
        with open(args.checkpoint, "rb") as f:
            cp = pickle.load(f)
    except FileNotFoundError:
        print(f"ERROR: {args.checkpoint} not found — optimizer hasn't completed a generation yet.")
        sys.exit(1)

    gen_reached = cp.get("generation", "?")
    print(f"Generations completed : {gen_reached}")

    lb = cp.get("logbook")
    if lb is None:
        print("ERROR: no logbook in checkpoint.")
        sys.exit(1)

    gen  = lb.select("gen")
    minv = np.array(lb.select("min"))

    print(f"Generations in logbook: {len(gen)}")
    print(f"Best (min) fitness    : {minv.min():.4f}  at gen {gen[np.argmin(minv)]}")

    # Best individual
    hof = cp.get("halloffame", [])
    best = hof[0] if hof else None
    # Simulated EPSP ratios for best individual (loaded early so we can print below)
    resdb = load_best_resdb(best) if best is not None else None

    if best:
        from plastyfire.evaluator_edges import EDGES_FIT_PARAM_NAMES, FITTED_TAU
        print("\nBest individual (hall of fame):")
        for name, val in zip(EDGES_FIT_PARAM_NAMES, best):
            print(f"  {name:<35} = {val:.6f}")
        print(f"  {'tau_effca_GB_GluSynapse':<35} = {FITTED_TAU}  (fixed)")
        print(f"  fitness values: {list(best.fitness.values)}")
        if resdb is not None:
            for pid in PROTOCOL_IDX:
                vals = resdb[resdb["protocol_id"] == pid]["epsp_ratio"].values
                n_valid = int(np.sum(~np.isnan(vals)))
                print(f"  {pid}: min={np.nanmin(vals):.4f}  max={np.nanmax(vals):.4f}  "
                      f"mean={np.nanmean(vals):.4f} ± {np.nanstd(vals):.4f}  (n={n_valid}/{len(vals)})")
        else:
            print("  (cache not found for hof individual)")

    # Best individual in last generation (min sum of fitness values)
    population = cp.get("population", [])
    if population:
        best_last = min(population, key=lambda ind: sum(ind.fitness.values))
        print(f"\nBest individual in last generation (gen {gen_reached}):")
        print(f"  fitness values: {list(best_last.fitness.values)}")

    # In vitro targets
    invitro = pd.read_csv(CSVF_NAME)
    invitro = invitro[invitro["protocol_id"].isin(PROTOCOL_IDX)].set_index("protocol_id")

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Left: min fitness over generations
    ax = axes[0]
    ax.plot(gen, minv, "g-o", ms=5, lw=2)
    ax.set_xlabel("Generation")
    ax.set_ylabel("Fitness (lower = better)")
    ax.set_title("Best fitness per generation")
    ax.grid(True, alpha=0.3)

    # Right: mean simulated EPSP ratio vs in vitro per protocol
    ax2 = axes[1]
    protocol_labels = []
    sim_means, sim_sems = [], []
    tgt_means, tgt_sems = [], []

    for pid in PROTOCOL_IDX:
        row = invitro.loc[pid]
        dt = int(row["dt_train"])
        label = f"{pid}\n(dt={dt:+d}ms)"
        protocol_labels.append(label)
        tgt_means.append(row["mean_epsp_ratio"])
        tgt_sems.append(row["sem_epsp_ratio"])

        if resdb is not None:
            vals = resdb[resdb["protocol_id"] == pid]["epsp_ratio"].values
            n_valid = int(np.sum(~np.isnan(vals)))
            sim_means.append(np.nanmean(vals))
            sim_sems.append(np.nanstd(vals) / np.sqrt(n_valid) if n_valid > 0 else 0)
            print(f"  {pid}: sim mean={np.nanmean(vals):.4f} ± {np.nanstd(vals):.4f}  "
                  f"target={row['mean_epsp_ratio']:.4f} ± {row['sem_epsp_ratio']:.4f}  (n={n_valid}/{len(vals)})")
        else:
            sim_means.append(None)
            sim_sems.append(None)

    x = np.arange(len(PROTOCOL_IDX))
    width = 0.35

    bars_tgt = ax2.bar(x - width/2, tgt_means, width, yerr=tgt_sems,
                       capsize=5, color="steelblue", alpha=0.8, label="In vitro target")
    if all(v is not None for v in sim_means):
        bars_sim = ax2.bar(x + width/2, sim_means, width, yerr=sim_sems,
                           capsize=5, color="tomato", alpha=0.8, label="Sim (best params)")

    ax2.axhline(1.0, color="gray", ls="--", lw=1, alpha=0.6)
    ax2.set_xticks(x)
    ax2.set_xticklabels(protocol_labels)
    ax2.set_ylabel("Mean EPSP ratio")
    ax2.set_title(f"Simulated vs in vitro EPSP ratio\n(best individual, gen {gen_reached})")
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    plt.savefig(args.out, dpi=150)
    print(f"\nSaved → {args.out}")


if __name__ == "__main__":
    main()
