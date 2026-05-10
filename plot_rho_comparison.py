#!/usr/bin/env python3
"""Plot rho traces from simulation_traces.pkl (CHINDEMI full) and rho.h5 (refitting) together."""
import pickle
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PKL_PATH = (
    "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/"
    "full_trace_results/CHINDEMI_PARAMS_custom_ratio_50_50/"
    "180164-197248/10Hz_-10ms/simulation_traces.pkl"
)
H5_PATH = (
    "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/"
    "refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/"
    "simulations/180164-197248/10Hz_-10ms/out/rho.h5"
)
OUT_PATH = "rho_comparison_180164-197248_10Hz_-10ms.png"


def load_pkl(path):
    with open(path, "rb") as f:
        d = pickle.load(f)
    t = np.array(d["t"]) / 1000.0          # ms -> s
    rho = np.array(d["rho_GB"], dtype=float)
    if rho.ndim == 1:
        rho = rho[:, None]
    # shape: (timepoints, n_synapses)
    syn_ids = list(d["synprop"]["synapseID"])
    return t, rho, syn_ids, d


def load_h5(path):
    with h5py.File(path, "r") as f:
        grp = f["report/S1nonbarrel_neurons"]
        t_info = grp["mapping/time"][:]       # [start_ms, stop_ms, dt_ms]
        elem_ids = grp["mapping/element_ids"][:]
        data = grp["data"][:]                 # (timepoints, n_synapses)
    t_start, t_stop, dt = t_info
    t = np.arange(data.shape[0]) * dt / 1000.0  # ms -> s
    return t, data, list(elem_ids.astype(int))


def main():
    t_pkl, rho_pkl, syn_ids_pkl, pkl_data = load_pkl(PKL_PATH)
    t_h5,  rho_h5,  syn_ids_h5            = load_h5(H5_PATH)

    print(f"PKL: {rho_pkl.shape}, synapses={syn_ids_pkl}, t=[{t_pkl[0]:.1f},{t_pkl[-1]:.1f}]s")
    print(f"H5 : {rho_h5.shape},  synapses={syn_ids_h5},  t=[{t_h5[0]:.1f},{t_h5[-1]:.1f}]s")

    # Match synapse columns
    common_ids = sorted(set(syn_ids_pkl) & set(syn_ids_h5))
    if not common_ids:
        raise RuntimeError("No common synapse IDs between pkl and h5.")

    n_syn = len(common_ids)
    fig, axes = plt.subplots(n_syn, 1, figsize=(14, 3.5 * n_syn), squeeze=False)
    fig.suptitle(
        "Rho traces comparison\n"
        "Blue: CHINDEMI_PARAMS_custom_ratio_50_50 (full sim, pkl)\n"
        "Red: refitting n100 seed19091997 (h5)\n"
        "Pair 180164→197248  |  10 Hz  -10 ms",
        fontsize=11,
    )

    pre_spikes  = np.array(pkl_data.get("prespikes",  [])) / 1000.0
    post_spikes = np.array(pkl_data.get("postspikes", [])) / 1000.0

    for row, sid in enumerate(common_ids):
        ax = axes[row, 0]
        col_pkl = syn_ids_pkl.index(float(sid))
        col_h5  = syn_ids_h5.index(sid)

        ax.plot(t_pkl, rho_pkl[:, col_pkl], color="steelblue", lw=1.2,
                label=f"pkl rho_GB  (syn {sid})", alpha=0.9)
        ax.plot(t_h5,  rho_h5[:, col_h5],  color="tomato",    lw=1.0,
                label=f"h5 rho  (syn {sid})", alpha=0.75, linestyle="--")

        ax.set_ylim(-0.05, 1.05)
        ax.axhline(0.5, color="gray", linestyle=":", alpha=0.5, lw=0.8)
        ax.set_ylabel(fr"$\rho$  (syn {sid})", fontsize=11)
        ax.legend(fontsize=9, loc="upper right")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        for sp in pre_spikes[:5]:   # first few to avoid clutter
            ax.axvline(sp, color="steelblue", alpha=0.15, lw=0.7)
        for sp in post_spikes[:5]:
            ax.axvline(sp, color="orange", alpha=0.15, lw=0.7)

    axes[-1, 0].set_xlabel("Time (s)", fontsize=11)
    plt.tight_layout()
    plt.savefig(OUT_PATH, dpi=200, bbox_inches="tight")
    print(f"Saved → {OUT_PATH}")


if __name__ == "__main__":
    main()
