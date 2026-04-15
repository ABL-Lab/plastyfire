"""
Compare full vs prefire simulation traces — first synapse, stacked, aligned to induction onset.

Usage:
    python compare_prefire_vs_full.py
    python compare_prefire_vs_full.py --zoom-start -0.5 --zoom-end 5.0
"""

import argparse
import pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── defaults ─────────────────────────────────────────────────────────────────
FULL_PKL = (
    "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/"
    "trace_results/CHINDEMI_PARAMS/180164-197248/10Hz_5ms/simulation_traces.pkl"
)
PREFIRE_PKL = (
    "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/"
    "refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/"
    "simulations/180164-197248/10Hz_5ms/simulation_traces_prefire.pkl"
)
OUT_PNG = (
    "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/"
    "compare_prefire_vs_full_10Hz_5ms.png"
)

TRACES  = ["rho_GB", "effcai_GB", "cai_CR", "shaft_cai"]
YLABELS = ["ρ (rho_GB)", "effcai_GB (a.u.)", "cai_CR (mM)", "shaft_cai (mM)"]
COLORS  = {"full": "#2176AE", "prefire": "#E84855"}


# ── helpers ───────────────────────────────────────────────────────────────────
def load(path):
    with open(path, "rb") as f:
        return pickle.load(f)


def make_time(data, n_pts):
    """Uniformly-spaced time axis in seconds, length n_pts.
    Uses data['t'][-1] as the simulation end time (always the CVode clock).
    """
    t = data.get("t")
    t_end_s = (np.asarray(t)[-1] / 1000.0) if (t is not None and len(t) >= 2) \
              else n_pts * 0.025 / 1000.0
    return np.linspace(0.0, t_end_s, n_pts)


def first_synapse(data, key):
    """Return the trace for the first synapse as a 1-D array."""
    val = data.get(key)
    if val is None:
        return None
    if isinstance(val, dict):
        arr = np.asarray(next(iter(val.values())))
    else:
        arr = np.asarray(val)
        if arr.ndim == 2:
            arr = arr[:, 0]
    return arr


# ── main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--full",          default=FULL_PKL)
    parser.add_argument("--prefire",       default=PREFIRE_PKL)
    parser.add_argument("--out",           default=OUT_PNG)
    parser.add_argument("--full-onset",    type=float, default=241.0,
                        help="Induction onset in the full sim (s, default: 241.0)")
    parser.add_argument("--prefire-onset", type=float, default=1.0,
                        help="Induction onset in the prefire sim (s, default: 1.0)")
    parser.add_argument("--zoom-start",    type=float, default=-0.5,
                        help="Zoom start relative to induction onset (s, default: -0.5)")
    parser.add_argument("--zoom-end",      type=float, default=5.0,
                        help="Zoom end relative to induction onset (s, default: 5.0)")
    args = parser.parse_args()

    full    = load(args.full)
    prefire = load(args.prefire)

    print("Full    keys:", sorted(full.keys()))
    print("Prefire keys:", sorted(prefire.keys()))

    # Sources: (onset_s, label, data_dict, color)
    sources = [
        (args.full_onset,    "full",    full,    COLORS["full"]),
        (args.prefire_onset, "prefire", prefire, COLORS["prefire"]),
    ]

    # Drop traces not present in either file
    active_traces  = [k for k in TRACES
                      if any(s[2].get(k) is not None for s in sources)]
    active_ylabels = [YLABELS[TRACES.index(k)] for k in active_traces]

    n = len(active_traces)
    fig, axes = plt.subplots(n, 1, figsize=(12, 3 * n), sharex=True)
    if n == 1:
        axes = [axes]

    fig.suptitle(
        f"Full (blue) vs Prefire (red) — first synapse, aligned to induction onset\n"
        f"Pair 180164-197248 | 10 Hz, +5 ms | CHINDEMI_PARAMS | "
        f"t=0 → first pulse  |  zoom {args.zoom_start}–{args.zoom_end} s",
        fontsize=11, fontweight="bold"
    )

    for ax, key, ylabel in zip(axes, active_traces, active_ylabels):
        ax.set_ylabel(ylabel, fontsize=9)
        for onset, tag, data, color in sources:
            arr = first_synapse(data, key)
            if arr is None:
                continue
            t = make_time(data, len(arr)) - onset   # shift so t=0 = first pulse
            mask = (t >= args.zoom_start) & (t <= args.zoom_end)
            ax.plot(t[mask], arr[mask], color=color, lw=1.2, label=tag)
        ax.axvline(0, color="k", lw=0.8, ls="--", alpha=0.5, label="onset")
        ax.set_title(key, fontsize=9, loc="left")
        ax.legend(fontsize=8, loc="upper right")
        ax.tick_params(labelsize=8)
        ax.set_xlim(args.zoom_start, args.zoom_end)

    axes[-1].set_xlabel("time relative to induction onset (s)", fontsize=9)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"\nSaved → {args.out}")

    # ── rho summary ────────────────────────────────────────────────────────
    rho_f = first_synapse(full,    "rho_GB")
    rho_p = first_synapse(prefire, "rho_GB")
    if rho_f is not None and rho_p is not None:
        print(f"\nFinal rho_GB  full   : {rho_f[-1]:.4f}")
        print(f"Final rho_GB  prefire: {rho_p[-1]:.4f}")
        print(f"Abs diff             : {abs(rho_f[-1] - rho_p[-1]):.4f}")


if __name__ == "__main__":
    main()
