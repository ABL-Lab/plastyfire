#!/usr/bin/env python3
"""Plot effcai traces for a single synapse comparing two protocols."""
import pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TRACE_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/trace_results/Chindemi_params"
PAIR = "180164-197248"

protos = ["10Hz_10ms", "10Hz_-10ms"]
colors = ["tab:red", "tab:blue"]
SYN_IDX = 0  # first synapse

fig, ax = plt.subplots(figsize=(14, 5))

for proto, color in zip(protos, colors):
    pkl_path = f"{TRACE_DIR}/{PAIR}/{proto}/simulation_traces.pkl"
    with open(pkl_path, "rb") as f:
        data = pickle.load(f)

    t = np.asarray(data["t"])
    effcai = np.asarray(data["cai_CR"])
    if effcai.ndim == 1:
        effcai = effcai.reshape(1, -1)
    elif effcai.shape[0] == len(t) and effcai.shape[1] != len(t):
        effcai = effcai.T

    print(f"{proto}: effcai shape={effcai.shape}, t range=[{t[0]:.1f}, {t[-1]:.1f}]ms, "
          f"n_syn={effcai.shape[0]}")
    print(f"  syn {SYN_IDX}: peak={effcai[SYN_IDX].max():.4f}, mean={effcai[SYN_IDX].mean():.6f}")

    ax.plot(t / 1000, effcai[SYN_IDX], label=proto, color=color, alpha=0.8, linewidth=0.5)

ax.set_xlabel("Time (s)")
ax.set_ylabel("effcai (effective calcium)")
ax.set_title(f"effcai trace comparison - Pair {PAIR}, Synapse {SYN_IDX}")
ax.legend()
ax.grid(True, alpha=0.3)
plt.axis([0,4,0,0.004])
out = "cai_compare_m_m.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"\nSaved to {out}")
