"""Figure for BASAL_NAKV.md: basal bAP vs distance, and effcai contrast ratios per variant."""
import sys, os
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import analyze as A
plt.style.use("/project/rrg-emuller/dhuruva/plastyfire/onerule.mplstyle")
df, fi, spk = A.load(); c = A.contrasts(df)
V = {"n1_k1_a1": "og-delta", "n1_k3_a1": "Kv x3", "n3_k1_a1": "Na x3", "n5_k1_a1": "Na x5", "n5_k1_aA": "Na x5 + Ka antic",
     "antic_delta": "antic-delta"}
COL = dict(zip(V, ["#222222", "#888888", "#1f77b4", "#d62728", "#ff7f0e", "#2ca02c"]))
m = spk.assign(ok=spk.n == spk.want).groupby("variant").ok.mean()
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw=dict(width_ratios=[1, 1.7]))
b = df[(df.proto == "s03_m120") & df.ok]
edges = np.array([0, 50, 100, 150, 200, 300]); mid = (edges[:-1] + edges[1:]) / 2
for v, lab in V.items():
    g = b[b.variant == v]; y = [g.bap[(g.dist >= lo) & (g.dist < hi)].median() for lo, hi in zip(edges[:-1], edges[1:])]
    ax[0].plot(mid, y, "o-", c=COL[v], label=f"{lab} ({m[v]:.0%} AP-matched)")
ax[0].axhline(30, ls=":", c="k"); ax[0].text(155, 24, "Nevian 2007 / Kampa:\n<=30 mV at >=200 um", fontsize=8)
ax[0].set_xlabel("basal path distance (um)"); ax[0].set_ylabel("bAP amplitude at synapse (mV)"); ax[0].legend(fontsize=7)
ax[0].set_title("A  single bAP at basal synapses (6 pairs)")
CN = ["sj0.1@+10 vs sj0.1@-10", "sj10@+10 vs sj10@-10", "sj20@+10 vs sj20@-10", "sj50@+10 vs sj50@-10",
      "sj50@-10 vs sj20@-10", "sj50@-10 vs sj0.1@-10", "s03_m25 vs s03_m120", "s03_b120 vs s03_m120"]
LB = ["+10/-10\n0.1 Hz", "+10/-10\n10 Hz", "+10/-10\n20 Hz", "+10/-10\n50 Hz", "-10:\n50/20 Hz", "-10:\n50/0.1 Hz", "-25/-120", "burst/single\n-120"]
w = 0.8 / len(V)
for i, (v, lab) in enumerate(V.items()):
    y = [c[(c.variant == v) & (c.contrast == k)].effcai_pk_ratio.squeeze() if ((c.variant == v) & (c.contrast == k)).any() else np.nan for k in CN]
    ax[1].bar(np.arange(len(CN)) + (i - len(V) / 2 + 0.5) * w, np.array(y, float), w, color=COL[v], label=lab)
ax[1].axhline(1, c="k", lw=0.8); ax[1].set_yscale("log"); ax[1].legend(fontsize=7); ax[1].set_xticks(range(len(CN)), LB, fontsize=7)
ax[1].set_ylabel("median per-synapse effcai peak ratio"); ax[1].set_title("B  protocol contrasts (spike-matched basal synapses)")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "..", "figs", "basal_nakv.png"), dpi=150)
