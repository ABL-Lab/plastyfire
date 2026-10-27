"""Unitary EPSP per pair, base SSCx (original emodels + edges) vs og-delta, L23->L5 and L5->L5; same readout
(pair_geometry.py / epsp_base_sscx.py: 30 spikes at 1 Hz, mean single-sweep peak). -> ebner/figs/epsp_base_vs_delta.png"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
F = "/home/dhuruva/ws/plastyfire/ebner/"
k = ["pregid", "postgid"]
sets = {"L23->L5": (pd.read_csv(F + "figs/epsp_base_l23l5.csv"), pd.read_csv(F + "pair_geometry_L23PC_L5TTPC.csv"), .42),
        "L5->L5": (pd.read_csv(F + "figs/epsp_base_l5l5.csv"), pd.read_csv(F + "figs/epsp_delta_l5l5.csv"), 1.3)}
fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
bins = np.linspace(0, 4, 81)
for i, (name, (b, d, lit)) in enumerate(sets.items()):
    m = b[k + ["epsp_mV"]].merge(d[k + ["epsp_mV"]], on=k, suffixes=("_base", "_delta")).dropna()
    for c, col, lab in [("epsp_mV_base", "tab:gray", "base SSCx"), ("epsp_mV_delta", "tab:blue", "og-delta")]:
        ax[i].hist(m[c].clip(upper=3.99), bins=bins, alpha=.5, color=col,
                   label="%s (n=%d, mean %.2f, median %.2f)" % (lab, len(m), m[c].mean(), m[c].median()))
    ax[i].axvline(lit, color="k", ls="--", label="literature %.2f" % lit)
    ax[i].set(xlabel="unitary EPSP per pair (mV)", ylabel="pairs", title=name); ax[i].legend(fontsize=8)
    ax[2].scatter(m.epsp_mV_base, m.epsp_mV_delta, s=10, label="%s r=%.2f" % (name, m.epsp_mV_base.corr(m.epsp_mV_delta)))
    print(name, "n=%d base mean %.3f median %.3f | delta mean %.3f median %.3f | delta/base median %.2f"
          % (len(m), m.epsp_mV_base.mean(), m.epsp_mV_base.median(), m.epsp_mV_delta.mean(), m.epsp_mV_delta.median(),
             (m.epsp_mV_delta / m.epsp_mV_base).median()))
ax[2].plot([0, 6], [0, 6], "k:"); ax[2].set(xlabel="base SSCx EPSP (mV)", ylabel="og-delta EPSP (mV)", title="per pair"); ax[2].legend()
fig.tight_layout(); fig.savefig(F + "figs/epsp_base_vs_delta.png", dpi=130)
