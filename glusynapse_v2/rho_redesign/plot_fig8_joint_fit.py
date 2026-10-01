"""Fig 8: joint L5 + L2/3->L5 fit (env FIT, default C1Ajd_s5), model vs data and per-target chi2."""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
FIT = os.environ.get("FIT", "C1Ajd_s5")
OUT = os.path.join(HERE, "figs", "fig8_joint_fit" + ("" if FIT == "C1Ajd_s5" else "_" + FIT))
C5, C23, CREF = "#4c78a8", "#e45756", "0.65"
# totals are recomputed from z below

d5 = pd.read_csv(os.path.join(RES, f"v4_{FIT}.csv"))
d23 = pd.read_csv(os.path.join(RES, f"v4_{FIT}_l23.csv"))
ref = pd.read_csv(os.path.join(RES, "v4_C1v2.csv"))
for d in (d5, d23, ref):
    d["key"] = d.target + "|" + d.condition
    d["chi2"] = d.z ** 2
refp = ref.set_index("key").pred
CH5, CH23, N5, N23 = d5.chi2.sum(), d23.chi2.sum(), len(d5), len(d23); TOT = CH5 + CH23

def is_val(k):
    k = k.lower()
    return "letzkus_3ap" in k and "-10" in k and "distal" in k
for d in (d5, d23):
    d["val"] = [is_val(k) for k in d.key]

def short(k):
    return (k.replace("sjostrom07_step200ms_", "S07 ").replace("sjostrom_", "")
             .replace("|control", "").replace("|", " "))

def panel(ax, d, color, title, overlay=False):
    x = np.arange(len(d))
    ax.errorbar(x, d.target_mean, yerr=d.target_sem, fmt="_", color="k", ms=12, capsize=3, lw=1.5, zorder=3, label="data (SEM)")
    if overlay:
        m = d.key.isin(refp.index)
        ax.plot(x[m.values] + 0.15, refp.loc[d.key[m]].values, marker="s", color=CREF, ls="none", ms=5,
                alpha=0.5, zorder=2, label="L5-only best fit (C1v2)")
    ok = ~d.val.values
    ax.plot(x[ok], d.pred[ok], marker="o", color=color, ls="none", ms=6, zorder=4, label=f"joint fit {FIT}")
    if (~ok).any():
        ax.plot(x[~ok], d.pred[~ok], marker="o", mfc="none", mec=color, mew=1.5, ls="none", ms=8, zorder=4,
                label="validation only (not fitted)")
    ax.axhline(1, color="0.7", lw=0.8, zorder=0)
    ax.set_xticks(x); ax.set_xticklabels([short(k) + (" [val]" if v else "") for k, v in zip(d.key, d.val)],
                                         rotation=60, ha="right", fontsize=6)
    ax.set_ylabel("weight ratio"); ax.set_title(title, fontsize=9)
    ax.legend(fontsize=6, loc="best")

fig = plt.figure(figsize=(14, 11))
gs = fig.add_gridspec(3, 6, height_ratios=[1.2, 1.0, 1.0], hspace=1.1, wspace=1.0)
axa = fig.add_subplot(gs[0, :]); axb = fig.add_subplot(gs[1, 0:3]); axc = fig.add_subplot(gs[2, :])
panel(axa, d5, C5, f"(A) L5->L5 targets: chi2 {CH5:.2f} over {N5} (joint total {TOT:.2f})", overlay=True)
panel(axb, d23, C23, f"(B) L2/3->L5 targets: chi2 {CH23:.2f} over {N23}")

allt = pd.concat([d5.assign(path="L5->L5"), d23.assign(path="L2/3->L5")]).sort_values("chi2", ascending=False).reset_index(drop=True)
x = np.arange(len(allt))
cols = [C5 if p == "L5->L5" else C23 for p in allt.path]
axc.bar(x, allt.chi2, color=cols, edgecolor=["k" if v else "none" for v in allt.val], hatch=None)
for i, v in enumerate(allt.val):
    if v: axc.bar(x[i], allt.chi2[i], color="none", edgecolor="k", hatch="//", lw=0.8)
axc.set_xticks(x); axc.set_xticklabels([short(k) + (" [val]" if v else "") for k, v in zip(allt.key, allt.val)],
                                       rotation=60, ha="right", fontsize=6)
axc.set_ylabel("chi2 contribution (z^2)")
axc.set_title(f"(C) per-target chi2, sorted: total {TOT:.2f} = {CH5:.2f} (L5, {N5}) + {CH23:.2f} (L2/3->L5, {N23})", fontsize=9)
axc.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=C5, label="L5->L5"),
                    plt.Rectangle((0, 0), 1, 1, color=C23, label="L2/3->L5"),
                    plt.Rectangle((0, 0), 1, 1, fc="none", ec="k", hatch="//", label="validation only")],
           fontsize=7, loc="upper right")
fig.subplots_adjust(top=0.95, bottom=0.08, left=0.06, right=0.98)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT + ".png", dpi=200); fig.savefig(OUT + ".pdf")
print("saved", OUT, "sum chi2 L5", d5.chi2.sum(), "L23", d23.chi2.sum(), "val rows", int(d5.val.sum() + d23.val.sum()))
