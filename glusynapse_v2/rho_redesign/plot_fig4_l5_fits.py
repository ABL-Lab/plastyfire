"""Fig 4: model vs data for the 29 L5 + Sjostrom 2007 targets, three fits."""
import re, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = "glusynapse_v2/"
FITS = [
    ("current rule", R + "results/reduced_gpu_subset_pl5sj07_td4_s2.csv", 39.45, 14, "#4c78a8", "o"),
    ("A0", R + "rho_redesign/results/v4_A0_s3.csv", 38.19, 14, "#f58518", "s"),
    ("B", R + "rho_redesign/results/v4_B_s2.csv", 35.23, 15, "#54a24b", "^"),
]
OUT = R + "rho_redesign/figs/fig4_l5_fits"
dfs = [pd.read_csv(f[1]) for f in FITS]
base = dfs[0].copy()
base["key"] = base.target + "|" + base.condition
for d in dfs:
    d["key"] = d.target + "|" + d.condition
    assert len(d) == 29, len(d)
keys = list(base.key)
P = np.array([d.set_index("key").loc[keys, "pred"].values for d in dfs])
Z = np.array([d.set_index("key").loc[keys, "z"].values for d in dfs])
kidx = {k: i for i, k in enumerate(keys)}
mean = base.target_mean.values; sem = base.target_sem.values
off = np.array([-0.12, 0, 0.12])

def draw(ax, ks, xs, labels=None, rot=0):
    xs = np.asarray(xs, float)
    ii = [kidx[k] for k in ks]
    ax.errorbar(xs, mean[ii], yerr=sem[ii], fmt="_", color="k", ms=14, capsize=3, lw=1.5, zorder=3, label="data (SEM)")
    for j, f in enumerate(FITS):
        ax.plot(xs + off[j] * 0.6, P[j, ii], marker=f[5], color=f[4], ls="none", ms=6, zorder=4)
    ax.axhline(1, color="0.7", lw=0.8, zorder=0)
    if labels is not None:
        ax.set_xticks(xs); ax.set_xticklabels(labels, rotation=rot, ha="right" if rot else "center", fontsize=7)
    ax.set_ylabel("weight ratio")

fig = plt.figure(figsize=(13, 10))
gs = fig.add_gridspec(3, 6, height_ratios=[1.2, 1.2, 0.9], hspace=0.75, wspace=0.9)
axa = fig.add_subplot(gs[0, 0:3]); axb = fig.add_subplot(gs[0, 3:6])
axc = fig.add_subplot(gs[1, 0:4]); axd = fig.add_subplot(gs[1, 4:6])
axz = fig.add_subplot(gs[2, :])

# (a) frequency
freqs = ["0.1", "10", "20", "40", "50"]
xpos = np.arange(len(freqs))
for sgn, ls, name in (("+", "-", "dt +10 ms"), ("-", "--", "dt -10 ms")):
    ks = [f"sjostrom_{f}hz_dt{sgn}10ms|control" for f in freqs]
    ii = [kidx[k] for k in ks]
    ax = axa
    ax.errorbar(xpos + (0.05 if sgn == "+" else -0.05), mean[ii], yerr=sem[ii], fmt="o" if sgn == "+" else "D",
                mfc="k" if sgn == "+" else "w", color="k", capsize=3, zorder=3, label=f"data {name}")
    for j, f in enumerate(FITS):
        ax.plot(xpos, P[j, ii], color=f[4], ls=ls, marker=f[5], ms=5, lw=1.2, alpha=0.9)
axa.set_xticks(xpos); axa.set_xticklabels(freqs); axa.set_xlabel("pairing frequency (Hz)")
axa.axhline(1, color="0.7", lw=0.8); axa.set_ylabel("weight ratio")
axa.set_title("(a) Sjostrom 2001 frequency (solid +10 ms, dashed -10 ms)", fontsize=9)
axa.legend(fontsize=6, loc="upper left")

# (b) timing + burst
tk = [f"sjostrom_0.1hz_dt{d}ms|control" for d in ("-200", "-120", "-25", "-10", "+10")]
bk = [f"sjostrom_burst5x20hz_r50_dt{d}ms|control" for d in ("-200", "-120")]
bl = ["0.1Hz -200", "0.1Hz -120", "0.1Hz -25", "0.1Hz -10", "0.1Hz +10", "burst r50 -200", "burst r50 -120"]
draw(axb, tk + bk, range(7), bl, rot=35)
axb.set_title("(b) 0.1 Hz timing and burst r50", fontsize=9)

# (c) pharmacology + 2007
ck = ["sjostrom_0.1hz_dt-10ms|mglu_block", "sjostrom_20hz_dt-10ms|mglu_block", "sjostrom_50hz_dt+10ms|mglu_block",
      "sjostrom_0.1hz_dt-10ms|nmdar_block", "sjostrom_20hz_dt-10ms|nmdar_block", "sjostrom_50hz_dt+10ms|nmdar_block",
      "sjostrom07_step200ms_pair|control", "sjostrom07_step200ms_pair|mglu_block", "sjostrom07_step200ms_pair|no_block",
      "sjostrom07_step200ms_pre_only|control", "sjostrom07_step200ms_post_only|control"]
ck = [k for k in ck if k in kidx]
cl = [k.replace("sjostrom07_step200ms_", "S07 ").replace("sjostrom_", "").replace("|", " ") for k in ck]
xc = np.arange(len(ck)); xc = xc + (xc >= 6) * 0.6
draw(axc, ck, xc, cl, rot=35)
axc.set_title("(c) pharmacology and Sjostrom 2007 step200ms", fontsize=9)

# (d) Markram 10 Hz
dk = ["10Hz_-10ms|control", "10Hz_5ms|control", "10Hz_10ms|control"]
draw(axd, dk, range(3), ["-10 ms", "+5 ms", "+10 ms"])
axd.set_title("(d) Markram / paired L5, 10 Hz", fontsize=9)

# z panel
x = np.arange(len(keys))
for j, f in enumerate(FITS):
    axz.plot(x + off[j], Z[j], marker=f[5], color=f[4], ls="none", ms=5,
             label=f"{f[0]}: chi2 {f[2]:.2f}, {f[3]} free")
axz.axhline(0, color="k", lw=0.6)
for v in (-2, 2): axz.axhline(v, color="r", lw=0.6, ls=":")
short = [k.replace("sjostrom07_step200ms_", "S07 ").replace("sjostrom_", "").replace("|control", "").replace("|", " ") for k in keys]
axz.set_xticks(x); axz.set_xticklabels(short, rotation=60, ha="right", fontsize=6)
axz.set_ylabel("z = (model - data)/SEM"); axz.set_title("per-target z-scores (all 29)", fontsize=9)
axz.legend(fontsize=8, ncol=3, loc="upper left")

h = [plt.Line2D([], [], color=f[4], marker=f[5], ls="none", label=f"{f[0]}: chi2 {f[2]:.2f}, {f[3]} free") for f in FITS]
h.append(plt.Line2D([], [], color="k", marker="_", ls="none", ms=12, label="data (SEM)"))
fig.legend(handles=h, loc="upper center", ncol=4, fontsize=9, frameon=False)
fig.subplots_adjust(top=0.92, bottom=0.12, left=0.06, right=0.98)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT + ".png", dpi=200); fig.savefig(OUT + ".pdf")
print("saved", OUT)
