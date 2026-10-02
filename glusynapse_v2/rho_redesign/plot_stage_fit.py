"""Standard result figure for the staged core-target fits (CORE_TARGETS.md).

One panel per pathway: data (mean +/- SEM) against the model prediction for every target. Fitted (core) targets come
first, then a dashed line, then the validation targets. A last panel shows the per-target z^2, with fitted and
validation targets separated the same way.

Env:
  FIT     csv prefix with predictions for all targets, e.g. s1C_s_val (reads FIT.csv, FIT_l23.csv, FIT_l23l23.csv)
  FITLBL  legend label for FIT (default FIT)
  STAGE   1, 2 or 3: which core targets were fitted (7, 11 or 18; default 1)
  OUT     output path without extension (default glusynapse_v2/figs/stage_fits/stage<STAGE>/<FIT>)
MEASURED 22286465: 3 plots 0:15 elapsed, 241 MB MaxRSS; request 300M 0:15:00, 1 CPU.
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = "/project/rrg-emuller/dhuruva/plastyfire/glusynapse_v2/rho_redesign"
FIT = os.environ["FIT"]
FITLBL = os.environ.get("FITLBL", FIT)
STAGE = int(os.environ.get("STAGE", "1"))
OUT = os.environ.get("OUT", f"/project/rrg-emuller/dhuruva/plastyfire/glusynapse_v2/figs/stage_fits/stage{STAGE}/{FIT}")
PATHS = [("", "L5→L5", "#4c78a8"), ("_l23", "L2/3→L5", "#e45756"), ("_l23l23", "L2/3→L2/3", "#54a24b")]

# CORE_TARGETS.md, in order; stage 1 = first 7, stage 2 = first 11, stage 3 = all 18
CORE = ["10Hz_10ms|control", "10Hz_-10ms|control", "sjostrom_0.1hz_dt+10ms|control", "sjostrom_0.1hz_dt-10ms|control",
        "sjostrom_20hz_dt+10ms|control", "sjostrom_40hz_dt+10ms|control", "sjostrom_50hz_dt-10ms|control",
        "sjostrom_50hz_dt+10ms|nmdar_block", "sjostrom_0.1hz_dt-10ms|mglu_block",
        "sjostrom07_step200ms_pair|control", "sjostrom07_step200ms_pair|no_block",
        "letzkus_1ap_dt+10ms|control", "letzkus_3ap_200hz_dt+10ms@proximal|control",
        "letzkus_3ap_200hz_dt+10ms@distal|control",
        "zilberter_1ap_dt+10ms|control", "zilberter_1ap_dt-10ms|control",
        "zilberter_train10_50hz_dt+4ms_last|control", "zilberter_train10_50hz_dt-10ms_last|control"]
FITTED = set(CORE[:{1: 7, 2: 11, 3: 18}[STAGE]])


def load(name, suf):
    d = pd.read_csv(os.path.join(HERE, f"{name}{suf}.csv"))
    d["key"] = d.target + "|" + d.condition
    d["chi2"] = d.z ** 2
    d["fitted"] = d.key.isin(FITTED)
    # fitted first (in CORE order), then validation in file order
    d["order"] = [CORE.index(k) if f else 100 + i for i, (k, f) in enumerate(zip(d.key, d.fitted))]
    return d.sort_values("order").reset_index(drop=True)


def short(k):
    k = (k.replace("sjostrom07_step200ms_", "S07 ").replace("sjostrom_", "").replace("zilberter_", "Z ")
          .replace("letzkus_", "Lz ").replace("egger1999_", "Egger ").replace("_last", ""))
    return (k.replace("|control", "").replace("|mglu_block", " mGluR bl.").replace("|nmdar_block", " NMDAR bl.")
             .replace("|no_block", " NO bl.").replace("|", " "))


fit = {suf: load(FIT, suf) for suf, _, _ in PATHS}


def sums(d):
    return d.chi2[d.fitted].sum(), d.chi2[~d.fitted].sum(), int(d.fitted.sum()), int((~d.fitted).sum())


def panel(ax, suf, label, color):
    d = fit[suf]
    x = np.arange(len(d))
    cf, cv, nf, nv = sums(d)
    ax.errorbar(x, d.target_mean, yerr=d.target_sem, fmt="_", color="k", ms=10, capsize=2.5, lw=1.2, zorder=3,
                label="data (mean ± SEM)")
    ax.plot(x[d.fitted], d.pred[d.fitted], marker="o", color=color, ls="none", ms=6, zorder=4, label=None if nf == 0 else
            f"{FITLBL}, fitted (χ² {cf:.1f}, n {nf})")
    ax.plot(x[~d.fitted], d.pred[~d.fitted], marker="o", mfc="white", mec=color, mew=1.4, ls="none", ms=5.5,
            zorder=4, label=f"{FITLBL}, validation (χ² {cv:.1f}, n {nv})")
    if 0 < nf < len(d):
        ax.axvline(nf - 0.5, color="0.3", ls="--", lw=1.0, zorder=1)
        ax.text(nf - 0.6, 1.0, "fitted ", transform=ax.get_xaxis_transform(), ha="right", va="top", fontsize=7)
        ax.text(nf - 0.4, 1.0, " validation", transform=ax.get_xaxis_transform(), ha="left", va="top", fontsize=7)
    ax.axhline(1, color="0.75", lw=0.8, zorder=0)
    ax.set_xticks(x)
    ax.set_xticklabels([short(k) for k in d.key], rotation=60, ha="right", fontsize=6)
    for t, f in zip(ax.get_xticklabels(), d.fitted):
        if f:
            t.set_fontweight("bold")
    ax.set_xlim(-0.7, len(d) - 0.3)
    ax.set_ylabel("EPSP ratio (after/before)", fontsize=8)
    ax.set_title(f"{label}: {nf} fitted, {nv} validation", fontsize=9, loc="left")
    ax.legend(fontsize=6, loc="upper left", bbox_to_anchor=(0, 0.93), frameon=False)


fig = plt.figure(figsize=(14, 13))
gs = fig.add_gridspec(3, 25, height_ratios=[1.1, 1.0, 1.0], hspace=1.15, wspace=0.0)
axa = fig.add_subplot(gs[0, :])
axb = fig.add_subplot(gs[1, 0:8])
axc = fig.add_subplot(gs[1, 10:25])
axd = fig.add_subplot(gs[2, :])
for ax, (suf, label, color), tag in zip((axa, axb, axc), PATHS, "ABC"):
    panel(ax, suf, f"({tag}) {label}", color)

rows = []
for suf, label, color in PATHS:
    d = fit[suf].copy()
    d["path"], d["color"] = label, color
    rows.append(d)
allt = pd.concat(rows)
allt = pd.concat([allt[allt.fitted].sort_values("chi2", ascending=False),
                  allt[~allt.fitted].sort_values("chi2", ascending=False)]).reset_index(drop=True)
x = np.arange(len(allt))
nf = int(allt.fitted.sum())
axd.bar(x[allt.fitted], allt.chi2[allt.fitted], color=allt.color[allt.fitted], width=0.8)
axd.bar(x[~allt.fitted], allt.chi2[~allt.fitted], color="white", edgecolor=allt.color[~allt.fitted], lw=1.0,
        width=0.8)
axd.axvline(nf - 0.5, color="0.3", ls="--", lw=1.0)
axd.set_xticks(x)
axd.set_xticklabels([short(k) for k in allt.key], rotation=70, ha="right", fontsize=5)
axd.set_xlim(-0.7, len(allt) - 0.3)
axd.set_ylabel("χ² contribution (z²)", fontsize=8)
CF, CV = allt.chi2[allt.fitted].sum(), allt.chi2[~allt.fitted].sum()
axd.set_title(f"(D) per-target z² of {FITLBL}: fitted {CF:.1f} (n {nf}, filled) | validation {CV:.1f} "
              f"(n {len(allt) - nf}, open)", fontsize=9, loc="left")
handles = [plt.Rectangle((0, 0), 1, 1, color=c, label=l) for _, l, c in PATHS]
axd.legend(handles=handles, fontsize=7, loc="upper right", frameon=False)
for ax in (axa, axb, axc, axd):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=7)
fig.subplots_adjust(top=0.97, bottom=0.08, left=0.06, right=0.99)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT + ".pdf")
fig.savefig(OUT + ".png", dpi=150)
print("saved", OUT, "fitted", round(CF, 2), nf, "validation", round(CV, 2), len(allt) - nf)
