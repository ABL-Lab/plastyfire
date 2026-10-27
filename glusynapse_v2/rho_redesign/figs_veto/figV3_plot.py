"""figV3_why_veto.png: Markram 10 Hz -10 vs +10 side by side (shared y axes), from figV3_traces.py + figV_outcomes.csv.

Rows: 1 spikes + spine Ca, first 3 pairings, one rho0 = 0 synapse (peaks annotated); 2 eCB pool / P1 at the first 5
arrivals vs the trigger k_E; 3 population control-lane release change d (mean +/- SEM over synapses), B4 vs veto off;
4 population rho by rho0 class, B4 vs Chindemi refit; 5 EPSP ratio: data, Chindemi refit, B4, B4 veto off.

  python glusynapse_v2/rho_redesign/figs_veto/figV3_plot.py
"""
import json, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

plt.style.use("/lustre09/project/6070394/dhuruva/plastyfire/onerule.mplstyle")
plt.rcParams["font.sans-serif"] = ["Helvetica Neue", "Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"]
plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
                     "legend.fontsize": 7.5, "axes.titlesize": 8.5})
BLUE, ORANGE, RED, GREY = "#68a8e0", "#f9a24b", "#ff4747", "0.55"
INK, POST = "0.2", "#c0392b"
WORK = "/scratch/dhuruva/figs_veto"
HERE = os.path.dirname(os.path.abspath(__file__))
KS = ("m10", "p10")
COLT = dict(m10="10 Hz trains, post 10 ms BEFORE pre (−10)", p10="10 Hz trains, pre 10 ms BEFORE post (+10)")
TITLES = ["1  The order is in the spine Ca",
          "2  eCB pool: above the trigger in both, cannot tell the order",
          "3  Release: −0.29 in both without the veto; the veto keeps +10 at 0",
          "4  ρ: nearly the same for −10 and +10",
          "5  Outcome (EPSP ratio, mean ± SEM)"]

z = np.load(os.path.join(WORK, "figV3.npz")); Z = {k: z[k] for k in z.files}
M = json.load(open(os.path.join(WORK, "figV3_meta.json")))
OC = pd.read_csv(os.path.join(WORK, "figV_outcomes.csv"))
EVF = M["ev_fields"]


def ev(k, f):
    return Z[f"{k}__EV"][:, EVF.index(f)]


def rel(k):
    return Z[f"{k}__t"] - Z[f"{k}__pre"][0]


def ticks(ax, x, y0, y1, col, lw=0.9):
    tr = ax.get_xaxis_transform()
    for v in x:
        ax.plot([v, v], [y0, y1], color=col, lw=lw, transform=tr, solid_capstyle="butt")


def mean_sem(A):
    A = np.asarray(A, float)
    return A.mean(0), A.std(0, ddof=1) / np.sqrt(len(A)) if len(A) > 1 else np.zeros(A.shape[1])


# ---------------------------------------------------------------- layout (inches, no layout engine)
W = 7.0
import textwrap
CAPTXT = ["1. The spine Ca knows the order: when the bAP follows glutamate (+10) the NMDA receptors are unblocked and "
          "the Ca peak is higher (row 1).",
          "2. The eCB pool (a 100 ms average of the bAP's VDCC Ca) and c* average this away: the pool is above the "
          "trigger in both trains (row 2), and ρ ends nearly the same (row 4).",
          "3. So the eCB step also fires in +10 trains and release drops to −0.29 (row 3, dashed); the veto is a patch "
          "that looks ahead for the bAP and cancels the step."]
CAPW = "\n".join(textwrap.fill(c, 118, subsequent_indent="    ") for c in CAPTXT)
COLH, GAP, BOT = 0.30, 0.50, 0.95
CAP = 0.14 + 0.142 * (CAPW.count("\n") + 1)
LEFT, RIGHT, CGAP = 0.78, 0.30, 0.40
HR = [1.15, 1.05, 0.95, 1.05, 0.95]
HEADS = [0.26, 0.44, 0.44, 0.44, 0.26]                   # heading band (+ legend line) above each row
cw = (W - LEFT - RIGHT - CGAP) / 2
H = CAP + COLH + sum(HR) + sum(HEADS) + (len(HR) - 1) * GAP + BOT
fig = plt.figure(figsize=(W, H), layout="none")
tops, y = [], H - CAP - COLH
for h, hd in zip(HR, HEADS):
    y -= hd; tops.append(y); y -= h + GAP
AX = [[fig.add_axes([(LEFT + j * (cw + CGAP)) / W, (tops[r] - HR[r]) / H, cw / W, HR[r] / H]) for j in range(2)]
      for r in range(len(HR))]
for j, k in enumerate(KS):
    fig.text((LEFT + j * (cw + CGAP) + cw / 2) / W, (H - CAP - 0.04) / H, COLT[k], ha="center", va="top",
             fontsize=8.5, fontweight="bold")
for r, t in enumerate(TITLES):
    fig.text(LEFT / W, (tops[r] + HEADS[r] - 0.19) / H, t, ha="left", va="bottom", fontsize=8.5, fontweight="bold",
             color="0.1")


def row_legend(r, handles, ncol):
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(LEFT / W, (tops[r] + 0.02) / H), ncol=ncol,
               frameon=False, handlelength=1.6, columnspacing=1.0, borderaxespad=0.0, fontsize=7.2)


# ---------------------------------------------------------------- row 1: spikes + spine Ca, first 3 pairings
peaks = {}
for j, k in enumerate(KS):
    ax = AX[0][j]; x = rel(k); pre = Z[f"{k}__pre"] - Z[f"{k}__pre"][0]; post = Z[f"{k}__post"] - Z[f"{k}__pre"][0]
    per = pre[1] - pre[0]
    lo, hi = -30.0, 2 * per + 80.0
    m = (x >= lo) & (x <= hi)
    cai = Z[f"{k}__cai"]
    ax.plot(x[m], cai[m], color=INK, lw=0.9)
    ticks(ax, pre[(pre >= lo) & (pre <= hi)], 0.90, 1.0, "k")
    ticks(ax, post[(post >= lo) & (post <= hi)], 0.77, 0.87, POST)
    pk = []
    for p in pre[:3]:
        mm = (x >= p - 15) & (x < p + 60)
        i = np.nanargmax(np.where(mm, cai, -np.inf)); pk.append(float(cai[i]))
        ax.annotate(f"{cai[i]:.2f} µM", (x[i], cai[i]), xytext=(4, 2), textcoords="offset points", fontsize=7.5,
                    color="0.1", ha="left", va="bottom")
    peaks[k] = pk
    ax.set_xlim(lo, hi)
    ax.set_xlabel("time from first pre spike (ms)")
    if j == 0:
        ax.set_ylabel("spine Ca (µM)")
        ax.text(0.01, 0.95, "pre", transform=ax.transAxes, fontsize=7, va="center", color="k")
        ax.text(0.01, 0.82, "post", transform=ax.transAxes, fontsize=7, va="center", color=POST)
top = max(max(v) for v in peaks.values()) * 1.6
for j in range(2):
    AX[0][j].set_ylim(0, top)

# ---------------------------------------------------------------- row 2: eCB pool at the first 5 arrivals
ymax = 0
for j, k in enumerate(KS):
    ax = AX[1][j]; x = rel(k); s = M["traces"][k]
    pre = Z[f"{k}__pre"] - Z[f"{k}__pre"][0]; post = Z[f"{k}__post"] - Z[f"{k}__pre"][0]
    per = pre[1] - pre[0]
    lo, hi = -30.0, 4 * per + 70.0
    m = (x >= lo) & (x <= hi)
    pool = Z[f"{k}__pool"] / s["uE"]
    ax.plot(x[m], pool[m], color=INK, lw=0.9)
    kE = s["thE"] / s["uE"]
    ax.axhline(kE, color="k", ls="--", lw=0.8)
    ta = ev(k, "t_a") + Z[f"{k}__t"][0] - Z[f"{k}__pre"][0]
    tv = ev(k, "trigval") / s["uE"]
    stp, vet = ev(k, "step") > 0.5, ev(k, "veto") > 0.5
    inz = (ta >= lo) & (ta <= hi)
    for a, v in zip(ta[inz], vet[inz]):
        if v:
            ax.axvspan(a + s["vt0"], a + s["vtv"], color=ORANGE, alpha=0.4, lw=0, zorder=0)
    ticks(ax, post[(post >= lo) & (post <= hi)], 0.88, 0.98, POST)
    ax.plot(ta[inz & stp], tv[inz & stp], "v", color=RED, ms=6, zorder=5)
    ax.plot(ta[inz & vet], tv[inz & vet], "o", mfc="white", mec="#1f4e79", mew=1.2, ms=6, zorder=5)
    ax.plot(ta[inz & ~stp & ~vet], tv[inz & ~stp & ~vet], ".", color=GREY, ms=6, zorder=5)
    ymax = max(ymax, np.nanmax(pool[m]), np.nanmax(tv[inz]))
    ax.set_xlim(lo, hi)
    ax.set_xlabel("time from first pre spike (ms)")
    if j == 0:
        ax.set_ylabel("eCB pool / P1")
    else:
        ax.text(1.0, kE, " trigger", transform=ax.get_yaxis_transform(), fontsize=7, va="center", ha="left")
for j in range(2):
    AX[1][j].set_ylim(0, ymax * 1.3)
row_legend(1, [Line2D([], [], ls="none", marker="v", color=RED, ms=6, label="release steps down"),
               Line2D([], [], ls="none", marker="o", mfc="white", mec="#1f4e79", ms=6, label="step vetoed"),
               Patch(color=ORANGE, alpha=0.4, lw=0, label="veto window")], 3)

# ---------------------------------------------------------------- row 3: population release change d
for j, k in enumerate(KS):
    ax = AX[2][j]; g = Z[f"{k}__grid"]
    ends = []
    for fk, col, ls, lab in (("B4", BLUE, "-", "B4 (with veto)"), ("KOV", RED, "--", "B4, veto off")):
        mu, se = mean_sem(Z[f"{k}__{fk}__D"])
        ax.fill_between(g, mu - se, mu + se, color=col, alpha=0.25, lw=0)
        ax.plot(g, mu, color=col, ls=ls, lw=1.3, label=lab)
        ends.append(f"{lab}: {mu[-1]:+.2f}")
    ax.axhline(0, color="0.85", lw=0.6, zorder=0)
    ax.axhline(-0.29, color=GREY, ls=":", lw=0.7)
    ax.text(0.97, 0.60 if k == "m10" else 0.24, "end of induction\n" + "\n".join(ends), transform=ax.transAxes, ha="right", va="center", fontsize=7.2, bbox=dict(fc="white", ec="0.8", lw=0.5, alpha=0.9, boxstyle="round,pad=0.3"))
    ax.set_ylim(-0.34, 0.06); ax.set_xlim(g[0], g[-1])
    ax.set_xlabel("induction time (s)")
    if j == 0:
        ax.set_ylabel("release change d\n(mean ± SEM, all synapses)")
row_legend(2, [Line2D([], [], color=BLUE, lw=1.3, label="B4 (with veto)"),
               Line2D([], [], color=RED, ls="--", lw=1.3, label="B4, veto off")], 2)

# ---------------------------------------------------------------- row 4: population rho by rho0 class
for j, k in enumerate(KS):
    ax = AX[3][j]; g = Z[f"{k}__grid"]; r0 = Z[f"{k}__rho0"]
    ends = []
    for fk, col, ls in (("B4", BLUE, "-"), ("CHR", ORANGE, ":")):
        Rm = Z[f"{k}__{fk}__R"]
        for c0 in (0.0, 1.0):
            sel = r0 == c0
            if not sel.any():
                continue
            mu, se = mean_sem(Rm[sel])
            ax.fill_between(g, mu - se, mu + se, color=col, alpha=0.18, lw=0)
            ax.plot(g, mu, color=col, ls=ls, lw=1.4)
            ends.append((c0, fk, mu[-1]))
    ax.text(0.02, 0.93, "synapses that start at ρ = 1", transform=ax.transAxes, fontsize=7, va="top", color="0.3")
    ax.text(0.02, 0.07, "synapses that start at ρ = 0", transform=ax.transAxes, fontsize=7, va="bottom", color="0.3")
    ED = {(c0, fk): v for c0, fk, v in ends}
    ax.text(0.97, 0.38, "\n".join(f"end ρ, start {c0:.0f}: B4 {ED[(c0, 'B4')]:.2f}, Chindemi {ED[(c0, 'CHR')]:.2f}" for c0 in (1.0, 0.0) if (c0, "B4") in ED), transform=ax.transAxes, ha="right", va="center", fontsize=7.0, bbox=dict(fc="white", ec="0.8", lw=0.5, alpha=0.9, boxstyle="round,pad=0.3"))
    ax.set_ylim(-0.05, 1.05); ax.set_xlim(g[0], g[-1])
    ax.set_xlabel("induction time (s)")
    if j == 0:
        ax.set_ylabel("ρ (mean ± SEM)")
row_legend(3, [Line2D([], [], color=BLUE, lw=1.4, label="B4"),
               Line2D([], [], color=ORANGE, ls=":", lw=1.4, label="Chindemi refit")], 2)

# ---------------------------------------------------------------- row 5: outcome
V = {}
for k in KS:
    b = OC[(OC.proto == k) & (OC.fit == "B4")].iloc[0]; c = OC[(OC.proto == k) & (OC.fit == "CHR")].iloc[0]
    v = OC[(OC.proto == k) & (OC.fit == "KOV")].iloc[0]
    V[k] = [("data", b.data_mean, b.data_sem, "k", True), ("Chindemi\nrefit", c.gpu_pred, c.gpu_sem, ORANGE, True),
            ("B4", b.gpu_pred, b.gpu_sem, BLUE, True), ("B4,\nveto off", v.cpu_mean, v.cpu_sem, RED, False)]
lo_ = min(mu - se for k in KS for _, mu, se, _, _ in V[k]); hi_ = max(mu + se for k in KS for _, mu, se, _, _ in V[k])
for j, k in enumerate(KS):
    ax = AX[4][j]
    for i, (lab, mu, se, col, filled) in enumerate(V[k]):
        ax.errorbar(i, mu, yerr=se, fmt="o", color=col, mfc=col if filled else "white", mec=col, ms=6, capsize=2,
                    lw=1.0, mew=1.2)
        ax.text(i + 0.14, mu, f"{mu:.2f}", fontsize=7.5, va="center", ha="left")
    ax.axhline(1.0, color="0.8", lw=0.6, zorder=0)
    ax.set_xticks(range(4)); ax.set_xticklabels([v[0] for v in V[k]])
    ax.set_xlim(-0.5, 3.6); ax.set_ylim(min(0.6, lo_ - 0.05), max(1.3, hi_ + 0.05))
    if j == 0:
        ax.set_ylabel("EPSP ratio")

# shared y axes: right column without tick labels
for r in range(len(HR)):
    AX[r][1].set_ylim(AX[r][0].get_ylim())
    AX[r][1].tick_params(labelleft=False)

cap = CAPTXT
fig.text(0.012, 1 - 0.07 / H, CAPW, ha="left", va="top", fontsize=7.6, linespacing=1.3)
ch = M["choice"]
fig.text(0.012, 0.05 / H, textwrap.fill(
    f"Rows 1–2: one L5→L5 synapse that starts at ρ = 0 (pair {ch['pair']}, synapse {ch['syn']}; of {ch['n_cand']} such "
    f"synapses, the one closest to the median end ρ at both protocols). Rows 3–4: all synapses of the simulated pairs. "
    f"B4 = best model (veto_peak_k 0.75, d_NO_max 0.14); veto off = same parameters, veto_peak_k 1e30; Chindemi refit = "
    f"ρ only. Data: mean ± SEM.", 150), ha="left", va="bottom", fontsize=6.3, color="0.3", linespacing=1.25)
out = os.path.join(HERE, "figV3_why_veto.png")
fig.savefig(out, dpi=300)
print(f"wrote {out} ({W} x {H:.2f} in)")
print("peaks (uM) first 3 pairings:", {k: [round(v, 3) for v in p] for k, p in peaks.items()})
for k in KS:
    print(k, "outcomes:", [(l.replace(chr(10), ' '), round(mu, 3), round(se, 3)) for l, mu, se, _, _ in V[k]])
print("caption:\n  " + "\n  ".join(cap))
