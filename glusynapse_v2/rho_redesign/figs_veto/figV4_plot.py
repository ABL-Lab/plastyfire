"""figV4_why_ecb.png: why the eCB step exists (0.1 Hz) and why it then needs the veto (10 Hz). From figV4_traces.py.

Columns: Sjostrom 0.1 Hz -10 | +10 | Markram 10 Hz -10 | +10 (L5 -> L5). Rows: 1 spikes + spine Ca of one rho0 = 0
synapse (B4); 2 eCB pool / P1 at the arrivals vs the trigger; 3 population control-lane d (B4, eCB off, veto off);
4 population rho by rho0 class (B4, Chindemi refit); 5 EPSP ratio (data, Chindemi refit, B4, eCB off, veto off).

  python glusynapse_v2/rho_redesign/figs_veto/figV4_plot.py
"""
import json, os, textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.style.use("/lustre09/project/6070394/dhuruva/plastyfire/onerule.mplstyle")
plt.rcParams["font.sans-serif"] = ["Helvetica Neue", "Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"]
plt.rcParams.update({"font.size": 7.5, "axes.labelsize": 7.5, "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "legend.fontsize": 7.2, "axes.titlesize": 8})
BLUE, ORANGE, RED, GREY, GREEN = "#68a8e0", "#f9a24b", "#ff4747", "0.55", "#3c9d5d"
INK, POST = "0.2", "#c0392b"
WORK = "/scratch/dhuruva/figs_veto"
HERE = os.path.dirname(os.path.abspath(__file__))
KS = ("s01m", "s01p", "m10", "p10")
COLT = dict(s01m="−10: post first", s01p="+10: pre first", m10="−10: post first", p10="+10: pre first")
GROUPS = [("0.1 Hz: one pairing every 10 s  →  why the eCB step exists", 0, 1),
          ("10 Hz trains  →  why the eCB step then needs the veto", 2, 3)]
TITLES = ["1  Spike timing and spine Ca (one synapse that starts at ρ = 0)",
          "2  eCB pool (bAP VDCC Ca, 100 ms memory) at each presynaptic arrival",
          "3  Release change d (all synapses, mean ± SEM)",
          "4  ρ (Chindemi's bistable efficacy), by start state",
          "5  Outcome: EPSP ratio (mean ± SEM over pairs)"]

z = np.load(os.path.join(WORK, "figV4.npz")); Z = {k: z[k] for k in z.files}
M = json.load(open(os.path.join(WORK, "figV4_meta.json")))
OC = pd.read_csv(os.path.join(WORK, "figV4_outcomes.csv"))
EVF = M["ev_fields"]


def oc(k, fk):
    return OC[(OC.proto == k) & (OC.fit == fk)].iloc[0]


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
    return A.mean(0), (A.std(0, ddof=1) / np.sqrt(len(A)) if len(A) > 1 else np.zeros(A.shape[1]))


def win(k, row):
    if k.startswith("s01"):
        return -40.0, 160.0
    return (-30.0, 280.0) if row == 0 else (-30.0, 480.0)


# ---------------------------------------------------------------- caption (numbers from the outcomes table)
c0, b0, e0 = oc("s01m", "CHR"), oc("s01m", "B4"), oc("s01m", "KOE")
cp, bp, vp = oc("p10", "CHR"), oc("p10", "B4"), oc("p10", "KOV")
CAPTXT = [
    f"1. Why eCB exists (0.1 Hz −10): one bAP 10 ms before glutamate is too little Ca for ρ to depress, so Chindemi's "
    f"rule gives no change ({c0['mean']:.2f}; data {c0.data_mean:.2f}). The data's LTD is presynaptic and CB1-dependent "
    f"(CB1 block abolishes it: Sjöström 2003 Fig 8B). The bAP's Ca is still in the pool when glutamate arrives, the eCB step fires, release drops and "
    f"B4 gives {b0['mean']:.2f} (eCB off: {e0['mean']:.2f}).",
    "2. At 0.1 Hz the pool reads the order correctly: at +10 the bAP comes after the arrival, the pool is empty when "
    "the arrival checks it, and no step fires. No veto is needed here.",
    f"3. In 10 Hz trains every arrival also sees the bAP of the previous pairing (90 ms earlier at +10), so the pool "
    f"cannot tell −10 from +10. Without the veto the step also fires at +10 ({vp['mean']:.2f}; data "
    f"{vp.data_mean:.2f}). The veto cancels a step when its own bAP follows within 17 ms (B4 {bp['mean']:.2f}). "
    f"Chindemi's ρ alone gets +10 ({cp['mean']:.2f}) but not 0.1 Hz −10."]
CAPW = "\n".join(textwrap.fill(c, 128, subsequent_indent="    ") for c in CAPTXT)

# ---------------------------------------------------------------- layout (inches, no layout engine)
W = 7.5
NC = len(KS)
GRPH, COLH, GAP, BOT = 0.22, 0.20, 0.42, 0.62
CAP = 0.14 + 0.135 * (CAPW.count("\n") + 1)
LEFT, RIGHT, CGAP, MGAP = 0.62, 0.42, 0.14, 0.34      # MGAP = extra gap between the 0.1 Hz and 10 Hz groups
HR = [0.95, 0.95, 0.85, 0.85, 0.95]
HEADS = [0.24, 0.40, 0.40, 0.40, 0.24]
cw = (W - LEFT - RIGHT - (NC - 1) * CGAP - MGAP) / NC
XL = [LEFT + j * (cw + CGAP) + (MGAP if j >= 2 else 0) for j in range(NC)]
H = CAP + GRPH + COLH + sum(HR) + sum(HEADS) + (len(HR) - 1) * GAP + BOT
fig = plt.figure(figsize=(W, H), layout="none")
tops, y = [], H - CAP - GRPH - COLH
for h, hd in zip(HR, HEADS):
    y -= hd; tops.append(y); y -= h + GAP
AX = [[fig.add_axes([XL[j] / W, (tops[r] - HR[r]) / H, cw / W, HR[r] / H]) for j in range(NC)]
      for r in range(len(HR))]
for txt, a, b in GROUPS:
    x0, x1 = XL[a], XL[b] + cw
    yg = H - CAP - 0.03
    fig.text((x0 + x1) / 2 / W, yg / H, txt, ha="center", va="top", fontsize=8.3, fontweight="bold")
    fig.add_artist(Line2D([x0 / W, x1 / W], [(yg - 0.19) / H] * 2, color="0.5", lw=0.6))
for j, k in enumerate(KS):
    fig.text((XL[j] + cw / 2) / W, (H - CAP - GRPH - 0.02) / H, COLT[k], ha="center", va="top", fontsize=7.8)
for r, t in enumerate(TITLES):
    fig.text(XL[0] / W, (tops[r] + HEADS[r] - 0.18) / H, t, ha="left", va="bottom", fontsize=8, fontweight="bold",
             color="0.1")


def row_legend(r, handles, ncol):
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(XL[0] / W, (tops[r] + 0.02) / H), ncol=ncol,
               frameon=False, handlelength=1.6, columnspacing=1.0, borderaxespad=0.0, fontsize=7.0)


# ---------------------------------------------------------------- row 1: spikes + spine Ca
top = 0
for j, k in enumerate(KS):
    ax = AX[0][j]; x = rel(k); p0 = Z[f"{k}__pre"][0]
    pre, post = Z[f"{k}__pre"] - p0, Z[f"{k}__post"] - p0
    lo, hi = win(k, 0)
    m = (x >= lo) & (x <= hi)
    cai = Z[f"{k}__cai"]
    ax.plot(x[m], cai[m], color=INK, lw=0.9)
    ticks(ax, pre[(pre >= lo) & (pre <= hi)], 0.90, 1.0, "k")
    ticks(ax, post[(post >= lo) & (post <= hi)], 0.77, 0.87, POST)
    for p in pre[(pre >= lo) & (pre <= hi)]:
        mm = (x >= p - 15) & (x < p + 60)
        i = np.nanargmax(np.where(mm, cai, -np.inf))
        ax.annotate(f"{cai[i]:.1f}", (x[i], cai[i]), xytext=(3, 1), textcoords="offset points", fontsize=6.8,
                    color="0.1", ha="left", va="bottom")
    top = max(top, np.nanmax(cai[m]))
    ax.set_xlim(lo, hi)
    ax.set_xlabel("ms from first pre spike")
    if j == 0:
        ax.set_ylabel("spine Ca (µM)")
        ax.text(0.02, 0.95, "pre", transform=ax.transAxes, fontsize=6.5, va="center", color="k")
        ax.text(0.02, 0.82, "post", transform=ax.transAxes, fontsize=6.5, va="center", color=POST)
for j in range(NC):
    AX[0][j].set_ylim(0, top * 1.6)

# ---------------------------------------------------------------- row 2: eCB pool at the arrivals
ymax = 0
for j, k in enumerate(KS):
    ax = AX[1][j]; x = rel(k); s = M["traces"][k]; p0 = Z[f"{k}__pre"][0]
    post = Z[f"{k}__post"] - p0
    lo, hi = win(k, 1)
    m = (x >= lo) & (x <= hi)
    pool = Z[f"{k}__pool"] / s["uE"]
    ax.plot(x[m], pool[m], color=INK, lw=0.9)
    kE = s["thE"] / s["uE"]
    ax.axhline(kE, color="k", ls="--", lw=0.7)
    ta = ev(k, "t_a") + Z[f"{k}__t"][0] - p0
    tv = ev(k, "trigval") / s["uE"]
    trig, stp, vet = ev(k, "trig") > 0.5, ev(k, "step") > 0.5, ev(k, "veto") > 0.5
    inz = (ta >= lo) & (ta <= hi)
    canc = trig & ~stp                                   # trigger crossed, step cancelled (veto / C-annulment)
    below = ~trig
    ticks(ax, post[(post >= lo) & (post <= hi)], 0.88, 0.98, POST)
    ax.plot(ta[inz & stp], tv[inz & stp], "v", color=RED, ms=5.5, zorder=5)
    ax.plot(ta[inz & canc], tv[inz & canc], "o", mfc="white", mec="#1f4e79", mew=1.1, ms=5.5, zorder=5)
    ax.plot(ta[inz & below], tv[inz & below], "o", color=GREY, ms=3.5, zorder=5)
    ymax = max(ymax, np.nanmax(pool[m]))
    ax.set_xlim(lo, hi)
    ax.set_xlabel("ms from first pre spike")
    if j == 0:
        ax.set_ylabel("eCB pool / P1")
    if j == NC - 1:
        ax.text(1.02, kE, "trigger", transform=ax.get_yaxis_transform(), fontsize=6.5, va="center", ha="left")
for j in range(NC):
    AX[1][j].set_ylim(0, max(ymax, 0.6) * 1.35)
row_legend(1, [Line2D([], [], ls="none", marker="v", color=RED, ms=5.5, label="pool above trigger → release steps down"),
               Line2D([], [], ls="none", marker="o", mfc="white", mec="#1f4e79", ms=5.5,
                      label="above trigger, step cancelled by the veto"),
               Line2D([], [], ls="none", marker="o", color=GREY, ms=3.5, label="below trigger, no step")], 3)

# ---------------------------------------------------------------- row 3: population release change d
DSET = (("B4", BLUE, "-", "B4"), ("KOE", GREEN, "-.", "B4, eCB off"), ("KOV", RED, "--", "B4, veto off"))
for j, k in enumerate(KS):
    ax = AX[2][j]; g = Z[f"{k}__grid"]
    for fk, col, ls, lab in DSET:
        mu, se = mean_sem(Z[f"{k}__{fk}__D"])
        ax.fill_between(g, mu - se, mu + se, color=col, alpha=0.22, lw=0)
        ax.plot(g, mu, color=col, ls=ls, lw=1.2)
    ax.axhline(0, color="0.85", lw=0.6, zorder=0)
    ax.axhline(-0.29, color=GREY, ls=":", lw=0.7)
    ax.set_ylim(-0.34, 0.18); ax.set_xlim(g[0], g[-1])
    ax.set_xlabel("induction time (s)")
    if j == 0:
        ax.set_ylabel("release change d")
        ax.text(0.98, -0.29, "eCB floor", transform=ax.get_yaxis_transform(), fontsize=6.3, va="bottom", ha="right",
                color="0.4")
row_legend(2, [Line2D([], [], color=c, ls=ls, lw=1.2, label=lab) for _, c, ls, lab in DSET], 3)

# ---------------------------------------------------------------- row 4: population rho by rho0 class
for j, k in enumerate(KS):
    ax = AX[3][j]; g = Z[f"{k}__grid"]; r0 = Z[f"{k}__rho0"]
    for fk, col, ls in (("B4", BLUE, "-"), ("CHR", ORANGE, ":")):
        Rm = Z[f"{k}__{fk}__R"]
        for cc in (0.0, 1.0):
            sel = r0 == cc
            if sel.any():
                mu, se = mean_sem(Rm[sel])
                ax.fill_between(g, mu - se, mu + se, color=col, alpha=0.18, lw=0)
                ax.plot(g, mu, color=col, ls=ls, lw=1.3)
    if j == 0:
        ax.text(0.03, 0.93, "start ρ = 1", transform=ax.transAxes, fontsize=6.5, va="top", color="0.3")
        ax.text(0.03, 0.07, "start ρ = 0", transform=ax.transAxes, fontsize=6.5, va="bottom", color="0.3")
        ax.set_ylabel("ρ")
    ax.set_ylim(-0.05, 1.05); ax.set_xlim(g[0], g[-1])
    ax.set_xlabel("induction time (s)")
row_legend(3, [Line2D([], [], color=BLUE, lw=1.3, label="B4"),
               Line2D([], [], color=ORANGE, ls=":", lw=1.3, label="Chindemi refit (ρ only, no eCB)")], 2)

# ---------------------------------------------------------------- row 5: outcome
OSET = (("data", "k", True), ("Chindemi", ORANGE, True), ("B4", BLUE, True), ("eCB\noff", GREEN, False),
        ("veto\noff", RED, False))
lo_, hi_ = 9, 0
for j, k in enumerate(KS):
    ax = AX[4][j]
    vals = [(oc(k, "B4").data_mean, oc(k, "B4").data_sem)] + [(oc(k, fk)["mean"], oc(k, fk)["sem"])
                                                              for fk in ("CHR", "B4", "KOE", "KOV")]
    for i, ((lab, col, filled), (mu, se)) in enumerate(zip(OSET, vals)):
        ax.errorbar(i, mu, yerr=se, fmt="o", color=col, mfc=col if filled else "white", mec=col, ms=4.5, capsize=1.5,
                    lw=0.9, mew=1.0)
        ax.text(i, mu + se + 0.02, f"{mu:.2f}", fontsize=6.0, va="bottom", ha="center")
        lo_, hi_ = min(lo_, mu - se), max(hi_, mu + se)
    ax.axhline(1.0, color="0.8", lw=0.6, zorder=0)
    ax.set_xticks(range(len(OSET))); ax.set_xticklabels([o[0] for o in OSET], fontsize=6.3, rotation=0)
    ax.set_xlim(-0.6, len(OSET) - 0.4)
    if j == 0:
        ax.set_ylabel("EPSP ratio")
for j in range(NC):
    AX[4][j].set_ylim(min(0.55, lo_ - 0.05), max(1.3, hi_ + 0.1))

for r in range(len(HR)):
    for j in range(1, NC):
        if r in (2, 3, 4) or (r in (0, 1)):
            AX[r][j].set_ylim(AX[r][0].get_ylim())
            AX[r][j].tick_params(labelleft=False)

fig.text(0.012, 1 - 0.07 / H, CAPW, ha="left", va="top", fontsize=7.2, linespacing=1.3)
ch = M["choice"]
fig.text(0.012, 0.04 / H, textwrap.fill(
    f"L5→L5. Rows 1–2: one synapse that starts at ρ = 0 (pair {ch['pair']}, synapse {ch['syn']}), under B4; 0.1 Hz "
    f"columns show the first pairing, 10 Hz columns the first pairings of the train. Rows 3–5: all synapses / pairs. "
    f"B4 = best model (k_V 0.75, d_NO_max 0.14); eCB off = B4 with A_eCB 0; veto off = B4 with veto_peak_k 1e30 (no "
    f"refit); Chindemi refit = ρ only. Data (mean ± SEM): 0.1 Hz −10 Sjöström 2001 Fig 7B; 0.1 Hz +10 Sjöström 2001 Fig 1D; 10 Hz ±10 Markram 1997.", 165),
    ha="left", va="bottom", fontsize=6.0, color="0.3", linespacing=1.25)
out = os.path.join(HERE, "figV4_why_ecb.png")
fig.savefig(out, dpi=300)
print(f"wrote {out} ({W} x {H:.2f} in)")
print("caption:\n  " + "\n  ".join(CAPTXT))
