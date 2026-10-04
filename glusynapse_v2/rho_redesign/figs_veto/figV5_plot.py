"""figV5_chindemi_amount.png: Chindemi's rho reads the amount of Ca (frequency), not the order. From figV5_traces.py.

Row 1: c* of one synapse, -10 (red) vs +10 (blue), with that synapse's theta_d / theta_p; 0.1 | 10 | 20 Hz.
Row 2 left: time per synapse in the depression zone and above theta_p over the induction (mean +/- SEM), -10 vs +10.
Row 2 right: EPSP ratio vs frequency: data (filled) and Chindemi refit (open), -10 red, +10 blue.
Data sources are read from ebner_targets.csv (source_paper, figure_or_table), not typed by hand.

  python glusynapse_v2/rho_redesign/figs_veto/figV5_plot.py
"""
import json, os, textwrap
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
RED, BLUE = "#e0473c", "#2f6fb5"
WORK = "/scratch/dhuruva/figs_veto"
HERE = os.path.dirname(os.path.abspath(__file__))
FR = ("0.1", "10", "20")
z = np.load(os.path.join(WORK, "figV5.npz")); Z = {k: z[k] for k in z.files}
M = json.load(open(os.path.join(WORK, "figV5_meta.json")))
ZN = pd.read_csv(os.path.join(WORK, "figV5_zones.csv"))
OC = pd.read_csv(os.path.join(WORK, "figV5_outcomes.csv")).set_index("proto")


def ticks(ax, x, y0, y1, col):
    tr = ax.get_xaxis_transform()
    for v in x:
        ax.plot([v, v], [y0, y1], color=col, lw=0.9, transform=tr, solid_capstyle="butt")


o = {k: OC.loc[k] for k in OC.index}
CAPTXT = [
    "Chindemi's ρ has one input, c* (the low-passed spine Ca). Between θd and θp ρ goes down, above θp ρ goes up. "
    "More pairings per second → more Ca. And pre-first (+10) always gives MORE Ca than post-first (−10): the NMDA "
    "receptors are already bound when the bAP arrives (row 2 left: blue ≥ red at every frequency).",
    f"So the rule can never depress −10 without depressing +10 as much or more. 0.1 Hz: little Ca in either order → "
    f"no change (−10: {o['0.1-']['mean']:.2f}, data {o['0.1-'].data_mean:.2f}). 10 Hz: both orders in the depression "
    f"zone → +10 also depressed ({o['10+']['mean']:.2f}, data {o['10+'].data_mean:.2f}). 20 Hz: both orders reach θp → "
    f"−10 potentiates ({o['20-']['mean']:.2f}, data {o['20-'].data_mean:.2f}).",
    "Bottom right: the data split by ORDER (−10 is LTD at every frequency, +10 LTP grows with frequency). The refit's "
    "two lines move together with FREQUENCY."]
CAPW = "\n".join(textwrap.fill(c, 125, subsequent_indent="   ") for c in CAPTXT)

W = 7.2
CAP = 0.16 + 0.14 * (CAPW.count("\n") + 1)
H = CAP + 0.85 + 1.55 + 0.75 + 1.75 + 0.75
fig = plt.figure(figsize=(W, H), layout="none")
L, R, G = 0.62, 0.18, 0.32
cw = (W - L - R - 2 * G) / 3
y1 = H - CAP - 0.85 - 1.55
A1 = [fig.add_axes([(L + j * (cw + G)) / W, y1 / H, cw / W, 1.55 / H]) for j in range(3)]
y2 = y1 - 0.75 - 1.75
w2 = (W - L - R - 0.75) / 2
A2 = fig.add_axes([L / W, y2 / H, w2 / W, 1.75 / H])
A3 = fig.add_axes([(L + w2 + 0.75) / W, y2 / H, w2 / W, 1.75 / H])
fig.text(0.012, 1 - 0.07 / H, CAPW, ha="left", va="top", fontsize=7.5, linespacing=1.3)
fig.text(L / W, (y1 + 1.55 + 0.60) / H, "1  c* of one synapse (Chindemi refit), first pairings",
         fontsize=8.5, fontweight="bold")
fig.text(L / W, (y2 + 1.75 + 0.30) / H, "2  Whole induction, all synapses", fontsize=8.5, fontweight="bold")

# ---- row 1: c* traces with theta_d / theta_p
ymax = 0
for j, f in enumerate(FR):
    ax = A1[j]
    lo, hi = (-40.0, 260.0) if f == "0.1" else (-30.0, 4 * 1000.0 / float(f) + 60.0)
    for s, col in (("-", RED), ("+", BLUE)):
        k = f + s
        if f"{k}__t" not in Z:
            continue
        t, c = Z[f"{k}__t"], Z[f"{k}__c"]
        m = (t >= lo) & (t <= hi)
        ax.plot(t[m], c[m], color=col, lw=1.0)
        ymax = max(ymax, float(np.nanmax(c[m])))
        post = Z[f"{k}__post"]; pre = Z[f"{k}__pre"]
        ticks(ax, post[(post >= lo) & (post <= hi)], 0.90 if s == "-" else 0.80, 0.97 if s == "-" else 0.87, col)
    ticks(ax, Z[f"{f}-__pre"][(Z[f'{f}-__pre'] >= lo) & (Z[f'{f}-__pre'] <= hi)], 0.70, 0.77, "k")
    s0 = M["traces"][f + "-"]
    ax.axhline(s0["td"], color="0.3", ls="--", lw=0.8); ax.axhline(s0["tp"], color="0.3", ls="-.", lw=0.8)
    ax.set_xlim(lo, hi); ax.set_title(f"{f} Hz", pad=3)
    ax.set_xlabel("ms from first pre spike")
    if j == 0:
        ax.set_ylabel("c* (spine Ca, low-passed)")
for j, f in enumerate(FR):
    ax = A1[j]; s0 = M["traces"][f + "-"]
    top = max(ymax, s0["tp"]) * 1.45
    ax.set_ylim(0, top)
    ax.axhspan(s0["td"], s0["tp"], color="#fbe1de", lw=0, zorder=0)
    ax.axhspan(s0["tp"], top, color="#dff0e4", lw=0, zorder=0)
    if j == 2:
        ax.text(1.02, s0["td"], "θd", transform=ax.get_yaxis_transform(), fontsize=7.5, va="center")
        ax.text(1.02, s0["tp"], "θp", transform=ax.get_yaxis_transform(), fontsize=7.5, va="center")
    if j == 0:
        ax.text(0.98, 0.70, "pre", transform=ax.transAxes, fontsize=6.5, ha="right", va="center")
fig.legend(handles=[Line2D([], [], color=RED, lw=1, label="−10 (post first)"),
                    Line2D([], [], color=BLUE, lw=1, label="+10 (pre first)"),
                    Patch(color="#fbe1de", label="depression zone"),
                    Patch(color="#dff0e4", label="potentiation zone")],
           loc="lower left", bbox_to_anchor=(L / W, (y1 + 1.55 + 0.22) / H), ncol=4, frameon=False, fontsize=7.2)

# ---- row 2 left: time in zones
x = np.arange(3); bw = 0.36
for s, col, dx in (("-", RED, -bw / 2), ("+", BLUE, bw / 2)):
    dep, pot, dse, pse = [], [], [], []
    for f in FR:
        q = ZN[ZN.proto == f + s]
        dep.append(q.t_dep.mean()); pot.append(q.t_pot.mean())
        dse.append(q.t_dep.std(ddof=1) / np.sqrt(len(q))); pse.append(q.t_pot.std(ddof=1) / np.sqrt(len(q)))
    A2.bar(x + dx, dep, bw, color=col, alpha=0.45, yerr=dse, ecolor="0.3", capsize=1.5, lw=0)
    A2.bar(x + dx, pot, bw, bottom=dep, color=col, alpha=1.0, yerr=pse, ecolor="0.3", capsize=1.5, lw=0,
           hatch="////", edgecolor="white")
A2.set_xticks(x); A2.set_xticklabels([f"{f} Hz" for f in FR])
A2.set_ylabel("time per synapse (s)")
A2.set_title("time c* spends in each zone", fontsize=8)
A2.legend(handles=[Patch(color="0.5", alpha=0.45, label="depression zone"),
                   Patch(facecolor="0.3", hatch="////", edgecolor="white", label="potentiation zone"),
                   Patch(color=RED, label="−10"), Patch(color=BLUE, label="+10")],
          loc="upper left", frameon=False, fontsize=7, ncol=2)

# ---- row 2 right: outcome vs frequency
xf = np.arange(3)
for s, col, lab in (("-", RED, "−10"), ("+", BLUE, "+10")):
    d = [o[f + s].data_mean for f in FR]; de = [o[f + s].data_sem for f in FR]
    m = [o[f + s]["mean"] for f in FR]; me = [o[f + s]["sem"] for f in FR]
    A3.errorbar(xf - 0.06, d, yerr=de, color=col, marker="o", ms=5, lw=1.4, capsize=2, label=f"data {lab}")
    A3.errorbar(xf + 0.06, m, yerr=me, color=col, marker="o", mfc="white", ms=5, lw=1.2, ls="--", capsize=2,
                label=f"Chindemi refit {lab}")
A3.axhline(1, color="0.8", lw=0.6, zorder=0)
A3.set_xticks(xf); A3.set_xticklabels([f"{f} Hz" for f in FR]); A3.set_xlim(-0.4, 2.4)
A3.set_ylabel("EPSP ratio"); A3.set_title("outcome: data vs Chindemi refit", fontsize=8)
A3.legend(loc="upper left", frameon=False, fontsize=6.8, ncol=2)
lo_ = min(min(o[k].data_mean - o[k].data_sem, o[k]["mean"] - o[k]["sem"]) for k in o)
hi_ = max(max(o[k].data_mean + o[k].data_sem, o[k]["mean"] + o[k]["sem"]) for k in o)
A3.set_ylim(lo_ - 0.05, hi_ + 0.25)

papers = sorted(set(s.split(", Fig")[0] for s in OC.source))           # source strings from ebner_targets.csv
figs = {s: ", ".join("Fig" + OC.loc[f + s, "source"].split(", Fig", 1)[1] for f in FR) for s in ("-", "+")}
foot = (f"L5→L5, Chindemi refit (CHR_s5). Row 1: pair {M['traces']['0.1-']['pair']}, synapse "
        f"{M['traces']['0.1-']['syn']} (start ρ = {M['traces']['0.1-']['rho0']:.0f}); θd/θp are its own thresholds. "
        f"Row 2: all synapses / pairs. Data: {'; '.join(papers)}. −10: {figs['-']}. +10: {figs['+']}.")
fig.text(0.012, 0.04 / H, textwrap.fill(foot, 170), ha="left", va="bottom", fontsize=5.8, color="0.3", linespacing=1.25)
out = os.path.join(HERE, "figV5_chindemi_amount.png")
fig.savefig(out, dpi=300)
print(f"wrote {out}")
print(OC[["mean", "sem", "data_mean", "data_sem", "source"]].to_string())
print(ZN.groupby("proto")[["t_dep", "t_pot"]].mean().to_string())
