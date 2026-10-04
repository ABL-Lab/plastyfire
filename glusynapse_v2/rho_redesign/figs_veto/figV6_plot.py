"""figV6_chindemi_all.png: the Chindemi refit on every L5 -> L5 protocol, in one figure.

Row 1: c* of one synapse, Sjostrom 2001 -10 (red) vs +10 (blue) at 0.1 / 10 / 20 / 40 / 50 Hz, with theta_d / theta_p.
Row 2: time per synapse in the depression / potentiation zones; EPSP ratio vs frequency (data vs refit).
Row 3: every other L5 control protocol: data vs refit (CHR_s5_val.csv).
Inputs: figV5_traces.py with FREQS=0.1,10,20,40,50 FIGTAG=figV6. Data values and sources come from the scored
tables (CHR_s5_val.csv; sources from ebner_targets.csv, else the source string in targets.py), never typed by hand.

  python glusynapse_v2/rho_redesign/figs_veto/figV6_plot.py
"""
import json, os, re, textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from glusynapse_v2 import targets as TG

plt.style.use("/lustre09/project/6070394/dhuruva/plastyfire/onerule.mplstyle")
plt.rcParams["font.sans-serif"] = ["Helvetica Neue", "Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"]
plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7.2, "ytick.labelsize": 7.2,
                     "legend.fontsize": 7.2, "axes.titlesize": 8.5})
RED, BLUE, ORANGE, DEP, POT = "#e0473c", "#2f6fb5", "#f39c32", "#fbe1de", "#dff0e4"
WORK = "/scratch/dhuruva/figs_veto"
HERE = os.path.dirname(os.path.abspath(__file__))
EBNER = "/lustre09/project/6070394/dhuruva/plastyfire/ebner/ebner_targets.csv"
VAL = "/scratch/dhuruva/latest_figs/variants/CHR_s5_val.csv"
TAG = "figV6"
FR = ("0.1", "10", "20", "40", "50")
z = np.load(os.path.join(WORK, f"{TAG}.npz")); Z = {k: z[k] for k in z.files}
M = json.load(open(os.path.join(WORK, f"{TAG}_meta.json")))
ZN = pd.read_csv(os.path.join(WORK, f"{TAG}_zones.csv"))
OC = pd.read_csv(os.path.join(WORK, f"{TAG}_outcomes.csv")).set_index("proto")
V = pd.read_csv(VAL); V = V[V.condition == "control"].set_index("target")
EB = pd.read_csv(EBNER)                                   # may hold several rows per protocol (other papers)

# other L5 protocols (row 3), grouped
GROUPS = [("Markram 10 Hz", ["10Hz_-10ms", "10Hz_10ms"]),
          ("Markram +5 ms, frequency", ["2Hz_5ms", "5Hz_5ms", "10Hz_5ms", "20Hz_5ms", "30Hz_5ms", "40Hz_5ms"]),
          ("Sjöström Δt 0, ±25", ["sjostrom_0.1hz_dt0ms", "sjostrom_20hz_dt-25ms", "sjostrom_20hz_dt0ms",
                                 "sjostrom_20hz_dt+25ms", "sjostrom_40hz_dt0ms"]),
          ("Sjöström 0.1 Hz timing", ["sjostrom_0.1hz_dt-25ms", "sjostrom_0.1hz_dt-120ms", "sjostrom_0.1hz_dt-200ms",
                                      "sjostrom_burst5x20hz_r50_dt-120ms", "sjostrom_burst5x20hz_r50_dt-200ms"]),
          ("Sjöström 2007 step", ["sjostrom07_step200ms_pair", "sjostrom07_step200ms_pre_only",
                                  "sjostrom07_step200ms_post_only"])]
LAB = {"10Hz_-10ms": "−10", "10Hz_10ms": "+10", "2Hz_5ms": "2", "5Hz_5ms": "5", "10Hz_5ms": "10", "20Hz_5ms": "20",
       "30Hz_5ms": "30", "40Hz_5ms": "40", "sjostrom_0.1hz_dt0ms": "0.1 Hz 0", "sjostrom_20hz_dt-25ms": "20 −25",
       "sjostrom_20hz_dt0ms": "20 0", "sjostrom_20hz_dt+25ms": "20 +25", "sjostrom_40hz_dt0ms": "40 0",
       "sjostrom_0.1hz_dt-25ms": "−25", "sjostrom_0.1hz_dt-120ms": "−120", "sjostrom_0.1hz_dt-200ms": "−200",
       "sjostrom_burst5x20hz_r50_dt-120ms": "burst −120", "sjostrom_burst5x20hz_r50_dt-200ms": "burst −200",
       "sjostrom07_step200ms_pair": "pair", "sjostrom07_step200ms_pre_only": "pre only",
       "sjostrom07_step200ms_post_only": "post only"}


def source(p):
    """(paper, paper + figure) for protocol p: ebner_targets.csv first, else the source string in targets.py."""
    rows = EB[EB.protocol_id == p]
    if p in V.index:                                       # the row whose value is the scored one
        rows = rows[(rows.mean_ratio - V.loc[p].target_mean).abs() < 0.006]
    if len(rows):
        r = rows.iloc[0]
        return r.source_paper, f"{r.source_paper}, {r.figure_or_table}"
    for v in vars(TG).values():
        if isinstance(v, dict):
            for key in ((p, "control"), p):
                t = v.get(key)
                if isinstance(t, tuple) and t and isinstance(t[-1], str):
                    s = t[-1]
                    if s.startswith("as "):                # pooled bar: "as -120 (pooled bar)" -> that row's source
                        tok = s[3:].split(" ")[0]
                        ref = [kk for kk in v if (kk[0] if isinstance(kk, tuple) else kk) != p and tok in
                               (kk[0] if isinstance(kk, tuple) else kk) and
                               (not isinstance(kk, tuple) or kk[1] == "control")]
                        if ref:
                            ref.sort(key=lambda kk: -len(os.path.commonprefix([kk[0] if isinstance(kk, tuple) else kk, p])))
                            s = v[ref[0]][-1] + " [pooled bar]"
                    m = re.match(r"^([A-Z][A-Za-z]+(?: et al\.)? \d{4})", s)
                    return (m.group(1) if m else s), s
    return "SOURCE NOT FOUND", "SOURCE NOT FOUND"


def ticks(ax, x, y0, y1, col):
    tr = ax.get_xaxis_transform()
    for v in x:
        ax.plot([v, v], [y0, y1], color=col, lw=0.8, transform=tr, solid_capstyle="butt")


for k in OC.index:                                          # data + source from the scored table, not the trace job
    q = OC.loc[k, "target"]
    OC.loc[k, "data_mean"], OC.loc[k, "data_sem"] = float(V.loc[q].target_mean), float(V.loc[q].target_sem)
    OC.loc[k, "source"] = source(q)[1]
OC["data_mean"] = OC.data_mean.astype(float); OC["data_sem"] = OC.data_sem.astype(float)
o = {k: OC.loc[k] for k in OC.index}
tot = {k: (ZN[ZN.proto == k].t_dep + ZN[ZN.proto == k].t_pot).mean() for k in OC.index}
more = [f for f in FR if tot[f + "+"] >= tot[f + "-"]]
CAPTXT = ["Chindemi's ρ reads one signal, c* (spine Ca): more pairings per second → more Ca. +10 gives at least as "
          f"much Ca as −10 at {'every frequency' if len(more) == len(FR) else ', '.join(more) + ' Hz'} (row 2 left).",
          "The data split by order; the refit follows frequency (row 2 right). Row 3: every other L5 protocol."]
CAPW = "\n".join(textwrap.fill(c, 140, subsequent_indent="   ") for c in CAPTXT)

# ---------------------------------------------------------------- layout (inches)
W = 9.0
CAP = 0.14 + 0.15 * (CAPW.count("\n") + 1)
L, R = 0.62, 0.30
H1, H2, H3 = 1.35, 1.75, 1.55
HD1, HD2, HD3, G = 0.72, 0.55, 0.40, 0.55
FOOTH = 0.62
H = CAP + HD1 + H1 + G + HD2 + H2 + G + HD3 + H3 + 0.55 + FOOTH
fig = plt.figure(figsize=(W, H), layout="none")
y1 = H - CAP - HD1 - H1
g1 = 0.22; cw = (W - L - R - 4 * g1) / 5
A1 = [fig.add_axes([(L + j * (cw + g1)) / W, y1 / H, cw / W, H1 / H]) for j in range(5)]
y2 = y1 - G - HD2 - H2
w2 = (W - L - R - 0.8) / 2
A2 = fig.add_axes([L / W, y2 / H, w2 / W, H2 / H])
A3 = fig.add_axes([(L + w2 + 0.8) / W, y2 / H, w2 / W, H2 / H])
y3 = y2 - G - HD3 - H3
A4 = fig.add_axes([L / W, y3 / H, (W - L - R) / W, H3 / H])
fig.text(0.012, 1 - 0.07 / H, CAPW, ha="left", va="top", fontsize=8, linespacing=1.3)


def head(y, txt):
    fig.text(L / W, y / H, txt, ha="left", va="bottom", fontsize=8.8, fontweight="bold")


head(y1 + H1 + HD1 - 0.22, "1  c* of one synapse, Sjöström ±10 trains (first pairings)")
head(y2 + H2 + HD2 - 0.22, "2  Whole induction, all synapses / pairs")
head(y3 + H3 + HD3 - 0.22, "3  Every other L5→L5 protocol")
fig.legend(handles=[Line2D([], [], color=RED, lw=1.1, label="−10 (post first)"),
                    Line2D([], [], color=BLUE, lw=1.1, label="+10 (pre first)"),
                    Patch(color=DEP, label="depression zone (θd–θp)"), Patch(color=POT, label="potentiation zone (>θp)")],
           loc="lower left", bbox_to_anchor=(L / W, (y1 + H1 + 0.24) / H), ncol=4, frameon=False)

# ---- row 1
ymax = max(float(np.nanmax(Z[f"{f}{s}__c"][(Z[f"{f}{s}__t"] >= -40) &
                                          (Z[f"{f}{s}__t"] <= (260 if f == "0.1" else 4000 / float(f) + 60))]))
           for f in FR for s in "-+" if f"{f}{s}__t" in Z)
top = max(ymax, max(M["traces"][f + "-"]["tp"] for f in FR)) * 1.5
for j, f in enumerate(FR):
    ax = A1[j]; s0 = M["traces"][f + "-"]
    lo, hi = (-40.0, 260.0) if f == "0.1" else (-30.0, 4 * 1000.0 / float(f) + 60.0)
    ax.axhspan(s0["td"], s0["tp"], color=DEP, lw=0, zorder=0); ax.axhspan(s0["tp"], top, color=POT, lw=0, zorder=0)
    for s, col, yy in (("-", RED, (0.90, 0.97)), ("+", BLUE, (0.80, 0.87))):
        k = f + s
        t, c = Z[f"{k}__t"], Z[f"{k}__c"]; m = (t >= lo) & (t <= hi)
        ax.plot(t[m], c[m], color=col, lw=1.0)
        post = Z[f"{k}__post"]; ticks(ax, post[(post >= lo) & (post <= hi)], *yy, col)
    pre = Z[f"{f}-__pre"]; ticks(ax, pre[(pre >= lo) & (pre <= hi)], 0.70, 0.77, "k")
    ax.axhline(s0["td"], color="0.35", ls="--", lw=0.7); ax.axhline(s0["tp"], color="0.35", ls="-.", lw=0.7)
    ax.set_xlim(lo, hi); ax.set_ylim(0, top); ax.set_title(f"{f} Hz", pad=3)
    ax.set_xlabel("ms")
    if j == 0:
        ax.set_ylabel("c*")
        ax.text(0.97, 0.735, "pre", transform=ax.transAxes, fontsize=6.3, ha="right", va="center")
    else:
        ax.tick_params(labelleft=False)
    if j == 4:
        ax.text(1.03, s0["td"], "θd", transform=ax.get_yaxis_transform(), fontsize=7.2, va="center")
        ax.text(1.03, s0["tp"], "θp", transform=ax.get_yaxis_transform(), fontsize=7.2, va="center")

# ---- row 2 left: time in zones
x = np.arange(len(FR)); bw = 0.36
for s, col, dx in (("-", RED, -bw / 2), ("+", BLUE, bw / 2)):
    q = [ZN[ZN.proto == f + s] for f in FR]
    dep = [a.t_dep.mean() for a in q]; pot = [a.t_pot.mean() for a in q]
    A2.bar(x + dx, dep, bw, color=col, alpha=0.45, lw=0)
    A2.bar(x + dx, pot, bw, bottom=dep, color=col, lw=0, hatch="////", edgecolor="white")
A2.set_xticks(x); A2.set_xticklabels([f"{f} Hz" for f in FR])
A2.set_ylabel("time per synapse (s)"); A2.set_title("time c* spends in each zone", fontsize=8)
A2.legend(handles=[Patch(color="0.5", alpha=0.45, label="depression zone"),
                   Patch(facecolor="0.3", hatch="////", edgecolor="white", label="potentiation zone")],
          loc="upper left", frameon=False, fontsize=7)

# ---- row 2 right: outcome vs frequency
for s, col, lab in (("-", RED, "−10"), ("+", BLUE, "+10")):
    d = [o[f + s].data_mean for f in FR]; de = [o[f + s].data_sem for f in FR]
    mm = [o[f + s]["mean"] for f in FR]; me = [o[f + s]["sem"] for f in FR]
    A3.errorbar(x - 0.07, d, yerr=de, color=col, marker="o", ms=4.5, lw=1.4, capsize=2, label=f"data {lab}")
    A3.errorbar(x + 0.07, mm, yerr=me, color=col, marker="o", mfc="white", ms=4.5, lw=1.1, ls="--", capsize=2,
                label=f"Chindemi refit {lab}")
A3.axhline(1, color="0.8", lw=0.6, zorder=0)
A3.set_xticks(x); A3.set_xticklabels([f"{f} Hz" for f in FR]); A3.set_xlim(-0.4, len(FR) - 0.6)
A3.set_ylabel("EPSP ratio"); A3.set_title("outcome, Sjöström ±10", fontsize=8)
A3.legend(loc="upper left", frameon=False, fontsize=6.8, ncol=2)
lo_ = min(min(o[k].data_mean - o[k].data_sem, o[k]["mean"] - o[k]["sem"]) for k in o)
hi_ = max(max(o[k].data_mean + o[k].data_sem, o[k]["mean"] + o[k]["sem"]) for k in o)
A3.set_ylim(lo_ - 0.05, hi_ + 0.35)

# ---- row 3: every other L5 protocol
xi, xt, xl, gx = 0, [], [], []
for gname, ps in GROUPS:
    x0 = xi
    for p in ps:
        if p not in V.index:
            continue
        r = V.loc[p]
        A4.plot([xi, xi], [r.target_mean, r.pred], color="0.75", lw=1.0, zorder=1)
        A4.errorbar(xi - 0.08, r.target_mean, yerr=r.target_sem, fmt="o", color="k", ms=4.2, capsize=1.5, lw=0.9,
                    zorder=3)
        A4.errorbar(xi + 0.08, r.pred, yerr=r.pred_sem, fmt="o", color=ORANGE, mfc="white", ms=4.2, capsize=1.5,
                    lw=0.9, mew=1.1, zorder=3)
        xt.append(xi); xl.append(LAB.get(p, p)); xi += 1
    gx.append((x0, xi - 1, gname)); xi += 0.8
for a, b, gname in gx:
    A4.text((a + b) / 2, 1.02, gname, transform=A4.get_xaxis_transform(), ha="center", va="bottom", fontsize=7.2,
            fontweight="bold")
    if a > 0:
        A4.axvline(a - 0.9, color="0.85", lw=0.6)
A4.axhline(1, color="0.8", lw=0.6, zorder=0)
A4.set_xticks(xt); A4.set_xticklabels(xl, fontsize=6.6, rotation=35, ha="right", rotation_mode="anchor")
A4.set_xlim(-0.6, xi - 0.6); A4.set_ylabel("EPSP ratio")
A4.legend(handles=[Line2D([], [], ls="none", marker="o", color="k", ms=4.2, label="data"),
                   Line2D([], [], ls="none", marker="o", mfc="white", mec=ORANGE, ms=4.2, label="Chindemi refit")],
          loc="upper left", frameon=False, fontsize=7, ncol=2, bbox_to_anchor=(0.0, 0.98))

# ---- sources (generated)
srcs = [(p, *source(p)) for p in [f"sjostrom_{f}hz_dt{s}10ms" for f in FR for s in "-+"] +
        [p for _, ps in GROUPS for p in ps]]
pd.DataFrame(srcs, columns=["protocol", "paper", "source"]).to_csv(os.path.join(HERE, f"{TAG}_sources.csv"),
                                                                   index=False)
papers = list(dict.fromkeys(s[1] for s in srcs))
foot = (f"L5→L5, Chindemi refit (CHR_s5), control. Row 1: pair {M['traces']['0.1-']['pair']}, synapse "
        f"{M['traces']['0.1-']['syn']}, its own θd/θp. Data sources per protocol (paper + figure): {TAG}_sources.csv "
        f"next to this figure. Papers: {'; '.join(papers)}.")
fig.text(0.012, 0.04 / H, textwrap.fill(foot, 190), ha="left", va="bottom", fontsize=6.0, color="0.3", linespacing=1.25)
out = os.path.join(HERE, f"{TAG}_chindemi_all.png")
fig.savefig(out, dpi=300)


# ---- every plotted number to JSON, so the figure can be restyled without rerunning the simulation
def fl(a, nd=5):
    return [round(float(v), nd) for v in np.asarray(a, float).ravel()]


src_of = {p: (pa, s) for p, pa, s in srcs}
D = dict(figure=os.path.basename(out), fit="Chindemi refit CHR_s5 (/scratch/dhuruva/latest_figs/variants/chr/CHR_s5.json)",
         pathway="L5->L5", condition="control", caption=CAPTXT, footnote=foot,
         units=dict(t="ms from first pre spike", c="c* (low-passed spine Ca, model units)", time_in_zone="s per synapse",
                    ratio="EPSP ratio (after/before)"),
         row1=dict(synapse=dict(pair=M["traces"]["0.1-"]["pair"], syn=M["traces"]["0.1-"]["syn"],
                                rho0=M["traces"]["0.1-"]["rho0"]), traces={}),
         row2=dict(time_in_zone={}, outcome={}), row3=[])
for f in FR:
    lo, hi = (-40.0, 260.0) if f == "0.1" else (-30.0, 4 * 1000.0 / float(f) + 60.0)
    for s in "-+":
        k = f + s
        t, c = Z[f"{k}__t"], Z[f"{k}__c"]; m = np.flatnonzero((t >= lo) & (t <= hi))
        m = m[::max(1, len(m) // 2000)]                                    # <= ~2000 points per trace
        pre, post = Z[f"{k}__pre"], Z[f"{k}__post"]
        D["row1"]["traces"][k] = dict(protocol=f"sjostrom_{f}hz_dt{s}10ms",
                                      freq_hz=float(f), dt_ms=(-10 if s == "-" else 10), window_ms=[lo, hi],
                                      t=fl(t[m], 3), c=fl(c[m]), pre=fl(pre[(pre >= lo) & (pre <= hi)], 3),
                                      post=fl(post[(post >= lo) & (post <= hi)], 3),
                                      theta_d=M["traces"][k]["td"], theta_p=M["traces"][k]["tp"])
        q = ZN[ZN.proto == k]
        D["row2"]["time_in_zone"][k] = dict(dep_mean=float(q.t_dep.mean()), dep_sem=float(q.t_dep.std(ddof=1) / np.sqrt(len(q))),
                                            pot_mean=float(q.t_pot.mean()), pot_sem=float(q.t_pot.std(ddof=1) / np.sqrt(len(q))),
                                            n_syn=int(len(q)))
        r = OC.loc[k]
        D["row2"]["outcome"][k] = dict(protocol=r.target, model_mean=float(r["mean"]), model_sem=float(r["sem"]),
                                       n_pairs=int(r.n), data_mean=float(r.data_mean), data_sem=float(r.data_sem),
                                       source=r.source)
for gname, ps in GROUPS:
    for p in ps:
        if p in V.index:
            r = V.loc[p]
            D["row3"].append(dict(group=gname, protocol=p, label=LAB.get(p, p), data_mean=float(r.target_mean),
                                  data_sem=float(r.target_sem), model_mean=float(r.pred), model_sem=float(r.pred_sem),
                                  n_pairs=int(r.n_pairs), paper=src_of[p][0], source=src_of[p][1]))
jout = os.path.join(HERE, f"{TAG}_data.json")
json.dump(D, open(jout, "w"), indent=1)
print(f"wrote {jout}")
print(f"wrote {out}")
print(pd.DataFrame(srcs, columns=["protocol", "paper", "source"]).to_string(index=False))
print({k: round(v, 3) for k, v in tot.items()})
