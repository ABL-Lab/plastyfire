"""figV1_chindemi_refit.png and figV2_best_B4_k075.png from figV_traces.py outputs (/scratch/dhuruva/figs_veto).

Columns: (a) Sjostrom 0.1 Hz dt -10, (b) Markram 10 Hz -10, (c) Markram 10 Hz +10; one L5->L5 synapse (figV_meta.json).
Rows: 1 pre / post spikes + spine Ca (cai recovered from effcai), 2 c* with theta_d / theta_p, 3 rho over the induction;
figV2 only: 4 eCB pool / P1 at each arrival vs k_E (step / vetoed / below trigger, veto windows shaded), 5 control-lane
release change d = clip(d_eCB + d_NO) with the veto (solid) and veto off (dashed); last row: EPSP ratio, model mean +/-
SEM over pairs vs data mean +/- SEM (B4 and CHR from their GPU val csvs; veto off from the CPU port, figV_outcomes.csv).

MEASURED 22479951: 0:12, 430 MB -> 600M, 0:15.
  python glusynapse_v2/rho_redesign/figs_veto/figV_plot.py [--work /scratch/dhuruva/figs_veto] [--figdir <this dir>]
"""
import argparse, json, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

plt.style.use("/lustre09/project/6070394/dhuruva/plastyfire/onerule.mplstyle")
plt.rcParams["font.sans-serif"] = ["Helvetica Neue", "Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"]
plt.rcParams.update({"font.size": 6.5, "axes.titlesize": 6.8, "axes.labelsize": 6.3, "xtick.labelsize": 5.8,
                     "ytick.labelsize": 5.8, "legend.fontsize": 5.6})
BLUE, ORANGE, RED, GREY = "#68a8e0", "#f9a24b", "#ff4747", "0.55"
INK, PRE, POST = "0.25", "k", "#c0392b"
TD, TP = "#3b6ea8", "#b03a2e"
HERE = os.path.dirname(os.path.abspath(__file__))
KS = ("s01", "m10", "p10")
COLT = dict(s01="a  Sjöström 0.1 Hz, dt −10 ms\n(post 10 ms before pre, single pairings)",
            m10="b  Markram 10 Hz, dt −10 ms\n(post before pre, trains)",
            p10="c  Markram 10 Hz, dt +10 ms\n(pre before post, trains)")


def zoom(pre):
    if len(pre) > 1 and pre[1] - pre[0] < 1000:
        per = pre[1] - pre[0]
        return pre[0] - 30.0, pre[0] + 5 * per + 25.0
    return pre[0] - 40.0, pre[0] + 140.0


class D:
    def __init__(self, work):
        z = np.load(os.path.join(work, "figV_traces.npz"))
        self.z = {k: z[k] for k in z.files}
        self.m = json.load(open(os.path.join(work, "figV_meta.json")))
        self.oc = pd.read_csv(os.path.join(work, "figV_outcomes.csv"))
        self.evf = self.m["ev_fields"]

    def a(self, k, f, fk=None):
        return self.z[f"{k}__{fk}__{f}" if fk else f"{k}__{f}"]

    def sc(self, k, fk):
        return self.m["traces"][f"{k}__{fk}"]

    def ev(self, k, fk, f):
        return self.a(k, "EV", fk)[:, self.evf.index(f)]

    def xz(self, k):
        t, pre = self.a(k, "t"), self.a(k, "pre")
        return t - pre[0], zoom(pre - pre[0])

    def ta(self, k, fk):
        t, pre = self.a(k, "t"), self.a(k, "pre")
        return self.ev(k, fk, "t_a") + t[0] - pre[0]


def ticks(ax, x, y0, y1, col, lw=0.7):
    tr = ax.get_xaxis_transform()
    for v in x:
        ax.plot([v, v], [y0, y1], color=col, lw=lw, transform=tr, solid_capstyle="butt", clip_on=True)


def row_spikes(ax, d, k, first):
    x, (lo, hi) = d.xz(k)
    m = (x >= lo) & (x <= hi)
    cai = d.a(k, "cai")
    ax.plot(x[m], cai[m], color=INK, lw=0.6)
    pre, post = d.a(k, "pre") - d.a(k, "pre")[0], d.a(k, "post") - d.a(k, "pre")[0]
    ticks(ax, pre[(pre >= lo) & (pre <= hi)], 0.90, 1.0, PRE)
    ticks(ax, post[(post >= lo) & (post <= hi)], 0.77, 0.87, POST)
    top = np.nanmax(cai[m])
    ax.set_ylim(0, top * 1.45)
    ax.set_xlim(lo, hi)
    if first:
        ax.set_ylabel("spine Ca (µM)")
        ax.text(0.01, 0.95, "pre", transform=ax.transAxes, fontsize=5.3, va="center", ha="left", color=PRE)
        ax.text(0.01, 0.82, "post", transform=ax.transAxes, fontsize=5.3, va="center", ha="left", color=POST)


def row_c(ax, d, k, fk, first, xl):
    x, (lo, hi) = d.xz(k)
    m = (x >= lo) & (x <= hi)
    c = d.a(k, "c", fk)
    sc = d.sc(k, fk)
    ax.plot(x[m], c[m], color=INK, lw=0.6)
    ax.axhline(sc["td"], color=TD, ls="--", lw=0.6)
    ax.axhline(sc["tp"], color=TP, ls="--", lw=0.6)
    ymax = max(np.nanmax(c[m]), sc["tp"]) * 1.12
    ax.set_ylim(0, ymax)
    ax.text(1.0, sc["td"], " θd", color=TD, fontsize=5.5, va="center", ha="left", transform=ax.get_yaxis_transform())
    ax.text(1.0, sc["tp"], " θp", color=TP, fontsize=5.5, va="center", ha="left", transform=ax.get_yaxis_transform())
    ax.set_xlim(lo, hi)
    if first:
        ax.set_ylabel("effective Ca c*")
    if xl:
        ax.set_xlabel("time from first pre spike (ms)")


def row_rho(ax, d, k, fk, col, first, xl):
    t, pre = d.a(k, "t"), d.a(k, "pre")
    xs = (t - pre[0]) / 1000.0
    r = d.a(k, "rho", fk)
    ax.plot(xs, r, color=col, lw=0.8)
    sc = d.sc(k, fk)
    ax.set_ylim(-0.06, 1.06)
    ax.set_xlim(xs[0], xs[-1])
    ax.text(0.98, 0.5, f"ρ {sc['rho0']:.0f} → {sc['end']['rho']:.2f}", transform=ax.transAxes, ha="right",
            va="center", fontsize=5.6, color="0.2")
    if first:
        ax.set_ylabel("ρ (synapse\nstate)")
    if xl:
        ax.set_xlabel("induction time (s)")


def row_pool(ax, d, k, first, xl):
    x, (lo, hi) = d.xz(k)
    m = (x >= lo) & (x <= hi)
    sc = d.sc(k, "B4")
    pool = d.a(k, "pool", "B4") / sc["uE"]
    ax.plot(x[m], pool[m], color=INK, lw=0.6)
    kE = sc["thE"] / sc["uE"]
    ax.axhline(kE, color="k", ls="--", lw=0.6)
    ax.text(1.0, kE, " k_E", fontsize=5.5, va="center", ha="left", transform=ax.get_yaxis_transform())
    ta = d.ta(k, "B4")
    tv = d.ev(k, "B4", "trigval") / sc["uE"]
    stp, vet, trg = (d.ev(k, "B4", f) > 0.5 for f in ("step", "veto", "trig"))
    inz = (ta >= lo) & (ta <= hi)
    for a, v in zip(ta[inz], vet[inz]):
        ax.axvspan(a + sc["vt0"], a + sc["vtv"], color=ORANGE if v else "0.85", alpha=0.45 if v else 0.5, lw=0, zorder=0)
    post = d.a(k, "post") - d.a(k, "pre")[0]
    ticks(ax, post[(post >= lo) & (post <= hi)], 0.88, 0.98, POST)
    s = inz & stp; v = inz & vet; n = inz & ~trg
    ax.plot(ta[s], tv[s], "v", color=RED, ms=3.6, zorder=5)
    ax.plot(ta[v], tv[v], "o", mfc="white", mec="#1f4e79", mew=0.8, ms=3.4, zorder=5)
    ax.plot(ta[n], tv[n], ".", color=GREY, ms=3.0, zorder=5)
    ymax = max(np.nanmax(pool[m]), np.nanmax(tv[inz]) if inz.any() else 0, kE) * 1.35
    ax.set_ylim(0, ymax)
    ax.set_xlim(lo, hi)
    N = len(ta)
    ax.text(1.0, 1.03, f"all {N} arrivals: {int(stp.sum())} step, {int(vet.sum())} vetoed, {int((~trg).sum())} below",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=5.0, color="0.15",
            bbox=dict(fc="0.96", ec="0.7", lw=0.4, pad=1.2, boxstyle="round,pad=0.25"))
    if first:
        ax.set_ylabel("eCB pool / P1")
    if xl:
        ax.set_xlabel("time from first pre spike (ms)")


def row_d(ax, d, k, first, xl):
    t, pre = d.a(k, "t"), d.a(k, "pre")
    xs = (t - pre[0]) / 1000.0
    ax.plot(xs, d.a(k, "dctl", "B4"), color=BLUE, lw=0.9, label="with veto")
    ax.plot(xs, d.a(k, "dctl", "KOV"), color=RED, lw=0.9, ls="--", label="veto off")
    dmin = d.sc(k, "B4")["dmin"]
    ax.axhline(dmin, color=GREY, ls=":", lw=0.6)
    ax.axhline(0, color="0.85", lw=0.5, zorder=0)
    ax.set_ylim(dmin - 0.05, 0.2)
    ax.set_xlim(xs[0], xs[-1])
    if first:
        ax.set_ylabel("release\nchange d")
        ax.legend(loc="upper right", frameon=False, handlelength=1.8, borderaxespad=0.1, fontsize=5.0, ncol=2)
        ax.text(0.98, dmin + 0.008, f"floor {dmin:g}", transform=ax.get_yaxis_transform(), ha="right",
                va="bottom", fontsize=5.0, color=GREY)
    if xl:
        ax.set_xlabel("induction time (s)")


def outcome_vals(d, k, which):
    oc = d.oc
    b = oc[(oc.proto == k) & (oc.fit == "B4")].iloc[0]
    v = dict(data=(b.data_mean, b.data_sem))
    for fk in which:
        r = oc[(oc.proto == k) & (oc.fit == fk)].iloc[0]
        if fk in ("B4", "CHR"):
            v[fk] = (r.gpu_pred, r.gpu_sem)
        else:
            v[fk] = (r.cpu_mean, r.cpu_sem)
    return v


def row_out(ax, d, k, which, first, ylim):
    v = outcome_vals(d, k, which)
    lab = dict(data="data", B4="B4", KOV="veto off", CHR="Chindemi\nrefit")
    col = dict(data="k", B4=BLUE, KOV=RED, CHR=ORANGE)
    names = ["data"] + list(which)
    for j, nm in enumerate(names):
        mu, se = v[nm]
        open_ = nm == "KOV"
        ax.errorbar(j, mu, yerr=se, fmt="o", color=col[nm], mfc="white" if open_ else col[nm], mec=col[nm], ms=3.6,
                    capsize=1.5, lw=0.8, mew=0.9)
        ax.text(j + 0.13, mu, f"{mu:.2f}", fontsize=5.4, va="center", ha="left", color="0.15")
    ax.axhline(1.0, color="0.8", lw=0.5, zorder=0)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels([lab[n] for n in names])
    ax.set_xlim(-0.5, len(names) - 0.3)
    ax.set_ylim(*ylim)
    if first:
        ax.set_ylabel("EPSP ratio\n(mean ± SEM)")


def make(d, kind, figdir):
    best = kind == "B4"
    rows = ["spk", "c", "rho"] + (["pool", "d"] if best else []) + ["out"]
    # explicit layout in inches (no layout engine): caption band | column titles | rows | footer band
    W = 7.0
    hr = dict(spk=0.85, c=0.8, rho=0.72, pool=0.95, d=0.72, out=0.85)       # axes heights (in)
    CAP, TIT, BOT = 0.58, 0.36, 0.62                                        # caption band, titles, bottom band (in)
    LEFT, RIGHT, CGAP = 0.62, 0.32, 0.42
    cw = (W - LEFT - RIGHT - 2 * CGAP) / 3
    gap = {}                                                                # gap below each row (in)
    for j, r in enumerate(rows[:-1]):
        gap[r] = 0.40
        if rows[j + 1] == "pool":
            gap[r] = 0.50                                                   # room for the arrival-count box
        if r == "pool":
            gap[r] = 0.62                                                   # room for the row-4 legend
    H = CAP + TIT + sum(hr[r] for r in rows) + sum(gap.values()) + BOT
    fig = plt.figure(figsize=(W, H), layout="none")
    ytop = {}
    y = H - CAP - TIT
    for r in rows:
        ytop[r] = y
        y -= hr[r] + gap.get(r, 0.0)

    def axpos(r, jc):
        return [(LEFT + jc * (cw + CGAP)) / W, (ytop[r] - hr[r]) / H, cw / W, hr[r] / H]
    fk = "B4" if best else "CHR"
    col = BLUE if best else ORANGE
    which = ("B4", "KOV") if best else ("CHR",)
    allv = [x for k in KS for (mu, se) in outcome_vals(d, k, which).values() for x in (mu - se, mu + se)]
    ylim = (min(0.5, min(allv) - 0.06), max(1.35, max(allv) + 0.08))
    for jc, k in enumerate(KS):
        first = jc == 0; mid = jc == 1
        for jr, rw in enumerate(rows):
            ax = fig.add_axes(axpos(rw, jc))
            if rw == "spk":
                row_spikes(ax, d, k, first)
                ax.set_title(COLT[k], loc="left", fontsize=6.6, pad=3)
            elif rw == "c":
                row_c(ax, d, k, fk, first, mid)
            elif rw == "rho":
                row_rho(ax, d, k, fk, col, first, mid)
            elif rw == "pool":
                row_pool(ax, d, k, first, mid)
            elif rw == "d":
                row_d(ax, d, k, first, mid)
            else:
                row_out(ax, d, k, which, first, ylim)
    # caption
    ci = d.m["choice"]
    sel = (f"Same L5→L5 synapse in every panel (pair {ci['pair']}, synapse {ci['syn']}): the median veto effect among "
           f"{ci['n_cand']} synapses that fit the story under B4 ({ci['mask']}).")
    V = {k: outcome_vals(d, k, which) for k in KS}
    if best:
        trg = d.ev("p10", "B4", "trig") > 0.5
        j1 = int(np.flatnonzero(trg)[0]) + 1 if trg.any() else None
        cap = [f"Best model (B4, veto_peak_k 0.75, d_NO_max 0.14). Row 4: the eCB pool reads only recent bAP Ca, so in the +10 train (c) "
               f"it is above the trigger from arrival {j1} on,",
               "although every bAP comes after its pre spike. The veto cancels those steps because the bAP's own VDCC "
               "current peaks 2–17 ms after the arrival (orange windows).",
               f"Without the veto (dashed, same parameters) release steps down and the +10 outcome drops from "
               f"{V['p10']['B4'][0]:.2f} to {V['p10']['KOV'][0]:.2f} (data {V['p10']['data'][0]:.2f})."]
        fn = "figV2_best_B4_k075.png"
    else:
        oc = pd.read_csv("/scratch/dhuruva/latest_figs/variants/CHR_s5_val.csv")
        mg = oc[(oc.target == "10Hz_-10ms") & (oc.condition == "mglu_block")].iloc[0]
        cap = ["Chindemi's rule refitted (ρ only: a00–a11, γd, γp fitted, apical tied to basal; no presynaptic release "
               "step). Its only route to LTD is ρ falling while c* sits between θd and θp.",
               f"The single 0.1 Hz pairing (a) gives {V['s01']['CHR'][0]:.2f} (data {V['s01']['data'][0]:.2f}): no LTD. "
               f"At 10 Hz −10 (b) it reaches {V['m10']['CHR'][0]:.2f} (data {V['m10']['data'][0]:.2f}) through ρ alone,",
               f"so blocking mGluR / eCB cannot remove it (model {mg.pred:.2f} vs data {mg.target_mean:.2f} under "
               f"mGluR block)."]
        fn = "figV1_chindemi_refit.png"
    if best:                                                                # row-4 legend, own band below row 4
        h = [Line2D([], [], ls="none", marker="v", color=RED, ms=3.6, label="step"),
             Line2D([], [], ls="none", marker="o", mfc="white", mec="#1f4e79", ms=3.4, label="vetoed"),
             Line2D([], [], ls="none", marker=".", color=GREY, ms=3.0, label="below trigger"),
             Line2D([], [], color="k", ls="--", lw=0.6, label="trigger k_E·P1"),
             Patch(color=ORANGE, alpha=0.45, lw=0, label="veto window, bAP inside"),
             Patch(color="0.85", alpha=0.5, lw=0, label="veto window, no bAP"),
             Line2D([], [], color=POST, lw=0.8, label="post spike (bAP)")]
        fig.legend(handles=h, loc="lower center", bbox_to_anchor=(0.5, (ytop["pool"] - hr["pool"] - 0.58) / H),
                   ncol=7, frameon=False, handlelength=1.2, columnspacing=1.0, fontsize=5.3)
    fig.text(0.012, 1 - 0.07 / H, "\n".join(cap), ha="left", va="top", fontsize=6.0, linespacing=1.3)
    import textwrap
    foot = textwrap.fill(sel + " Model: mean ± SEM over simulated pairs; data: mean ± SEM.", 175)
    fig.text(0.012, 0.04 / H, foot, ha="left", va="bottom", fontsize=5.0, color="0.3", linespacing=1.25)
    out = os.path.join(figdir, fn)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"wrote {out} ({7.0:.1f} x {H:.2f} in)", flush=True)
    print("  caption:\n    " + "\n    ".join(cap), flush=True)


def report(d):
    print("\nOUTCOMES (as plotted):")
    for k in KS:
        v = outcome_vals(d, k, ("B4", "KOV", "CHR"))
        print(f"  {k}: " + ", ".join(f"{n} {mu:.3f} ± {se:.3f}" for n, (mu, se) in v.items()))
    oc = d.oc
    for _, r in oc.iterrows():
        if r.fit in ("B4", "CHR") and abs(r.cpu_mean - r.gpu_pred) > 0.01:
            print(f"  WARNING {r.proto} {r.fit}: CPU port mean {r.cpu_mean:.3f} vs GPU val {r.gpu_pred:.3f}")
    for k in KS:
        sc = d.sc(k, "B4"); t, pre = d.a(k, "t"), d.a(k, "pre")
        ta = d.ta(k, "B4"); tv = d.ev(k, "B4", "trigval") / sc["uE"]
        qa = d.a(k, "qa", "B4") / sc["thP"]
        f = {n: d.ev(k, "B4", n) > 0.5 for n in ("trig", "veto", "step", "annul")}
        post = d.a(k, "post") - pre[0]
        nshow = 12 if k != "s01" else 3
        print(f"\n{k} B4 per arrival (k_E {sc['thE'] / sc['uE']:.3f}; veto when VDCC peak in (t_a+{sc['vt0']:g}, "
              f"t_a+{sc['vtv']:g}] > {sc['k_peak']:g} I1):")
        for j in range(min(nshow, len(ta))):
            nb = post[(post > ta[j] - 100) & (post < ta[j] + 20)] - ta[j]
            st = "STEP" if f["step"][j] else ("VETO" if f["veto"][j] else ("annul" if f["annul"][j] else
                                                                            ("trig?" if f["trig"][j] else "below")))
            print(f"  {j + 1:2d}  t_a {ta[j]:8.1f} ms  pool/P1 {tv[j]:6.3f}  {st:5s}  VDCC peak/thP {qa[j]:5.2f}  "
                  f"bAPs at t_a{''.join(f' {x:+.1f}' for x in nb)} ms")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default="/scratch/dhuruva/figs_veto")
    ap.add_argument("--figdir", default=HERE)
    A = ap.parse_args()
    d = D(A.work)
    report(d)
    make(d, "CHR", A.figdir)
    make(d, "B4", A.figdir)
