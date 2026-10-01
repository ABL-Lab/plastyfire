"""Figures of traces_v5.py (v5c rule, fit v5_S1_V5c_s6_39, emodel delta-split1): one figure per group
(results/traces_S1_V5c_<group>.png / .pdf) and an overview (results/traces_S1_V5c_overview.png / .pdf).

Per record (one column): (a) soma V, (b) spine Ca, (c) effective spine Ca / theta_p, (d) spine VDCC Ca pool with the
gate threshold theta_V, (e) eCB trigger pool W with theta_eCB and the eCB steps, all over one induction repetition;
(f) rho over the whole induction; (g) weight ratio: data +/- SEM, offline fit (all the fit's pairs, fit csv pred),
offline kernel and BCL live on this pair. BCL live solid, offline kernel dashed. Synapses: the two whose peak effective
Ca is nearest theta_p and the one farthest from it (chosen on the group's first record, same synapses in every column).

    python glusynapse_v2/bcl_validation/plot_traces_v5.py [--out /scratch/dhuruva/bcl_traces_S1_V5c]
"""
import argparse, glob, json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
OUT = "/scratch/dhuruva/bcl_traces_S1_V5c"
TAG = "S1_V5c"
SYN_COL = ("#0072B2", "#D55E00", "#009E73")              # Okabe-Ito blue, vermillion, green
GROUP_ORDER = ("markram10hz", "sj01freq", "sj03r50", "sj07step", "letzkus")
GROUP_SHORT = {"markram10hz": "Markram 1997", "sj01freq": "Sjöström 2001", "sj03r50": "Sjöström 2003",
               "sj07step": "Sjöström 2007", "letzkus": "Letzkus 2006"}
YLAB = ("soma V (mV)", "spine Ca ($\\mu$M)", "eff. spine Ca\n/ $\\theta_p$", "spine VDCC\nCa pool (a.u.)",
        "eCB trigger\npool W (a.u.)", "$\\rho$", "weight ratio")
LIVE_LS, OFF_LS = "-", (0, (3, 1.6))


def style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
        "legend.fontsize": 6, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5, "lines.linewidth": 0.9, "axes.spines.top": False,
        "axes.spines.right": False, "hatch.linewidth": 0.5, "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.dpi": 300, "mathtext.default": "regular"})
    return plt


def load(out):
    recs = []
    for f in sorted(glob.glob(os.path.join(out, "*__*.npz"))):
        z = np.load(f)
        d = {k: z[k] for k in z.files if k != "meta"}
        d["meta"] = json.loads(str(z["meta"]))
        recs.append(d)
    recs.sort(key=lambda d: (GROUP_ORDER.index(d["meta"]["group"]) if d["meta"]["group"] in GROUP_ORDER else 99,
                             d["meta"]["member"]))
    return recs


def pick_syn(d):
    """-> [(gid, tag, peak / theta_p)]: 2 nearest theta_p (|log peak/theta_p|), 1 farthest."""
    tp = np.asarray(d["tp"], float); mx = np.asarray(d["off_maxE"], float)
    ok = tp > 0
    m = np.full(len(tp), np.inf)
    m[ok] = np.abs(np.log(np.maximum(mx[ok], 1e-30) / tp[ok]))
    order = [int(i) for i in np.argsort(m) if np.isfinite(m[i])]
    near = order[:2]; rest = [i for i in order if i not in near]
    sel = [(i, "near $\\theta_p$") for i in near] + ([(rest[-1], "far from $\\theta_p$")] if rest else [])
    return [(int(d["syn"][i]), tag, float(mx[i] / tp[i])) for i, tag in sel]


def row_of(d, gid):
    w = np.flatnonzero(np.asarray(d["syn"]) == gid)
    return int(w[0]) if len(w) else None


def live_events(t, necb):
    """times where the eCB step counter increases."""
    if necb.size < 2:
        return np.zeros(0)
    return t[1:][np.diff(necb) > 0.5]


def off_events(d, i):
    s = np.asarray(d["off_ev_syn"]); return np.asarray(d["off_ev_t"])[s == i]


def letter(ax, s):
    ax.text(-0.30, 1.04, s, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom", ha="left")


def bars(ax, m, show_ylabel):
    vals = [m["exp"], m["pred"], m["ratio_off"], m.get("ratio_live", np.nan)]
    lab = ["data", "offline\nall pairs", "offline\nthis pair", "BCL\nthis pair"]
    st = [dict(fc="0.80", ec="0.30"), dict(fc="white", ec="0.20"), dict(fc="white", ec="0.20", hatch="////"),
          dict(fc="0.20", ec="0.20")]
    x = np.arange(4)
    for k in range(4):
        ax.bar(x[k], vals[k], width=0.7, lw=0.6, **st[k])
    ax.errorbar(x[0], m["exp"], yerr=m["sem"], fmt="none", ecolor="k", elinewidth=0.8, capsize=2)
    top = np.nanmax([m["exp"] + m["sem"]] + vals[1:])
    for k in range(4):
        if np.isfinite(vals[k]):
            yy = (m["exp"] + m["sem"]) if k == 0 else vals[k]
            ax.text(x[k], yy + 0.03 * top, f"{vals[k]:.2f}", ha="center", va="bottom", fontsize=5.5)
    ax.axhline(1.0, color="0.5", lw=0.5, ls=":")
    ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=5.5)
    ax.set_xlim(-0.6, 3.6)
    return top


def group_figure(plt, recs, out_base):
    from matplotlib.lines import Line2D
    from matplotlib.transforms import blended_transform_factory as btf
    m0 = recs[0]["meta"]; ncol = len(recs)
    sel = pick_syn(recs[0])
    fig = plt.figure(figsize=(max(4.2, 1.0 + 2.3 * ncol), 9.6))
    gs = fig.add_gridspec(7, ncol, height_ratios=[1.0, 0.75, 0.75, 0.75, 0.75, 1.0, 1.05], hspace=0.42, wspace=0.28,
                          left=0.13 if ncol == 1 else 0.09, right=0.98, top=0.85 if ncol > 1 else 0.82, bottom=0.05)
    axes = [[None] * ncol for _ in range(7)]
    top = {r: 0.0 for r in range(7)}
    for j, d in enumerate(recs):
        m = d["meta"]; tr = m["t_rep0"]
        for r in range(7):
            sx = axes[0][j] if (0 < r < 5) else None
            axes[r][j] = fig.add_subplot(gs[r, j], sharex=sx)
        A = [axes[r][j] for r in range(7)]
        tf = np.asarray(d["live_tf"]) - tr
        has_off = "off_t_f" in d
        for ax in A[:5]:
            for ta in m["rep_arrivals"]:
                ax.axvline(ta - tr, color="0.65", lw=0.5, ls=":", zorder=0)
        A[0].plot(tf, d["live_v_soma_f"], color="0.15", lw=0.7)
        for c, (gid, tag, q) in enumerate(sel):
            i = row_of(d, gid)
            if i is None:
                continue
            col = SYN_COL[c]; tp = float(d["tp"][i])
            A[1].plot(tf, 1e3 * d["live_cai_CR_f"][i], color=col, lw=0.8)
            top[1] = max(top[1], float(np.nanmax(1e3 * d["live_cai_CR_f"][i])))
            if tp > 0:
                y = d["live_effcai_GB_f"][i] / tp; A[2].plot(tf, y, color=col, ls=LIVE_LS)
                top[2] = max(top[2], float(np.nanmax(y)))
                if has_off:
                    A[2].plot(np.asarray(d["off_t_in"]) - tr, d["off_E_f"][i] / tp, color=col, ls=OFF_LS)
            A[3].plot(tf, d["live_Vg_GB_f"][i], color=col, ls=LIVE_LS)
            A[4].plot(tf, d["live_W_GB_f"][i], color=col, ls=LIVE_LS)
            top[3] = max(top[3], float(np.nanmax(d["live_Vg_GB_f"][i]))); top[4] = max(top[4], float(np.nanmax(d["live_W_GB_f"][i])))
            if has_off:
                to = np.asarray(d["off_t_f"]) - tr
                A[3].plot(to, d["off_V_f"][i], color=col, ls=OFF_LS); A[4].plot(to, d["off_W_f"][i], color=col, ls=OFF_LS)
            tw = btf(A[4].transData, A[4].transAxes)
            el = live_events(np.asarray(d["live_tf"]), d["live_necb_GB_f"][i]) - tr
            A[4].plot(el, np.full(len(el), 0.97 - 0.06 * c), "v", ms=4, color=col, transform=tw, clip_on=False)
            eo = off_events(d, i) - tr; eo = eo[(eo >= tf[0]) & (eo <= tf[-1])] if len(tf) else eo
            A[4].plot(eo, np.full(len(eo), 0.97 - 0.06 * c), "v", ms=6, mfc="none", mec=col, mew=0.7, transform=tw,
                      clip_on=False)
            tcl = np.asarray(d["live_tc"]) / 1e3
            A[5].plot(tcl, d["live_rho_GB_c"][i], color=col, ls=LIVE_LS)
            A[5].plot(np.asarray(d["tc_off"]) / 1e3, d["off_rho_c"][i], color=col, ls=OFF_LS)
            t5 = btf(A[5].transData, A[5].transAxes)
            ev = off_events(d, i) / 1e3
            A[5].plot(ev, np.full(len(ev), 1.0 - 0.05 * c), "|", ms=4, color=col, mew=0.7, transform=t5, clip_on=False)
        A[2].axhline(1.0, color="k", lw=0.6, ls=(0, (1, 1.5)))
        A[2].text(1.0, 1.0, "$\\theta_p$ ", transform=btf(A[2].transAxes, A[2].transData), ha="right", va="bottom",
                  fontsize=6)
        A[3].axhline(m["theta_V"], color="k", lw=0.6, ls=(0, (1, 1.5)))
        A[3].text(1.0, m["theta_V"], f"$\\theta_V$ = {m['theta_V']:.2f} ", transform=btf(A[3].transAxes, A[3].transData),
                  ha="right", va="bottom", fontsize=6)
        A[4].axhline(m["theta_eCB"], color="k", lw=0.6, ls=(0, (1, 1.5)))
        A[4].text(1.0, m["theta_eCB"], f"$\\theta_{{eCB}}$ = {m['theta_eCB']:.1f} ",
                  transform=btf(A[4].transAxes, A[4].transData), ha="right", va="bottom", fontsize=6)
        A[5].axhline(m["rho_star"], color="0.6", lw=0.5, ls=":")
        A[5].axvspan(m["w0"] / 1e3, m["w1"] / 1e3, color="0.88", lw=0, zorder=0)
        A[5].set_ylim(-0.04, 1.04); A[5].set_xlim(0, m["tstop"] / 1e3)
        A[5].set_xlabel("induction time (s)")
        A[4].set_xlabel("time from first pre arrival (ms)")
        if len(tf):
            A[0].set_xlim(tf[0], tf[-1])
        for ax in A[:4]:
            plt.setp(ax.get_xticklabels(), visible=False)
        top[6] = max(top[6], bars(A[6], m, j == 0))
        ttl = m["label"] + ("\n(offline: NMDAR-block lane, $\\rho$ frozen)" if m.get("off_frozen") else "")
        A[0].set_title(ttl, fontsize=7, pad=4)
    top[3] = max(top[3] * 1.08, m0["theta_V"] * 1.35); top[4] = max(top[4] * 1.08, m0["theta_eCB"] * 1.35)
    top[2] = max(top[2] * 1.08, 1.35); top[1] = top[1] * 1.08; top[6] = top[6] * 1.25
    for j in range(ncol):
        for r in (1, 2, 3, 4, 6):
            axes[r][j].set_ylim(-0.03 * top[r] if r != 6 else 0.0, top[r])
        if j:
            for r in range(7):
                plt.setp(axes[r][j].get_yticklabels(), visible=False)
    if ncol > 1:
        lo = min(axes[0][j].get_ylim()[0] for j in range(ncol)); hi = max(axes[0][j].get_ylim()[1] for j in range(ncol))
        for j in range(ncol):
            axes[0][j].set_ylim(lo, hi)
    for r in range(7):
        axes[r][0].set_ylabel(YLAB[r]); letter(axes[r][0], "abcdefg"[r])
    h = [Line2D([], [], color=SYN_COL[c], lw=1.2, label=f"syn {gid}: {tag} (peak {q:.2f} $\\theta_p$)")
         for c, (gid, tag, q) in enumerate(sel)]
    h += [Line2D([], [], color="0.2", ls=LIVE_LS, label="BCL live (GluSynapseV5)"),
          Line2D([], [], color="0.2", ls=OFF_LS, label="offline fit kernel"),
          Line2D([], [], color="0.65", ls=":", label="own pre arrival"),
          Line2D([], [], color="0.2", marker="v", ls="none", ms=4, label="eCB step (BCL)"),
          Line2D([], [], color="0.2", marker="v", ls="none", ms=6, mfc="none", label="eCB step (offline)")]
    fig.legend(handles=h, loc="upper center", ncol=3 if ncol > 2 else 2, frameon=False, bbox_to_anchor=(0.5, 0.955),
               fontsize=6, handlelength=2.2, columnspacing=1.2)
    fig.suptitle(f"{m0['group_title']}   (pair {m0['pair']}, emodel delta-split1, v5c fit S1_V5c_s6_39)",
                 fontsize=8, y=0.99)
    fig.savefig(out_base + ".png"); fig.savefig(out_base + ".pdf")
    plt.close(fig)
    print("wrote", out_base + ".png/.pdf", flush=True)


def overview(plt, recs, out_base):
    from matplotlib.lines import Line2D
    n = len(recs); nc = 3; nr = int(np.ceil((n + 1) / nc))
    fig = plt.figure(figsize=(7.2, 1.45 * nr + 2.6))
    gs = fig.add_gridspec(nr + 1, nc, height_ratios=[1.0] * nr + [1.6], hspace=0.75, wspace=0.25, left=0.08,
                          right=0.98, top=0.95, bottom=0.12)
    selg = {}
    for d in recs:
        selg.setdefault(d["meta"]["group"], pick_syn(d))
    for k, d in enumerate(recs):
        m = d["meta"]; ax = fig.add_subplot(gs[k // nc, k % nc])
        for c, (gid, tag, q) in enumerate(selg[m["group"]]):
            i = row_of(d, gid)
            if i is None:
                continue
            ax.plot(np.asarray(d["live_tc"]) / 1e3, d["live_rho_GB_c"][i], color=SYN_COL[c], ls=LIVE_LS, lw=0.8)
            ax.plot(np.asarray(d["tc_off"]) / 1e3, d["off_rho_c"][i], color=SYN_COL[c], ls=OFF_LS, lw=0.8)
        ax.axhline(m["rho_star"], color="0.6", lw=0.5, ls=":")
        ax.set_ylim(-0.04, 1.04); ax.set_xlim(0, m["tstop"] / 1e3)
        ax.set_title(f"{GROUP_SHORT.get(m['group'], m['group'])}: {m['label']}", fontsize=6, pad=2)
        if k % nc == 0:
            ax.set_ylabel("$\\rho$")
        if k // nc == nr - 1 or k + nc >= n:
            ax.set_xlabel("induction time (s)", fontsize=6)
    axl = fig.add_subplot(gs[(n) // nc, n % nc]); axl.axis("off")
    axl.legend(handles=[Line2D([], [], color="0.2", ls=LIVE_LS, label="BCL live"),
                        Line2D([], [], color="0.2", ls=OFF_LS, label="offline fit kernel"),
                        Line2D([], [], color=SYN_COL[0], label="synapses near $\\theta_p$"),
                        Line2D([], [], color=SYN_COL[2], label="synapse far from $\\theta_p$")],
               loc="center", frameon=False, fontsize=6)
    ax = fig.add_subplot(gs[nr, :]); x = np.arange(n)
    M = [d["meta"] for d in recs]
    ax.errorbar(x, [m["exp"] for m in M], yerr=[m["sem"] for m in M], fmt="_", color="k", ms=10, capsize=2, lw=0.9,
                label="data $\\pm$ SEM", zorder=3)
    ax.plot(x - 0.18, [m["pred"] for m in M], "o", mfc="white", mec="0.2", ms=4.5, mew=0.8, label="offline fit (all pairs)")
    ax.plot(x + 0.0, [m["ratio_off"] for m in M], "s", mfc="none", mec="#0072B2", ms=4, mew=0.8,
            label="offline kernel (this pair)")
    ax.plot(x + 0.18, [m.get("ratio_live", np.nan) for m in M], "D", color="#D55E00", ms=3.5,
            label="BCL live (this pair)")
    ax.axhline(1.0, color="0.6", lw=0.5, ls=":")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{GROUP_SHORT.get(m['group'], m['group'])}\n{m['label']}" for m in M], rotation=40, ha="right",
                       fontsize=5.5)
    ax.set_ylabel("weight ratio"); ax.set_xlim(-0.6, n - 0.4)
    ax.legend(ncol=4, frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.22), fontsize=6)
    fig.savefig(out_base + ".png"); fig.savefig(out_base + ".pdf")
    plt.close(fig)
    print("wrote", out_base + ".png/.pdf", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--res", default=RES)
    a = ap.parse_args()
    plt = style()
    recs = [d for d in load(a.out) if d["meta"].get("ok", True)]
    assert recs, f"no records in {a.out}"
    print(f"{'record':<62} {'data':>11} {'pred':>6} {'off':>6} {'BCL':>6} {'|drho|':>7} {'|dd|':>6} {'eCB l/o':>8} "
          f"{'t_off':>6}")
    for d in recs:
        m = d["meta"]
        print(f"{m['group'] + ' ' + m['pair'] + ' ' + m['proto'] + ' ' + m['cond']:<62} {m['exp']:5.2f}+/-{m['sem']:.2f} "
              f"{m['pred']:6.3f} {m['ratio_off']:6.3f} {m.get('ratio_live', np.nan):6.3f} "
              f"{m.get('rho_end_maxdiff_live_off', np.nan):7.3f} {m.get('d_end_maxdiff_live_off', np.nan):6.3f} "
              f"{m.get('n_live_ecb', -1):>3}/{m.get('n_off_ecb', -1):<4} {m.get('t_peak_offset_ms', np.nan):6.2f}")
    os.makedirs(a.res, exist_ok=True)
    for g in dict.fromkeys(d["meta"]["group"] for d in recs):
        group_figure(plt, [d for d in recs if d["meta"]["group"] == g], os.path.join(a.res, f"traces_{TAG}_{g}"))
    overview(plt, recs, os.path.join(a.res, f"traces_{TAG}_overview"))


if __name__ == "__main__":
    main()
