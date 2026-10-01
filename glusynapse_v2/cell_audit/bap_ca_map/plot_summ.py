"""Bin, tabulate and plot bAP Ca map results (out_*.json)."""
import glob, json, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

HERE = "/project/rrg-emuller/dhuruva/plastyfire/glusynapse_v2/cell_audit/bap_ca_map"
BINS = {"basal": [0, 40, 80, 120, 160, 250, 400], "apical": [0, 100, 250, 450, 700, 1000]}
PROTOS = ["1ap", "3ap_50hz", "5ap_20hz"]
LS = {"1ap": "-", "3ap_50hz": "--", "5ap_20hz": ":"}
COL = {"delta": "C0", "antic_delta": "C3"}
rows = []            # (group, cell, variant, proto, kind, dist, vpk, shaft_pk, spine_pk, vdcc)
for f in sorted(glob.glob(f"{HERE}/out_*.json")):
    d = json.load(open(f)); grp = "L23" if d["l23"] else "L5"
    NP = {"1ap": 1, "3ap_50hz": 3, "5ap_20hz": 5}
    for r in d["results"]:
        if len(r["soma_spikes"]) != NP[r["proto"]]:
            print("EXCLUDED (soma fires %d of %d):" % (len(r["soma_spikes"]), NP[r["proto"]]), f, r["proto"]); continue
        for s in r["sites"]:
            rows.append((grp, d["cell"], d["variant"], r["proto"], s["kind"], s["dist"], s["vpk"], s["shaft_cai_pk_nospine"] * 1e3,
                         s["cacr_pk"] * 1e3, s["vdcc_int"], s["cai_pk"] * 1e3, d["cell"] + s["sec"].split("]")[-2][-1] + s["sec"].split(".")[-1] + str(s["x"])))
R = np.array(rows, dtype=object)
def sel(grp, var, proto, kind):
    m = (R[:, 0] == grp) & (R[:, 2] == var) & (R[:, 3] == proto) & (R[:, 4] == kind)
    return R[m]
def binned(x, y, edges):
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (x >= lo) & (x < hi)
        out.append(np.median(y[m]) if m.sum() else np.nan)
    return np.array(out)
lines = []
def pr(a, b):
    ka = {k: v for k, v in zip(a[:, 11], a[:, 8].astype(float))}
    r = [v / ka[k] for k, v in zip(b[:, 11], b[:, 8].astype(float)) if k in ka and ka[k] > 0]
    return np.median(r) if r else np.nan
for grp in ("L5", "L23"):
    for kind in ("basal", "apical"):
        E = BINS[kind]; ctr = [(a + b) / 2 for a, b in zip(E[:-1], E[1:])]
        if not len(sel(grp, "delta", "1ap", kind)): continue
        fig, ax = plt.subplots(2, 4, figsize=(16, 7))
        for var in ("delta", "antic_delta"):
            for p in PROTOS:
                S = sel(grp, var, p, kind)
                if not len(S): continue
                x = S[:, 5].astype(float)
                for j, (col, lab) in enumerate([(6, "local bAP amplitude (mV)"), (7, "shaft cai peak (uM, no spines)"), (8, "spine cai_CR peak (uM)"), (9, "spine VDCC charge (a.u.)")]):
                    y = S[:, col].astype(float)
                    ax[0, j].plot(x, y, ".", c=COL[var], ms=2, alpha=0.25)
                    b = binned(x, y, E); ax[0, j].plot(ctr, b, LS[p], c=COL[var], marker="o", ms=3, label=f"{var} {p}")
                    ax[0, j].set_ylabel(lab); ax[0, j].set_xlabel("path distance (um)")
                    if j > 0: ax[0, j].set_yscale("log")
                    if col in (7, 8, 6):
                        ax[1, {6: 0, 7: 1, 8: 2}[col]].plot(ctr, b / b[0] if b[0] > 0 else b, LS[p], c=COL[var], marker="o", ms=3)
                # spine/shaft ratio and spine vs V
                ax[1, 3].plot(ctr, binned(x, S[:, 8].astype(float) / np.maximum(S[:, 7].astype(float), 1e-9), E), LS[p], c=COL[var], marker="o", ms=3)
        ax[0, 0].axhspan(0, 30, xmin=0, xmax=0, color="k")
        if kind == "basal":
            ax[0, 0].plot([200, 400], [30, 30], "k_", ms=20); ax[0, 0].annotate("Nevian 2007: <=30 mV\nat >=200 um [text bound]", (200, 32), fontsize=7)
            ax[0, 2].errorbar([30], [1.4], yerr=[0.6], fmt="ks", label="Chindemi/Sabatini bAP <60 um (1.4+-0.6)")
        ax[0, 0].legend(fontsize=6); ax[0, 2].legend(fontsize=6)
        for j, t in enumerate(["bAP amp / proximal", "shaft Ca / proximal", "spine Ca / proximal", "spine Ca / shaft Ca (uM/uM)"]):
            ax[1, j].set_title(t, fontsize=9); ax[1, j].set_xlabel("path distance (um)")
            if j < 3: ax[1, j].set_yscale("log")
        ax[1, 3].set_yscale("log")
        fig.suptitle(f"{grp} {kind}: medians over sites/cells per bin (dots = single sites). Nevian 2007 bAP-Ca curves not digitised (full text unavailable)", fontsize=9)
        fig.tight_layout(); fig.savefig(f"{HERE}/figs/bap_ca_{grp}_{kind}.png", dpi=110); plt.close(fig)
        # tables
        lines.append(f"\n### {grp} {kind}: median over sites (both cells), spine/shaft in uM, dV in mV")
        lines.append("| bin um | var | n | bAP mV | shaft 1AP | spine 1AP | spine 3AP50 | spine 5AP20 | 3AP50/1AP spine (per-site) | 5AP20/1AP spine (per-site) | spine/shaft 1AP |")
        lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
        for var in ("delta", "antic_delta"):
            S1, S3, S5 = (sel(grp, var, p, kind) for p in PROTOS)
            if not len(S1): continue
            for lo, hi in zip(E[:-1], E[1:]):
                m = lambda S: (S[:, 5].astype(float) >= lo) & (S[:, 5].astype(float) < hi)
                a, b3, b5 = S1[m(S1)], S3[m(S3)], S5[m(S5)]
                if not len(a): continue
                f = lambda A, c: np.median(A[:, c].astype(float)) if len(A) else np.nan
                g = lambda x: f"{x:.3g}"
                lines.append(f"| {lo}-{hi} | {var} | {len(a)} | {f(a,6):.1f} | {g(f(a,7))} | {g(f(a,8))} | {g(f(b3,8))} | {g(f(b5,8))} | "
                             f"{g(pr(a, b3))} | {g(pr(a, b5))} | {g(np.median(a[:,8].astype(float)/np.maximum(a[:,7].astype(float),1e-12)))} |")
        # does spine Ca follow shaft Ca or V?  Spearman over sites, 1AP, og
        for var in ("delta", "antic_delta"):
            S = sel(grp, var, "1ap", kind)
            if len(S) > 5:
                sp, sh, v, x = (S[:, c].astype(float) for c in (8, 7, 6, 5))
                r1 = spearmanr(sp, sh)[0]; r2 = spearmanr(sp, v)[0]; r3 = spearmanr(sh, v)[0]
                # log-slope of spine vs shaft where both >1e-4 uM
                ok = (sp > 1e-4) & (sh > 1e-4)
                sl = np.polyfit(np.log10(sh[ok]), np.log10(sp[ok]), 1)[0] if ok.sum() > 5 else np.nan
                lines.append(f"- {grp} {kind} {var} 1AP Spearman: spine~shaft {r1:.2f}, spine~V {r2:.2f}, shaft~V {r3:.2f}; log-log slope spine vs shaft {sl:.2f} (n={len(S)})")
open(f"{HERE}/tables.md", "w").write("\n".join(lines)); print("\n".join(lines))
