"""Model spine/shaft Ca (delta vs ljp25_g031) against the experimental data (CA_VS_DATA.md).

Inputs: out/bapmap_*.json (bap_ca_map_glob.py, job 22092211), ../ljp25_validation/rows_ljp25.json (diag_burst, 68 L5->L5 synapses).
Output: out/summary.txt (all numbers used in CA_VS_DATA.md) and figs/ca_vs_data.png (one panel per quantity).
"""
import glob, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = os.path.dirname(os.path.abspath(__file__))
VARS = ["delta", "ljp25_g031"]
COL = {"delta": "#1f77b4", "ljp25_g031": "#d62728"}
LAB = {"delta": "delta (ljp 0)", "ljp25_g031": "ljp25_g031"}
OUT = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s); OUT.append(s)

# ---------------- bAP map (distance, spine/shaft, tau, bursts) ----------------
maps = {}
for f in sorted(glob.glob(f"{D}/out/bapmap_*.json")):
    d = json.load(open(f)); g = "L23" if d["l23"] else "L5"
    maps.setdefault((g, d["tag"]), []).append(d)

NP = {"1ap": 1, "3ap_50hz": 3, "5ap_20hz": 5, "3ap_100hz": 3, "3ap_200hz": 2}   # 200 Hz: both L5 cells fire only 2 of 3 (Letzkus pulse); kept as 2AP@200Hz
def sites(g, tag, proto, kind):
    """per-site dicts for cells whose soma fired all pulses of this proto; key (cell, sec, x)."""
    out = {}
    for d in maps.get((g, tag), []):
        for r in d["results"]:
            if r["proto"] != proto: continue
            if len(r["soma_spikes"]) < NP[proto]:
                P(f"  excluded {g} {d['cell']} {tag} {proto}: {len(r['soma_spikes'])}/{NP[proto]} APs"); continue
            for s in r["sites"]:
                if s["kind"] == kind: out[(d["cell"], s["sec"].split("].")[-1], s["x"])] = s   # template hash differs per run
    return out

BINS = {"basal": [0, 40, 80, 120, 160, 250, 400], "apical": [0, 100, 250, 450, 700, 1000]}
def binned(vals, dist, edges):
    vals, dist = np.asarray(vals, float), np.asarray(dist, float)
    r = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (dist >= a) & (dist < b) & np.isfinite(vals)
        r.append((0.5 * (a + b), np.median(vals[m]) if m.any() else np.nan,
                  np.percentile(vals[m], 25) if m.any() else np.nan, np.percentile(vals[m], 75) if m.any() else np.nan, int(m.sum())))
    return np.array(r)

M = {}   # (g, tag, kind) -> dict of binned arrays
for g in ("L5", "L23"):
    for tag in VARS:
        for kind in ("basal", "apical"):
            one = sites(g, tag, "1ap", kind)
            if not one: continue
            k = list(one); dist = [one[i]["dist"] for i in k]
            e = BINS[kind] if g == "L5" else [0, 40, 80]
            res = dict(dist=dist,
                       spine=binned([one[i]["cacr_pk"] * 1e3 for i in k], dist, e),
                       shaft=binned([one[i]["shaft_cai_pk_nospine"] * 1e3 for i in k], dist, e),
                       ratio=binned([one[i]["cacr_pk"] / one[i]["shaft_cai_pk_nospine"] for i in k], dist, e),
                       vpk=binned([one[i]["vpk"] for i in k], dist, e),
                       tau_sp=binned([one[i]["tau_cacr"] for i in k], dist, e),
                       tau_sh=binned([one[i]["tau_cai_nospine"] for i in k], dist, e))
            for bp in ("3ap_50hz", "5ap_20hz", "3ap_100hz", "3ap_200hz"):
                b = sites(g, tag, bp, kind); kk = [i for i in k if i in b]
                if not kk: continue
                dd = [one[i]["dist"] for i in kk]
                res[bp + "_pk"] = binned([b[i]["cacr_pk"] / one[i]["cacr_pk"] for i in kk], dd, e)
                res[bp + "_int"] = binned([b[i]["vdcc_int"] / one[i]["vdcc_int"] for i in kk], dd, e)
            M[(g, tag, kind)] = res
            P(f"\n### {g} {kind} {tag} (median per bin [IQR], n)")
            for q in [x for x in res if x != "dist"]:
                P(f"{q:14s} " + " | ".join(f"{c:.0f}: {m:.3g} [{lo:.2g}-{hi:.2g}] n{n:.0f}" for c, m, lo, hi, n in res[q]))

# ---------------- synaptic protocols (EPSP, pairings) from rows_ljp25.json ----------------
rows = json.load(open(f"{D}/../ljp25_validation/rows_ljp25.json"))
R = {}
for r in rows:
    R.setdefault((r["var"], r["pair"], r["sid"]), {})[r["proto"]] = r
def syn_stats(var, sel, cond_release=True):
    """per-synapse peak Ca (uM) and total Ca charge (nmda+vdcc) for epsp / 1ap / +10 / -10; release in all EPSP-containing runs."""
    out = []
    for (v, pair, sid), d in R.items():
        if v != var or not all(p in d for p in ("epsp", "1ap", "1ap@+10", "1ap@-10")): continue
        if not sel(d["1ap"]): continue
        if cond_release and not all(d[p]["nmda"] > 0 for p in ("epsp", "1ap@+10", "1ap@-10")): continue
        pk = {p: d[p]["ca"] for p in d}; it = {p: d[p]["nmda"] + d[p]["vdcc"] for p in d}
        out.append((pk, it))
    return out
SEL = {"basal<60": lambda r: r["kind"] == "basal" and r["dist"] < 60,
       "basal50-150": lambda r: r["kind"] == "basal" and 50 <= r["dist"] < 150,
       "all": lambda r: True}
S = {}
P("\n### synaptic protocols (release in epsp, +10, -10 runs); ratios per synapse, median [IQR]")
for var in VARS:
    for lab, sel in SEL.items():
        q = syn_stats(var, sel)
        if not q: continue
        f = lambda xs: (np.median(xs), np.percentile(xs, 25), np.percentile(xs, 75))
        res = {}
        for mode, idx in (("pk", 0), ("int", 1)):
            X = [x[idx] for x in q]
            res[f"epsp_{mode}"] = f([x["epsp"] for x in X])
            res[f"ap_{mode}"] = f([x["1ap"] for x in X])
            res[f"epsp/ap_{mode}"] = f([x["epsp"] / x["1ap"] for x in X])
            res[f"nl+10_{mode}"] = f([x["1ap@+10"] / (x["epsp"] + x["1ap"]) for x in X])
            res[f"nl-10_{mode}"] = f([x["1ap@-10"] / (x["epsp"] + x["1ap"]) for x in X])
        S[(var, lab)] = res
        P(f"{var:11s} {lab:12s} n={len(q):2d} " + "; ".join(f"{k} {m:.3g} [{lo:.2g}-{hi:.2g}]" for k, (m, lo, hi) in res.items()))
    # unconditional EPSP given release only in the EPSP run (for the absolute EPSP Ca vs Sabatini)
    for lab, sel in SEL.items():
        e = [d["epsp"]["ca"] for (v, _, _), d in R.items() if v == var and "epsp" in d and sel(d["epsp"]) and d["epsp"]["nmda"] > 0]
        if e: P(f"{var:11s} {lab:12s} EPSP Ca given release (epsp run): mean {np.mean(e):.3g} +- {np.std(e):.2g} median {np.median(e):.3g} n={len(e)}")
    for lab, sel in SEL.items():
        a = [d["1ap"]["ca"] for (v, _, _), d in R.items() if v == var and "1ap" in d and sel(d["1ap"])]
        if a: P(f"{var:11s} {lab:12s} bAP spine Ca (1ap run): mean {np.mean(a):.3g} +- {np.std(a):.2g} median {np.median(a):.3g} n={len(a)}")

open(f"{D}/out/summary.txt", "w").write("\n".join(OUT) + "\n")

# ---------------- figure ----------------
fig, ax = plt.subplots(2, 4, figsize=(19, 9)); ax = ax.ravel()
DATA = dict(color="k", mfc="w", ms=8, capsize=3, lw=1.3, zorder=5)
off = {"delta": -0.12, "ljp25_g031": 0.12}

# A: bAP spine Ca amplitude (uM)
a = ax[0]
for var in VARS:
    x0 = off[var]
    pts = []
    s = [d["1ap"]["ca"] for (v, _, _), d in R.items() if v == var and "1ap" in d and SEL["basal<60"](d["1ap"])]
    pts.append((0 + x0, np.mean(s), np.std(s)))
    b = M[("L5", var, "basal")]["spine"]; i = 2   # 80-120 bin
    a.errorbar(1 + x0, b[i, 1], yerr=[[b[i, 1] - b[i, 2]], [b[i, 3] - b[i, 1]]], fmt="o", color=COL[var], capsize=3)
    ap_ = M[("L5", var, "apical")]["spine"]; j = 0
    a.errorbar(2 + x0, ap_[j, 1], yerr=[[ap_[j, 1] - ap_[j, 2]], [ap_[j, 3] - ap_[j, 1]]], fmt="o", color=COL[var], capsize=3)
    a.errorbar(pts[0][0], pts[0][1], yerr=pts[0][2], fmt="o", color=COL[var], capsize=3, label=LAB[var])
a.errorbar(-0.3, 1.7, yerr=0.6, fmt="s", **DATA, label="Sabatini 2002 (CA1, dye-free; via Chindemi)")
a.errorbar(0.3, 1.4, yerr=0.6, fmt="^", color="gray", capsize=3, label="Chindemi 2022 target (model)")
for x in (1.35, 2.35):
    a.errorbar(x, 1.05, yerr=[[1.05 - 0.59], [7.98 - 1.05]], fmt="D", **DATA)
a.plot([1.5, 2.5], [0.7, 0.7], "k--", lw=1)
a.text(1.5, 0.75, "Cornelisse dye-free model", fontsize=7)
a.errorbar([], [], fmt="D", **DATA, label="Cornelisse 2007 Tab1 (L5, ~100 um, 0-dye extrap.)")
a.set_yscale("log"); a.set_ylim(0.005, 12)
a.set_xticks([0, 1, 2]); a.set_xticklabels(["basal <60 um", "basal 80-120 um", "apical 0-100 um"], fontsize=8)
a.set_ylabel("spine d[Ca] per bAP (uM)"); a.set_title("A  bAP spine Ca amplitude"); a.legend(fontsize=6.5, loc="lower left")

# B: spine/shaft ratio vs distance (basal)
a = ax[1]
for var in VARS:
    for kind, mk in (("basal", "o-"), ("apical", "s:")):
        r = M[("L5", var, kind)]["ratio"]
        a.plot(r[:, 0], r[:, 1], mk, color=COL[var], label=f"{LAB[var]} {kind}")
a.axhspan(1.75, 2.7, color="k", alpha=0.12, label="Cornelisse: 0-dye 2.7 (Tab1), model 1.75")
a.axhline(1.2, color="k", ls="--", lw=1, label="Cornelisse dF/F 1st AP ~1.2 [digitised Fig1E]")
a.set_yscale("log"); a.set_xlabel("path distance (um)"); a.set_ylabel("spine / shaft d[Ca] (1 AP)")
a.set_title("B  spine / shaft ratio"); a.legend(fontsize=6.5)

# C: decay tau
a = ax[2]
for var in VARS:
    b = M[("L5", var, "basal")]
    for xi, key in ((0, "tau_sp"), (1, "tau_sh")):
        r = b[key][:2]      # 0-80 um bins, where the signal is large
        a.errorbar(xi + off[var], np.nanmedian(r[:, 1]), yerr=[[np.nanmedian(r[:, 1]) - np.nanmin(r[:, 2])], [np.nanmax(r[:, 3]) - np.nanmedian(r[:, 1])]],
                   fmt="o", color=COL[var], capsize=3, label=LAB[var] if xi == 0 else None)
a.plot([-0.35], [12], "s", **{k: v for k, v in DATA.items() if k != "capsize"}, label="Sabatini 2002 0-dye (via Chindemi)")
a.plot([-0.3, 0.7], [26.8, 81.7], "D", **{k: v for k, v in DATA.items() if k != "capsize"}, label="Cornelisse 0-dye (1-comp, Tab1 kappas)")
a.plot([-0.2], [15], "v", color="k", mfc="w", ms=8, label="Cornelisse 0-dye multi-comp model ~15 ms [digitised Fig7A]")
a.plot([-0.25, 0.75], [91.2, 200.9], "x", color="k", ms=8, label="Cornelisse raw, 100 uM OGB-1")
a.set_yscale("log"); a.set_xticks([0, 1]); a.set_xticklabels(["spine", "shaft"]); a.set_ylabel("1/e decay (ms)")
a.set_title("C  decay tau (1 AP)"); a.legend(fontsize=6.5)

# D: basal distance dependence (normalised to the 0-40 bin)
a = ax[3]
for var in VARS:
    r = M[("L5", var, "basal")]["spine"]; v = M[("L5", var, "basal")]["vpk"]
    a.plot(r[:, 0], r[:, 1] / r[0, 1], "o-", color=COL[var], label=f"{LAB[var]} spine Ca")
    a.plot(v[:, 0], v[:, 1] / v[0, 1], "--", color=COL[var], lw=1, label=f"{LAB[var]} bAP mV")
a.fill_between([0, 80], 0.5, 2, color="k", alpha=0.12, label="Koester&Sakmann 98: AP spine Ca\n'comparable' up to 80 um (text)")
a.axvspan(50, 150, ymin=0, ymax=0.04, color="g", alpha=0.5, label="NS06 spines 50-150 um respond to 1 AP")
a.set_yscale("log"); a.set_ylim(1e-4, 3); a.set_xlabel("basal path distance (um)"); a.set_ylabel("relative to 0-40 um")
a.set_title("D  distance dependence, L5 basal"); a.legend(fontsize=6.5, loc="lower left")

# E: apical distance (absolute)
a = ax[4]
for var in VARS:
    r = M[("L5", var, "apical")]["spine"]; s = M[("L5", var, "apical")]["shaft"]
    a.plot(r[:, 0], r[:, 1], "o-", color=COL[var], label=f"{LAB[var]} spine")
    a.plot(s[:, 0], s[:, 1], "s:", color=COL[var], label=f"{LAB[var]} shaft")
a.errorbar(100, 1.05, yerr=[[0.46], [6.93]], fmt="D", **DATA, label="Cornelisse spine (~100 um, 0-dye)")
a.errorbar(110, 0.383, yerr=[[0.153], [0.917]], fmt="v", **DATA, label="Cornelisse dendrite")
a.set_yscale("log"); a.set_xlabel("apical path distance (um)"); a.set_ylabel("d[Ca] per bAP (uM)")
a.set_title("E  apical (trunk+obliques)"); a.legend(fontsize=6.5)

# F: burst / single
a = ax[5]
labs = ["3AP50", "3AP100", "2AP@200*"]; keys = ["3ap_50hz", "3ap_100hz", "3ap_200hz"]
for var in VARS:
    b = M[("L5", var, "basal")]
    for xi, k in enumerate(keys):
        for mode, mk, dx in (("pk", "o", -0.08), ("int", "^", 0.08)):
            if k + "_" + mode not in b: continue
            r = b[k + "_" + mode]
            prox, dist_ = r[0, 1], np.nanmedian(r[2:4, 1])
            a.plot([xi + off[var] + dx] * 2, [prox, dist_], "-", color=COL[var], lw=0.8)
            a.plot(xi + off[var] + dx, prox, mk, color=COL[var], mfc=COL[var], label=f"{LAB[var]} {'peak' if mode=='pk' else 'Ca charge'} (filled <40 um, open 80-160)" if xi == 0 else None)
            a.plot(xi + off[var] + dx, dist_, mk, color=COL[var], mfc="w")
a.axhspan(1.22, 1.70, xmin=0, xmax=0.36, color="k", alpha=0.12, label="0-dye peak, 3AP50: 1.22-1.70 (tau 12-27 ms)")
a.plot([-0.3, 0.7], [2.3, 3.2], "D", **{k: v for k, v in DATA.items() if k != "capsize"}, label="NS06 dG/R (500 uM OGB-6F) [digitised Fig5E,F]")
a.plot(-0.35, 1.15, "v", color="k", mfc="w", ms=8, label="Cornelisse 0-dye model, 50 Hz spine ~1.15 [digitised Fig7A]")
a.axhline(3, color="k", ls=":", lw=1, label="linear sum of charge (3)")
a.annotate("Kampa&Stuart 06: distal basal\nsupralinear at 200 Hz (text);\nmodel fires 2 of 3 APs", (2, 2.6), fontsize=7, ha="center"); a.set_ylim(0.9, 6.5)
a.set_xticks(range(3)); a.set_xticklabels(labs); a.set_ylabel("burst / 1 AP (spine)"); a.set_title("F  burst / single AP")
a.legend(fontsize=5.5, loc="upper left")

# G: EPSP spine Ca (ratio to bAP)
a = ax[6]
for var in VARS:
    for xi, lab in enumerate(("basal<60", "basal50-150")):
        if (var, lab) not in S: continue
        for mode, mk, dx in (("pk", "o", -0.06), ("int", "^", 0.06)):
            m, lo, hi = S[(var, lab)][f"epsp/ap_{mode}"]
            a.errorbar(xi + off[var] + dx, m, yerr=[[m - lo], [hi - m]], fmt=mk, color=COL[var], capsize=3,
                       label=f"{LAB[var]} {'peak' if mode=='pk' else 'Ca charge'}" if xi == 0 else None)
a.plot(1.35, 0.045 / 0.033, "D", **{k: v for k, v in DATA.items() if k != "capsize"}, label="NS06 EPSP/1AP dG/R [digitised Fig5E]")
a.plot(-0.35, 0.7 / 1.7, "s", **{k: v for k, v in DATA.items() if k != "capsize"}, label="Sabatini 0.7/1.7 uM (0-dye peaks)")
a.axhspan(0.7, 1.4, xmin=0.5, xmax=1, color="k", alpha=0.08, label="Koester&Sakmann 98 'comparable' (<80 um)")
a.set_yscale("log"); a.set_xticks([0, 1]); a.set_xticklabels(["basal <60", "basal 50-150"]); a.set_ylabel("EPSP / 1 AP spine Ca")
a.set_title("G  EPSP spine Ca (given release)"); a.legend(fontsize=6.5)

# H: supralinearity
a = ax[7]
for var in VARS:
    for xi, lab in enumerate(("basal<60", "basal50-150")):
        if (var, lab) not in S: continue
        for sgn, dx0 in (("+10", 0), ("-10", 2)):
            for mode, mk, dx in (("pk", "o", -0.06), ("int", "^", 0.06)):
                m, lo, hi = S[(var, lab)][f"nl{sgn}_{mode}"]
                a.errorbar(dx0 + xi * 0.8 + off[var] + dx, m, yerr=[[m - lo], [hi - m]], fmt=mk, color=COL[var], capsize=2,
                           label=f"{LAB[var]} {'peak' if mode=='pk' else 'Ca charge'}" if (xi == 0 and sgn == "+10") else None)
a.errorbar([0.8 + 0.35], [1.73], yerr=0.28, fmt="D", **DATA, label="NS06 1AP [digitised Fig5E]; 3AP +10: 1.8+-0.1 (text)")
a.errorbar([2.8 + 0.35], [1.03], yerr=0.18, fmt="D", **DATA)
a.axhline(1, color="k", ls=":", lw=1)
a.set_xticks([0, 0.8, 2, 2.8]); a.set_xticklabels(["+10\n<60", "+10\n50-150", "-10\n<60", "-10\n50-150"], fontsize=8)
a.set_ylabel("pairing / (EPSP + 1 AP)"); a.set_title("H  pairing supralinearity (1 AP)"); a.legend(fontsize=6.5)

fig.suptitle("Spine Ca: model (L5 TTPC, og-delta emodel; delta = ljp_VDCC 0, ljp25_g031) vs data.  Open black = data; dye-loaded data compared with model Ca charge (triangles)", fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(f"{D}/figs/ca_vs_data.png", dpi=130)
print("saved")
