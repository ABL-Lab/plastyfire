"""Figure: cell_audit/figs/nmda_vdcc_calibration.png (sp1 + cal1 diag_burst, clamp_vdcc)."""
import json, glob
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.style.use("/project/rrg-emuller/dhuruva/plastyfire/onerule.mplstyle")
D = "/project/rrg-emuller/dhuruva/plastyfire/glusynapse_v2/cell_audit/"
R = "/project/rrg-emuller/dhuruva/plastyfire/glusynapse_v2/results/"
rows = []
for tag in ("sp1", "cal1"):
    for f in sorted(glob.glob(R + f"diag_burst_*_{tag}.json")):
        d = json.load(open(f))
        for r in d["results"]:
            if r["proto"] not in ("epsp", "1ap", "1ap@+10", "1ap@-10"): continue
            for s in r["syn"]:
                rows.append(dict(pair=d["pair"], var=r["variant"], proto=r["proto"], kind=s["kind"], dist=s["dist"],
                                 ca=(s["cacr_pk"] - 70e-6) * 1e3, nmda=s["nmda_int"], vdcc=s["vdcc_int"],
                                 vpk=s.get("vpk", max(e["vpk"] for e in s["per_spike"]))))
json.dump(rows, open(D + "nmda_vdcc/rows_all.json", "w"))
sel = lambda var, proto: [r for r in rows if r["var"] == var and r["proto"] == proto]
C = {"og": "#1f4e79", "s5": "#2a9d8f", "s10": "#e9a03b", "g2": "#8e5ea2", "s10g3": "#c0392b", "og_basal": "#7f7f7f"}
LAB = {"og": "og-delta, mod as calibrated", "s5": "ljp_VDCC +5 mV", "s10": "ljp_VDCC +10 mV", "g2": "gca_bar ×2",
       "s10g3": "s10g3 (delta-sv)", "og_basal": "og (passive) basal"}
fig, ax = plt.subplots(2, 2, figsize=(7.0, 5.2))
# A: bAP spine Ca vs distance
a = ax[0, 0]
for var in ("og", "s5", "g2", "og_basal"):
    q = sel(var, "1ap")
    if not q: continue
    for kind, mk in (("basal", "o"), ("apical", "^")):
        qq = [r for r in q if r["kind"] == kind]
        a.scatter([r["dist"] for r in qq], [max(r["ca"], 1e-3) for r in qq], s=7 if var == "og" else 4, marker=mk,
                  color=C[var], alpha=0.9 if var == "og" else 0.5, lw=0, label=LAB[var] if kind == "basal" else None)
a.axvspan(0, 60, color="0.9", zorder=-5)
a.errorbar([30], [1.4], [[0.6], [0.6]], fmt="s", color="k", ms=3, capsize=2, label="Chindemi 2022 in silico (basal <60 µm)")
a.errorbar([45], [1.7], [[0.6], [0.6]], fmt="D", mfc="w", color="k", ms=3, capsize=2, label="Sabatini 2002 in vitro")
a.set_yscale("log"); a.set_ylim(1e-3, 20); a.set_xlabel("Path distance from soma (µm)"); a.set_ylabel("Spine Ca per bAP (µM)")
a.set_title("A  Single bAP: spine Ca vs distance (VDCC only)", loc="left"); a.legend(fontsize=5, frameon=False, loc="lower left")
# B: local v peak vs distance
b = ax[0, 1]
for var in ("og", "og_basal"):
    q = sel(var, "1ap")
    for kind, mk in (("basal", "o"), ("apical", "^")):
        qq = [r for r in q if r["kind"] == kind]
        b.scatter([r["dist"] for r in qq], [r["vpk"] for r in qq], s=6, marker=mk, color=C[var], lw=0,
                  alpha=0.9 if var == "og" else 0.5, label=f"{LAB[var]} ({kind})")
b.axhline(-5.9, color="k", ls=":", lw=0.6); b.text(250, -4.5, "VDCC V½ (m) −5.9 mV", fontsize=5, ha="right")
b.axvspan(0, 60, color="0.9", zorder=-5)
b.set_xlabel("Path distance from soma (µm)"); b.set_ylabel("Peak local v during bAP (mV)")
b.set_title("B  bAP peak at the synapse", loc="left"); b.legend(fontsize=5, frameon=False)
# C: clamp transfer
c = ax[1, 0]
cl = json.load(open(D + "nmda_vdcc/clamp_vdcc.json"))
for var, ls in (("mod", "-"), ("g2", "--"), ("s10g3", ":")):
    for td, lw in ((1.0, 0.6), (2.0, 1.2)):
        q = sorted([r for r in cl if r["var"] == var and r["td"] == td and r["vol"] == 0.153], key=lambda r: r["vp"])
        c.plot([r["vp"] for r in q], [r["ca"] for r in q], ls=ls, lw=lw, color={"mod": C["og"], "g2": C["g2"], "s10g3": C["s10g3"]}[var],
               label=f"clamp {var}, half-width {q[0]['hw']:.1f} ms")
q = sel("og", "1ap")
c.scatter([r["vpk"] for r in q], [max(r["ca"], 1e-3) for r in q], s=5, color=C["og"], lw=0, alpha=0.7, label="cell synapses (og-delta)")
c.axhspan(0.8, 2.0, color="0.9", zorder=-5); c.text(-48, 1.0, "Chindemi 1.4±0.6 µM", fontsize=5)
c.set_yscale("log"); c.set_ylim(1e-3, 20); c.set_xlim(-50, 30)
c.set_xlabel("Peak bAP voltage at spine (mV)"); c.set_ylabel("Spine Ca per bAP (µM)")
c.set_title("C  Mod alone: clamped bAP → VDCC Ca (vol 0.153 µm³)", loc="left"); c.legend(fontsize=5, frameon=False, loc="lower right")
# D: VDCC share
d = ax[1, 1]
vars_ = [v for v in ("og", "og_basal", "s5", "g2", "s10", "s10g3") if sel(v, "1ap@+10")]
protos = ("epsp", "1ap@+10", "1ap@-10")
w = 0.8 / len(vars_)
for i, var in enumerate(vars_):
    ys = []
    for p in protos:
        q = sel(var, p); ys.append(np.sum([r["vdcc"] for r in q]) / np.sum([r["vdcc"] + r["nmda"] for r in q]))
    d.bar(np.arange(len(protos)) + (i - (len(vars_) - 1) / 2) * w, ys, w, color=C[var], label=LAB[var])
d.set_xticks(range(len(protos))); d.set_xticklabels(["EPSP only", "pre→post +10 ms", "post→pre −10 ms"])
d.set_ylabel("VDCC share of spine Ca charge"); d.set_ylim(0, 1)
d.set_title("D  VDCC share, all 68 synapses pooled", loc="left"); d.legend(fontsize=5, frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(D + "figs/nmda_vdcc_calibration.png", dpi=200)
print("ok")
