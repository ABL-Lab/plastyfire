"""Chindemi Fig2/Supp A.1-style summary from diag_burst sp1 jsons: spine Ca (cai_CR peak - rest) per case, NMDA/VDCC charge share."""
import json, glob, sys
import numpy as np
R = "/project/rrg-emuller/dhuruva/plastyfire/glusynapse_v2/cell_audit/ljp25_validation/out/"
tag = sys.argv[1] if len(sys.argv) > 1 else "sp1"
rows = []
for f in sorted(glob.glob(R + f"diag_burst_*_{tag}.json")):
    d = json.load(open(f))
    for r in d["results"]:
        if r["proto"] not in ("epsp", "1ap", "1ap@+10", "1ap@-10"): continue
        for s in r["syn"]:
            rows.append(dict(pair=d["pair"], var=r["variant"], proto=r["proto"], sid=s["sid"], kind=s["kind"], dist=s["dist"],
                             ca=(s["cacr_pk"] - 70e-6) * 1e3, nmda=s["nmda_int"], vdcc=s["vdcc_int"], vpk=s.get("vpk", max(e["vpk"] for e in s["per_spike"])),
                             nsp=len(r["soma_spikes"])))
json.dump(rows, open(f"/project/rrg-emuller/dhuruva/plastyfire/glusynapse_v2/cell_audit/ljp25_validation/rows_{tag}.json", "w"))
def st(x): x = np.asarray(x); return f"{np.median(x):.3g} [{np.percentile(x,25):.2g}-{np.percentile(x,75):.2g}] mean {x.mean():.3g}±{x.std():.2g} n={len(x)}"
for var in sorted(set(r["var"] for r in rows)):
    print("=====", var)
    for proto in ("epsp", "1ap", "1ap@+10", "1ap@-10"):
        for lab, sel in (("all", lambda r: True), ("basal<60", lambda r: r["kind"] == "basal" and r["dist"] < 60),
                         ("basal60-150", lambda r: r["kind"] == "basal" and 60 <= r["dist"] < 150),
                         ("basal>=150", lambda r: r["kind"] == "basal" and r["dist"] >= 150), ("apical", lambda r: r["kind"] == "apical")):
            q = [r for r in rows if r["var"] == var and r["proto"] == proto and sel(r)]
            if not q: continue
            sh = np.sum([r["vdcc"] for r in q]) / np.sum([r["vdcc"] + r["nmda"] for r in q])
            shs = np.median([r["vdcc"] / (r["vdcc"] + r["nmda"]) for r in q])
            print(f"{proto:8s} {lab:11s} Ca uM {st([r['ca'] for r in q])} | vdcc share pooled {sh:.2f} med {shs:.2f} | vpk med {np.median([r['vpk'] for r in q]):.1f}")
