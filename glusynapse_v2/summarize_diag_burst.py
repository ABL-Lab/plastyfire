"""Table of diag_burst results: per pair x variant x protocol, soma spikes, synaptic VDCC charge of the burst
relative to the single bAP (same variant), per-spike charge relative to spike 1, local bAP peak, kBK state.

    python glusynapse_v2/summarize_diag_burst.py glusynapse_v2/results/diag_burst_1*.json
"""
import json, sys
import numpy as np, pandas as pd

rows = []
for f in sys.argv[1:]:
    D = json.load(open(f))
    one = {r["variant"]: np.array([s["ica_total"] for s in r["syn"]]) for r in D["results"] if r["proto"] == "1ap"}
    for r in D["results"]:
        s = r["syn"]
        ic = np.array([[e["ica_int"] for e in x["per_spike"]] for x in s])
        vp = np.array([[e["vpk"] for e in x["per_spike"]] for x in s])
        tot = np.array([x["ica_total"] for x in s])
        row = dict(pair=D["pair"], variant=r["variant"], proto=r["proto"], n_soma=len(r["soma_spikes"]),
                   burst_over_1ap=float(np.median(tot / one[r["variant"]])) if r["variant"] in one else np.nan,
                   per_spike=np.round(np.median(ic / ic[:, :1], 0), 2).tolist(),
                   vpk=np.round(np.median(vp, 0), 1).tolist(),
                   q1ap=float(np.median(one.get(r["variant"], [np.nan]))))
        if "bk_p_start" in s[0]["per_spike"][0]:
            row["bk_p_start"] = np.round(np.median([[e["bk_p_start"] for e in x["per_spike"]] for x in s], 0), 3).tolist()
            row["bk_p_pk"] = np.round(np.median([[e["bk_p_pk"] for e in x["per_spike"]] for x in s], 0), 3).tolist()
        rows.append(row)
df = pd.DataFrame(rows)
pd.set_option("display.width", 250, "display.max_rows", 500, "display.max_colwidth", 40)
print(df.to_string(index=False))
# the headline: burst/1AP charge per variant (median over pairs) for the protocols that matter
print("\nburst / single-bAP synaptic VDCC charge, median over pairs (linear sum of 3 bAPs = 3.0)")
print(df.pivot_table(index="variant", columns="proto", values="burst_over_1ap", aggfunc="median").round(2).to_string())
print("\nsoma spikes, min over pairs")
print(df.pivot_table(index="variant", columns="proto", values="n_soma", aggfunc="min").to_string())
print("\nsingle-bAP charge relative to delta (median over pairs)")
q = df[df.proto == "1ap"].pivot_table(index="pair", columns="variant", values="q1ap")
print((q.div(q["delta"], axis=0)).median().round(2).to_string())

# pairing runs (diag_burst --protos with "@dt" and "epsp"): NMDA Ca, peak cai_CR, peak effcai per synapse
pr = []
for f in sys.argv[1:]:
    D = json.load(open(f))
    R = {(r["variant"], r["proto"]): r for r in D["results"]}
    if not any(p == "epsp" for _, p in R):
        continue
    for (v, p), r in R.items():
        if "@" not in p:
            continue
        base = p.split("@")[0]
        rel = np.array([x["nmda_int"] > 0 for x in R[(v, "epsp")]["syn"]])     # released on the EPSP-only run (same
        get = lambda q, k: np.array([x[k] for x in R[(v, q)]["syn"]])[rel]     # random streams in every run)
        row = dict(pair=D["pair"], variant=v, proto=p, n_soma=len(r["soma_spikes"]), n_rel=int(rel.sum()),
                   effcai_epsp=float(np.median(get("epsp", "effcai_pk"))), effcai=float(np.median(get(p, "effcai_pk"))))
        for k in ("nmda_int", "cacr_pk", "effcai_pk"):
            # supralinearity: pairing / (EPSP alone + bAPs alone), per synapse, median
            row[k + "_supra"] = float(np.median(get(p, k) / (get("epsp", k) + get(base, k))))
        row["effcai_vs_1ap"] = float(np.median(get(p, "effcai_pk") / get("1ap" + p[len(base):], "effcai_pk")))
        row["nmda_vs_epsp"] = float(np.median(get(p, "nmda_int") / get("epsp", "nmda_int")))
        row["vdcc_share"] = float(np.median(get(p, "vdcc_int") / (get(p, "vdcc_int") + get(p, "nmda_int"))))
        pr.append(row)
if pr:
    P = pd.DataFrame(pr)
    print("\n=== pairings (median over synapses; supra = pair / (epsp + bAPs); 1 = linear)")
    print(P.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    print("\neffcai peak, pairing / same-dt 1AP pairing (Nevian target ratio 3AP50+10 / 1AP+10 = 2.01/1.04), median over pairs")
    print(P.pivot_table(index="variant", columns="proto", values="effcai_vs_1ap", aggfunc="median").round(2).to_string())
    print("\neffcai supralinearity, median over pairs")
    print(P.pivot_table(index="variant", columns="proto", values="effcai_pk_supra", aggfunc="median").round(2).to_string())
    print("\nNMDA Ca charge, pairing / EPSP alone (bAP relief of the Mg block), median over pairs")
    print(P.pivot_table(index="variant", columns="proto", values="nmda_vs_epsp", aggfunc="median").round(2).to_string())
    print("\nVDCC share of synaptic Ca charge in the pairing, median over pairs")
    print(P.pivot_table(index="variant", columns="proto", values="vdcc_share", aggfunc="median").round(3).to_string())
