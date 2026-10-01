"""Synapse-level pairing test of delta-bd candidates vs og-delta (diag_burst.py --bd-cands runs, TAG _syn1).

Per pair x variant x pairing, over synapses that released on the EPSP-only run: median effcai peak gain of the burst
pairing over the same-dt 1AP pairing, the burst pairing's effcai relative to og on the same synapse, and the -10/+10
order. Only runs where the soma fires the full AP count (1 or 3) enter, for both the candidate and og.

    python glusynapse_v2/emodel_bd/compare_syn.py 'glusynapse_v2/results/diag_burst_*_syn1.json' [--out csv]
"""
import argparse, glob, json
import numpy as np, pandas as pd


def rows(f):
    D = json.load(open(f))
    R = {(r["variant"], r["proto"]): r for r in D["results"]}
    out = []
    og_rel = None
    for (v, p), r in R.items():
        if "@" not in p or not p.startswith("3ap"):
            continue
        base, _, dt = p.partition("@")
        one = R.get((v, "1ap@" + dt)); ep = R.get((v, "epsp")); ogr = R.get(("og", p))
        if one is None or ep is None or ogr is None:
            continue
        rel = np.array([x["nmda_int"] > 0 for x in R[("og", "epsp")]["syn"]])
        g = lambda rr, k: np.array([x[k] for x in rr["syn"]])[rel]
        out.append(dict(pair=D["pair"], variant=v, proto=p, n3=len(r["soma_spikes"]), n1=len(one["soma_spikes"]),
                        og_n3=len(ogr["soma_spikes"]), n_rel=int(rel.sum()),
                        gain=float(np.median(g(r, "effcai_pk") / g(one, "effcai_pk"))),
                        vs_og=float(np.median(g(r, "effcai_pk") / g(ogr, "effcai_pk"))),
                        one_vs_og=float(np.median(g(one, "effcai_pk") / g(R[("og", "1ap@" + dt)], "effcai_pk"))),
                        vdcc_share=float(np.median(g(r, "vdcc_int") / (g(r, "vdcc_int") + g(r, "nmda_int")))),
                        effcai=float(np.median(g(r, "effcai_pk")))))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("glob"); ap.add_argument("--out")
    a = ap.parse_args()
    df = pd.DataFrame([r for f in sorted(glob.glob(a.glob)) for r in rows(f)])
    df["full"] = (df.n3 >= 3) & (df.n1 >= 1) & (df.og_n3 >= 3)
    pd.set_option("display.width", 250, "display.max_rows", 500)
    m = df[df.full]
    print("pairs with full AP count (candidate and og), per pairing:")
    print(df.pivot_table(index="variant", columns="proto", values="full", aggfunc="sum").to_string())
    for k, t in (("gain", "effcai peak, 3AP pairing / same-dt 1AP pairing (Nevian 50 Hz +10 target 1.93)"),
                 ("vs_og", "3AP pairing effcai peak relative to og on the same synapse"),
                 ("one_vs_og", "1AP pairing effcai peak relative to og on the same synapse"),
                 ("vdcc_share", "VDCC share of synaptic Ca in the 3AP pairing")):
        print(f"\n{t}, median over spike-matched pairs")
        print(m.pivot_table(index="variant", columns="proto", values=k, aggfunc="median").round(3).to_string())
    # -10 vs +10: Nevian 3AP 50 Hz LTP at +10, LTD at -10; want effcai(+10) > effcai(-10)
    w = m.pivot_table(index=["pair", "variant"], columns="proto", values="effcai")
    for f in ("3ap_50hz", "3ap_100hz", "3ap_200hz"):
        if f"{f}@+10" in w and f"{f}@-10" in w:
            print(f"\n{f}: effcai(+10)/effcai(-10), median over pairs")
            print((w[f"{f}@+10"] / w[f"{f}@-10"]).groupby("variant").median().round(3).to_string())
    if a.out:
        df.to_csv(a.out, index=False)


if __name__ == "__main__":
    main()
