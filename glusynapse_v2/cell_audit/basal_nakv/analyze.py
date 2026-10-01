"""Discriminability of per-synapse effcai / VDCC between fitted-protocol contrasts, per basal variant.

    python glusynapse_v2/cell_audit/basal_nakv/analyze.py
Per synapse (basal only, spike-matched runs), paired over the same synapses:
  ratio = median over synapses of A/B;  AUC = P(A_syn > B_syn) (paired sign fraction, 0.5 = no separation)
"""
import glob, json, os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CONTR = [("sj0.1@+10", "sj0.1@-10"), ("sj10@+10", "sj10@-10"), ("sj20@+10", "sj20@-10"), ("sj50@+10", "sj50@-10"),
         ("sj50@-10", "sj0.1@-10"), ("sj50@-10", "sj20@-10"), ("sj50@+10", "sj0.1@+10"), ("s03_m25", "s03_m120"), ("s03_b120", "s03_m120")]
NPOST = {"sj0.1": 1, "sj10": 5, "sj20": 5, "sj50": 5, "s03_m25": 1, "s03_m120": 1, "s03_b120": 5}
SIG = ("effcai_pk", "effcai_int", "vdcc_int", "vpk")


def load():
    rows, fi, spk = [], [], []
    for f in sorted(glob.glob(os.path.join(HERE, "out", "*.json"))):
        J = json.load(open(f))
        for r in J["results"]:
            ns = len(r["soma_spikes"])
            if r["proto"] == "fi":
                fi.append(dict(pair=J["pair"], variant=r["variant"], n=ns)); continue
            spk.append(dict(pair=J["pair"], variant=r["variant"], proto=r["proto"], n=ns, want=NPOST[r["proto"].split("@")[0]]))
            for s in r["syn"]:
                if s["kind"] != "basal":
                    continue
                rows.append(dict(pair=J["pair"], variant=r["variant"], proto=r["proto"], sid=s["sid"], dist=s["dist"],
                                 ok=ns == NPOST[r["proto"].split("@")[0]], nmda_int=s["nmda_int"], bap=s["per_spike"][0]["vpk"] - s["v_rest"], **{k: s[k] for k in SIG}))
    return pd.DataFrame(rows), pd.DataFrame(fi), pd.DataFrame(spk)


def contrasts(df):
    out = []
    for v, g in df.groupby("variant"):
        w = {k: g.pivot_table(index=["pair", "sid"], columns="proto", values=k) for k in SIG}
        okp = g.pivot_table(index=["pair", "sid"], columns="proto", values="ok")
        for a, b in CONTR:
            if a not in okp or b not in okp:
                continue
            m = (okp[a] == 1) & (okp[b] == 1)
            e = dict(variant=v, contrast=f"{a} vs {b}", n=int(m.sum()))
            for k in ("effcai_pk", "effcai_int", "vdcc_int"):
                A, B = w[k][a][m], w[k][b][m]
                e[f"{k}_ratio"] = float(np.median(A / B)); e[f"{k}_auc"] = float(np.mean(A > B))
            out.append(e)
    return pd.DataFrame(out)


if __name__ == "__main__":
    df, fi, spk = load()
    print("synapses:", df.groupby("variant").sid.count().to_dict())
    c = contrasts(df)
    c.to_csv(os.path.join(HERE, "contrasts.csv"), index=False)
    pd.set_option("display.width", 250)
    print(c.pivot(index="variant", columns="contrast", values="effcai_pk_ratio").round(2).to_string())
    print(c.pivot(index="variant", columns="contrast", values="effcai_pk_auc").round(2).to_string())
    print(c.pivot(index="variant", columns="contrast", values="n").to_string())
    print(c.pivot(index="variant", columns="contrast", values="vdcc_int_ratio").round(2).to_string())
    og = fi[fi.variant == "n1_k1_a1"].set_index("pair").n
    fi["rel"] = fi.apply(lambda r: r.n / og.get(r.pair, np.nan), axis=1)
    print("fI (0.5 nA, 600 ms) spikes:", fi.groupby("variant").n.median().to_dict(), "rel og:", fi.groupby("variant").rel.median().round(2).to_dict())
    spk["match"] = spk.n == spk.want
    print("post AP count matched:", spk.groupby("variant").match.mean().round(2).to_dict())
    print(spk[~spk.match].groupby(["variant", "proto"]).n.agg(list).to_string())
    # bAP at the synapse (1-AP protocol s03_m120 pre after; vpk dominated by bAP) vs distance
    b = df[(df.proto == "s03_m120") & df.ok].copy(); b["bin"] = pd.cut(b.dist, [0, 50, 100, 150, 400])
    print(b.pivot_table(index="variant", columns="bin", values="bap", aggfunc="median", observed=False).round(1).to_string())
    print("NMDA share of synaptic Ca (sj20@+10):", (df[df.proto == "sj20@+10"].groupby("variant").apply(
        lambda g: np.median(g.nmda_int / (g.nmda_int + g.vdcc_int)))).round(3).to_dict())
