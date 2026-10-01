"""Spine Ca (cai_CR - min) recovered from effcai: u_k = (eff[k+1] - a_k eff[k]) / b_k (exact inverse of
extract.effcai_from_cai, = step-mean of cai_CR - min over [t_k, t_k+1)). Per protocol: per-arrival peak and
charge above thresholds in [arrival, arrival + 50 ms]. Which threshold separates LTP from non-LTP?"""
import sys, numpy as np
sys.path.insert(0, "glusynapse_v2")
from batch_v2 import BatchV2
TAU = 200.0
LTP = {"nevian_3ap_50hz_dt+10ms": 2.01, "nevian_3ap_100hz_dt+10ms": 2.29, "nevian_2ap_50hz_dt+10ms": 1.95,
       "nevian_3ap_50hz_dt-10ms": 1.42, "sjostrom_50hz_dt+10ms": 1.57, "sjostrom_50hz_dt-10ms": 1.70,
       "sjostrom_40hz_dt+10ms": 1.53, "sjostrom_40hz_dt-10ms": 1.51, "sjostrom_20hz_dt+10ms": 1.31}
NON = {"nevian_1ap_dt+10ms": 1.04, "nevian_3ap_20hz_dt+10ms": 1.09, "sjostrom_0.1hz_dt+10ms": 0.97,
       "sjostrom_10hz_dt+10ms": 1.16, "nevian_3ap_50hz_dt-30ms": 0.98, "nevian_3ap_50hz_dt-50ms": 0.68,
       "nevian_3ap_50hz_dt+50ms": 0.92, "nevian_1ap_dt-10ms": 0.80, "sjostrom_10hz_dt-10ms": 0.57,
       "10Hz_10ms": 1.20, "10Hz_-10ms": 0.79}
P = list(LTP) + list(NON)
args = [x for x in sys.argv[1:] if not x.startswith("--")]
pairs = set(args[0].split(",")) if args else {"180351-198084", "181015-184976"}
B = BatchV2(["glusynapse_v2/extracted/ebner_delta-prefire", "glusynapse_v2/extracted/markram_delta-prefire-tr"],
            protocols=P, pairs=pairs, fast=False, signals=())
TH = (0.0, 5e-4, 1e-3)
NORM = "--norm" in sys.argv
print(f"{'protocol':26s} {'target':>6s} {'pk50':>8s} " + " ".join(f"Q>{t:g}" for t in TH) + "   (median over syn x arrivals; Q = mM ms)")
for r in sorted(B.recs, key=lambda r: (r["proto"] not in LTP, P.index(r["proto"]))):
    t = r["t"]; e = r["effcai"].astype(np.float64); h = np.diff(t)
    a = np.exp(-h / TAU); b = TAU * (1 - a)
    u = (e[:, 1:] - a * e[:, :-1]) / b                      # (n_syn, n-1), step means of cai_CR - min
    pk, Q = [], {x: [] for x in TH}
    cn = r["c_pre"] + r["c_post"]                            # per-synapse single pre + single post effcai peak
    for i in range(u.shape[0]):
        for k in np.unique(r["arr"][i])[:20]:
            k1 = np.searchsorted(t, t[k] + 50.0)
            seg = u[i, k:k1]; hh = h[k:k1]
            pk.append(seg.max())
            for x in TH:
                Q[x].append(np.sum(np.maximum(seg - x, 0) * hh) / (cn[i] if NORM else 1.0))
    tgt = LTP.get(r["proto"], NON.get(r["proto"]))
    print(f"{r['pair'][:6]} {r['proto']:26s} {tgt:6.2f} {np.median(pk):8.2e} " + " ".join(f"{np.median(Q[x]):8.1e}" for x in TH), flush=True)
