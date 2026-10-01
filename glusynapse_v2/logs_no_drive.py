"""What separates burst pairings (Nevian 3AP: LTP) from single pairings (Nevian 1AP, Sjostrom 0.1 Hz: none)?
Per presynaptic arrival: VDCC charge and peak in windows after the arrival, and effcai peak. One pair, CPU."""
import sys, numpy as np
sys.path.insert(0, "glusynapse_v2")
from batch_v2 import BatchV2, _grid
P = ["nevian_3ap_50hz_dt+10ms", "nevian_3ap_100hz_dt+10ms", "nevian_2ap_50hz_dt+10ms", "nevian_3ap_50hz_dt-10ms",
     "nevian_1ap_dt+10ms", "nevian_3ap_20hz_dt+10ms", "sjostrom_0.1hz_dt+10ms", "sjostrom_10hz_dt+10ms",
     "sjostrom_50hz_dt+10ms", "sjostrom_50hz_dt-10ms", "nevian_3ap_50hz_dt-50ms", "10Hz_10ms"]
pair = sys.argv[1] if len(sys.argv) > 1 else "180351-198084"
B = BatchV2(["glusynapse_v2/extracted/ebner_delta-prefire", "glusynapse_v2/extracted/markram_delta-prefire-tr"],
            protocols=P, pairs={pair}, fast=False, signals=("vdcc",))
print(f"{'protocol':28s} {'Q10':>9s} {'Q30':>9s} {'Q60':>9s} {'pk30':>9s} {'dEff150':>7s}  (medians over synapses x arrivals; Q in nA*ms)")
for r in sorted(B.recs, key=lambda r: P.index(r["proto"])):
    t = r["t"]
    v = r["vdcc"]; arr = r["arr"]; e = r["effcai"]
    h = np.diff(t, append=t[-1])
    Q = {w: [] for w in (10, 30, 60)}; pk = []; ep = []
    for i in range(v.shape[0]):
        for a in np.unique(arr[i])[:20]:
            for w in Q:
                k = np.searchsorted(t, t[a] + w); Q[w].append(np.sum(v[i, a:k] * h[a:k]))
            k = np.searchsorted(t, t[a] + 30); pk.append(v[i, a:k].max())
            k = np.searchsorted(t, t[a] + 150); ep.append(e[i, a:k].max() - e[i, max(a - 1, 0)])
    print(f"{r['proto']:28s} " + " ".join(f"{np.median(Q[w]):9.2e}" for w in Q) + f" {np.median(pk):9.2e} {np.median(ep):7.3f}")
