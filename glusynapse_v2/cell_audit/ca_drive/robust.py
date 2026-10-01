"""Robustness of the recommended own-K Hill-onset drive (CA_DRIVE_DESIGN.md): k and n scans, a history that includes the
own EPSP influx peaks (K_i = k x max(bAP, EPSP peak)), the hand veto and no veto.  Login node.  -> out/robust.json"""
import json, os
import numpy as np
import common as C
import cands as K
from extra import hon
A = C.A
S, t, dist, pid = C.build()
out = {}
for lj in (0, 25):
    s = S[lj]; qb = s["ap1"]["q"].max(1); qe = s["epsp"]["q"].max(1)
    for hist, ref in (("bAP", qb), ("max(bAP,EPSP)", np.maximum(qb, qe))):
        for n in (4.0, 8.0):
            for k in (0.1, 0.2, 0.32, 0.5, 0.7):
                for veto in (("nmda", "hand", "none") if (hist == "bAP" and n == 8.0 and k == 0.32) else ("nmda",)):
                    r = K.evaluate(s, t, dist, pid, lambda p, Kv=k * ref[:, None], n=n: hon(s[p]["q"], Kv, n), veto=veto)
                    key = f"ljp{lj}|{hist}|n{n:.0f}|k{k}|{veto}"
                    out[key] = {kk: r[kk] for kk in ("f_syn", "bins", "conn_major", "conn_meanT", "te", "tauT", "trains")}
                    print(key, out[key], flush=True)
json.dump(out, open(os.path.join(C.HERE, "out", "robust.json"), "w"), indent=0)
