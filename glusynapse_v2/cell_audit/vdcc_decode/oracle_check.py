"""Oracle check: VDCC-only spine Ca integral normalised per synapse to its own bAP area (so each bAP gives S ~ 1),
fine theta_Te grid; plus within-bin spread of the per-bAP drive.  Login node.  -> printed, oracle_check.json"""
import json, os
import numpy as np
import replay as RP
import tolerance as TL
A = RP.A
A.MARGIN = 1.0
recs = A.load("og"); t = recs[0]["t"]
dist = np.concatenate([r["dist"] for r in recs])
gk = np.load(os.path.join(RP.HERE, "gK.npz")); g, K = gk["g"], gk["K"]
S = {}
for p in ("ap1", "burst", "epsp"):
    v = RP.stack(recs, p, "v")[0]
    S[p] = dict(v=v, cvd0=RP.spine(v, np.zeros_like(v), g, K)[1], cvd25=RP.spine(v, np.zeros_like(v), g, K, ljp=25.0)[1])
out = {"spread_IQR_fold": {}, "oracle": {}}
bins = {"<60": dist < 60, "60-150": (dist >= 60) & (dist < 150), ">150": dist >= 150}
for key in ("v", "cvd0", "cvd25"):
    x = S["ap1"][key]; area = np.maximum(x - x[:, :1], 0).sum(1)
    for b, m in bins.items():
        q = np.percentile(area[m], [25, 75]); out["spread_IQR_fold"][f"{key}|{b}"] = float(q[1] / q[0])
for key in ("cvd0", "cvd25"):
    x0 = S["ap1"][key]; rest = x0[:, :1]; area = (np.maximum(x0 - rest, 0).sum(1) * A.DT)[:, None]
    best = (0, None)
    for te in np.arange(0.5, 1.001, 0.01):
        for tt in (20.0, 40.0):
            for w in (0.0, 15.0):
                T = {}
                for p in S:
                    u = np.maximum(S[p][key] - rest, 0) / area
                    if w:
                        u = u * A.veto_mask(t[p], A.pre_times(p) + 0.1, w)[None, :]
                    T[p] = A.lp(np.maximum(A.lp(u, A.TAU_E1) - te, 0), tt)
                at = lambda p, xx: T[p][:, np.searchsorted(t[p], xx + 0.1)]
                L = np.stack([at(p, xx) for p, xx in A.LTD], 1); N = np.stack([at(p, xx) for p, xx in A.NOLTD], 1)
                f = A.score(L, N, at("epsp", 50.0), strict=True)[1]
                if f > best[0]:
                    best = (float(f), dict(te=round(float(te), 2), tauT=tt, w=w))
    out["oracle"][key] = best
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(RP.HERE, "oracle_check.json"), "w"), indent=1)
