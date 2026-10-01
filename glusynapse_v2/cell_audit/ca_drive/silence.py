"""Own-K Hill-onset drive (recommended, CA_DRIVE_DESIGN.md): (a) history built from the vetoed influx (bAP and own EPSP peaks
weighted by 1 - B_NMDA), (b) after long silence K_i has decayed to the floor K0: does the own EPSP (pre-only) then gate a
2nd pre 50 ms later, with theta_Tg from the fitted parameter set?  Login node.  -> out/silence.json"""
import json, os
import numpy as np
import common as C
import cands as K
from extra import hon
A = C.A
S, t, dist, pid = C.build()
out = {}
for lj in (0, 25):
    s = S[lj]; qb = s["ap1"]["q"].max(1)
    wE = C.veto_weight("nmda", t["epsp"], C.arrivals("epsp"))[None, :]
    qe_v = (s["epsp"]["q"] * wE).max(1)
    ref = np.maximum(qb, qe_v)
    r = K.evaluate(s, t, dist, pid, lambda p: hon(s[p]["q"], 0.32 * ref[:, None], 8.0), veto="nmda")
    out[f"ljp{lj}|vetoed history"] = {k: r[k] for k in ("f_syn", "bins", "conn_major", "conn_meanT", "te", "tauT")}
    out[f"ljp{lj}|frac own vetoed EPSP peak > bAP peak"] = round(float(np.mean(qe_v > qb)), 3)
    # fitted set with bAP history
    te, tt = 0.9, 20.0
    inc = {p: hon(s[p]["q"], 0.32 * qb[:, None], 8.0) * C.veto_weight("nmda", t[p], C.arrivals(p))[None, :] for p in C.PROTOS}
    L, N, E, TR = C.rows({p: C.rule(inc[p], te, tt) for p in C.PROTOS}, t)
    th = A.score(L, N, E)[2]
    for k0 in (1e-4, 1e-3, 1e-2):
        K0 = k0 * np.median(qb)
        e = C.rule(hon(s["epsp"]["q"], K0, 8.0) * wE, te, tt)[:, np.searchsorted(t["epsp"], 50.1)]
        e_nov = C.rule(hon(s["epsp"]["q"], K0, 8.0), te, tt)[:, np.searchsorted(t["epsp"], 50.1)]
        out[f"ljp{lj}|silence K0={k0}xmedian"] = dict(frac_E_le_thTg_nmda_veto=round(float(np.mean(e <= th)), 3),
                                                     frac_E_le_thTg_no_veto=round(float(np.mean(e_nov <= th)), 3),
                                                     frac_bAP_below_K0=round(float(np.mean(qb < 3 * K0)), 3))
for k, v in out.items():
    print(k, v)
json.dump(out, open(os.path.join(C.HERE, "out", "silence.json"), "w"), indent=0)
