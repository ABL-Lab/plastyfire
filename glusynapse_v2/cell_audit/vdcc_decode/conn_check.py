"""Connection-level scores (T summed or averaged over the contacts of each pair) for the event drives, and the
oracle-normalised Ca integral with and without the -10 ms row.  Login node.  -> conn_check.json"""
import json, os
import numpy as np
import replay as RP
import counterfactual as CF
A = RP.A
A.MARGIN = 1.0
recs = A.load("og"); t = recs[0]["t"]
dist = np.concatenate([r["dist"] for r in recs])
pid = np.concatenate([np.full(len(r["dist"]), i) for i, r in enumerate(recs)])
gk = np.load(os.path.join(RP.HERE, "gK.npz")); g, K = gk["g"], gk["K"]
S = {}
for p in CF.PROTOS:
    v = RP.stack(recs, p, "v")[0]
    S[p] = dict(v=v, nvdcc=-RP.stack(recs, p, "ica_vdcc")[0], ivd25=-RP.spine(v, np.zeros_like(v), g, K, ljp=25.0)[0])
    S[p]["vrel"] = v
    if p in ("ap1", "burst", "epsp"):
        S[p]["cvd25"] = RP.spine(v, np.zeros_like(v), g, K, ljp=25.0)[1]
med = lambda k: np.median((S["ap1"][k] - S["ap1"][k][:, :1]).max(1))


def d1(x, rest, p, arr, tt):   # D1 exactly as check_d1 (th 3 mV rel, no hyst as check_hp, w 15, te .85, tauT 40)
    return A.c_event(x, rest, dict(th=3.0, te=0.85, tauT=40.0, rel=1, hyst=0.0, w=15.0), arr, tt)


out = {}
cands = {"D1 V event + veto (t_drive 3)": ("vrel", d1, {}),
         "C1v VDCC current event + veto, og": ("nvdcc", CF.event, dict(th=med("nvdcc") * 0.003, te=0.85, tauT=40.0, w=15.0)),
         "C1v VDCC current event + veto, R -25 mV": ("ivd25", CF.event, dict(th=med("ivd25") * 0.003, te=0.85, tauT=40.0, w=15.0))}
for name, (key, fn, p) in cands.items():
    L, N, E = CF.rows(S, key, fn, p, t)
    th = A.score(L, N, E)[2]
    ok = (L.min(1) >= 2 * th) & (np.maximum(N.max(1), E) <= th)
    mean = lambda X: CF.conn(X, pid) / np.bincount(pid)[:, None] if X.ndim > 1 else CF.conn(X, pid) / np.bincount(pid)
    out[name] = dict(f_syn=float(ok.mean()), by_dist=[float(ok[(dist >= a) & (dist < b)].mean()) for a, b in ((0, 60), (60, 150), (150, 250), (250, 1e9))],
                     f_conn_sum=A.score(CF.conn(L, pid), CF.conn(N, pid), CF.conn(E, pid))[1],
                     f_conn_mean=A.score(mean(L), mean(N), mean(E))[1],
                     f_conn_all_contacts_pass=float(np.mean([ok[pid == i].all() for i in np.unique(pid)])),
                     f_conn_majority_pass=float(np.mean([ok[pid == i].mean() >= 0.5 for i in np.unique(pid)])))
    print(name, out[name], flush=True)
# oracle Ca integral, -10 row dropped
x0 = S["ap1"]["cvd25"]; rest = x0[:, :1]; area = (np.maximum(x0 - rest, 0).sum(1) * A.DT)[:, None]
best = {}
for te in np.arange(0.3, 1.001, 0.02):
    for tt in (20.0, 40.0, 80.0):
        T = {p: A.lp(np.maximum(A.lp(np.maximum(S[p]["cvd25"] - rest, 0) / area, A.TAU_E1) - te, 0), tt) for p in ("ap1", "burst", "epsp")}
        at = lambda p, xx: T[p][:, np.searchsorted(t[p], xx + 0.1)]
        L = np.stack([at(p, xx) for p, xx in A.LTD], 1); N = np.stack([at(p, xx) for p, xx in A.NOLTD], 1); E = at("epsp", 50.0)
        for lab, cols in (("all rows", slice(None)), ("without -10", slice(1, None))):
            f = A.score(L[:, cols], N, E)[1]
            if f > best.get(lab, (0,))[0]:
                best[lab] = (float(f), dict(te=round(float(te), 2), tauT=tt))
out["oracle cvd25 integral"] = best
print(best)
json.dump(out, open(os.path.join(RP.HERE, "conn_check.json"), "w"), indent=1)
