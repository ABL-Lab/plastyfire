"""Counterfactual Ca drives for the Sjostrom 2003 window (score as local_t/check_hp.py: MARGIN 1 = LTD rows >= 2 theta_Tg,
no-LTD rows and the own-EPSP read 50 ms later <= theta_Tg, one common parameter set).  Login node, numpy/scipy.
Signals are replayed on the recorded local V (replay.py; exact for the spine VDCC, open-loop for the shaft channels).
Drives (volume-averaged Ca, integrated, per ECB_TRIGGER_LIT.md Q5):
  lin : T = lp(pos(x - rest - th)/a, tauT)                         (linear above threshold)
  two : S = lp(pos(x - rest - th)/a, 100); T = lp(pos(S - te), tauT)  (t_drive 2 style supralinearity)
  ev  : crossing count (C1), for reference
  +v  : own-pre veto of the input for 15 ms after each own pre arrival
Scores: per synapse (191, and by distance) and per connection (T summed over the contacts of each of the 24 pairs).
    python cell_audit/vdcc_decode/counterfactual.py  -> counterfactual.json, figs/counterfactual.png"""
import itertools, json, os
import numpy as np
from scipy.signal import lfilter
import replay as RP
A, C = RP.A, RP.C
HERE = RP.HERE
PROTOS = ["ap1", "burst", "epsp"] + A.TRAINS
A.MARGIN = 1.0
F = RP.F
KE, G0, CAINF = 62.0, 0.24, 6.5e-5
DIAMS = np.logspace(np.log10(0.2), np.log10(8.0), 80)


def cad(ica, diam):
    """CaDynamics 'cad' (modified_mechanisms/CaDynamics.mod) on a replayed ica (mA/cm2), per-row diam (um)."""
    drive = np.maximum(-1e4 * ica * (4 / diam[:, None]) / (2 * F), 0) / (1 + KE)
    tau = (1 + KE) / (G0 * 4 / diam)
    out = np.empty_like(ica)
    for i in range(ica.shape[0]):
        a = np.exp(-A.DT / tau[i])
        out[i] = CAINF + lfilter([tau[i] * (1 - a)], [1, -a], drive[i])
    return out


def fit_diam(ica, cai):
    err = []
    for d in DIAMS:
        err.append(((cad(ica, np.full(len(ica), d)) - cai) ** 2).sum(1))
    err = np.array(err)
    k = err.argmin(0)
    return DIAMS[k], 1 - err[k, np.arange(len(k))] / ((cai - cai.mean(1, keepdims=True)) ** 2).sum(1)


def build(recs):
    """stacked signals per protocol: dict p -> dict key -> (N, T)."""
    S = {p: {k: RP.stack(recs, p, k)[0] for k in ("v", "cai", "cacr", "ica_vdcc", "ica_nmda")} for p in PROTOS}
    gk = np.load(os.path.join(HERE, "gK.npz")); g, K = gk["g"], gk["K"]
    ih, il = RP.dend(S["ap1"]["v"])
    diam, r2 = fit_diam(ih + il, S["ap1"]["cai"])
    info = dict(diam_median=float(np.median(diam)), diam_IQR=np.percentile(diam, [25, 75]).round(2).tolist(),
                cai_fit_R2_median=float(np.median(r2)), cai_fit_R2_p10=float(np.percentile(r2, 10)))
    zero = lambda p: np.zeros_like(S[p]["v"])
    for p in PROTOS:
        s = S[p]; v = s["v"]
        s["nvdcc"] = -s["ica_vdcc"]
        for ljp in (0, 5, 10, 15, 20, 25, 30):
            I, Cv = RP.spine(v, zero(p), g, K, ljp=float(ljp))      # VDCC-only spine Ca (NMDAR share excluded)
            s[f"cvd{ljp}"] = Cv; s[f"ivd{ljp}"] = -I
        for name, hs, ls in (("used", 12.3, 10.0), ("hsrc", 0.0, 10.0), ("lsrc", 12.3, 0.0), ("src", 0.0, 0.0)):
            ih, il = RP.dend(v, hshift=hs, lshift=ls)
            s[f"cai_{name}"] = cad(ih + il, diam)
    # validation of the shaft replay on the other protocols
    info["cai_replay_vs_rec_peak_ratio"] = {p: float(np.median((S[p]["cai_used"].max(1) - CAINF) /
                                                              np.maximum(S[p]["cai"].max(1) - CAINF, 1e-12)))
                                            for p in ("ap1", "burst", "epsp", "sj20@-10")}
    return S, info


def drive(x, rest, p, arr, t):
    u = np.maximum(x - rest - p["th"], 0) / p["a"]
    if p.get("w"):
        u = u * A.veto_mask(t, arr, p["w"])[None, :]
    if p["kind"] == "lin":
        return A.lp(u, p["tauT"])
    S = A.lp(u, A.TAU_E1)
    return A.lp(np.maximum(S - p["te"], 0), p["tauT"])


def event(x, rest, p, arr, t):
    c = A.crossings(x - rest, p["th"])
    if p.get("w"):
        c = c * A.veto_mask(t, arr, p["w"])[None, :]
    return A.lp(np.maximum(A.imp_decay(c, A.TAU_E1) - p["te"], 0), p["tauT"])


def rows(S, key, fn, p, t):
    rest = S["ap1"][key][:, :1]
    T = {pp: fn(S[pp][key], rest, p, A.pre_times(pp) + 0.1, t[pp]) for pp in PROTOS}
    at = lambda pp, tt: T[pp][:, np.searchsorted(t[pp], tt + 0.1)]
    L = np.stack([at(pp, tt) for pp, tt in A.LTD], 1); N = np.stack([at(pp, tt) for pp, tt in A.NOLTD], 1)
    return L, N, at("epsp", 50.0)


def conn(X, pid):
    u = np.unique(pid)
    return np.stack([X[pid == i].sum(0) for i in u])


def grids(S, key):
    pk = (S["ap1"][key] - S["ap1"][key][:, :1]).max(1)
    med = np.median(pk[pk > 0])
    area = np.median(np.maximum(S["ap1"][key] - S["ap1"][key][:, :1], 0).sum(1) * A.DT)
    out = []
    for q, tt, w in itertools.product((0.0, 0.01, 0.03, 0.1, 0.3), (20.0, 40.0, 80.0), (0.0, 15.0)):
        out.append((drive, dict(kind="lin", th=med * q, a=area, tauT=tt, w=w)))
        for te in (0.1, 0.3, 1.0):
            out.append((drive, dict(kind="two", th=med * q, a=area, tauT=tt, te=te, w=w)))
    for q, te, tt, w in itertools.product((0.003, 0.01, 0.03, 0.1, 0.3), (0.6, 0.85), (20.0, 40.0), (0.0, 15.0)):
        out.append((event, dict(kind="ev", th=med * q, te=te, tauT=tt, w=w)))
    return out


def main():
    recs = A.load("og")
    t = recs[0]["t"]
    dist = np.concatenate([r["dist"] for r in recs])
    pid = np.concatenate([np.full(len(r["dist"]), i) for i, r in enumerate(recs)])
    S, info = build(recs)
    print(json.dumps(info))
    keys = (["cacr", "nvdcc"] + [f"cvd{j}" for j in (0, 5, 10, 15, 20, 25, 30)] + [f"ivd{j}" for j in (0, 25)]
            + ["cai", "cai_used", "cai_hsrc", "cai_lsrc", "cai_src"])
    res = dict(info=info, cands={})
    bins = ((0, 60), (60, 150), (150, 1e9))
    for key in keys:
        best = {}
        for fn, p in grids(S, key):
            L, N, E = rows(S, key, fn, p, t)
            fs = A.score(L, N, E, strict=True)
            fc = A.score(conn(L, pid), conn(N, pid), conn(E, pid), strict=True)
            k = p["kind"] + ("+v" if p.get("w") else "")
            if k not in best or fs[1] > best[k]["f_syn"]:
                th = fs[2]
                ok = (L.min(1) >= 2 * th) & (np.maximum(N.max(1), E) <= th) if np.isfinite(th) else np.zeros(len(L), bool)
                best[k] = dict(f_syn=fs[1], by_dist=[float(ok[(dist >= a) & (dist < b)].mean()) for a, b in bins],
                               f_conn_same_params=fc[1], params={kk: (round(vv, 12) if isinstance(vv, float) else vv) for kk, vv in p.items()})
            kc = k + "|conn"
            if kc not in best or fc[1] > best[kc]["f_conn"]:
                best[kc] = dict(f_conn=fc[1], params={kk: (round(vv, 12) if isinstance(vv, float) else vv) for kk, vv in p.items()})
        res["cands"][key] = best
        print(key, " ".join(f"{k}:{v['f_syn']:.2f}[{'/'.join(f'{x:.2f}' for x in v['by_dist'])}] c{v['f_conn_same_params']:.2f}"
                            for k, v in best.items() if "|" not in k),
              "| best conn " + " ".join(f"{k[:-5]}:{v['f_conn']:.2f}" for k, v in best.items() if "|" in k), flush=True)
    json.dump(res, open(os.path.join(HERE, "counterfactual.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
