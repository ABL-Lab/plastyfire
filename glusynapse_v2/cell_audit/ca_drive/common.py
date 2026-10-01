"""Shared signals + scoring for the Ca-drive design (CA_DRIVE_DESIGN.md).  Login node, numpy/scipy only.
Signals per ljp (0, 25): VDCC-only spine influx q (nA, >0 inward, rest-subtracted) and VDCC-only spine Ca c (mM,
rest-subtracted), replayed on the recorded local V of local_t/out with the per-synapse g, K of vdcc_decode/gK.npz.
Score = local_t/check_hp.py (MARGIN 1, strict: LTD rows >= 2 th_Tg; no-LTD rows and own-EPSP read 50 ms later <= th_Tg,
one common parameter set), per synapse by distance bin and per connection (majority of contacts; mean T over contacts)
as vdcc_decode/conn_check.py.  The replay cache lives in the scratch dir given by $CA_CACHE (default: this dir)."""
import os, sys
import numpy as np
from scipy.signal import lfilter
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "vdcc_decode"))
import replay as RP          # noqa: E402  (read-only use)
A = RP.A
A.MARGIN = 1.0
PROTOS = ["ap1", "burst", "epsp"] + A.TRAINS
BINS = ((0, 60), (60, 150), (150, 250), (250, 1e9))
CACHE = os.environ.get("CA_CACHE", HERE)
TAU_D_NMDA = 70.0                            # GluSynapse_v2.mod tau_d_NMDA (existing parameter)


def build():
    recs = A.load("og")
    t = recs[0]["t"]
    dist = np.concatenate([r["dist"] for r in recs])
    pid = np.concatenate([np.full(len(r["dist"]), i) for i, r in enumerate(recs)])
    f = os.path.join(CACHE, "ca_drive_signals.npz")
    if os.path.exists(f):
        z = np.load(f)
        S = {lj: {p: {k: z[f"{lj}|{p}|{k}"].astype(np.float64) for k in ("q", "c", "v")} for p in PROTOS} for lj in (0, 25)}
        return S, t, dist, pid
    gk = np.load(os.path.join(RP.HERE, "gK.npz")); g, K = gk["g"], gk["K"]
    S, save = {}, {}
    for lj in (0, 25):
        S[lj] = {}
        for p in PROTOS:
            v = RP.stack(recs, p, "v")[0]
            I, Cc = RP.spine(v, np.zeros_like(v), g, K, ljp=float(lj))
            S[lj][p] = dict(q=-I, c=Cc, v=v)
        for p in PROTOS:                      # rest = pre-stimulus value of ap1
            for k in ("q", "c"):
                S[lj][p][k] = S[lj][p][k] - S[lj]["ap1"][k][:, :1]
            for k in ("q", "c", "v"):
                save[f"{lj}|{p}|{k}"] = S[lj][p][k].astype(np.float32)
    np.savez(f, **save)
    return S, t, dist, pid


def arrivals(p):
    return A.pre_times(p) + 0.1


def nmda_bound(tt, arr):
    """own glutamate-bound NMDAR state = the mod's B_NMDA (jumps at the own release, decays with tau_d_NMDA 70 ms),
    normalised to 1 per release, summed over own arrivals, capped at 1.  No new parameter."""
    b = np.zeros_like(tt)
    for a in arr:
        b += np.where(tt >= a, np.exp(-np.maximum(tt - a, 0.0) / TAU_D_NMDA), 0.0)
    return np.minimum(b, 1.0)


def veto_weight(kind, tt, arr, w=15.0, thb=0.8):
    """weight on the eCB priming input: none | hand (0 for w ms after own pre) | nmda (1 - b) | nmda_th (0 while b > thb)."""
    if kind == "none" or len(arr) == 0:
        return np.ones_like(tt)
    if kind == "hand":
        return A.veto_mask(tt, arr, w)
    b = nmda_bound(tt, arr)
    if kind == "nmda":
        return 1.0 - b
    if kind == "nmda_th":
        return (b <= thb).astype(float)
    raise ValueError(kind)


def rule(inc, te, tauT, arr=(), tt=None, reset=False):
    """S += inc (per step), S' = -S/tau_E1; optional reset of S to 0 at each own pre arrival (Ca-then-glutamate read
    consumes the primed state); T = lp(pos(S - te), tauT)."""
    a = np.exp(-A.DT / A.TAU_E1)
    if not reset or len(arr) == 0:
        S = lfilter([1.0], [1.0, -a], inc, axis=-1)
    else:
        idx = np.searchsorted(tt, arr)
        S = np.empty_like(inc); st = 0; s0 = np.zeros(inc.shape[0])
        for k in list(idx) + [inc.shape[1]]:
            if k > st:
                seg, zi = lfilter([1.0], [1.0, -a], inc[:, st:k], axis=-1, zi=(a * s0)[:, None])
                S[:, st:k] = seg; s0 = seg[:, -1]
            s0 = np.zeros(inc.shape[0]); st = k       # reset at the arrival step (read happens just before)
    return A.lp(np.maximum(S - te, 0.0), tauT)


def rows(Tp, t):
    at = lambda p, x: Tp[p][:, np.searchsorted(t[p], x + 0.1)]
    L = np.stack([at(p, x) for p, x in A.LTD], 1)
    N = np.stack([at(p, x) for p, x in A.NOLTD], 1)
    E = at("epsp", 50.0)
    TR = {}
    for k in A.TRAINS:
        f = float(k[2:4].strip("@")); pre = np.arange(5) * 1000.0 / f - float(k.split("@")[1])
        TR[k] = np.stack([at(k, x) for x in pre], 1)
    return L, N, E, TR


def conn_mean(X, pid):
    return np.stack([X[pid == i].mean(0) for i in np.unique(pid)])


def full_score(L, N, E, TR, dist, pid):
    f, fw, th = A.score(L, N, E, strict=True)
    ok = (L.min(1) >= 2 * th) & (np.maximum(N.max(1), E) <= th) if np.isfinite(th) else np.zeros(len(L), bool)
    u = np.unique(pid)
    gate = lambda k: round(float(np.median(np.tanh(np.maximum(TR[k] - th, 0) / th).mean(1))), 2) if np.isfinite(th) else None
    return dict(f_syn=round(float(fw), 3),
                bins=[round(float(ok[(dist >= a) & (dist < b)].mean()), 2) for a, b in BINS],
                conn_major=round(float(np.mean([ok[pid == i].mean() >= 0.5 for i in u])), 2),
                conn_meanT=round(float(A.score(conn_mean(L, pid), conn_mean(N, pid), conn_mean(E, pid))[1]), 2),
                th_Tg=float(th), trains={k: gate(k) for k in ("sj20@-10", "sj20@+10", "sj50@-10")},
                medL=np.round(np.median(L, 0), 2).tolist(), medN=np.round(np.median(N, 0), 3).tolist(),
                medE=round(float(np.median(E)), 3))
