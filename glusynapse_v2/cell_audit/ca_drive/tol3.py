"""Candidate 3: can a different T nonlinearity make the window insensitive to the per-bAP drive lambda?
Ideal impulses at the somatic AP times (ap1, burst), S += lambda, S' = -S/tau_E1, T = lp(g(S), tau_T); widest contiguous
lambda range (fold) that one theta_Tg passes (LTD rows >= 2 theta_Tg, no-LTD rows <= theta_Tg), as vdcc_decode/tolerance.py.
Then the best g on real data (linear VDCC influx, population scale) via common.py.  Login node.  -> out/tol3.json"""
import itertools, json, os
import numpy as np
from scipy.signal import lfilter
import common as C
A = C.A
LAMS = np.exp(np.linspace(np.log(1e-3), np.log(1e3), 1201))


def G(kind, S, Sslow=None):
    if kind == "lin":   return np.maximum(S - 1, 0)                     # current rule (theta_Te = 1 unit)
    if kind == "sqrt":  return np.sqrt(np.maximum(S - 1, 0))
    if kind == "step":  return (S > 1).astype(float)                     # Heaviside: duration code
    if kind == "log":   return np.maximum(np.log(np.maximum(S, 1e-30)), 0)
    if kind == "hill4": return S ** 4 / (1 + S ** 4)                    # divisive / saturating
    if kind == "ratio": return np.maximum(S / np.maximum(Sslow, 1e-30) - 0.5, 0)   # S_fast / S_slow (lambda-free)
    raise ValueError(kind)


def decay(imp, tau):
    return lfilter([1.0], [1.0, -np.exp(-A.DT / tau)], imp, axis=-1)


def fold(kind, tauE, tauT, t, spk):
    T = {}
    for p in ("ap1", "burst"):
        imp = np.zeros((len(LAMS), len(t[p])))
        for sp in spk[p]:
            imp[:, np.searchsorted(t[p], sp)] = LAMS
        S = decay(imp, tauE); Ss = decay(imp, 2000.0)
        T[p] = A.lp(G(kind, S, Ss), tauT)
    at = lambda p, x: T[p][:, np.searchsorted(t[p], x + 0.1)]
    lo = np.min([at(p, x) for p, x in A.LTD], 0); hi = np.max([at(p, x) for p, x in A.NOLTD], 0)
    best = (1.0, None)
    for th in np.unique(np.concatenate([hi, lo / 2]))[::3]:
        if th <= 0:
            continue
        ok = (lo >= 2 * th) & (hi <= th)
        if ok.any():
            idx = np.where(ok)[0]; r = max(np.split(idx, np.where(np.diff(idx) > 1)[0] + 1), key=len)
            f = LAMS[r[-1]] / LAMS[r[0]]
            if f > best[0]:
                best = (float(f), [float(LAMS[r[0]]), float(LAMS[r[-1]])])
    return best


def main():
    S, t, dist, pid = C.build()
    recs = A.load("og"); spk = {p: recs[0]["spk"][p] for p in ("ap1", "burst")}
    res = {"fold": {}}
    for kind, tauE, tauT in itertools.product(("lin", "sqrt", "step", "log", "hill4", "ratio"), (25.0, 50.0, 100.0, 200.0),
                                              (10.0, 20.0, 40.0)):
        f = fold(kind, tauE, tauT, t, spk)
        res["fold"][f"{kind}|tauE{tauE:.0f}|tauT{tauT:.0f}"] = f
    for kind in ("lin", "sqrt", "step", "log", "hill4", "ratio"):
        b = max(((k, v) for k, v in res["fold"].items() if k.startswith(kind + "|")), key=lambda kv: kv[1][0])
        print(f"{kind:6s} best fold {b[1][0]:8.2f} at {b[0]}  lambda {b[1][1]}", flush=True)
    json.dump(res, open(os.path.join(C.HERE, "out", "tol3.json"), "w"), indent=0)


if __name__ == "__main__":
    main()
