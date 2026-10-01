"""Candidate synapse-local Ca drives for the eCB-LTD trace (CA_DRIVE_DESIGN.md), scored per synapse / bin / connection.
Login node.   CA_CACHE=<scratch> python cell_audit/ca_drive/cands.py [family ...]  -> out/cands_<family>.json
Every drive gives an increment per step inc(t) to S (S' = -S/tau_E1 + inc), T = lp(pos(S - te), tauT), gate at own pre.
  EV     : C1 reference, upward crossings of q above 0.3% of the median bAP peak
  LINQ   : linear Ca entry, inc = q dt / (population-median bAP charge)          (literature linear sensor, no norm.)
  HRc/HRq: Hill rate, inc = Hill(x; K, n) dt / w0 (x = c or q)                   (saturating sensor, time-above-K)
  HOc/HOq: Hill onset, inc = pos(d Hill(x; K, n))                                  (saturating sensor, bound-fraction rise)
  NRM    : own-history normalisation, inc = q dt / M_i, M_i = own (vetoed) per-event charge (steady state, see hist)
  NRMc   : same on the Ca integral (c dt / own bAP c area): the VDCC_DECODE oracle, for contrast
Veto of the input: none | hand (15 ms) | nmda (x (1 - b_NMDA(t))) | nmda_th (0 while b_NMDA > 0.8); '+R' also resets S
at the own pre arrival (Ca-then-glutamate read consumes the primed state)."""
import itertools, json, os, sys
import numpy as np
import common as C
A = C.A
TE = (0.3, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 1.0, 1.2, 1.5)
TAUT = (20.0, 40.0)
KQ = 10.0 ** np.arange(-4.0, 0.51, 0.5)
NH = (1.0, 2.0, 4.0, 8.0)
VETOS = ("none", "hand", "nmda", "nmda_th", "nmda+R")


def hill(x, K, n):
    x = np.maximum(x, 0.0) / K
    return x ** n / (1.0 + x ** n)


def evaluate(S, t, dist, pid, incf, grid_te=TE, veto="none"):
    """incf(p) -> increments (N, T) before veto.  Best (f_syn, conn) over te x tauT."""
    kind, reset = (veto[:-2], True) if veto.endswith("+R") else (veto, False)
    inc = {}
    for p in C.PROTOS:
        arr = C.arrivals(p)
        inc[p] = incf(p) * C.veto_weight(kind, t[p], arr)[None, :]
    best = None
    for te, tt in itertools.product(grid_te, TAUT):
        Tp = {p: C.rule(inc[p], te, tt, C.arrivals(p), t[p], reset) for p in C.PROTOS}
        L, N, E, TR = C.rows(Tp, t)
        f = A.score(L, N, E, strict=True)[1]
        if best is None or f > best[0]:
            best = (f, te, tt, (L, N, E, TR))
    f, te, tt, (L, N, E, TR) = best
    r = C.full_score(L, N, E, TR, dist, pid); r.update(te=te, tauT=tt, veto=veto)
    return r


def charge(x):
    return np.maximum(x, 0).sum(1) * A.DT


def run_family(fam, S, t, dist, pid, lj, vetos=VETOS):
    s = S[lj]; out = []
    qmed = np.median(s["ap1"]["q"].max(1)); cmed = np.median(s["ap1"]["c"].max(1))
    Amed = np.median(charge(s["ap1"]["q"]))
    if fam == "EV":
        th = 0.003 * qmed
        for v in vetos:
            out.append(dict(evaluate(S[lj], t, dist, pid, lambda p: A.crossings(s[p]["q"], th), veto=v), th_rel=0.003))
    elif fam == "LINQ":
        for v in vetos:
            out.append(evaluate(S[lj], t, dist, pid, lambda p: np.maximum(s[p]["q"], 0) * A.DT / Amed,
                                grid_te=tuple(sorted({x * f for x in TE for f in (0.01, 0.1, 1.0)})), veto=v))
    elif fam in ("HRc", "HRq", "HOc", "HOq"):
        key = "c" if fam.endswith("c") else "q"; med = cmed if key == "c" else qmed
        for kq, n in itertools.product(KQ, NH):
            K = kq * med
            if fam.startswith("HR"):
                w0 = np.median(hill(s["ap1"][key], K, n).sum(1) * A.DT)
                incf = lambda p, K=K, n=n, w0=w0: hill(s[p][key], K, n) * A.DT / w0
            else:
                def incf(p, K=K, n=n):
                    h = hill(s[p][key], K, n); d = np.zeros_like(h); d[:, 1:] = np.maximum(np.diff(h, axis=1), 0)
                    return d
            for v in ("hand", "nmda"):
                out.append(dict(evaluate(S[lj], t, dist, pid, incf, veto=v), K_rel=float(kq), n=n))
    elif fam in ("NRM", "NRMc"):
        key = "q" if fam == "NRM" else "c"
        for v in vetos:
            kind = v[:-2] if v.endswith("+R") else v
            wE = C.veto_weight(kind, t["epsp"], C.arrivals("epsp"))[None, :]
            Ab = charge(s["ap1"][key]); Ae = charge(s["epsp"][key] * wE)
            for hist, rho, q0 in (("pair", 0.0, 0.0), ("mix rho1", 1.0, 0.0), ("mix rho3", 3.0, 0.0),
                                  ("pair, floor p10", 0.0, 10), ("pair, floor p25", 0.0, 25),
                                  ("pre-only E, floor p10", -1, 10), ("pre-only E, floor p25", -1, 25),
                                  ("pre-only E, no floor", -1, 0)):
                M0 = np.percentile(Ab, q0) if q0 else 0.0
                M = np.maximum((Ab + max(rho, 0) * Ae) / (1 + max(rho, 0)), M0)
                ME = np.maximum(Ae, M0) if rho < 0 else M        # history for the E row (pre-only baseline)
                incf = lambda p, M=M, ME=ME: np.maximum(s[p][key], 0) * A.DT / (ME if p == "epsp" else M)[:, None]
                out.append(dict(evaluate(S[lj], t, dist, pid, incf, veto=v), hist=hist))
    return out


def main():
    fams = sys.argv[1:] or ["EV", "LINQ", "HRc", "HRq", "HOc", "HOq", "NRM", "NRMc"]
    S, t, dist, pid = C.build()
    for fam in fams:
        fn = os.path.join(C.HERE, "out", f"cands_{fam}.json"); os.makedirs(os.path.dirname(fn), exist_ok=True); res = {}
        for lj in (0, 25):
            rs = run_family(fam, S, t, dist, pid, lj)
            rs.sort(key=lambda r: (-r["f_syn"], -r["conn_meanT"]))
            res[f"{fam}|ljp{lj}"] = rs
            for r in rs[:4]:
                print(fam, lj, {k: r[k] for k in r if k not in ("medL", "medN")}, flush=True)
            json.dump(res, open(fn, "w"), indent=0)


if __name__ == "__main__":
    main()
