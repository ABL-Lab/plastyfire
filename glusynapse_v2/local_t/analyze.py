"""Offline test of synapse-local T drives on local_t/out/*.npz (sims.py). Login node, numpy only.

Every candidate maps one local trace x(t) of a synapse to T(t); the eCB gate at a pre arrival is tanh(pos(T - theta_Tg)).
Read points (T at the pre arrival = pre + 0.1 ms):
  LTD    : ap1 at 10/25/50 ms (Sjostrom 2003 -10/-25/-50), burst at 120/200 ms after the last AP
  no LTD : ap1 at 100/120/200 ms. (+10/+25 single pairings: pre arrives before the post AP, so T = 0 for every
           post-driven drive; checked with the epsp run: T at any time in the epsp run must stay <= theta_Tg.)
  trains : sj{10,20,50}@{+-10}: mean gate over the 5 arrivals (reported, not scored; 20 Hz -10 wants LTD).
Scores per candidate over all synapses of all pairs (global parameters, one common theta_Tg):
  f_order: fraction of synapses with min(LTD) > max(no LTD) (per-synapse ordering, any theta)
  f_win  : fraction with min(LTD) >= (1 + m) theta_Tg and max(no LTD) <= theta_Tg for the best common theta_Tg
           (m = 0.25 margin, so the weakest LTD row still gates >= tanh(0.25 theta_Tg))
    python glusynapse_v2/local_t/analyze.py [--variant og]
"""
import argparse, glob, itertools, os
import numpy as np
from scipy.signal import lfilter

HERE = os.path.dirname(os.path.abspath(__file__))
DT = 0.1
TAU_E1 = 100.0
LTD = [("ap1", 10), ("ap1", 25), ("ap1", 50), ("burst", 320), ("burst", 400)]
NOLTD = [("ap1", 100), ("ap1", 120), ("ap1", 200)]
TRAINS = [f"sj{f}@{d}" for f in (10, 20, 50) for d in ("+10", "-10")]
MARGIN = 0.25


def lp(u, tau):                    # T' = -T/tau + u (exact for piecewise-constant u), per row
    a = np.exp(-DT / tau)
    return lfilter([tau * (1 - a)], [1.0, -a], u, axis=-1)


def imp_decay(imp, tau):           # S += imp, S' = -S/tau
    return lfilter([1.0], [1.0, -np.exp(-DT / tau)], imp, axis=-1)


_XC = {}
HYST = 2.0                         # V detectors re-arm only after V falls HYST mV below th


def crossings(x, th, hyst=0.0):    # 1 at each upward crossing of th (re-armed below th - hyst)
    c = np.zeros_like(x)
    if hyst <= 0:
        c[:, 1:] = (x[:, 1:] > th) & (x[:, :-1] <= th)
        return c
    armed = np.ones(x.shape[0], bool)
    for k in range(1, x.shape[1]):
        up = armed & (x[:, k] > th); c[up, k] = 1; armed &= ~up; armed |= x[:, k] < th - hyst
    return c


def veto_mask(t, arr, w):          # 0 within w ms after an own pre arrival (the own EPSP / glutamate-bound NMDARs)
    m = np.ones_like(t)
    for a in arr:
        m[(t >= a) & (t < a + w)] = 0.0
    return m


# candidates: name -> (signal key, fn(x, rest, p) -> T, param grid). rest: per-synapse pre-T0 value of the signal
def c_event(x, rest, p, arr=(), t=None):   # A1/B1/C1: crossing count S, supralinear pos(S - theta_Te), low-pass tau_T
    ck = (id(x), p["th"], p.get("hyst", 0.0), bool(p.get("rel")))
    if ck not in _XC:
        _XC[ck] = crossings(x - (rest if p.get("rel") else 0), p["th"], p.get("hyst", 0.0))
    c = _XC[ck]
    if p.get("w") and len(arr):                                       # D: own-pre veto
        c = c * veto_mask(t, arr, p["w"])[None, :]
    S = imp_decay(c, TAU_E1)
    return lp(np.maximum(S - p["te"], 0), p["tauT"])


def c_time(x, rest, p):            # A2: time above threshold (ms) instead of a count
    S = lp((x - (rest if p.get("rel") else 0) > p["th"]).astype(float), TAU_E1) / 1.0
    return lp(np.maximum(S - p["te"], 0), p["tauT"])


def c_twostage(x, rest, p):        # B2/C2: pos(x - rest - th)/scale -> S (tau_E1) -> pos(S - te) -> T
    S = lp(np.maximum(x - rest - p["th"], 0), TAU_E1)
    return lp(np.maximum(S - p["te"], 0), p["tauT"])


def c_linear(x, rest, p):          # A3: one low-pass of pos(x - th)
    return lp(np.maximum(x - (rest if p.get("rel") else 0) - p["th"], 0), p["tauT"])


def c_soma(x, rest, p):            # A0 reference = t_drive 2 (somatic AP impulses; cell-wide, NOT allowed)
    return c_event(x, rest, dict(p, th=-10.0))


TE = (0.6, 0.8, 0.9, 1.0, 1.2, 1.5)
TAUT = (20.0, 40.0, 80.0)
CANDS = {
    "A0 REFERENCE t_drive 2 (soma V, not allowed)": ("vsoma", c_soma, [dict(te=te, tauT=tt) for te in TE for tt in TAUT]),
    "A1 V event (abs th)": ("v", c_event, [dict(th=th, te=te, tauT=tt) for th in (-60, -55, -50, -45, -40)
                                           for te in TE for tt in TAUT]),
    "A1r V event (th rel. rest)": ("v", c_event, [dict(th=th, te=te, tauT=tt, rel=1) for th in (10, 15, 20, 25, 30)
                                                  for te in TE for tt in TAUT]),
    "A1h V event (rel., 2 mV hysteresis)": ("v", c_event, [dict(th=th, te=te, tauT=tt, rel=1, hyst=HYST)
                                                           for th in (3, 4, 5, 6, 8, 10) for te in TE for tt in TAUT]),
    "D1 V event (rel.) + own-pre veto": ("v", c_event, [dict(th=th, te=te, tauT=tt, rel=1, hyst=HYST, w=w)
                                                        for th in (4, 5, 6, 8, 10) for te in (0.6, 0.8, 0.9)
                                                        for tt in (20.0, 40.0) for w in (10.0, 15.0, 20.0, 30.0)]),
    "A2 V time-above": ("v", c_time, [dict(th=th, te=te, tauT=tt) for th in (-55, -50, -45)
                                      for te in (0.5, 1, 2, 3, 5) for tt in TAUT]),
    "A3 V linear": ("v", c_linear, [dict(th=th, tauT=tt, rel=1) for th in (5, 10, 20, 30) for tt in (20, 50, 100, 200)]),
    "B1 shaft Ca event": ("cai", c_event, [dict(th=th, te=te, tauT=tt, rel=1) for th in (1e-5, 3e-5, 1e-4, 3e-4)
                                           for te in TE for tt in TAUT]),
    "B2 shaft Ca two-stage": ("cai", c_twostage, None),
    "Bs spine Ca (cacr) event": ("cacr", c_event, [dict(th=th, te=te, tauT=tt, rel=1) for th in (1e-5, 3e-5, 1e-4, 3e-4)
                                                   for te in TE for tt in TAUT]),
    "C1 VDCC event": ("nvdcc", c_event, None),
    "C2 VDCC two-stage": ("nvdcc", c_twostage, None),
}


def pre_times(p):
    if p == "epsp":
        return np.array([0.0])
    if p.startswith("sj"):
        f = float(p[2:4].strip("@"))
        return np.arange(5) * 1000.0 / f - float(p.split("@")[1])
    return np.array([])


def load(variant):
    recs = []
    for f in sorted(glob.glob(os.path.join(HERE, "out", "*.npz"))):
        z = np.load(f)
        g = lambda p, k: z[f"{variant}|{p}|{k}"]
        r = dict(pair=os.path.basename(f)[:-4], dist=g("ap1", "dist"), kind=g("ap1", "kind"), sig={}, t={}, spk={}, arr={})
        for p in ["ap1", "burst", "epsp"] + TRAINS:
            if f"{variant}|{p}|t" not in z.files:
                continue
            r["t"][p] = g(p, "t").astype(np.float64); r["spk"][p] = g(p, "soma_spikes")
            r["arr"][p] = pre_times(p) + 0.1
            r["sig"][p] = {k: g(p, k).astype(np.float64) for k in ("v", "cai", "cacr", "effcai", "ica_vdcc", "ica_nmda")}
            r["sig"][p]["nvdcc"] = -r["sig"][p]["ica_vdcc"]
            imp = np.zeros_like(r["sig"][p]["v"][:1])                      # soma AP train as a -10 mV crossing
            for sp in g(p, "soma_spikes"):
                imp[0, np.searchsorted(r["t"][p], sp)] = 1.0
            r["sig"][p]["vsoma"] = np.repeat(-80.0 + 90.0 * imp, len(r["dist"]), 0)
        recs.append(r)
    return recs


def scale_grids(recs):             # data-scaled grids for continuous Ca / current candidates
    for key, name in (("cai", "B2 shaft Ca two-stage"), ("nvdcc", "C2 VDCC two-stage"), ("nvdcc", "C1 VDCC event")):
        pk = np.concatenate([(r["sig"]["ap1"][key] - r["sig"]["ap1"][key][:, :1]).max(1) for r in recs])
        med = np.median(pk[pk > 0]) if np.any(pk > 0) else 1.0
        if name.startswith("C1"):
            CANDS[name] = (key, c_event, [dict(th=med * q, te=te, tauT=tt, rel=1) for q in (0.01, 0.03, 0.1, 0.3)
                                          for te in TE for tt in TAUT])
        else:      # integral of pos(x - th) over one bAP ~ pk*w; normalise so one median bAP gives S ~ 1
            area = np.median(np.concatenate([np.maximum(r["sig"]["ap1"][key] - r["sig"]["ap1"][key][:, :1], 0).sum(1) * DT
                                             for r in recs]))
            CANDS[name] = (key, lambda x, rest, p, a=area: c_twostage(x / a, rest / a, p),
                           [dict(th=med * q / area, te=te, tauT=tt) for q in (0.0, 0.03, 0.1, 0.3)
                            for te in (0.2, 0.5, 1.0, 2.0) for tt in TAUT])


def evaluate(recs, key, fn, p):
    L, N, E, TR, D = [], [], [], {k: [] for k in TRAINS}, []
    for r in recs:
        rest = r["sig"]["ap1"][key][:, :1]
        T = {pp: (fn(r["sig"][pp][key], rest, p, r["arr"][pp], r["t"][pp]) if fn is c_event
                  else fn(r["sig"][pp][key], rest, p)) for pp in r["sig"]}
        at = lambda pp, tt: T[pp][:, np.searchsorted(r["t"][pp], tt + 0.1)]
        L.append(np.stack([at(pp, tt) for pp, tt in LTD], 1)); N.append(np.stack([at(pp, tt) for pp, tt in NOLTD], 1))
        E.append(at("epsp", 50.0) if "epsp" in T else np.zeros(len(rest)))   # what a 2nd pre 50 ms later would see
        for k in TRAINS:
            if k in T:
                f = float(k[2:4].strip("@")); post = np.arange(5) * 1000.0 / f; pre = post - float(k.split("@")[1])
                TR[k].append(np.stack([at(k, x) for x in pre], 1))
        D.append(r["dist"])
    L, N, E, D = np.concatenate(L), np.concatenate(N), np.concatenate(E), np.concatenate(D)
    return L, N, E, {k: np.concatenate(v) for k, v in TR.items() if v}, D


def score(L, N, E, strict=True):
    """strict: the own EPSP alone must not gate a 2nd pre spike 50 ms later (pre-only trains give no LTD)."""
    lo, hi = L.min(1), N.max(1)
    if strict:
        hi = np.maximum(hi, E)
    f_order = np.mean(lo > hi)
    best = (0.0, np.nan)
    for th in np.unique(np.concatenate([hi, lo / (1 + MARGIN)])):
        if th <= 0:
            continue
        f = np.mean((lo >= (1 + MARGIN) * th) & (hi <= th))
        if f > best[0]:
            best = (f, th)
    return f_order, best[0], best[1]


def fidelity(recs, key, th, rel):
    """fraction of synapses whose local signal crosses th exactly once per post AP (ap1: 1, burst: 5) and never
    for the own EPSP (epsp: 0)."""
    ok = []
    for r in recs:
        n = {}
        for p in ("ap1", "burst", "epsp"):
            x = r["sig"][p][key]
            n[p] = crossings(x - (r["sig"]["ap1"][key][:, :1] if rel else 0), th).sum(1)
        ok.append((n["ap1"] == 1) & (n["burst"] == 5) & (n["epsp"] == 0))
    return np.mean(np.concatenate(ok))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--variant", default="og")
    ap.add_argument("--loose", dest="strict", action="store_false"); a = ap.parse_args()
    recs = load(a.variant)
    print(f"{a.variant}: {len(recs)} pairs, {sum(len(r['dist']) for r in recs)} synapses")
    scale_grids(recs)
    for key, rel, ths in (("v", 0, (-65, -60, -55, -50, -45, -40, -30)), ("v", 1, (5, 10, 15, 20, 25, 30))):
        print(f"fidelity {key} {'rel' if rel else 'abs'}: " +
              " ".join(f"{th}:{fidelity(recs, key, th, rel):.2f}" for th in ths))
    med = np.median(np.concatenate([(-r["sig"]["ap1"]["ica_vdcc"]).max(1) for r in recs]))
    print("fidelity VDCC (x median bAP peak): " + " ".join(f"{q}:{fidelity(recs, 'nvdcc', med * q, 1):.2f}"
                                                         for q in (0.001, 0.003, 0.01, 0.03, 0.1)))
    best_all = {}
    for name, (key, fn, grid) in CANDS.items():
        best = None
        for p in grid:
            L, N, E, TR, D = evaluate(recs, key, fn, p)
            s = score(L, N, E, strict=a.strict)
            if best is None or (s[1], s[0]) > (best[0][1], best[0][0]):
                best = (s, p, L, N, E, TR, D)
        (fo, fw, th), p, L, N, E, TR, D = best
        best_all[name] = best
        med = lambda X: np.round(np.median(X, 0), 3).tolist()
        print(f"\n{name}: best {p}\n  f_order {fo:.2f}  f_win {fw:.2f} (theta_Tg {th:.3g})"
              f"\n  median T  LTD {med(L)}  noLTD {med(N)}  epsp->pre+50 {np.median(E):.3g} (frac > theta_Tg {np.mean(E > th) if np.isfinite(th) else np.nan:.2f})"
              f"\n  median ratio min(LTD)/max(noLTD) {np.median(L.min(1) / np.maximum(N.max(1), 1e-12)):.2f}")
        for lo_, hi_ in ((0, 60), (60, 150), (150, 1e9)):
            m = (D >= lo_) & (D < hi_)
            if m.any():
                lo, hi = L[m].min(1), N[m].max(1)
                ok = (lo >= (1 + MARGIN) * th) & (hi <= th) if np.isfinite(th) else np.zeros(m.sum(), bool)
                print(f"    dist {lo_}-{hi_}: n {m.sum()}  f_order {np.mean(lo > hi):.2f}  f_win {ok.mean():.2f}")
        if np.isfinite(th):
            gate = {k: np.round(np.median(np.tanh(np.maximum(v - th, 0)).mean(1)), 3) for k, v in TR.items()}
            print(f"  trains, median mean gate over 5 arrivals: {gate}")
    np.save(os.path.join(HERE, f"best_{a.variant}{'' if a.strict else '_loose'}.npy"), {k: (v[0], v[1]) for k, v in best_all.items()}, allow_pickle=True)


if __name__ == "__main__":
    main()
