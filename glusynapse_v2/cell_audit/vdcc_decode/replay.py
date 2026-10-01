"""Open-loop replay of Ca channel gating on the recorded local V of local_t/out (0.1 ms, 4 linear-interpolated
substeps, exponential Euler).  The spine VDCC carries ~0.1 nS at most, so its gating does not feed back on V:
replaying it with another ljp is exact up to sampling.  Dendritic Ca_HVA2/Ca_LVAst replays are open-loop estimates.
Module (used by decode.py / counterfactual.py)."""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "local_t"))
sys.path.insert(0, HERE)
import analyze as A          # noqa: E402  (local_t/analyze.py, read-only use)
import channels as C         # noqa: E402
R, F, CEL = 8.315, 96485.309, 34.0
NSUB = 4
ECA_SHAFT = 1000 * R * (CEL + 273.15) / (2 * F) * np.log(2.0 / 6.5e-5)   # ~137 mV (cao 2 mM, cai 65 nM)


def nernst(ci, co=2.0):
    return 1000 * R * (CEL + 273.15) / (2 * F) * np.log(co / np.maximum(ci, 1e-9))


def gate(v, fn, **kw):
    """m^2 h(t) for a channel fn(v) -> (minf, mtau, hinf, htau), replayed on v (n, T) at 0.1 ms."""
    mi, mt, hi, ht = fn(v[:, 0], **kw)
    m, h = mi.copy(), hi.copy()
    out = np.empty_like(v)
    out[:, 0] = m * m * h
    dt = A.DT / NSUB
    for k in range(1, v.shape[1]):
        for s in range(1, NSUB + 1):
            vv = v[:, k - 1] + (v[:, k] - v[:, k - 1]) * s / NSUB
            mi, mt, hi, ht = fn(vv, **kw)
            m += (1 - np.exp(-dt / mt)) * (mi - m); h += (1 - np.exp(-dt / ht)) * (hi - h)
        out[:, k] = m * m * h
    return out


def spine(v, inmda, g, K, ljp=0.0, mtau=1.0, cmin=70e-6, tau=12.0):
    """spine VDCC current (nA, GluSynapse sign) and cai_CR (mM), with cai_CR feeding Eca_syn as in the mod."""
    fn = lambda x, ljp=ljp: (lambda r: (r[0], r[1] * mtau, r[2], r[3]))(C.rtype(x, ljp))
    mi, mt, hi, ht = fn(v[:, 0])
    m, h, c = mi.copy(), hi.copy(), np.full(v.shape[0], cmin)
    I, Cc = np.empty_like(v), np.empty_like(v)
    dt = A.DT / NSUB
    for k in range(v.shape[1]):
        for s in range(NSUB if k else 1):
            vv = v[:, max(k - 1, 0)] + (v[:, k] - v[:, max(k - 1, 0)]) * (s + 1) / NSUB if k else v[:, 0]
            mi, mt, hi, ht = fn(vv)
            if k:
                m += (1 - np.exp(-dt / mt)) * (mi - m); h += (1 - np.exp(-dt / ht)) * (hi - h)
            i = g * m * m * h * (vv - nernst(c))
            if k:
                c += dt * (-K * (inmda[:, k] + i) - (c - cmin) / tau)
        I[:, k], Cc[:, k] = i, c
    return I, Cc


def fit_gK(v, ivdcc, inmda, cacr):
    """per-synapse g (ica_VDCC = g m^2 h (v - Eca(cacr))) and K (cacr' = -K (ica_nmda + ica_vdcc) - (c - cmin)/12)
    from the recorded og traces (list of (n, T) arrays over protocols)."""
    num = den = kn = kd = 0.0
    for vv, iv, im, cc in zip(v, ivdcc, inmda, cacr):
        x = gate(vv, C.rtype) * (vv - nernst(cc))
        num = num + (iv * x).sum(1); den = den + (x * x).sum(1)
        lhs = np.diff(cc, axis=1) / A.DT + (cc[:, :-1] - 70e-6) / 12.0
        rhs = -(im + iv)[:, :-1]
        kn = kn + (lhs * rhs).sum(1); kd = kd + (rhs * rhs).sum(1)
    return num / np.maximum(den, 1e-30), kn / np.maximum(kd, 1e-30)


def dend(v, hshift=12.3, lshift=10.0, g_hva=0.0023897659621395684, g_lva=0.0023897659621395684, vshift=0.0):
    """dendritic Ca_HVA2 + Ca_LVAst current density (mA/cm2, inward < 0) replayed on v; vshift moves both curves."""
    gh = gate(v + vshift, C.hva, hshift=hshift); gl = gate(v + vshift, C.lva, shift=lshift)
    return g_hva * gh * (v - ECA_SHAFT), g_lva * gl * (v - ECA_SHAFT)


def stack(recs, p, key):
    """concatenate one signal of protocol p over all pairs (same t grid) -> (N, T), and the per-rec slices."""
    arrs = [r["sig"][p][key] for r in recs]
    idx = np.cumsum([0] + [a.shape[0] for a in arrs])
    return np.concatenate(arrs), [slice(idx[i], idx[i + 1]) for i in range(len(arrs))]
