"""Candidate D (RHO_REDESIGN.md): potentiation needs a recent own spine VDCC Ca event. CPU reference for gpu_v4_rho
(vgate) and CPU scoring. Option key vgate (0/1) in the fit filters.

    G += 1 at each upward crossing of -ica_VDCC through K_ca at the synapse (cev, the t_drive 4 events, unweighted),
    G' = -G / tau_E1 (100 ms, the t_drive 4 constant), pot = [x_p > theta_p] * [G > 1/2]  (window tau_E1 ln 2 = 69 ms).
No new parameter. A synapse whose spine never sees a depolarisation-gated Ca event (no bAP Ca, no strong local EPSP)
cannot potentiate; its depression band is unchanged.
"""
import math
import numpy as np
from numba import njit

import model_v2 as MV
import rho_v4
from batch_v2 import RHO_STAR_GB, GAMMA_D_GB, GAMMA_P_GB, TAU_IND_GB, TAU_EFFCA


def gate_events(B):
    """-> ptr (N+1,), step index (M,), count (M,) of the unweighted cev events per synapse, record order."""
    st, ct, ptr = [], [], [0]
    for r in B.recs:
        t = r["t"]
        for i in range(len(r["syn"])):
            ev = np.asarray(r["cev"][i], np.float64); ev = ev[np.isfinite(ev)]
            u, c = np.unique(np.clip(np.searchsorted(t, ev - 1e-9), 0, len(t) - 1), return_counts=True)
            st.append(u); ct.append(c.astype(np.float64)); ptr.append(ptr[-1] + len(u))
    cat = lambda xs, dt: np.concatenate(xs).astype(dt) if sum(len(x) for x in xs) else np.zeros(1, dt)
    return np.asarray(ptr, np.int64), cat(st, np.int64), cat(ct, np.float64)


@njit(cache=True)
def rho_loop_v4d(E, H, rec, L, td, tp, rho0, gd, gp, rs, tau_fast, tau_eff, kt, vgate, tau_g, gptr, gstep, gcnt):
    """rho_v4.rho_loop_v4 plus the gate, same step order as the GPU kernel (no skipping)."""
    out = np.empty(rho0.shape[0])
    for i in range(rho0.shape[0]):
        r = rho0[i]; a = td[i]; b = tp[i]; j = rec[i]; F = 0.0; G = 0.0
        kg = gptr[i]; kge = gptr[i + 1]
        ng = gstep[kg] if kg < kge else -1
        for k in range(L[i]):
            x = E[i, k]
            xp = F if tau_fast > 0.0 else x
            if vgate and k == ng:
                G += gcnt[kg]; kg += 1
                ng = gstep[kg] if kg < kge else -1
            pot = 1.0 if xp > b else 0.0
            if vgate and G <= 0.5:
                pot = 0.0
            dep = 1.0 if x > a else 0.0
            r += H[j, k] * (-r * (1 - r) * (rs - r) + pot * gp * (1 - r) - dep * (1 - pot) * gd * r)
            r = min(max(r, 0.0), 1.0)
            hm = H[j, k] * kt
            if vgate:
                G = G * math.exp(-hm / tau_g)
            if tau_fast > 0.0 and hm > 0.0:
                aE = math.exp(-hm / tau_eff); bE = tau_eff * (1.0 - aE)
                u = (E[i, k + 1] - aE * x) / bE
                aF = math.exp(-hm / tau_fast)
                F = aF * F + tau_fast * (1.0 - aF) * u
        out[i] = r
    return out


def rho_all(B, a, P, gev=None):
    """rho_f for every synapse (record order) with options rho_gamma, tau_fast, vgate, gamma_d/gamma_p."""
    gamma, tau_fast, _ = rho_v4.opts(P)
    vg = bool({**MV.DEFAULTS, **P}.get("vgate", 0))
    E, H = B._stack()
    th = [rho_v4.thetas(r["c_pre"], r["c_post"], a, gamma) for r in B.recs]
    td = np.concatenate([t[0] for t in th]); tp = np.concatenate([t[1] for t in th])
    if gev is None:
        gev = gate_events(B) if vg else (np.zeros(len(B._rho0) + 1, np.int64), np.zeros(1, np.int64), np.zeros(1))
    gd, gp = rho_v4.rates(P)
    return rho_loop_v4d(E, H, B._rec, B._len, td, tp, B._rho0, gd, gp, float(RHO_STAR_GB),
                        float(tau_fast), float(TAU_EFFCA), 1e3 * TAU_IND_GB, vg, float({**MV.DEFAULTS, **P}["tau_E1"]),
                        *gev)
