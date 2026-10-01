"""Post rho rule v4 candidates (RHO_REDESIGN.md). CPU reference (numba) for gpu_v4_rho.py and for CPU scoring.

Every option defaults to the v3 rule (batch_v2.thetas + batch_v2._rho_loop + binary readout), so v4 with the defaults
is v3. Options are carried in the filters / pre dict of a fit (saved in its json):

  rho_gamma  (A) threshold de-normalisation. theta_d = a00 c_pre + a01 cq, theta_p = a10 c_pre + a11 cq with
             cq = C_REF**(1 - gamma) * c_post**gamma. gamma 1 = v3 (cq = c_post exactly). gamma 0: the thresholds no
             longer scale with the synapse's own bAP Ca, so a weak-bAP (distal) synapse needs the same absolute
             pairing Ca as a proximal one. C_REF only sets units (absorbed by a01, a11).
  tau_fast   (B) potentiation reads a fast filter F of spine Ca (ms), depression keeps effcai (tau_effca 278 ms):
             F' = -F / tau_fast + cacr, cacr = spine Ca above rest recovered exactly from effcai per grid step
             (batch_v2.cacr_from_effcai). pot = F > theta_p, dep = effcai > theta_d. 0 = off (v3);
             tau_fast = tau_effca reproduces v3 up to rounding.
  rho_sigma  (C) probabilistic readout: synapse potentiated with probability Phi((rho_f - 0.5) / sigma) instead
             of the step rho_f >= 0.5. 0 = off (v3).
  vamp_mode, theta_V  ROUND2 C1 (1) / C2 (2): pot also needs the own filtered VDCC current V > theta_V
             (scan_vgate_amp.py is the CPU reference; gpu_v4_rho vamp kernel). 0 = off.
  gamma_d, gamma_p  rho depression / potentiation rates (uniform, one value each). Default = GAMMA_D_GB, GAMMA_P_GB
             (v3); free in fit_v4 --fit-gamma, bounds (50, 200), (150, 300) as analytical_method/fit.py.
"""
import math
import numpy as np
from numba import njit

import batch_v2
from batch_v2 import RHO_STAR_GB, GAMMA_D_GB, GAMMA_P_GB, TAU_IND_GB, TAU_EFFCA, RHO_BINARY_THRESHOLD

C_REF = 0.01          # mM, unit of cq when gamma < 1 (fixed; any value is absorbed by a01, a11)
V4_DEFAULTS = dict(rho_gamma=1.0, tau_fast=0.0, rho_sigma=0.0)
RATE_DEFAULTS = dict(gamma_d=float(GAMMA_D_GB), gamma_p=float(GAMMA_P_GB))
VAMP_DEFAULTS = dict(vamp_mode=0, theta_V=0.0)


def vamp(P):
    """(vamp_mode, theta_V) of a pre / filters dict (ROUND2 C1 / C2, scan_vgate_amp.py): 0 = off (A0 / v3)."""
    Q = {**VAMP_DEFAULTS, **(P or {})}
    return int(Q["vamp_mode"]), float(Q["theta_V"])


def opts(P):
    Q = {**V4_DEFAULTS, **(P or {})}
    return float(Q["rho_gamma"]), float(Q["tau_fast"]), float(Q["rho_sigma"])


def rates(P):
    """(gamma_d, gamma_p) of a pre / filters dict; absent = the v3 constants."""
    Q = {**RATE_DEFAULTS, **(P or {})}
    return float(Q["gamma_d"]), float(Q["gamma_p"])


def cq(c_post, gamma):
    """c_post (gamma 1, the same array: bitwise v3) or C_REF^(1-gamma) c_post^gamma."""
    if gamma == 1.0:
        return c_post
    return C_REF ** (1.0 - gamma) * np.maximum(c_post, 0.0) ** gamma


def thetas(c_pre, c_post, a, gamma=1.0):
    """batch_v2.BatchV2.thetas (basal a's, BAP_GATE honoured) with c_post -> cq."""
    c = cq(c_post, gamma)
    td, tp = a["a00"] * c_pre + a["a01"] * c, a["a10"] * c_pre + a["a11"] * c
    if batch_v2.BAP_GATE is not None:
        off = c_post < batch_v2.BAP_GATE
        td, tp = np.where(off, np.inf, td), np.where(off, np.inf, tp)
    return td, tp


@njit(cache=True)
def rho_loop_v4(E, H, rec, L, td, tp, rho0, gd, gp, rs, tau_fast, tau_eff, kt):
    """batch_v2._rho_loop with pot read from F when tau_fast > 0 (same step order as gpu_v4_rho). H in units of
    s / tau_ind (B._H); kt = 1e3 tau_ind converts to ms."""
    out = np.empty(rho0.shape[0])
    for i in range(rho0.shape[0]):
        r = rho0[i]; a = td[i]; b = tp[i]; j = rec[i]; F = 0.0
        for k in range(L[i]):
            x = E[i, k]
            xp = F if tau_fast > 0.0 else x
            if not (x <= a and x <= b and xp <= b and (r == 0.0 or r == 1.0)):
                pot = 1.0 if xp > b else 0.0
                dep = 1.0 if x > a else 0.0
                r += H[j, k] * (-r * (1 - r) * (rs - r) + pot * gp * (1 - r) - dep * (1 - pot) * gd * r)
                r = min(max(r, 0.0), 1.0)
            if tau_fast > 0.0:
                hm = H[j, k] * kt
                if hm > 0.0:
                    aE = math.exp(-hm / tau_eff); bE = tau_eff * (1.0 - aE)
                    u = (E[i, k + 1] - aE * x) / bE
                    aF = math.exp(-hm / tau_fast)
                    F = aF * F + tau_fast * (1.0 - aF) * u
        out[i] = r
    return out


def rho_all(B, a, P):
    """rho_f for every synapse of BatchV2 B (record order), v4 rule."""
    gamma, tau_fast, _ = opts(P)
    E, H = B._stack()
    td = np.concatenate([thetas(r["c_pre"], r["c_post"], a, gamma)[0] for r in B.recs])
    tp = np.concatenate([thetas(r["c_pre"], r["c_post"], a, gamma)[1] for r in B.recs])
    gd, gp = rates(P)
    return rho_loop_v4(E, H, B._rec, B._len, td, tp, B._rho0, gd, gp, float(RHO_STAR_GB),
                       float(tau_fast), float(TAU_EFFCA), 1e3 * TAU_IND_GB), td, tp


def rb_of(rho, sigma):
    """Readout state: step at 0.5 (sigma 0, v3) or Phi((rho - 0.5) / sigma)."""
    rho = np.asarray(rho, float)
    if sigma <= 0.0:
        return (rho >= RHO_BINARY_THRESHOLD).astype(float)
    from scipy.special import ndtr
    return ndtr((rho - RHO_BINARY_THRESHOLD) / sigma)


def ratio(b, rho0, rho_f, dpre_f, sigma=0.0, dpre0=0.0):
    """PairBasis.ratio with the v4 readout (baseline rho0 is 0/1, so its readout is unchanged)."""
    def epsp(rb, dpre):
        u0 = b.Use_d + rb * (b.Use_p - b.Use_d)
        scale = np.minimum(1.0, u0 * (1.0 + dpre)) / u0
        mean = float(np.sum((b.e0 * b.w + rb * b.delta) * scale))
        var = float(np.sum(rb * b.svar) + (1 - rb.sum()) ** 2 * b.v0)
        return mean, np.sqrt(max(var, 0.0))
    b_m, b_s = epsp(rb_of(rho0, 0.0), np.full(len(rho0), dpre0))
    a_m, _ = epsp(rb_of(rho_f, sigma), dpre_f)
    if b_m == 0:
        return np.nan
    return (a_m / b_m) * (1.0 + min((b_s / b_m) ** 2, 0.25))
