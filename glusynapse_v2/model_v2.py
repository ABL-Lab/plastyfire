"""Offline GluSynapse_v2 presynaptic model (tsk T15) and EPSP readout (T16, Option B).

The post rule (rho) is unchanged, so rho still comes from analytical_method (model.integrate_rho or
batch.Batch). This module adds dpre, which in prefire traces never feeds back into the calcium,
so it can be computed from shaft_cai and the spike times alone.

Mod equations (glusynapse_v2/mod/GluSynapse_v2.mod, z_het = 0):
    T' = -T/tau_T  + u_T,   N' = -N/tau_NO + u_N,  with the drive u (pre_drive):
        0: u_T = pos(cai - rest - theta_T)/ca_scale,  u_N = pos(cai - rest - theta_NO)/ca_scale
        1: u_T = pos(-ica_VDCC - theta_Ti)/i_scale,   u_N = pos(-ica_VDCC - theta_NOi)/i_scale
    t_drive 1 (T25, L5 eCB pathway) replaces u_T by u_T = pos(effcai/c_post - theta_Te)/Te_scale: eCB release
    needs postsynaptic Ca from any source (bAPs, subthreshold NMDAR/VDCC, Sjostrom 2003/2004), and effcai/c_post
    is 1 at the peak of a single bAP at that synapse. batch_v2._sig passes effcai/c_post as sig.
    t_drive 2 (T25): cell-level eCB driven by postsynaptic APs, the same for all synapses of the cell:
        S' = -S/tau_E1,  S += 1 at each post AP;   u_T = pos(S - theta_Te)/Te_scale   (supralinear: bursts)
    sig is then the post-AP impulse train (batch_v2._sig). Local Ca can't make the Sjostrom window: at most of
    our synapses the attenuated bAP barely lifts spine Ca, and a linear filter can't give LTD for a burst ending
    200 ms before the pre spike but none for one AP 120 ms before (DECISIONS 2026-09-29).
    All drives: the gate at a pre spike is tanh(pos(T - theta_Tg)); theta_Tg = 0 is the old tanh(T) (T >= 0).
    The feature functions take sig = shaft_cai (drive 0) or -ica_VDCC (drive 1, extract_v2 'vdcc').
    Z' = -Z/tau_Z,                     Z += 1 at each spike arrival
    dpre' = A_NO tanh(N) tanh(pos(Z - theta_Z)) (dmax - dpre) / (1e3 tau_ind)   (theta_Z = 0 before v2.2)
    at each spike arrival: dpre -= A_mglu tanh(pos(T - theta_Tg)) (dpre - dmin)

A_NO and A_mglu enter only through two per-synapse quantities, so the solution is exact event by event:
    K_j   = integral of tanh(N) tanh(pos(Z - theta_Z)) dt over the interval before arrival j   (ms)
    tT_j  = tanh(T) at arrival j
    between arrivals: dmax - dpre  *= exp(-A_NO K_j / (1e3 tau_ind))
    at arrival j:     dpre - dmin  *= 1 - A_mglu tT_j
Filters (theta_*, tau_*) cost one lfilter pass over the trace. The amplitudes cost a recurrence over
the spike count, which is what DE sweeps.
"""
import os
import numpy as np
from scipy.signal import lfilter

HERE = os.path.dirname(os.path.abspath(__file__))
EDGES = os.environ.get("PLASTYFIRE_EDGES", "/project/rrg-emuller/dhuruva/plastyfire/data/dhuruva_modified_edges.h5")
EDGE_POP = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"
TAU_IND = 70.0          # s, same as GluSynapse tau_ind_GB

DEFAULTS = dict(ca_sh_rest=6.5e-5, ca_scale=1e-3, theta_T=1e-4, tau_T=50.0,
                theta_NO=5e-4, tau_NO=5.0, tau_Z=30.0, theta_Z=0.0, A_mglu=0.0, A_NO=0.0,
                dpre_min=-0.8, dpre_max=1.0, dpre0=0.0,
                pre_drive=0, i_scale=1e-6, theta_Ti=1e-7, theta_NOi=1e-7,
                no_drive=0, theta_NOe=0.05, e_scale=0.01,
                theta_NOc=0.0, c_scale=1e-3, theta_N=0.0,
                t_drive=0, theta_Te=0.5, Te_scale=1.0, tau_E1=100.0, theta_Tg=0.0,
                theta_V=3.0, w_V=15.0, tau_b=300.0,    # t_drive 3 (D1): fixed constants, not fitted
                # t_drive 4 (Ca drive, CA_DRIVE_DESIGN.md): fixed constants, not fitted.
                # K_ca (nA) = 0.2% of the median per-bAP peak of -ica_VDCC over the 191 L5 subset24 synapses
                # (local_t/out/*.npz, og variant, ap1: median peak 1.2344e-5 nA). Events at K_ca * K_mult are
                # read from the extracted cev (1), cev_lo (0.5) or cev_hi (2). tau_d_NMDA is the mod value.
                K_ca=2.469e-8, K_mult=1.0, tau_d_NMDA=70.0)


# ---------------------------------------------------------------- edge parameters
def edge_params(syn_ids, cache=os.path.join(HERE, "edge_params.npz")):
    """Use_d, Use_p, gmax_d, gmax_p and delay per edge id. Cached locally so the shared h5 is read once."""
    syn_ids = np.asarray(syn_ids, dtype=np.int64)
    have = {}
    if os.path.isfile(cache):
        c = np.load(cache)
        have = {k: dict(zip(c["syn"], c[k])) for k in c.files if k != "syn"}
    missing = np.setdiff1d(syn_ids, np.fromiter(have.get("delay", {}).keys(), np.int64))
    if missing.size:
        import h5py
        with h5py.File(EDGES, "r") as f:
            p = f[f"edges/{EDGE_POP}/0"]
            idx = np.sort(missing)
            new = {k: p[src][idx] for k, src in (("Use_d", "Use_d_TM"), ("Use_p", "Use_p_TM"),
                                                 ("gmax_d", "gmax_d_AMPA"), ("gmax_p", "gmax_p_AMPA"),
                                                 ("delay", "delay"))}
        for k, v in new.items():
            have.setdefault(k, {}).update(zip(idx, v))
        allsyn = np.fromiter(have["delay"].keys(), np.int64)
        np.savez(cache, syn=allsyn, **{k: np.array([have[k][s] for s in allsyn]) for k in have})
    return {k: np.array([have[k][s] for s in syn_ids], dtype=np.float64) for k in have}


# ---------------------------------------------------------------- filters
# The grid is t_ms (n_t,), uniform or the Ebner windowed grid (0.25 ms around spikes, 5 ms in the gaps).
# Fine steps use NEURON's Euler coefficients; the 5 ms gap step uses the exact exponential step, which
# stays stable for tau_Z ~ 10 ms (tests/test_offline_gate.py checks both grids against the mod).
def _runs(h):
    """[(s, e)] step ranges k in [s, e) of constant step h[k]."""
    br = np.flatnonzero(np.abs(np.diff(h)) > 1e-9) + 1
    ed = np.concatenate([[0], br, [len(h)]])
    return list(zip(ed[:-1], ed[1:]))


def _coef(h, tau):
    """NEURON's Euler step where it is accurate (h/tau < 0.1), the exact exponential step otherwise."""
    if h / tau < 0.1:
        return 1.0 - h / tau, h
    a = np.exp(-h / tau)
    return a, tau * (1 - a)


def _lowpass_grid(u, t_ms, tau):
    """x[k+1] = a_k x[k] + b_k u[k], x[0] = 0 (see _coef)."""
    h = np.diff(t_ms); y = np.zeros_like(u)
    for s, e in _runs(h):
        a, b = _coef(h[s], tau)
        y[:, s + 1:e + 1], _ = lfilter([b], [1.0, -a], u[:, s:e], axis=-1, zi=a * y[:, s:s + 1])
    return y


def _impulse_grid(imp, t_ms, tau):
    """y[k] = a_{k-1} y[k-1] + imp[k]: jumps on the arrival sample, then decays."""
    h = np.diff(t_ms); y = np.zeros_like(imp); y[:, 0] = imp[:, 0]
    for s, e in _runs(h):
        a, _b = _coef(h[s], tau)
        y[:, s + 1:e + 1], _ = lfilter([1.0], [1.0, -a], imp[:, s + 1:e + 1], axis=-1, zi=a * y[:, s:s + 1])
    return y


V_HYST = 2.0    # mV: a crossing re-arms only after v - vbar falls this far below theta_V (local_t/analyze.HYST)


def v_events(v, dt, t0=0.0, theta_V=None, tau_b=None, hyst=V_HYST):
    """t_drive 3 (D1): times (ms) of upward crossings of (v - vbar) > theta_V, v the segment voltage on a uniform
    grid of step dt (ms) starting at t0, vbar its slow EMA (tau_b, starts at v[0]). Full-resolution trace in,
    event times out; one crossing per excursion (re-armed below theta_V - hyst)."""
    theta_V = DEFAULTS["theta_V"] if theta_V is None else theta_V
    tau_b = DEFAULTS["tau_b"] if tau_b is None else tau_b
    v = np.asarray(v, np.float64); a = np.exp(-dt / tau_b)
    vb, _ = lfilter([1.0 - a], [1.0, -a], v, zi=[a * v[0]])
    x = v - vb
    hi = x > theta_V; lo = x < theta_V - hyst
    # hysteresis without a Python loop: keep samples that are hi or lo, one event per hi run
    idx = np.flatnonzero(hi | lo); m = hi[idx]
    first = np.ones(len(idx), bool); first[1:] = m[1:] != m[:-1]
    out = idx[first & m]
    return t0 + dt * np.asarray(out, np.float64)


CA_HYST = 0.5     # t_drive 4: a crossing of K re-arms only after -ica_VDCC falls below K * CA_HYST


def ca_events(x, dt, t0=0.0, K=None, hyst=CA_HYST, block=2048):
    """t_drive 4: times (ms) of upward crossings of x (= -ica_VDCC, nA, uniform grid of step dt starting at t0)
    through the fixed threshold K (default DEFAULTS['K_ca']); one event per excursion (re-armed below K * hyst).
    Speed: blocks whose maximum is below K * hyst are all 'lo' and are represented by their first sample
    (repeated lo samples change nothing), so only active blocks are scanned in full (identical result)."""
    K = DEFAULTS["K_ca"] if K is None else K
    x = np.asarray(x, np.float64); n = len(x); nb = n // block
    if nb >= 4:
        bmax = x[:nb * block].reshape(nb, block).max(axis=1)
        act = np.repeat(bmax >= K * hyst, block)
        keep = np.ones(n, bool); keep[:nb * block] = act
        keep[np.arange(nb)[~(bmax >= K * hyst)] * block] = True          # first sample of each quiet block
        pos = np.flatnonzero(keep); xs = x[pos]
    else:
        pos = None; xs = x
    hi = xs > K; lo = xs < K * hyst
    idx = np.flatnonzero(hi | lo); m = hi[idx]
    first = np.ones(len(idx), bool); first[1:] = m[1:] != m[:-1]
    out = idx[first & m]
    if pos is not None:
        out = pos[out]
    return t0 + dt * np.asarray(out, np.float64)


def glu_weight(ev, arr_t, tau_d=None):
    """t_drive 4: weight 1 - b(t_event) per event. b = the synapse's own glutamate-bound NMDA state: +1 at each own
    release (here: each own pre arrival), decay tau_d_NMDA, capped at 1; arrivals at or before the event count.
    ev: event times (NaN padded row ok) -> (finite events, weights)."""
    tau_d = DEFAULTS["tau_d_NMDA"] if tau_d is None else tau_d
    ev = np.asarray(ev, np.float64); ev = ev[np.isfinite(ev)]
    arr_t = np.sort(np.asarray(arr_t, np.float64))
    b = np.zeros(len(ev))
    for a in arr_t:
        m = ev >= a - 1e-9                        # same sample convention as the impulse grid
        b[m] += np.exp(-(ev[m] - a) / tau_d)
    return ev, 1.0 - np.minimum(b, 1.0)


def veto(ev, arr_t, w):
    """Drop events (padded NaN row ok) within w ms after any of one synapse's own pre arrival times arr_t."""
    ev = np.asarray(ev, np.float64); ev = ev[np.isfinite(ev)]
    arr_t = np.sort(np.asarray(arr_t, np.float64))
    if not len(ev) or not len(arr_t):
        return ev
    j = np.searchsorted(arr_t, ev, side="right") - 1          # last arrival at or before each event
    dtl = np.where(j >= 0, ev - arr_t[np.maximum(j, 0)], np.inf)
    return ev[dtl >= w]


def grid(n_t, dt=None, t0=0.0, t_ms=None):
    return np.asarray(t_ms, np.float64) if t_ms is not None else t0 + dt * np.arange(n_t)


def arrival_index(prespikes, delay, t_ms):
    """First grid index at or after each spike arrival (prespike + NetCon delay), per synapse: (n_syn, n_spk)."""
    ta = np.asarray(prespikes)[None, :] + np.asarray(delay)[:, None]
    return np.clip(np.searchsorted(t_ms, ta - 1e-9), 0, len(t_ms) - 1)


T_KEYS = ("pre_drive", "ca_sh_rest", "ca_scale", "theta_T", "i_scale", "theta_Ti", "tau_T", "t_drive", "theta_Te",
          "Te_scale", "tau_E1", "theta_Tg")
NO_KEYS = ("pre_drive", "ca_sh_rest", "ca_scale", "theta_NO", "i_scale", "theta_NOi", "tau_NO", "tau_Z", "theta_Z",
           "no_drive", "theta_NOe", "e_scale", "theta_NOc", "c_scale", "theta_N")


def drive(sig, P, which):
    """u_T (which='T') or u_N (which='NO'), see the module docstring. sig: shaft_cai or -ica_VDCC, or
    effcai for which='NO' with no_drive 1 (v2.3: NO needs NMDAR-driven spine Ca, nNOS sits on PSD-95/NMDAR),
    or spine Ca cai_CR - min_ca (cacr, recovered from effcai) with no_drive 2 (v2.4)."""
    x = sig.astype(np.float64)
    if which == "T" and P["t_drive"] == 1:
        return np.maximum(x - P["theta_Te"], 0.0) / P["Te_scale"]
    if which == "NO" and P["no_drive"] == 1:
        return np.maximum(x - P["theta_NOe"], 0.0) / P["e_scale"]
    if which == "NO" and P["no_drive"] == 2:
        return np.maximum(x - P["theta_NOc"], 0.0) / P["c_scale"]
    if P["pre_drive"]:
        return np.maximum(x - P["theta_" + which + "i"], 0.0) / P["i_scale"]
    return np.maximum(x - P["ca_sh_rest"] - P["theta_" + which], 0.0) / P["ca_scale"]


def feature_T(sig, t_ms, arr, P):
    """tanh(T) at each arrival (n_syn, n_spk). Depends only on T_KEYS."""
    P = {**DEFAULTS, **P}
    if P["t_drive"] in (2, 3, 4):
        S = _impulse_grid(np.asarray(sig, np.float64), t_ms, P["tau_E1"])
        u = np.maximum(S - P["theta_Te"], 0.0) / P["Te_scale"]
    else:
        u = drive(sig, P, "T")
    T = _lowpass_grid(u, t_ms, P["tau_T"])
    return np.tanh(np.maximum(T[np.arange(len(sig))[:, None], arr] - P["theta_Tg"], 0.0))


def feature_K(sig, t_ms, arr, P):
    """Integral (ms) of tanh(N) tanh(Z) before each arrival plus the tail after the last (n_syn, n_spk + 1).
    Depends only on NO_KEYS."""
    P = {**DEFAULTS, **P}
    n_syn, n_t = sig.shape
    rows = np.arange(n_syn)[:, None]
    N = _lowpass_grid(drive(sig, P, "NO"), t_ms, P["tau_NO"])
    imp = np.zeros_like(N)
    np.add.at(imp, (np.broadcast_to(rows, arr.shape), arr), 1.0)
    Z = _impulse_grid(imp, t_ms, P["tau_Z"])
    g = np.tanh(np.maximum(N - P["theta_N"], 0.0)) * np.tanh(np.maximum(Z - P["theta_Z"], 0.0)); del N, Z, imp
    h = np.diff(t_ms)
    cum = np.concatenate([np.zeros((n_syn, 1)), np.cumsum(g[:, :-1] * h, axis=1)], axis=1)   # left point
    idx = np.concatenate([np.zeros((n_syn, 1), np.int64), arr, np.full((n_syn, 1), n_t - 1)], axis=1)
    return np.diff(np.take_along_axis(cum, idx, axis=1), axis=1)


def pre_features(sig, t_ms, arr, P):
    """-> (tT, K), see feature_T and feature_K. arr is from arrival_index."""
    return feature_T(sig, t_ms, arr, P), feature_K(sig, t_ms, arr, P)


def dpre_final(tT, K, A_mglu, A_NO, dpre_min=-0.8, dpre_max=1.0, dpre0=0.0, tau_ind=TAU_IND):
    """Exact event-driven recurrence; vectorised over synapses."""
    n_syn, n_spk = tT.shape
    d = np.full(n_syn, float(dpre0))
    r = A_NO / (1e3 * tau_ind)
    for j in range(n_spk):
        d = dpre_max - (dpre_max - d) * np.exp(-r * K[:, j])
        d = d - A_mglu * tT[:, j] * (d - dpre_min)
    return dpre_max - (dpre_max - d) * np.exp(-r * K[:, n_spk])


# ---------------------------------------------------------------- EPSP readout, Option B
def use_target(rho_bin, dpre, Use_d, Use_p):
    return np.minimum(1.0, (Use_d + rho_bin * (Use_p - Use_d)) * (1.0 + dpre))


def epsp_v2(basis_df, rho_bin, dpre, ep):
    """Basis EPSP with each synapse's contribution scaled by its Use ratio (v2/v1).

    The depressed EPSP e0 is split over synapses in proportion to gmax_d*Use_d. Synapse i contributes
    e0_i + rho_i*Delta_i (as in v1), times Use_tgt(rho_i, dpre_i) / Use_tgt(rho_i, 0).
    With dpre = 0 this is exactly model.epsp_from_basis (same mean, same variance)."""
    n = len(rho_bin)
    cfg = dict(zip(basis_df["config"], zip(basis_df["mean"], basis_df["std"])))
    e0, s0 = cfg[",".join(["0"] * n)]
    w = ep["gmax_d"] * ep["Use_d"]; w = w / w.sum()
    delta = np.zeros(n); svar = np.zeros(n)
    for i in range(n):
        s = ["0"] * n; s[i] = "1"
        m, sd = cfg[",".join(s)]
        delta[i], svar[i] = m - e0, sd ** 2
    rb = (np.asarray(rho_bin) >= 0.5).astype(float)
    scale = use_target(rb, dpre, ep["Use_d"], ep["Use_p"]) / use_target(rb, 0.0, ep["Use_d"], ep["Use_p"])
    mean = float(np.sum((e0 * w + rb * delta) * scale))
    k = rb.sum()
    var = float(np.sum(rb * svar) + (1 - k) ** 2 * s0 ** 2)
    return mean, np.sqrt(max(var, 0.0))


def epsp_ratio_v2(basis_df, rho0, rho_f, dpre_f, ep, dpre0=0.0):
    """Same definition as model.epsp_ratio (cv2 term included)."""
    b_m, b_s = epsp_v2(basis_df, rho0, np.full(len(rho0), dpre0), ep)
    a_m, _ = epsp_v2(basis_df, rho_f, dpre_f, ep)
    if b_m == 0:
        return None
    cv2 = min((b_s / b_m) ** 2, 0.25)
    return (a_m / b_m) * (1.0 + cv2)
