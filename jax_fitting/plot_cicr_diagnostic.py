#!/usr/bin/env python3
"""
Diagnostic plotting utilities for CICR parameter fitting.
Automatically scales abstract state variables vs. physical Ca2+ concentrations.
"""

import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

try:
    from numba import njit as _njit
except ImportError:
    def _njit(*a, **kw):
        def _wrap(fn):
            return fn
        return _wrap

# ══════════════════════════════════════════════════════════════════
# Numba-compiled single-synapse ODE helpers
# ══════════════════════════════════════════════════════════════════

@_njit(cache=True)
def _compute_effcai(cai, t, tau_effca, min_ca):
    """Leaky-integrator effcai for one synapse trace."""
    n = len(cai)
    effcai = np.zeros(n)
    e = 0.0
    for i in range(n - 1):
        dt = t[i + 1] - t[i]
        decay = np.exp(-dt / tau_effca)
        e = e * decay + (cai[i] - min_ca) * tau_effca * (1.0 - decay)
        if e < 0.0:
            e = 0.0
        effcai[i + 1] = e
    return effcai

@_njit(cache=True)
def _compute_rho(effcai, t, theta_d, theta_p, gamma_d, gamma_p, rho0):
    """Euler-integrate the Graupner-Brunel rho ODE, return full trace."""
    n = len(t)
    rho = np.zeros(n)
    rho[0] = rho0
    inv_tau = 1.0 / (70.0 * 1000.0)
    rho_star = 0.5
    r = rho0
    for i in range(n - 1):
        dt = t[i + 1] - t[i]
        
        # CaMKII VETO LOGIC:
        pot = 1.0 if effcai[i] > theta_p else 0.0
        dep = 1.0 if (effcai[i] > theta_d and pot == 0.0) else 0.0
        
        drho = (
            -r * (1.0 - r) * (rho_star - r)
            + pot * gamma_p * (1.0 - r)
            - dep * gamma_d * r
        ) * inv_tau
        r = min(1.0, max(0.0, r + dt * drho))
        rho[i + 1] = r
    return rho

@_njit(cache=True)
def _apply_primed_cicr_debug(
    cai_1syn, t, c_pre, c_post, tau_IP3, IP3_threshold,
    a_evt, b_evt, a_er, b_er, g_serca, g_release,
    sat_cap, sigmoid_slope, K_act,
):
    """Fallback Single-synapse Primed CICR"""
    K_serca = 0.1
    fc, fe = 0.83, 0.17
    caAvg = 0.0017
    MAX_DT_SUB = 0.1

    n = len(cai_1syn)
    cai_total = np.copy(cai_1syn)
    priming = np.zeros(n)
    ca_er_out = np.zeros(n)
    ca_cicr_out = np.zeros(n)

    evt_thresh = a_evt * c_pre + b_evt * c_post
    gamma_er   = a_er  * c_pre + b_er  * c_post

    ca_er = max(0.0, (caAvg - fc * cai_1syn[0]) / fe)
    P_val, ca_cicr, was_above, evt_peak = 0.0, 0.0, 0, 0.0

    priming[0], ca_er_out[0], ca_cicr_out[0] = P_val, ca_er, ca_cicr

    for i in range(n - 1):
        dt_full = t[i + 1] - t[i]
        if dt_full <= 0.0:
            cai_total[i + 1], priming[i + 1], ca_er_out[i + 1], ca_cicr_out[i + 1] = cai_1syn[i + 1] + ca_cicr, P_val, ca_er, ca_cicr
            continue

        is_above  = 1 if cai_1syn[i] > evt_thresh else 0
        evt_amp   = 0.0
        if is_above == 1 and was_above == 0:
            evt_peak = cai_1syn[i]
        elif is_above == 1:
            if cai_1syn[i] > evt_peak: evt_peak = cai_1syn[i]
        elif was_above == 1:
            raw = evt_peak - evt_thresh
            evt_amp  = (raw / sat_cap) if raw < sat_cap else 1.0
            evt_peak = 0.0
        was_above = is_above

        P_val = P_val * np.exp(-dt_full / tau_IP3) + evt_amp
        arg = -sigmoid_slope * (P_val - IP3_threshold)
        P_gate = 0.0 if arg > 500.0 else (1.0 if arg < -500.0 else 1.0 / (1.0 + np.exp(arg)))

        n_sub = max(1, int(np.ceil(dt_full / MAX_DT_SUB)))
        dt_sub = dt_full / n_sub
        for _ in range(n_sub):
            ca_cyt = cai_1syn[i] + ca_cicr
            ca_cicr_um = 1000.0 * ca_cicr
            J_serca = min(g_serca * ca_cicr_um**2 / (K_serca**2 + ca_cicr_um**2), ca_cicr / dt_sub) if ca_cicr > gamma_er else 0.0
            J_rel = g_release * P_gate * (ca_cyt / (ca_cyt + K_act)) * (ca_er - ca_cyt) if ca_er > ca_cyt else 0.0
            J_net = J_rel - J_serca
            ca_er = max(0.0, ca_er - (fc / fe) * J_net * dt_sub)
            ca_cicr = max(0.0, ca_cicr + J_net * dt_sub)

        cai_total[i + 1], priming[i + 1], ca_er_out[i + 1], ca_cicr_out[i + 1] = cai_1syn[i + 1] + ca_cicr, P_val, ca_er, ca_cicr

    return cai_total, priming, ca_er_out, ca_cicr_out

# ══════════════════════════════════════════════════════════════════
# Parameter loading & Trace Computation
# ══════════════════════════════════════════════════════════════════

def params_from_json(json_path, param_names=None, default_params=None):
    with open(json_path) as f: data = json.load(f)
    dp = dict(default_params) if default_params else {}
    if "best_parameters" in data:
        for name, entry in data["best_parameters"].items():
            dp[name] = entry["value"] if isinstance(entry, dict) else float(entry)
        return dp
    raise ValueError(f"No best_parameters found in {json_path}")

def pick_random_pair_syn(protocol_data, protocols, seed=None):
    rng = np.random.default_rng(seed)
    ref_pairs = protocol_data[protocols[0]]
    valid = [i for i in range(len(ref_pairs)) if all(i < len(protocol_data[p]) for p in protocols)]
    pair_idx = int(rng.choice(valid))
    peaks = ref_pairs[pair_idx]["cai"].max(axis=1)
    good = np.where(peaks > 0.001)[0]
    syn_idx = int(rng.choice(good)) if len(good) > 0 else int(rng.choice(np.arange(len(peaks))))
    return pair_idx, syn_idx

def compute_single(pd_item, syn_idx, dp, debug_sim_fn=None):
    s = syn_idx
    t = pd_item["t"]
    cp, cq, is_apical = pd_item["c_pre"][s], pd_item["c_post"][s], pd_item["is_apical"][s]
    gamma_d = float(dp.get("gamma_d", dp.get("gamma_d_GB_GluSynapse", 150.0)))
    gamma_p = float(dp.get("gamma_p", dp.get("gamma_p_GB_GluSynapse", 200.0)))

    # Rescale c_pre/c_post from recording tau to the model's current tau_eff
    TRACE_TAU_EFFCA = 200.0
    tau_eff = float(dp.get("tau_eff", TRACE_TAU_EFFCA))
    scale = tau_eff / TRACE_TAU_EFFCA
    cp_s, cq_s = cp * scale, cq * scale

    if is_apical:
        theta_d = dp["a20"] * cp_s + dp["a21"] * cq_s
        theta_p = dp["a30"] * cp_s + dp["a31"] * cq_s
    else:
        theta_d = dp["a00"] * cp_s + dp["a01"] * cq_s
        theta_p = dp["a10"] * cp_s + dp["a11"] * cq_s

    cai_raw = pd_item["cai"][s].copy()

    if debug_sim_fn is not None:
        dp = dict(dp)
        dp["rho0"] = pd_item["rho0"][s]
        sim_out = debug_sim_fn(cai_raw, t, cp, cq, is_apical, dp)
        cai_total, priming, ca_er, ca_cicr = sim_out["cai_total"], sim_out["priming"], sim_out["ca_er"], sim_out["ca_cicr"]
    else:
        # Fallback to standard Primed-CICR
        cai_total, priming, ca_er, ca_cicr = _apply_primed_cicr_debug(
            cai_raw, t, cp, cq, float(dp.get("tau_IP3", 500.0)), float(dp.get("IP3_threshold", 3.0)),
            float(dp.get("a_evt", 0.003)), float(dp.get("b_evt", 0.003)), float(dp.get("a_er", 0.005)), float(dp.get("b_er", 0.005)),
            float(dp.get("g_serca", 1.9565)), float(dp.get("g_release", 0.5)), float(dp.get("sat_cap", 0.05)),
            float(dp.get("sigmoid_slope", 5.0)), float(dp.get("K_act", 0.0005))
        )
        sim_out = {}

    TAU_EFFCA = float(dp.get("tau_eff", 200.0))
    MIN_CA = 70e-6

    effcai_no = _compute_effcai(cai_raw, t, TAU_EFFCA, MIN_CA)
    effcai_ci = sim_out.get("effcai", _compute_effcai(cai_total, t, TAU_EFFCA, MIN_CA))

    rho_no = _compute_rho(effcai_no, t, theta_d, theta_p, gamma_d, gamma_p, pd_item["rho0"][s])
    rho_ci = sim_out.get("rho", _compute_rho(effcai_ci, t, theta_d, theta_p, gamma_d, gamma_p, pd_item["rho0"][s]))

    return {
        "t": t, "cai_raw": cai_raw, "cai_total": cai_total, 
        "effcai_no": effcai_no, "effcai_ci": effcai_ci,
        "rho_no": rho_no, "rho_ci": rho_ci,
        "priming": priming, "ca_er": ca_er, "ca_cicr": ca_cicr,
        "theta_d": theta_d, "theta_p": theta_p, "cp": cp, "cq": cq,
    }

# ══════════════════════════════════════════════════════════════════
# Plotting
# ══════════════════════════════════════════════════════════════════

_PROTO_COLORS = [
    dict(line='#d62728', mGluR='#9467bd', er='#17becf', cicr_ca='#d62728'),
    dict(line='#e377c2', mGluR='#8c564b', er='#bcbd22', cicr_ca='#e377c2'),
    dict(line='#8c564b', mGluR='#7f7f7f', er='#aec7e8', cicr_ca='#8c564b'),
]
_SHARED = dict(thresh='#937860', theta_d='#DD8452', theta_p='#55A868', cicr_ref='#C44E52')

def _xlim_from_cai(t_s, cai_raw):
    above = cai_raw > cai_raw[0] * 2.0
    return float(min(t_s[np.where(above)[0][-1]] * 1.5, t_s[-1])) if np.any(above) else float(t_s[-1])

def plot_diagnostic(results, protocols=None, pair_idx=0, syn_idx=0, param_label="", output="cicr_diagnostic.png"):
    if protocols is None: protocols = list(results.keys())
    protocols = [p for p in protocols if p in results]
    if not protocols: return

    plt.rcParams.update({'font.size': 14, 'axes.titlesize': 16, 'axes.labelsize': 14, 'figure.dpi': 150})
    fig, axes = plt.subplots(5, 1, figsize=(10, 16), sharex=True)
    linestyles = ['-', '--', ':']

    t_s0 = results[protocols[0]]["t"] / 1000.0
    xlim = _xlim_from_cai(t_s0, results[protocols[0]]["cai_raw"])

    has_abstract_er = any(np.max(results[p]["ca_er"]) > 0.1 for p in protocols)
    ax4_twin = axes[4].twinx() if has_abstract_er else None

    for i, proto in enumerate(protocols):
        d = results[proto]
        t_s = d["t"] / 1000.0
        ls, pc, sfx = linestyles[i % len(linestyles)], _PROTO_COLORS[i % len(_PROTO_COLORS)], f" ({proto})"
        
        axes[0].plot(t_s, d["cai_total"] * 1000, color=pc['line'], ls=ls, alpha=0.85, label=proto)
        axes[1].plot(t_s, d["effcai_ci"], color=pc['line'], ls=ls, alpha=0.85, label=proto)
        axes[2].plot(t_s, d["rho_ci"], color=pc['line'], ls=ls, alpha=0.9, linewidth=1.5, label=proto)
        axes[3].plot(t_s, d["priming"], color=pc['mGluR'], ls=ls, alpha=0.9, label='Latent Trace / Priming' + sfx)
        
        if has_abstract_er:
            ax4_twin.plot(t_s, d["ca_er"], color=pc['er'], ls=ls, alpha=0.9, label='ER State' + sfx)
        else:
            axes[4].plot(t_s, d["ca_er"] * 1000.0, color=pc['er'], ls=ls, alpha=0.9, label=r'ca$_\mathrm{ER}$' + sfx)
        
        axes[4].plot(t_s, d["ca_cicr"] * 1000, color=pc['cicr_ca'], ls=ls, alpha=0.9, label=r'ca$_\mathrm{CICR}$' + sfx)

        if i == 0:
            axes[0].set_title("CICR Protocol Comparison", fontsize=18, fontweight='bold')
            axes[0].set_ylabel(r'Ca$^{2+}$ (µM)')
            axes[1].axhline(d["theta_d"], color=_SHARED['theta_d'], ls=':', label=fr'$\theta_d$ = {d["theta_d"]:.3f}')
            axes[1].axhline(d["theta_p"], color=_SHARED['theta_p'], ls=':', label=fr'$\theta_p$ = {d["theta_p"]:.3f}')
            axes[1].set_ylabel('effcai')
            axes[2].axhline(0.5, color='gray', ls=':', alpha=0.5, label=r'$\rho^*$ = 0.5')
            axes[2].set_ylim(-0.05, 1.05); axes[2].set_ylabel(r'$\rho$')
            axes[3].set_ylabel('Eligibility / State')
            axes[4].set_ylabel(r'ca$_\mathrm{CICR}$ (µM)' if has_abstract_er else r'ER / CICR Ca$^{2+}$ (µM)')
            axes[4].set_xlabel('Time (s)', fontsize=14)
            if has_abstract_er:
                ax4_twin.set_ylabel('ER State', color=_PROTO_COLORS[0]['er'])

    for ax in axes[:4]:
        ax.set_xlim(0, xlim)
        ax.legend(loc='upper right', fontsize=8, ncol=2)

    axes[4].set_xlim(0, xlim)
    
    # Combine legends for axes[4] and ax4_twin
    lines, labels = axes[4].get_legend_handles_labels()
    if has_abstract_er:
        lines2, labels2 = ax4_twin.get_legend_handles_labels()
        axes[4].legend(lines + lines2, labels + labels2, loc='upper right', fontsize=8, ncol=2)
    else:
        axes[4].legend(loc='upper right', fontsize=8, ncol=2)

    plt.tight_layout(h_pad=0.8)
    fig.suptitle(f'Diagnostic — {param_label}\nPair #{pair_idx}, Synapse #{syn_idx}', fontsize=18, fontweight='bold', y=1.02)
    plt.savefig(output, dpi=200, bbox_inches='tight', facecolor='white')
    print(f"\nSaved → {output}")
    plt.close(fig)

def print_stats(label, d):
    pass # Overridden by cicr_common.py