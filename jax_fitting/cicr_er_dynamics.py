#!/usr/bin/env python3
"""
CICR Model V5 — ER Store Dynamics with IP3R/RyR Competition

Equations:
  Ca_e     = max(0, cai_raw - cai_rest)
  dIP3/dt  = delta_IP3 * Ca_e/(K_PLC + Ca_e) - IP3/tau_IP3
  dP/dt    = phi_serca * Hill(Ca_e) * (1-P) - J_RyR - J_IP3R - k_leak*(P - P_ss)
  dCa_RyR/dt = eta_RyR * J_RyR - Ca_RyR/tau_RyR
  drho/dt  = GB(effcai, gamma_p_eff, gamma_d)

  J_RyR  = V_RyR * P * Ca_e^n_T / (K_T^n_T + Ca_e^n_T)
  J_IP3R = V_IP3R * P * IP3/(K_IP3+IP3) * Ca_e/(K_act+Ca_e)
  Hill(x) = x^n_H / (K_H^n_H + x^n_H)
  gamma_p_eff = gamma_p + g_RyR * Ca_RyR

Fixed constants: n_H=4, K_H=0.0004, K_act=0.0003, K_IP3=1.0, K_PLC=0.0005,
                 P_ss=0.3, tau_RyR=100, cai_rest=70e-6
"""

import numpy as np
import jax.numpy as jnp
from cicr_common import CICRModel, compute_effcai_piecewise_linear_jax

N_H = 4
K_H = 0.0004
K_ACT = 0.0003
K_IP3 = 1.0
K_PLC = 0.0005
P_SS = 0.3
TAU_RYR = 100.0
CAI_REST = 70e-6
TRACE_TAU_REF = 200.0


def _hill(x, K, n):
    xn = x**n
    return xn / (K**n + xn + 1e-30)


def _v5_debug(cai_1syn, t, cp, cq, is_apical, dp):
    delta_IP3 = float(dp["delta_IP3"])
    tau_IP3   = float(dp["tau_IP3"])
    phi_serca = float(dp["phi_serca"])
    k_leak    = float(dp["k_leak"])
    V_RyR     = float(dp["V_RyR"])
    K_T       = float(dp["K_T"])
    n_T       = float(dp["n_T"])
    V_IP3R    = float(dp["V_IP3R"])
    g_RyR     = float(dp["g_RyR"])
    tau_eff   = float(dp["tau_eff"])
    gamma_d   = float(dp.get("gamma_d_GB_GluSynapse", 150.0))
    gamma_p   = float(dp.get("gamma_p_GB_GluSynapse", 150.0))

    scale = tau_eff / TRACE_TAU_REF
    cp_s, cq_s = cp * scale, cq * scale
    if is_apical:
        theta_d = dp["a20"] * cp_s + dp["a21"] * cq_s
        theta_p = dp["a30"] * cp_s + dp["a31"] * cq_s
    else:
        theta_d = dp["a00"] * cp_s + dp["a01"] * cq_s
        theta_p = dp["a10"] * cp_s + dp["a11"] * cq_s
    theta_p = max(theta_p, theta_d + 0.01)

    n = len(cai_1syn)
    cai_total = np.copy(cai_1syn)
    ip3_out   = np.zeros(n)
    p_out     = np.zeros(n)
    ca_ryr_out = np.zeros(n)
    effcai_out = np.zeros(n)
    rho_out   = np.zeros(n)

    IP3, P, ca_ryr, effcai = 0.0, P_SS, 0.0, 0.0
    rho = float(dp.get("rho0", 0.5))

    ip3_out[0] = IP3; p_out[0] = P; ca_ryr_out[0] = ca_ryr
    effcai_out[0] = effcai; rho_out[0] = rho

    for i in range(n - 1):
        dt = t[i+1] - t[i]
        if dt <= 0:
            ip3_out[i+1] = IP3; p_out[i+1] = P; ca_ryr_out[i+1] = ca_ryr
            effcai_out[i+1] = effcai; rho_out[i+1] = rho
            cai_total[i+1] = cai_1syn[i+1] + ca_ryr
            continue

        ca_e = max(0.0, cai_1syn[i] - CAI_REST)

        dIP3 = delta_IP3 * ca_e / (K_PLC + ca_e + 1e-12) - IP3 / tau_IP3
        IP3 = max(0.0, IP3 + dt * dIP3)

        T = _hill(ca_e, K_T, n_T)
        S = _hill(ca_e, K_H, N_H)
        J_ryr  = V_RyR * P * T
        J_ip3r = V_IP3R * P * (IP3 / (K_IP3 + IP3)) * (ca_e / (K_ACT + ca_e + 1e-12))

        dP = phi_serca * S * (1.0 - P) - J_ryr - J_ip3r - k_leak * (P - P_SS)
        P = min(1.0, max(0.0, P + dt * dP))

        decay_ryr = np.exp(-dt / TAU_RYR)
        ca_ryr = ca_ryr * decay_ryr + J_ryr * TAU_RYR * (1.0 - decay_ryr)

        eff_decay = np.exp(-dt / tau_eff)
        ca_total_i = ca_e + ca_ryr
        effcai = effcai * eff_decay + ca_total_i * tau_eff * (1.0 - eff_decay)

        gp_eff = gamma_p + g_RyR * ca_ryr
        pot = 1.0 if effcai > theta_p else 0.0
        dep = 1.0 if effcai > theta_d else 0.0
        drho = (-rho*(1-rho)*(0.5-rho) + gp_eff*(1-rho)*pot - gamma_d*rho*dep) / 70000.0
        rho = min(1.0, max(0.0, rho + dt * drho))

        cai_total[i+1]  = cai_1syn[i+1] + ca_ryr
        ip3_out[i+1]    = IP3
        p_out[i+1]      = P
        ca_ryr_out[i+1] = ca_ryr
        effcai_out[i+1] = effcai
        rho_out[i+1]    = rho

    return {
        "cai_total": cai_total,
        "priming":   ip3_out,
        "ca_er":     p_out,
        "ca_cicr":   ca_ryr_out,
        "effcai":    effcai_out,
        "rho":       rho_out,
    }


class ERDynamicsCICRModel(CICRModel):
    DESCRIPTION = "V5 ER Dynamics (IP3R/RyR Competition)"

    FIT_PARAMS = [
        ("gamma_d_GB_GluSynapse", 50.0, 250.0),
        ("gamma_p_GB_GluSynapse", 50.0, 300.0),
        ("a00", 1.0, 2.5),
        ("a01", 1.0, 5.0),
        ("a10", 1.0, 5.0),
        ("a11", 1.0, 5.0),
        ("a20", 1.0, 2.5),
        ("a21", 1.0, 5.0),
        ("a30", 1.0, 10.0),
        ("a31", 1.0, 5.0),
        ("delta_IP3", 0.01, 5.0),
        ("tau_IP3", 500.0, 5000.0),
        ("phi_serca", 0.0001, 0.05),
        ("k_leak", 1e-6, 0.001),
        ("V_RyR", 0.0001, 0.05),
        ("K_T", 0.001, 0.005),
        ("n_T", 4.0, 15.0),
        ("V_IP3R", 0.0001, 0.01),
        ("g_RyR", 10.0, 50000.0),
        ("tau_eff", 50.0, 500.0),
    ]

    DEFAULT_PARAMS = {
        "gamma_d_GB_GluSynapse": 69.61242611028226, "gamma_p_GB_GluSynapse": 166.96578019613747,
        "a00": 1.0003490121729528, "a01": 1.0877438002708641, "a10": 1.3917207001380245, "a11": 4.129505808304705,
        "a20": 2.3079274077684175, "a21": 1.9619384872398133, "a30": 4.456560780131431, "a31": 4.80485148600644,
        "delta_IP3": 2.3564846028724626, "tau_IP3": 4198.8339980611045,
        "phi_serca": 0.00010926586342370828, "k_leak": 0.00038983484712417705,
        "V_RyR": 0.0008831061403601626, "K_T": 0.0016997377463884833, "n_T": 9.677747356489755,
        "V_IP3R": 0.009883023779226525, "g_RyR": 14060.428640146583, "tau_eff": 306.6936809019515,
    }

    def unpack_params(self, x):
        gb = {
            'gamma_d': x[0], 'gamma_p': x[1],
            'a00': x[2], 'a01': x[3], 'a10': x[4], 'a11': x[5],
            'a20': x[6], 'a21': x[7], 'a30': x[8], 'a31': x[9],
        }
        cicr = {
            'delta_IP3': x[10], 'tau_IP3': x[11],
            'phi_serca': x[12], 'k_leak': x[13],
            'V_RyR': x[14], 'K_T': x[15], 'n_T': x[16],
            'V_IP3R': x[17], 'g_RyR': x[18], 'tau_eff': x[19],
        }
        return gb, cicr

    def get_init_fn(self):
        def init_fn(cai_first, rho0):
            return (0.0, P_SS, 0.0, 0.0, rho0)
        return init_fn

    def get_debug_sim_fn(self):
        return _v5_debug

    def get_step_factory(self):
        def scan_factory(params, syn_params):
            gb, cicr = params
            c_pre_opt, c_post_opt, is_apical, cai_pre, t_pre, cai_post, t_post = syn_params

            c_pre = jnp.max(compute_effcai_piecewise_linear_jax(cai_pre, t_pre, tau_effca=cicr['tau_eff']))
            c_post = jnp.max(compute_effcai_piecewise_linear_jax(cai_post, t_post, tau_effca=cicr['tau_eff']))

            theta_d = jnp.where(is_apical,
                                gb['a20']*c_pre + gb['a21']*c_post,
                                gb['a00']*c_pre + gb['a01']*c_post)
            raw_theta_p = jnp.where(is_apical,
                                    gb['a30']*c_pre + gb['a31']*c_post,
                                    gb['a10']*c_pre + gb['a11']*c_post)
            theta_p = jnp.maximum(raw_theta_p, theta_d + 0.01)

            n_H = N_H
            K_H_c = K_H
            K_act_c = K_ACT
            K_IP3_c = K_IP3
            K_PLC_c = K_PLC
            P_ss_c = P_SS
            tau_ryr_c = TAU_RYR
            cai_rest = CAI_REST

            def scan_step(carry, inputs):
                IP3, P, ca_ryr, effcai, rho = carry
                cai_raw, dt = inputs

                ca_e = jnp.maximum(0.0, cai_raw - cai_rest)

                dIP3 = cicr['delta_IP3'] * ca_e / (K_PLC_c + ca_e + 1e-12) - IP3 / cicr['tau_IP3']
                IP3_new = jnp.maximum(0.0, IP3 + dt * dIP3)

                ca_e_nT = ca_e ** cicr['n_T']
                K_T_nT  = cicr['K_T'] ** cicr['n_T']
                T = ca_e_nT / (K_T_nT + ca_e_nT + 1e-30)

                ca_e_nH = ca_e ** n_H
                K_H_nH  = K_H_c ** n_H
                S = ca_e_nH / (K_H_nH + ca_e_nH + 1e-30)

                J_ryr  = cicr['V_RyR'] * P * T
                J_ip3r = cicr['V_IP3R'] * P * (IP3_new / (K_IP3_c + IP3_new + 1e-12)) * (ca_e / (K_act_c + ca_e + 1e-12))

                dP = cicr['phi_serca'] * S * (1.0 - P) - J_ryr - J_ip3r - cicr['k_leak'] * (P - P_ss_c)
                P_new = jnp.clip(P + dt * dP, 0.0, 1.0)

                decay_ryr = jnp.exp(-dt / tau_ryr_c)
                ca_ryr_new = ca_ryr * decay_ryr + J_ryr * tau_ryr_c * (1.0 - decay_ryr)

                eff_decay = jnp.exp(-dt / cicr['tau_eff'])
                ca_total = ca_e + ca_ryr_new
                effcai_new = effcai * eff_decay + ca_total * cicr['tau_eff'] * (1.0 - eff_decay)

                gp_eff = gb['gamma_p'] + cicr['g_RyR'] * ca_ryr_new
                pot = jnp.where(effcai_new > theta_p, 1.0, 0.0)
                dep = jnp.where(effcai_new > theta_d, 1.0, 0.0)

                drho = (-rho*(1.0-rho)*(0.5-rho)
                        + gp_eff*(1.0-rho)*pot
                        - gb['gamma_d']*rho*dep) / 70000.0
                rho_new = jnp.where(dt > 0, jnp.clip(rho + dt * drho, 0.0, 1.0), rho)

                return (IP3_new, P_new, ca_ryr_new, effcai_new, rho_new), None

            return scan_step
        return scan_factory


if __name__ == "__main__":
    ERDynamicsCICRModel().run()