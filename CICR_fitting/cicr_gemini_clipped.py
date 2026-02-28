#!/usr/bin/env python3
"""
Clipped-Latent ER-CICR Plasticity Model.

Grounded in Neymotin et al. (2015) ER calcium wave dynamics.
1. SERCA loading of ER stores is soft-clipped so narrow high-amplitude
   transients (bAPs) cannot instantly fill the ER — only broad, moderate
   events (EPSPs) accumulate significant ER Ca²⁺.
2. IP3R activation gate (steep Hill function) requires high cytosolic Ca²⁺
   to open, so only bAP-scale events trigger CICR release from primed ER.
3. Total cytosolic Ca²⁺ (direct + CICR) is low-pass filtered to drive
   Graupner-Brunel bistable plasticity (Chindemi et al. 2022).
"""

import numpy as np
import jax.numpy as jnp
from cicr_common import CICRModel

def _apply_cicr_debug(cai_1syn, t, cp, cq, is_apical, dp):
    tau_er_leak = float(dp["tau_er_leak"])
    k_serca    = float(dp["k_serca"])
    K_clip     = float(dp["K_clip"])
    K_ip3r     = float(dp["K_ip3r"])
    n_H        = float(dp["n_H"])
    g_ip3r     = float(dp["g_ip3r"])
    tau_cicr   = float(dp["tau_cicr"])

    n = len(cai_1syn)
    cai_rest = 70e-6

    cai_total   = np.copy(cai_1syn)
    ca_er_out   = np.zeros(n)
    ca_cicr_out = np.zeros(n)

    ca_er   = 0.0
    ca_cicr = 0.0

    ca_er_out[0]   = ca_er
    ca_cicr_out[0] = ca_cicr

    for i in range(n - 1):
        dt = t[i + 1] - t[i]
        if dt <= 0.0:
            ca_er_out[i + 1]   = ca_er
            ca_cicr_out[i + 1] = ca_cicr
            continue

        cai_evoked = max(0.0, cai_1syn[i] - cai_rest)

        cai_clip = K_clip * np.tanh(cai_evoked / (K_clip + 1e-12))

        rate_serca = k_serca * cai_clip
        inv_tau_leak = 1.0 / tau_er_leak
        tau_eff_er = 1.0 / (rate_serca + inv_tau_leak + 1e-12)
        er_inf = rate_serca * tau_eff_er
        ca_er = er_inf + (ca_er - er_inf) * np.exp(-dt / tau_eff_er)

        ip3r_gate = (cai_evoked**n_H) / (K_ip3r**n_H + cai_evoked**n_H + 1e-12)
        J_ip3r = g_ip3r * ca_er * ip3r_gate

        decay_cicr = np.exp(-dt / tau_cicr)
        ca_cicr = ca_cicr * decay_cicr + J_ip3r * tau_cicr * (1.0 - decay_cicr)

        cai_total[i + 1]   = cai_1syn[i + 1] + ca_cicr
        ca_er_out[i + 1]   = ca_er
        ca_cicr_out[i + 1] = ca_cicr

    return {
        "cai_total": cai_total,
        "priming":   ca_er_out,
        "ca_er":     ca_er_out,
        "ca_cicr":   ca_cicr_out,
    }


class ClippedLatentCICRModel(CICRModel):
    DESCRIPTION = "Clipped-Latent ER-CICR (JAX)"

    FIT_PARAMS = [
        ("gamma_d", 50.0, 300.0),
        ("gamma_p", 150.0, 600.0),
        ("a00", 0.01, 5.0),
        ("a01", 0.01, 5.0),
        ("a10", 0.01, 5.0),
        ("a11", 0.01, 5.0),
        ("a20", 0.01, 5.0),
        ("a21", 0.01, 5.0),
        ("a30", 0.01, 5.0),
        ("a31", 0.01, 5.0),
        ("tau_er_leak", 100.0, 2000.0),
        ("k_serca", 100.0, 50000.0),
        ("K_clip", 0.0001, 0.001),
        ("K_ip3r", 0.001, 0.005),
        ("n_H", 2.0, 12.0),
        ("g_ip3r", 0.00001, 0.01),
        ("tau_cicr", 10.0, 500.0),
        ("tau_eff", 50.0, 500.0),
    ]

    DEFAULT_PARAMS = {
        "gamma_d": 88.7017, "gamma_p": 268.3207,
        "a00": 0.4778, "a01": 0.6934, "a10": 1.0587, "a11": 3.4651,
        "a20": 0.7323, "a21": 0.1475, "a30": 1.8814, "a31": 4.4472,
        "tau_er_leak": 1271.4256, "k_serca": 27269.7497, "K_clip": 0.0008,
        "K_ip3r": 0.0037, "n_H": 3.4732,
        "g_ip3r": 0.0010, "tau_cicr": 192.5917,
        "tau_eff": 111.4966,
    }

    def unpack_params(self, x):
        gb = {'gamma_d': x[0], 'gamma_p': x[1],
              'a00': x[2], 'a01': x[3], 'a10': x[4], 'a11': x[5],
              'a20': x[6], 'a21': x[7], 'a30': x[8], 'a31': x[9]}
        cicr = {'tau_er_leak': x[10], 'k_serca': x[11], 'K_clip': x[12],
                'K_ip3r': x[13], 'n_H': x[14], 'g_ip3r': x[15],
                'tau_cicr': x[16], 'tau_eff': x[17]}
        return gb, cicr

    def get_init_fn(self):
        def init_fn(cai_first, rho0):
            return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, rho0)
        return init_fn

    def get_debug_sim_fn(self):
        return _apply_cicr_debug

    def get_step_factory(self):
        def scan_factory(params, syn_params):
            gb, cicr = params
            c_pre, c_post, is_apical = syn_params

            theta_d = jnp.where(is_apical,
                                gb['a20']*c_pre + gb['a21']*c_post,
                                gb['a00']*c_pre + gb['a01']*c_post)
            raw_theta_p = jnp.where(is_apical,
                                    gb['a30']*c_pre + gb['a31']*c_post,
                                    gb['a10']*c_pre + gb['a11']*c_post)
            theta_p = jnp.maximum(raw_theta_p, theta_d + 0.01)

            cai_rest = 70e-6

            def scan_step(carry, inputs):
                ca_er, ca_cicr, _d1, _d2, _d3, effcai, rho = carry
                cai_raw, dt = inputs

                cai_evoked = jnp.maximum(0.0, cai_raw - cai_rest)

                cai_clip = cicr['K_clip'] * jnp.tanh(cai_evoked / (cicr['K_clip'] + 1e-12))

                rate_serca = cicr['k_serca'] * cai_clip
                inv_tau_leak = 1.0 / cicr['tau_er_leak']
                tau_eff_er = 1.0 / (rate_serca + inv_tau_leak + 1e-12)
                er_inf = rate_serca * tau_eff_er
                ca_er_new = er_inf + (ca_er - er_inf) * jnp.exp(-dt / tau_eff_er)

                ip3r_gate = (cai_evoked**cicr['n_H']) / (cicr['K_ip3r']**cicr['n_H'] + cai_evoked**cicr['n_H'] + 1e-12)
                J_ip3r = cicr['g_ip3r'] * ca_er_new * ip3r_gate

                ca_cicr_decay = jnp.exp(-dt / cicr['tau_cicr'])
                ca_cicr_new = ca_cicr * ca_cicr_decay + J_ip3r * cicr['tau_cicr'] * (1.0 - ca_cicr_decay)

                ca_total = cai_evoked + ca_cicr_new
                eff_decay = jnp.exp(-dt / cicr['tau_eff'])
                effcai_new = effcai * eff_decay + ca_total * cicr['tau_eff'] * (1.0 - eff_decay)

                pot = jnp.where(effcai_new > theta_p, 1.0, 0.0)
                dep = jnp.where(effcai_new > theta_d, 1.0, 0.0)

                drho = (-rho*(1.0-rho)*(0.5-rho)
                        + pot*gb['gamma_p']*(1.0-rho)
                        - dep*gb['gamma_d']*rho) / 70000.0

                rho_new = jnp.clip(rho + dt * drho, 0.0, 1.0)

                return (ca_er_new, ca_cicr_new, 0.0, 0.0, 0.0, effcai_new, rho_new), None
            return scan_step
        return scan_factory


if __name__ == "__main__":
    ClippedLatentCICRModel().run()