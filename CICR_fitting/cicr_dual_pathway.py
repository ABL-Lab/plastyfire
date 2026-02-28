#!/usr/bin/env python3
"""
Release-First CICR with dynamic theta_p shift and tunable tau_eff.

Fixes over previous version:
  1. tau_eff is a fit parameter (was hardcoded at 200ms)
  2. Debug function computes REAL effcai and rho traces including
     the dynamic theta_p shift, so diagnostic plots are accurate.
"""

import numpy as np
import jax.numpy as jnp
from cicr_common import CICRModel


def _apply_release_first_debug(cai_1syn, t, cp, cq, is_apical, dp):
    k_e        = float(dp["k_e"])
    tau_ER     = float(dp["tau_ER"])
    A_ER       = float(dp["A_ER"])
    tau_CICR   = float(dp["tau_CICR"])
    cicr_shift = float(dp["cicr_shift"])
    K_shift    = float(dp["K_shift"])
    tau_eff    = float(dp["tau_eff"])
    gamma_d    = float(dp["gamma_d"])
    gamma_p    = float(dp["gamma_p"])

    if is_apical:
        theta_d = dp["a20"] * cp + dp["a21"] * cq
        raw_theta_p = dp["a30"] * cp + dp["a31"] * cq
    else:
        theta_d = dp["a00"] * cp + dp["a01"] * cq
        raw_theta_p = dp["a10"] * cp + dp["a11"] * cq
    theta_p_base = max(raw_theta_p, theta_d + 0.01)

    cai_rest = 70e-6
    K_trig = 1e-4
    n = len(cai_1syn)
    e_ER_out      = np.zeros(n)
    ca_cicr_out   = np.zeros(n)
    effcai_out    = np.zeros(n)
    rho_out       = np.zeros(n)
    theta_p_out   = np.zeros(n)

    e_ER    = 0.0
    ca_cicr = 0.0
    effcai  = 0.0
    rho     = 0.0

    theta_p_out[0] = theta_p_base

    for i in range(n - 1):
        dt = t[i + 1] - t[i]
        if dt <= 0.0:
            e_ER_out[i + 1]    = e_ER
            ca_cicr_out[i + 1] = ca_cicr
            effcai_out[i + 1]  = effcai
            rho_out[i + 1]     = rho
            theta_p_out[i + 1] = theta_p_base
            continue

        cai_raw    = cai_1syn[i]
        cai_evoked = max(0.0, cai_raw - cai_rest)

        e_decayed = e_ER * np.exp(-dt / tau_ER)
        trigger   = cai_evoked / (K_trig + cai_evoked + 1e-12)
        flux      = A_ER * e_decayed * trigger
        released  = min(flux * dt, max(0.0, e_decayed))
        e_ER      = max(0.0, e_decayed - released) + k_e * cai_evoked * dt
        ca_cicr   = ca_cicr * np.exp(-dt / tau_CICR) + released

        decay_eff = np.exp(-dt / tau_eff)
        effcai    = effcai * decay_eff + cai_evoked * tau_eff * (1.0 - decay_eff)

        shift      = cicr_shift * ca_cicr / (K_shift + ca_cicr + 1e-12)
        theta_p_eff = theta_p_base + shift

        pot = 1.0 if effcai > theta_p_eff else 0.0
        dep = 1.0 if effcai > theta_d else 0.0

        drho = (-rho*(1.0-rho)*(0.5-rho)
                + pot*gamma_p*(1.0-rho)
                - dep*gamma_d*rho) / 70000.0
        rho = min(1.0, max(0.0, rho + dt * drho))

        e_ER_out[i + 1]    = e_ER
        ca_cicr_out[i + 1] = ca_cicr
        effcai_out[i + 1]  = effcai
        rho_out[i + 1]     = rho
        theta_p_out[i + 1] = theta_p_eff

    return {
        "cai_total":  cai_1syn,
        "priming":    e_ER_out,
        "ca_er":      e_ER_out,
        "ca_cicr":    ca_cicr_out,
        "effcai":     effcai_out,
        "rho":        rho_out,
        "theta_p_eff": theta_p_out,
    }


class DualPathwayCICRModel(CICRModel):
    DESCRIPTION = "Dual-Pathway CICR (JAX Batched)"

    FIT_PARAMS = [
        ("gamma_d", 50.0, 200.0), ("gamma_p", 150.0, 300.0),
        ("a00", 0.01, 5.0), ("a01", 0.01, 5.0),
        ("a10", 0.01, 5.0), ("a11", 0.01, 5.0),
        ("a20", 0.01, 5.0), ("a21", 0.01, 5.0),
        ("a30", 0.01, 5.0), ("a31", 0.01, 5.0),
        ("tau_eff", 30.0, 300.0),
        ("k_e", 0.1, 50.0),
        ("tau_ER", 20.0, 200.0),
        ("A_ER", 0.1, 100.0),
        ("tau_CICR", 50.0, 500.0),
        ("cicr_shift", 0.001, 0.15),
        ("K_shift", 0.01, 5.0),
    ]

    DEFAULT_PARAMS = {
        "gamma_d": 100.0, "gamma_p": 250.0,
        "a00": 0.4, "a01": 2.5, "a10": 0.9, "a11": 2.7,
        "a20": 0.7, "a21": 2.3, "a30": 1.4, "a31": 1.7,
        "tau_eff": 200.0,
        "k_e": 10.0, "tau_ER": 80.0,
        "A_ER": 20.0, "tau_CICR": 300.0,
        "cicr_shift": 0.03, "K_shift": 1.0,
    }

    def unpack_params(self, x):
        gb = {'gamma_d': x[0], 'gamma_p': x[1],
              'a00': x[2], 'a01': x[3], 'a10': x[4], 'a11': x[5],
              'a20': x[6], 'a21': x[7], 'a30': x[8], 'a31': x[9]}
        cicr = {'tau_eff': x[10],
                'k_e': x[11], 'tau_ER': x[12],
                'A_ER': x[13], 'tau_CICR': x[14],
                'cicr_shift': x[15], 'K_shift': x[16]}
        return gb, cicr

    def get_init_fn(self):
        def init_fn(cai_first, rho0):
            return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, rho0)
        return init_fn

    def get_debug_sim_fn(self):
        return _apply_release_first_debug

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
            theta_p_base = jnp.maximum(raw_theta_p, theta_d + 0.01)

            cai_rest = 70e-6
            K_trig = 1e-4

            def scan_step(carry, inputs):
                e_ER, ca_cicr, _d1, _d2, _d3, effcai, rho = carry
                cai_raw, dt = inputs

                cai_evoked = jnp.maximum(0.0, cai_raw - cai_rest)

                e_decayed = e_ER * jnp.exp(-dt / cicr['tau_ER'])
                trigger   = cai_evoked / (K_trig + cai_evoked + 1e-12)
                flux      = cicr['A_ER'] * e_decayed * trigger
                released  = jnp.minimum(flux * dt, jnp.maximum(0.0, e_decayed))
                e_new     = jnp.maximum(0.0, e_decayed - released) + cicr['k_e'] * cai_evoked * dt

                ca_cicr_new = ca_cicr * jnp.exp(-dt / cicr['tau_CICR']) + released

                decay_eff  = jnp.exp(-dt / cicr['tau_eff'])
                effcai_new = effcai * decay_eff + cai_evoked * cicr['tau_eff'] * (1.0 - decay_eff)

                theta_p_shift = cicr['cicr_shift'] * ca_cicr_new / (cicr['K_shift'] + ca_cicr_new + 1e-12)
                theta_p_eff = theta_p_base + theta_p_shift

                pot = jnp.where(effcai_new > theta_p_eff, 1.0, 0.0)
                dep = jnp.where(effcai_new > theta_d, 1.0, 0.0)

                drho = (-rho*(1.0-rho)*(0.5-rho)
                        + pot*gb['gamma_p']*(1.0-rho)
                        - dep*gb['gamma_d']*rho) / 70000.0
                rho_new = jnp.clip(rho + dt * drho, 0.0, 1.0)

                return (e_new, ca_cicr_new, 0.0, 0.0, 0.0, effcai_new, rho_new), None
            return scan_step
        return scan_factory


if __name__ == "__main__":
    DualPathwayCICRModel().run()