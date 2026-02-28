#!/usr/bin/env python3
"""
Primed-Junctional CICR Model.

A functionally abstracted, biophysically grounded model of ER-mediated Calcium-Induced Calcium Release.

1. ER Priming: Wide, moderate-amplitude events (EPSPs) slowly prime the ER. 
   High-amplitude events (bAPs) are soft-clipped (K_prime_limit) to prevent artificial self-priming.
2. Junctional Trigger: A steep, high-threshold Ca2+ sensor (K_trig, n_trig) strictly isolates 
   massive bAP-driven Ca2+ influxes at ER-PM junctions.
3. Synergistic CICR: Massive ER Ca2+ release only occurs when the ER is both Primed AND Triggered.
4. Total Calcium is low-pass filtered (tau_eff) to drive Graupner-Brunel bistable plasticity.
"""

import numpy as np
import jax.numpy as jnp
from cicr_common import CICRModel

def _apply_primed_junctional_debug(cai_1syn, t, cp, cq, is_apical, dp):
    """Plain-NumPy exact mirror of the JAX scan loop for diagnostic plotting.

    If dp contains 'interp_dt', interpolates trace to finer resolution,
    then resamples output back to original time grid for plotting.
    """
    tau_prime     = float(dp["tau_prime"])
    k_prime       = float(dp["k_prime"])
    K_prime_limit = float(dp["K_prime_limit"])
    K_trig        = float(dp["K_trig"])
    n_trig        = float(dp["n_trig"])
    Vmax_CICR     = float(dp["Vmax_CICR"])
    tau_CICR      = float(dp["tau_CICR"])

    cai_rest = 70e-6
    interp_dt = float(dp.get("interp_dt", 0))
    t_orig = t

    # Optionally interpolate to finer grid
    if interp_dt > 0:
        t_fine = np.arange(t[0], t[-1], interp_dt)
        cai_1syn = np.interp(t_fine, t, cai_1syn)
        t = t_fine

    n = len(cai_1syn)
    cai_total   = np.copy(cai_1syn)
    prime_out   = np.zeros(n)
    ca_cicr_out = np.zeros(n)

    prime_state = 0.0
    ca_cicr = 0.0

    prime_out[0]   = prime_state
    ca_cicr_out[0] = ca_cicr

    for i in range(n - 1):
        dt = t[i + 1] - t[i]
        if dt <= 0.0:
            prime_out[i + 1]   = prime_state
            ca_cicr_out[i + 1] = ca_cicr
            continue

        cai_evoked = max(0.0, cai_1syn[i] - cai_rest)

        # 1. Soft-Clipped ER Priming
        cai_clip = K_prime_limit * np.tanh(cai_evoked / (K_prime_limit + 1e-12))
        rate_load = k_prime * cai_clip
        inv_tau_p = 1.0 / tau_prime
        tau_eff_p = 1.0 / (rate_load + inv_tau_p + 1e-12)
        P_inf = rate_load * tau_eff_p
        prime_state = P_inf + (prime_state - P_inf) * np.exp(-dt / tau_eff_p)

        # 2. Junctional Trigger
        trigger = (cai_evoked**n_trig) / (K_trig**n_trig + cai_evoked**n_trig + 1e-12)

        # 3. Synergistic CICR Release
        J_rate  = Vmax_CICR * prime_state * trigger
        decay_cicr = np.exp(-dt / tau_CICR)
        ca_cicr = ca_cicr * decay_cicr + J_rate * tau_CICR * (1.0 - decay_cicr)

        cai_total[i + 1]   = cai_1syn[i + 1] + ca_cicr
        prime_out[i + 1]   = prime_state
        ca_cicr_out[i + 1] = ca_cicr

    # Resample back to original time grid if interpolated
    if interp_dt > 0:
        cai_total   = np.interp(t_orig, t, cai_total)
        prime_out   = np.interp(t_orig, t, prime_out)
        ca_cicr_out = np.interp(t_orig, t, ca_cicr_out)

    return {
        "cai_total": cai_total,
        "priming":   prime_out,
        "ca_er":     prime_out,
        "ca_cicr":   ca_cicr_out,
    }


class PrimedJunctionalCICRModel(CICRModel):
    DESCRIPTION = "Primed-Junctional CICR Synergy (JAX)"

    FIT_PARAMS = [
        ("gamma_d", 50.0, 300.0), 
        ("gamma_p", 150.0, 600.0),
        ("a00", 1.0, 5.0), 
        ("a01", 1.0, 5.0), 
        ("a10", 1.0, 5.0), 
        ("a11", 1.0, 5.0),
        ("a20", 1.0, 5.0), 
        ("a21", 1.0, 5.0), 
        ("a30", 1.0, 5.0), 
        ("a31", 1.0, 5.0),
        ("tau_prime", 100.0, 2000.0),     # ER Priming decay time
        ("k_prime", 100.0, 50000.0),      # ER Priming rate
        ("K_prime_limit", 0.0001, 0.001), # Soft-clip limit isolating EPSPs
        ("K_trig", 0.001, 0.005),         # Trigger Ca2+ threshold
        ("n_trig", 2.0, 12.0),            # Hill steepness of the trigger
        ("Vmax_CICR", 0.00001, 0.01),     # Max Ca2+ efflux from ER
        ("tau_CICR", 10.0, 800.0),        # Decay time of the released Ca2+
        ("tau_eff", 50.0, 500.0)          # Chindemi integrator time constant
    ]

    
    # Pre-loaded with the optimized parameters you generated
    DEFAULT_PARAMS = {
        "gamma_d": 68.3033, "gamma_p": 261.5785,
        "a00": 1.0133, "a01": 1.5540, "a10": 1.5453, "a11": 4.7975,
        "a20": 1.0579, "a21": 3.0020, "a30": 4.7076, "a31": 4.3526,
        "tau_prime": 731.6300, "k_prime": 41485.0180, "K_prime_limit": 0.0006,
        "K_trig": 0.0019, "n_trig": 4.7157, 
        "Vmax_CICR": 0.0004, "tau_CICR": 79.9375, 
        "tau_eff": 408.6018
    }

    def unpack_params(self, x):
        gb = {'gamma_d': x[0], 'gamma_p': x[1],
              'a00': x[2], 'a01': x[3], 'a10': x[4], 'a11': x[5],
              'a20': x[6], 'a21': x[7], 'a30': x[8], 'a31': x[9]}
        cicr = {'tau_prime': x[10], 'k_prime': x[11], 'K_prime_limit': x[12],
                'K_trig': x[13], 'n_trig': x[14], 'Vmax_CICR': x[15],
                'tau_CICR': x[16], 'tau_eff': x[17]}
        return gb, cicr

    def get_init_fn(self):
        def init_fn(cai_first, rho0):
            # State vector length = 7 (matches framework assumptions in cicr_common logic)
            # (prime_state, ca_cicr, dummy, dummy, dummy, effcai, rho)
            return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, rho0)
        return init_fn

    def get_debug_sim_fn(self):
        return _apply_primed_junctional_debug

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
                prime_state, ca_cicr, _d1, _d2, _d3, effcai, rho = carry
                cai_raw, dt = inputs

                cai_evoked = jnp.maximum(0.0, cai_raw - cai_rest)

                # 1. Soft-Clipped ER Priming
                cai_clip = cicr['K_prime_limit'] * jnp.tanh(cai_evoked / (cicr['K_prime_limit'] + 1e-12))
                
                rate_load = cicr['k_prime'] * cai_clip
                inv_tau_p = 1.0 / cicr['tau_prime']
                tau_eff_p = 1.0 / (rate_load + inv_tau_p + 1e-12)
                P_inf = rate_load * tau_eff_p
                P_new = P_inf + (prime_state - P_inf) * jnp.exp(-dt / tau_eff_p)

                # 2. Junctional Trigger
                trigger = (cai_evoked**cicr['n_trig']) / (cicr['K_trig']**cicr['n_trig'] + cai_evoked**cicr['n_trig'] + 1e-12)
                
                # 3. Synergistic CICR Release
                J_rate = cicr['Vmax_CICR'] * P_new * trigger

                ca_cicr_decay = jnp.exp(-dt / cicr['tau_CICR'])
                ca_cicr_new = ca_cicr * ca_cicr_decay + J_rate * cicr['tau_CICR'] * (1.0 - ca_cicr_decay)

                # 4. Chindemi c* Leaky Integrator
                ca_total = cai_evoked + ca_cicr_new
                eff_decay = jnp.exp(-dt / cicr['tau_eff'])
                effcai_new = effcai * eff_decay + ca_total * cicr['tau_eff'] * (1.0 - eff_decay) 

                # 5. Graupner-Brunel Plasticity Core
                pot = jnp.where(effcai_new > theta_p, 1.0, 0.0)
                dep = jnp.where(effcai_new > theta_d, 1.0, 0.0)

                drho = (-rho*(1.0-rho)*(0.5-rho)
                        + pot*gb['gamma_p']*(1.0-rho)
                        - dep*gb['gamma_d']*rho) / 70000.0
                
                rho_new = jnp.clip(rho + dt * drho, 0.0, 1.0)

                return (P_new, ca_cicr_new, 0.0, 0.0, 0.0, effcai_new, rho_new), None
            
            return scan_step
        return scan_factory

if __name__ == "__main__":
    PrimedJunctionalCICRModel().run()