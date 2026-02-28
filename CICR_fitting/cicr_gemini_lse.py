#!/usr/bin/env python3
"""
Dynamic-Tau STAPCD Eligibility Model.

Inspired by Caya-Bissonnette et al. (STAPCD) and Chindemi et al. (2022).
Instead of just adding CICR amplitude, this model captures the core experimental
observation: ER priming dramatically LENGTHENS the decay time constant of the 
calcium transient when hit by a subsequent instructive bAP burst.

1. EPSP Ca2+ loads the ER state (E_state) via exact integration.
2. A steep Hill trigger isolates high-amplitude events (bAPs).
3. The conjunctive STAPCD Activation (S) dynamically stretches the effcai 
   low-pass filter tau, creating massive integrals for +10ms and fast decay for -10ms.
"""

import numpy as np
import jax.numpy as jnp
from cicr_common import CICRModel

def _apply_dynamic_tau_debug(cai_1syn, t, cp, cq, is_apical, dp):
    """Plain-NumPy exact mirror of the JAX scan loop for diagnostic plotting."""
    alpha_E     = float(dp["alpha_E"])
    tau_E       = float(dp["tau_E"])
    K_trig      = float(dp["K_trig"])
    n_H         = float(dp["n_H"])
    tau_base    = float(dp["tau_base"])
    tau_stretch = float(dp["tau_stretch"])
    A_CICR      = float(dp["A_CICR"])
    
    n = len(cai_1syn)
    cai_rest = 70e-6
    
    cai_total   = np.copy(cai_1syn)
    E_trace_out = np.zeros(n)
    S_out       = np.zeros(n)
    
    E_state = 0.0
    
    E_trace_out[0] = E_state
    S_out[0]       = 0.0
    
    for i in range(n - 1):
        dt = t[i + 1] - t[i]
        if dt <= 0.0:
            E_trace_out[i + 1] = E_state
            S_out[i + 1]       = 0.0
            continue
            
        cai_evoked = max(0.0, cai_1syn[i] - cai_rest)
        
        # 1. Exact integration of ER Loading State (Eligibility Trace)
        # Bounded unconditionally between 0 and 1
        rate_load = alpha_E * cai_evoked
        inv_tau_E = 1.0 / tau_E
        tau_eff_E = 1.0 / (rate_load + inv_tau_E + 1e-12)
        E_inf = rate_load * tau_eff_E
        E_state = E_inf + (E_state - E_inf) * np.exp(-dt / tau_eff_E)
        
        # 2. Steep Coincidence Trigger (The Instructive Readout)
        trigger = (cai_evoked**n_H) / (K_trig**n_H + cai_evoked**n_H + 1e-12)
        
        # 3. Conjunctive STAPCD Activation 
        S = E_state * trigger
        
        # 4. Amplitude Boost (Optional CICR addition)
        ca_boost = cai_evoked + A_CICR * S
        
        cai_total[i + 1]   = ca_boost + cai_rest
        E_trace_out[i + 1] = E_state
        S_out[i + 1]       = S  # Plotted in the 'ca_er' row to visualize trigger
        
    return {
        "cai_total": cai_total,
        "priming":   E_trace_out,         # Visualizing the slowly decaying ER state
        "ca_er":     S_out,               # Visualizing the instantaneous STAPCD trigger
        "ca_cicr":   cai_total - cai_1syn # Visualizing amplitude boost
    }


class DynamicTauEligibilityModel(CICRModel):
    DESCRIPTION = "Dynamic-Tau STAPCD Eligibility (JAX)"

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
        ("alpha_E", 100.0, 50000.0),    # ER loading rate multiplier
        ("tau_E", 200.0, 2000.0),       # ER trace decay time (eligibility window)
        ("K_trig", 0.0005, 0.01),       # Ca2+ threshold to trigger STAPCD
        ("n_H", 1.0, 8.0),              # Hill coefficient for trigger steepness
        ("tau_base", 10.0, 200.0),      # Fast Ca2+ decay without STAPCD
        ("tau_stretch", 100.0, 3000.0), # Massive stretching of tau when STAPCD triggers
        ("A_CICR", 0.0, 0.05)           # Direct amplitude boost from CICR
    ]

    DEFAULT_PARAMS = {
        "gamma_d": 109.9853, "gamma_p": 155.8440,
        "a00": 3.6260, "a01": 0.3992, "a10": 4.6528, "a11": 4.6151,
        "a20": 1.7521, "a21": 2.6680, "a30": 4.2831, "a31": 4.5037,
        "alpha_E": 33191.7142, "tau_E": 1397.1478,
        "K_trig": 0.0010, "n_H": 4.6538, 
        "tau_base": 163.7045, "tau_stretch": 692.0080, 
        "A_CICR": 0.0482 
    }

    def unpack_params(self, x):
        gb = {'gamma_d': x[0], 'gamma_p': x[1],
              'a00': x[2], 'a01': x[3], 'a10': x[4], 'a11': x[5],
              'a20': x[6], 'a21': x[7], 'a30': x[8], 'a31': x[9]}
        st = {'alpha_E': x[10], 'tau_E': x[11], 'K_trig': x[12],
              'n_H': x[13], 'tau_base': x[14], 'tau_stretch': x[15],
              'A_CICR': x[16]}
        return gb, st

    def get_init_fn(self):
        def init_fn(cai_first, rho0):
            # State vector length = 7 to perfectly match JAX unpacking inside cicr_common.py
            # (E_state, dummy, dummy, dummy, dummy, effcai, rho)
            return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, rho0)
        return init_fn

    def get_debug_sim_fn(self):
        return _apply_dynamic_tau_debug

    def get_step_factory(self):
        def scan_factory(params, syn_params):
            gb, st = params
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
                E_state, _d1, _d2, _d3, _d4, effcai, rho = carry
                cai_raw, dt = inputs

                cai_evoked = jnp.maximum(0.0, cai_raw - cai_rest)

                # 1. Exact integration of ER Loading State (Eligibility Trace)
                rate_load = st['alpha_E'] * cai_evoked
                inv_tau_E = 1.0 / st['tau_E']
                tau_eff_E = 1.0 / (rate_load + inv_tau_E + 1e-12)
                E_inf = rate_load * tau_eff_E
                E_new = E_inf + (E_state - E_inf) * jnp.exp(-dt / tau_eff_E)

                # 2. Instructive Signal Trigger (Readout of the Trace)
                trigger = (cai_evoked**st['n_H']) / (st['K_trig']**st['n_H'] + cai_evoked**st['n_H'] + 1e-12)
                
                # 3. Conjunctive STAPCD Activation 
                S = E_new * trigger

                # 4. Dynamic Tau stretching and Amplitude Boost
                tau_dyn = st['tau_base'] + st['tau_stretch'] * S
                ca_boost = cai_evoked + st['A_CICR'] * S

                # 5. EffCai leaky integrator (Scale boost by 1000 to match threshold scales 0.1 - 5.0)
                eff_decay = jnp.exp(-dt / tau_dyn)
                effcai_new = effcai * eff_decay + (ca_boost * 1000.0) * (1.0 - eff_decay) 

                # 6. Graupner-Brunel Plasticity Core
                pot = jnp.where(effcai_new > theta_p, 1.0, 0.0)
                dep = jnp.where(effcai_new > theta_d, 1.0, 0.0)

                drho = (-rho*(1.0-rho)*(0.5-rho)
                        + pot*gb['gamma_p']*(1.0-rho)
                        - dep*gb['gamma_d']*rho) / 70000.0
                
                rho_new = jnp.clip(rho + dt * drho, 0.0, 1.0)

                return (E_new, 0.0, 0.0, 0.0, 0.0, effcai_new, rho_new), None
            return scan_step
        return scan_factory


if __name__ == "__main__":
    DynamicTauEligibilityModel().run()