#!/usr/bin/env python3
"""
Synergistic Coincidence-Gated CICR Graupner-Brunel Model.

CICR adds calcium to the cytosolic pool, driving it over the potentiation
threshold. To perfectly separate +10ms and -10ms, CICR acts as a steep 
coincidence detector using a high-order Hill equation on Ca_raw.
"""

import numpy as np
import jax.numpy as jnp
from cicr_common import CICRModel

def _apply_synergistic_cicr_debug(cai_1syn, t, cp, cq, is_apical, dp):
    """Plain-NumPy exact mirror of the JAX scan loop for diagnostic plotting."""
    tau_IP3  = float(dp["tau_IP3"])
    k_IP3    = float(dp["k_IP3"])
    K_IP3    = float(dp["K_IP3"])
    K_act    = float(dp["K_act"])
    n_H      = float(dp["n_H"])
    A_ER     = float(dp["A_ER"])
    tau_ER   = float(dp["tau_ER"])
    tau_CICR = float(dp["tau_CICR"])
    
    cai_rest = 70e-6
    n = len(cai_1syn)
    
    cai_total   = np.copy(cai_1syn)
    IP3_out     = np.zeros(n)
    R_ER_out    = np.zeros(n)
    ca_cicr_out = np.zeros(n)
    
    IP3 = 0.0
    R_ER = 1.0
    ca_cicr = 0.0
    
    IP3_out[0] = IP3
    R_ER_out[0] = R_ER
    ca_cicr_out[0] = ca_cicr
    
    for i in range(n - 1):
        dt = t[i + 1] - t[i]
        if dt <= 0.0:
            IP3_out[i + 1]     = IP3
            R_ER_out[i + 1]    = R_ER
            ca_cicr_out[i + 1] = ca_cicr
            continue
            
        cai_evoked = max(0.0, cai_1syn[i] - cai_rest)
        
        # 1. IP3 Integration (Stimulus Memory)
        decay_ip3 = np.exp(-dt / tau_IP3)
        IP3 = IP3 * decay_ip3 + k_IP3 * cai_evoked * tau_IP3 * (1.0 - decay_ip3)
        
        # 2. Steep Coincidence Detection (IP3R Gating)
        ip3_term = (IP3**2) / (K_IP3**2 + IP3**2 + 1e-12)
        ca_term  = (cai_evoked**n_H) / (K_act**n_H + cai_evoked**n_H + 1e-12)
        J_rate   = A_ER * ip3_term * ca_term
        
        # 3. ER Dynamics (Unconditionally stable exact integration)
        tau_eff_R = 1.0 / (1.0 / tau_ER + J_rate + 1e-12)
        R_inf = (1.0 / tau_ER) * tau_eff_R
        R_ER_new = R_inf + (R_ER - R_inf) * np.exp(-dt / tau_eff_R)
        
        J_cicr = J_rate * R_ER
        R_ER = R_ER_new
        
        # 4. CICR Calcium Integration
        decay_cicr = np.exp(-dt / tau_CICR)
        ca_cicr = ca_cicr * decay_cicr + J_cicr * tau_CICR * (1.0 - decay_cicr)
        
        cai_total[i + 1]   = cai_1syn[i + 1] + ca_cicr
        IP3_out[i + 1]     = IP3
        R_ER_out[i + 1]    = R_ER # Plotted as a dimensionless fraction
        ca_cicr_out[i + 1] = ca_cicr
        
    return {
        "cai_total": cai_total,
        "priming":   IP3_out,
        "ca_er":     R_ER_out, 
        "ca_cicr":   ca_cicr_out,
    }

class SynergisticCICRModel(CICRModel):
    DESCRIPTION = "Synergistic Coincidence-Gated CICR (JAX)"

    # Defined required upper and lower bounds for differential evolution
    FIT_PARAMS = [
        ("gamma_d", 50.0, 300.0), 
        ("gamma_p", 150.0, 500.0),
        ("a00", 0.01, 5.0), 
        ("a01", 0.01, 5.0), 
        ("a10", 0.01, 5.0), 
        ("a11", 0.01, 5.0),
        ("a20", 0.01, 5.0), 
        ("a21", 0.01, 5.0), 
        ("a30", 0.01, 5.0), 
        ("a31", 0.01, 5.0),
        ("tau_IP3", 100.0, 2000.0), 
        ("k_IP3", 0.1, 50.0), 
        ("K_IP3", 0.01, 5.0),
        ("K_act", 0.0001, 0.05), 
        ("n_H", 1.0, 8.0), 
        ("A_ER", 0.1, 50.0),
        ("tau_ER", 10.0, 1000.0), 
        ("tau_CICR", 10.0, 500.0)
    ]

    DEFAULT_PARAMS = {
        "gamma_d": 100.0, "gamma_p": 250.0,
        "a00": 0.4, "a01": 2.5, "a10": 0.9, "a11": 2.7,
        "a20": 0.7, "a21": 2.3, "a30": 1.4, "a31": 1.7,
        "tau_IP3": 500.0, "k_IP3": 10.0, "K_IP3": 0.5,
        "K_act": 0.001, "n_H": 3.0, "A_ER": 5.0,
        "tau_ER": 200.0, "tau_CICR": 50.0,
    }

    def unpack_params(self, x):
        gb = {'gamma_d': x[0], 'gamma_p': x[1],
              'a00': x[2], 'a01': x[3], 'a10': x[4], 'a11': x[5],
              'a20': x[6], 'a21': x[7], 'a30': x[8], 'a31': x[9]}
        cicr = {'tau_IP3': x[10], 'k_IP3': x[11], 'K_IP3': x[12],
                'K_act': x[13], 'n_H': x[14], 'A_ER': x[15],
                'tau_ER': x[16], 'tau_CICR': x[17]}
        return gb, cicr

    def get_init_fn(self):
        def init_fn(cai_first, rho0):
            # State vector length = 7 (matches framework assumptions in cicr_common logic)
            # (IP3, R_ER, ca_cicr, dummy, dummy, effcai, rho)
            return (0.0, 1.0, 0.0, 0.0, 0.0, 0.0, rho0)
        return init_fn

    def get_debug_sim_fn(self):
        return _apply_synergistic_cicr_debug

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
                IP3, R_ER, ca_cicr, _d1, _d2, effcai, rho = carry
                cai_raw, dt = inputs

                cai_evoked = jnp.maximum(0.0, cai_raw - cai_rest)

                # 1. IP3 Integration
                IP3_decay = jnp.exp(-dt / cicr['tau_IP3'])
                IP3_new = IP3 * IP3_decay + cicr['k_IP3'] * cai_evoked * cicr['tau_IP3'] * (1.0 - IP3_decay)

                # 2. Steep Coincidence Detection (The Discriminator)
                ip3_term = (IP3_new**2) / (cicr['K_IP3']**2 + IP3_new**2 + 1e-12)
                ca_term  = (cai_evoked**cicr['n_H']) / (cicr['K_act']**cicr['n_H'] + cai_evoked**cicr['n_H'] + 1e-12)
                
                J_rate = cicr['A_ER'] * ip3_term * ca_term

                # 3. ER Dynamics (Unconditionally Stable)
                tau_eff_R = 1.0 / (1.0 / cicr['tau_ER'] + J_rate + 1e-12)
                R_inf = (1.0 / cicr['tau_ER']) * tau_eff_R
                R_ER_new = R_inf + (R_ER - R_inf) * jnp.exp(-dt / tau_eff_R)
                
                J_cicr = J_rate * R_ER

                # 4. Ca_CICR Release integration
                ca_cicr_decay = jnp.exp(-dt / cicr['tau_CICR'])
                ca_cicr_new = ca_cicr * ca_cicr_decay + J_cicr * cicr['tau_CICR'] * (1.0 - ca_cicr_decay)

                # 5. EffCai integration (Synergistically summing Raw + CICR)
                ca_total = cai_evoked + ca_cicr_new
                eff_decay = jnp.exp(-dt / 200.0) # standard tau_Ca = 200ms
                effcai_new = effcai * eff_decay + ca_total * 200.0 * (1.0 - eff_decay)

                # 6. Graupner-Brunel Plasticity Core
                pot = jnp.where(effcai_new > theta_p, 1.0, 0.0)
                dep = jnp.where(effcai_new > theta_d, 1.0, 0.0)

                drho = (-rho*(1.0-rho)*(0.5-rho)
                        + pot*gb['gamma_p']*(1.0-rho)
                        - dep*gb['gamma_d']*rho) / 70000.0
                
                rho_new = jnp.clip(rho + dt * drho, 0.0, 1.0)

                return (IP3_new, R_ER_new, ca_cicr_new, 0.0, 0.0, effcai_new, rho_new), None
            return scan_step
        return scan_factory

if __name__ == "__main__":
    SynergisticCICRModel().run()