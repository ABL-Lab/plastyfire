#!/usr/bin/env python3
"""
Biphasic IP3R CICR variant.
Uses smooth, continuous kinetics with Calcium-Dependent Inactivation (bell-shaped gating)
to differentiate spike shapes without discrete event counting.
"""

import jax.numpy as jnp
from cicr_common import CICRModel

class BiphasicCICRModel(CICRModel):
    DESCRIPTION = "Biphasic IP3R CICR (Continuous JAX)"
    
    FIT_PARAMS = [
        # GB Params
        ("gamma_d", 50.0, 200.0), ("gamma_p", 150.0, 300.0),
        ("a00", 1.0, 20.0), ("a01", 1.0, 10.0),
        ("a10", 1.0, 20.0), ("a11", 1.0, 10.0),
        ("a20", 1.0, 20.0), ("a21", 1.0, 10.0),
        ("a30", 1.0, 20.0), ("a31", 1.0, 10.0),
        
        # Biphasic CICR Params
        ("tau_IP3", 100.0, 3000.0),       # IP3 memory decay
        ("v_prod", 0.001, 10.0),          # Max IP3 production rate
        ("K_prod", 0.0001, 0.005),        # Ca threshold to produce IP3
        ("K_P", 0.1, 10.0),               # IP3 threshold to open ER
        ("K_act", 0.0001, 0.003),         # Ca threshold to ACTIVATE ER
        ("K_inh", 0.001, 0.01),           # Ca threshold to INHIBIT ER (The new key)
        ("a_er", 0.0, 0.05),              # SERCA baseline pre
        ("b_er", 0.0, 0.05),              # SERCA baseline post
        ("g_serca", 0.1, 50.0),           # SERCA pump rate
        ("g_release", 0.1, 50.0),         # ER release rate
    ]

    DEFAULT_PARAMS = {
        "gamma_d": 181.4, "gamma_p": 209.9,
        "a00": 1.03, "a01": 1.94, "a10": 1.92, "a11": 3.56,
        "a20": 3.16, "a21": 2.69, "a30": 7.73, "a31": 2.74,
        "tau_IP3": 1500.0, "v_prod": 1.0, "K_prod": 0.001,
        "K_P": 2.0, "K_act": 0.0005, "K_inh": 0.003, 
        "a_er": 0.005, "b_er": 0.005,
        "g_serca": 25.0, "g_release": 5.0,
    }

    def unpack_params(self, x):
        gb = {'gamma_d': x[0], 'gamma_p': x[1], 'a00': x[2], 'a01': x[3],
              'a10': x[4], 'a11': x[5], 'a20': x[6], 'a21': x[7], 'a30': x[8], 'a31': x[9]}
        cicr = {'tau_IP3': x[10], 'v_prod': x[11], 'K_prod': x[12],
                'K_P': x[13], 'K_act': x[14], 'K_inh': x[15], 
                'a_er': x[16], 'b_er': x[17], 'g_serca': x[18], 'g_release': x[19]}
        return gb, cicr

    def get_step_factory(self):
        def scan_factory(params, syn_params):
            gb, cicr = params
            c_pre, c_post, is_apical = syn_params
            
            gamma_er = cicr['a_er'] * c_pre + cicr['b_er'] * c_post
            theta_d = jnp.where(is_apical, gb['a20']*c_pre + gb['a21']*c_post, gb['a00']*c_pre + gb['a01']*c_post)
            theta_p = jnp.where(is_apical, gb['a30']*c_pre + gb['a31']*c_post, gb['a10']*c_pre + gb['a11']*c_post)
            
            def scan_step(carry, inputs):
                # Unpack all 7 carry elements to match init in cicr_common.py
                # (was_above and event_peak are unused in this variant)
                P, ca_er, ca_cicr, _was_above, _event_peak, effcai, rho = carry
                cai_raw, dt = inputs
                
                # 1. Smooth IP3 Production
                prod_rate = cicr['v_prod'] * (cai_raw**2) / (cicr['K_prod']**2 + cai_raw**2 + 1e-12)
                P = P * jnp.exp(-dt / cicr['tau_IP3']) + prod_rate * dt
                
                # 2. Biphasic IP3R/RyR Gate (The Bell Curve)
                ca_cyt = cai_raw + ca_cicr
                
                gate_IP3 = (P**3) / (cicr['K_P']**3 + P**3 + 1e-12)
                gate_act = (ca_cyt**2) / (cicr['K_act']**2 + ca_cyt**2 + 1e-12)
                gate_inh = (cicr['K_inh']**4) / (cicr['K_inh']**4 + ca_cyt**4 + 1e-12) # Steeper drop-off for inhibition
                
                P_open = gate_IP3 * gate_act * gate_inh
                
                # 3. ER Dynamics
                ca_cicr_um = 1000.0 * ca_cicr
                J_serca_raw = cicr['g_serca'] * (ca_cicr_um**2) / (0.1**2 + ca_cicr_um**2)
                J_serca = jnp.where(ca_cicr > gamma_er, jnp.minimum(J_serca_raw, ca_cicr / dt), 0.0)
                
                J_release = jnp.where(ca_er > ca_cyt, cicr['g_release'] * P_open * (ca_er - ca_cyt), 0.0)
                J_net = J_release - J_serca
                
                ca_er = jnp.maximum(0.0, ca_er - (0.83 / 0.17) * J_net * dt)
                ca_cicr = jnp.maximum(0.0, ca_cicr + J_net * dt)
                
                # 4. Graupner-Brunel Plasticity
                decay = jnp.exp(-dt / 200.0)
                effcai = effcai * decay + (cai_raw + ca_cicr - 70e-6) * 200.0 * (1.0 - decay)
                
                dep = jnp.where(effcai > theta_d, 1.0, 0.0)
                pot = jnp.where(effcai > theta_p, 1.0, 0.0)
                drho = (-rho*(1.0-rho)*(0.5-rho) + pot*gb['gamma_p']*(1.0-rho) - dep*gb['gamma_d']*rho) / 70000.0
                rho = jnp.clip(rho + dt * drho, 0.0, 1.0)
                
                return (P, ca_er, ca_cicr, 0.0, 0.0, effcai, rho), None
                
            return scan_step
        return scan_factory

if __name__ == "__main__":
    BiphasicCICRModel().run()