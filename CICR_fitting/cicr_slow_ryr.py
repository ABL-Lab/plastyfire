#!/usr/bin/env python3
"""
Slow-RyR Accumulator CICR variant.
Models ER release as a slow, continuous buildup (>2 seconds) rather than fast spikes,
aligning with experimental observations of prolonged calcium extrusion.
"""

import jax.numpy as jnp
from cicr_common import CICRModel

class SlowRyRCICRModel(CICRModel):
    DESCRIPTION = "Slow-RyR Accumulator CICR (JAX Batched)"
    
    FIT_PARAMS = [
        # GB Params
        ("gamma_d", 50.0, 200.0), ("gamma_p", 150.0, 300.0),
        ("a00", 1.0, 20.0), ("a01", 1.0, 10.0),
        ("a10", 1.0, 20.0), ("a11", 1.0, 10.0),
        ("a20", 1.0, 20.0), ("a21", 1.0, 10.0),
        ("a30", 1.0, 20.0), ("a31", 1.0, 10.0),
        
        # Slow-RyR CICR Params
        ("tau_IP3", 50.0, 1000.0),        # IP3 memory decay (fast spark)
        ("v_prod", 0.001, 10.0),          # Max IP3 production rate
        ("K_prod", 0.0001, 0.005),        # Ca threshold to produce IP3
        ("tau_RyR", 500.0, 5000.0),       # RyR accumulator time constant (The slow tide)
        ("K_P", 0.1, 10.0),               # IP3 threshold to open RyR
        ("n_prime", 1.0, 8.0),            # IP3 Hill steepness
        ("K_act", 0.0001, 0.005),         # Ca threshold to activate RyR
        ("m_ca", 1.0, 8.0),               # Ca Hill steepness
        ("a_er", 0.0, 0.05),              # SERCA baseline pre
        ("b_er", 0.0, 0.05),              # SERCA baseline post
        ("g_serca", 0.1, 50.0),           # SERCA pump rate
        ("g_release", 0.1, 50.0),         # ER release rate
    ]

    DEFAULT_PARAMS = {
        "gamma_d": 181.4, "gamma_p": 209.9,
        "a00": 1.03, "a01": 1.94, "a10": 1.92, "a11": 3.56,
        "a20": 3.16, "a21": 2.69, "a30": 7.73, "a31": 2.74,
        "tau_IP3": 200.0, "v_prod": 1.0, "K_prod": 0.001,
        "tau_RyR": 2000.0, "K_P": 2.0, "n_prime": 3.0, 
        "K_act": 0.0005, "m_ca": 3.0,
        "a_er": 0.005, "b_er": 0.005,
        "g_serca": 25.0, "g_release": 5.0,
    }

    def unpack_params(self, x):
        gb = {'gamma_d': x[0], 'gamma_p': x[1], 'a00': x[2], 'a01': x[3],
              'a10': x[4], 'a11': x[5], 'a20': x[6], 'a21': x[7], 'a30': x[8], 'a31': x[9]}
        cicr = {'tau_IP3': x[10], 'v_prod': x[11], 'K_prod': x[12],
                'tau_RyR': x[13], 'K_P': x[14], 'n_prime': x[15],
                'K_act': x[16], 'm_ca': x[17], 'a_er': x[18], 'b_er': x[19],
                'g_serca': x[20], 'g_release': x[21]}
        return gb, cicr

    def get_step_factory(self):
        def scan_factory(params, syn_params):
            gb, cicr = params
            c_pre, c_post, is_apical = syn_params
            
            gamma_er = cicr['a_er'] * c_pre + cicr['b_er'] * c_post
            
            # Calculate raw thresholds
            theta_d = jnp.where(is_apical, gb['a20']*c_pre + gb['a21']*c_post, gb['a00']*c_pre + gb['a01']*c_post)
            raw_theta_p = jnp.where(is_apical, gb['a30']*c_pre + gb['a31']*c_post, gb['a10']*c_pre + gb['a11']*c_post)
            
            # FLAW 2 FIX: Mathematically force theta_p to be at least 0.01 above theta_d
            theta_p = jnp.maximum(raw_theta_p, theta_d + 0.01)
            
            def scan_step(carry, inputs):
                P, ca_er, ca_cicr, RyR_open, _dummy, effcai, rho = carry
                cai_raw, dt = inputs
                
                # FLAW 1 FIX: Subtract resting baseline (70 nM) before CICR evaluation
                cai_evoked = jnp.maximum(0.0, cai_raw - 70e-6)
                
                # 1. Fast IP3 Production (Driven ONLY by evoked Ca)
                f_ca = (cai_evoked**2) / (cicr['K_prod']**2 + cai_evoked**2 + 1e-12)
                prod_rate = cicr['v_prod'] * f_ca
                P_new = P * jnp.exp(-dt / cicr['tau_IP3']) + prod_rate * dt
                
                # 2. Slow RyR Coincidence Target 
                ca_cyt_evoked = cai_evoked + ca_cicr
                P_safe = jnp.maximum(0.0, P_new)
                
                gate_IP3 = (P_safe**cicr['n_prime']) / (cicr['K_P']**cicr['n_prime'] + P_safe**cicr['n_prime'] + 1e-12)
                gate_act = (ca_cyt_evoked**cicr['m_ca']) / (cicr['K_act']**cicr['m_ca'] + ca_cyt_evoked**cicr['m_ca'] + 1e-12)
                RyR_inf = gate_IP3 * gate_act
                
                # 3. Slow RyR Accumulator Integration (The Tide)
                decay_RyR = jnp.exp(-dt / cicr['tau_RyR'])
                RyR_open_new = RyR_open * decay_RyR + RyR_inf * (1.0 - decay_RyR)
                
                # 4. ER Release & SERCA Dynamics (SERCA still pumps total CICR load)
                ca_cicr_um = 1000.0 * ca_cicr
                J_serca_raw = cicr['g_serca'] * (ca_cicr_um**2) / (0.1**2 + ca_cicr_um**2)
                J_serca = jnp.where(ca_cicr > gamma_er, jnp.minimum(J_serca_raw, ca_cicr / dt), 0.0)
                
                ca_cyt_total = cai_raw + ca_cicr
                J_release = jnp.where(ca_er > ca_cyt_total, cicr['g_release'] * RyR_open_new * (ca_er - ca_cyt_total), 0.0)
                J_net = J_release - J_serca
                
                ca_er_new = jnp.maximum(0.0, ca_er - (0.83 / 0.17) * J_net * dt)
                ca_cicr_new = jnp.maximum(0.0, ca_cicr + J_net * dt)
                
                # 5. Graupner-Brunel Plasticity
                decay_eff = jnp.exp(-dt / 200.0)
                effcai_new = effcai * decay_eff + (cai_raw + ca_cicr_new - 70e-6) * 200.0 * (1.0 - decay_eff)
                
                dep = jnp.where(effcai_new > theta_d, 1.0, 0.0)
                pot = jnp.where(effcai_new > theta_p, 1.0, 0.0)
                drho = (-rho*(1.0-rho)*(0.5-rho) + pot*gb['gamma_p']*(1.0-rho) - dep*gb['gamma_d']*rho) / 70000.0
                rho_new = jnp.clip(rho + dt * drho, 0.0, 1.0)
                
                return (P_new, ca_er_new, ca_cicr_new, RyR_open_new, 0.0, effcai_new, rho_new), None
                
            return scan_step
        return scan_factory

if __name__ == "__main__":
    SlowRyRCICRModel().run()