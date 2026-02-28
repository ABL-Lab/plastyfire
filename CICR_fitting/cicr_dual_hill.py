#!/usr/bin/env python3
"""
Dual-Hill CICR variant (analog calcium amplitude tracking).
Rewritten entirely in JAX for massive parallel vectorization.
"""

import jax.numpy as jnp
from cicr_common import CICRModel

class DualHillCICRModel(CICRModel):
    DESCRIPTION = "Dual-Hill CICR (JAX Batched)"
    
    FIT_PARAMS = [
        ("gamma_d", 50.0, 200.0), ("gamma_p", 150.0, 300.0),
        ("a00", 1.0, 10.0), ("a01", 1.0, 10.0),      # Increased to 10
        ("a10", 1.0, 10.0), ("a11", 1.0, 10.0),      # Increased to 10
        ("a20", 1.0, 20.0), ("a21", 1.0, 10.0),      # Increased to 20/10
        ("a30", 1.0, 20.0), ("a31", 1.0, 10.0),      # Increased to 20/10
        ("tau_IP3", 100.0, 5000.0),
        ("K_P", 0.1, 10.0),
        ("n_prime", 3.0, 8.0),
        ("a_evt", 0.0, 0.05), ("b_evt", 0.0, 0.05),
        ("a_er", 0.0, 0.05), ("b_er", 0.0, 0.05),
        ("g_serca", 0.01, 50.0),
        ("g_release", 0.001, 20.0),
        ("sat_cap", 0.0001, 0.5),                    # Lowered min bound to 0.0001
        ("K_Ca", 0.00005, 0.005),
        ("m_ca", 3.0, 8.0),
    ]

    DEFAULT_PARAMS = {
        "gamma_d": 181.4450, "gamma_p": 209.9664,
        "a00": 1.0344, "a01": 1.9424, "a10": 1.9207, "a11": 3.5667,
        "a20": 3.1613, "a21": 2.6934, "a30": 7.7314, "a31": 2.7481,
        "tau_IP3": 2357.0203, "K_P": 8.0160, "n_prime": 3.0173,
        "a_evt": 0.0138, "b_evt": 0.0490, "a_er": 0.0390, "b_er": 0.0104,
        "g_serca": 34.5415, "g_release": 1.1009, "sat_cap": 0.0812,
        "K_Ca": 0.0048, "m_ca": 3.7324,
    }

    def unpack_params(self, x):
        # 1. Start with everything locked to default
        p = dict(self.DEFAULT_PARAMS)
        
        # 2. Overwrite only the parameters the optimizer is currently fitting
        for name, val in zip(self.PARAM_NAMES, x):
            p[name] = val
            
        # 3. Route them to their specific dictionaries
        gb = {k: p[k] for k in ['gamma_d', 'gamma_p', 'a00', 'a01', 'a10', 
                                'a11', 'a20', 'a21', 'a30', 'a31']}
        cicr = {k: p[k] for k in p if k not in gb}
        
        return gb, cicr

    def get_step_factory(self):
        def scan_factory(params, syn_params):
            gb, cicr = params
            c_pre, c_post, is_apical = syn_params
            
            thresh = cicr['a_evt'] * c_pre + cicr['b_evt'] * c_post
            gamma_er = cicr['a_er'] * c_pre + cicr['b_er'] * c_post
            theta_d = jnp.where(is_apical, gb['a20']*c_pre + gb['a21']*c_post, gb['a00']*c_pre + gb['a01']*c_post)
            theta_p = jnp.where(is_apical, gb['a30']*c_pre + gb['a31']*c_post, gb['a10']*c_pre + gb['a11']*c_post)
            
            def scan_step(carry, inputs):
                Priming, ca_er, ca_cicr, was_above, event_peak, effcai, rho = carry
                cai_raw, dt = inputs
                
                is_above = jnp.where(cai_raw > thresh, 1.0, 0.0)
                
                rising = is_above * (1.0 - was_above)
                maint = is_above * was_above
                falling = (1.0 - is_above) * was_above
                
                event_peak = jnp.where(rising, cai_raw, event_peak)
                event_peak = jnp.where(maint, jnp.maximum(event_peak, cai_raw), event_peak)
                event_amp = jnp.where(falling, jnp.minimum((event_peak - thresh)/cicr['sat_cap'], 1.0), 0.0)
                event_peak = jnp.where(falling, 0.0, event_peak)
                
                Priming = Priming * jnp.exp(-dt / cicr['tau_IP3']) + event_amp
                P_safe = jnp.maximum(0.0, Priming)
                P_gate = (P_safe**cicr['n_prime']) / (cicr['K_P']**cicr['n_prime'] + P_safe**cicr['n_prime'] + 1e-12)
                
                ca_raw_safe = jnp.maximum(0.0, cai_raw)
                Ca_gate = (ca_raw_safe**cicr['m_ca']) / (cicr['K_Ca']**cicr['m_ca'] + ca_raw_safe**cicr['m_ca'] + 1e-12)
                
                ca_cyt = cai_raw + ca_cicr
                ca_cicr_um = 1000.0 * ca_cicr
                J_serca_raw = cicr['g_serca'] * (ca_cicr_um**2) / (0.1**2 + ca_cicr_um**2)
                J_serca = jnp.where(ca_cicr > gamma_er, jnp.minimum(J_serca_raw, ca_cicr / dt), 0.0)
                J_release = jnp.where(ca_er > ca_cyt, cicr['g_release'] * P_gate * Ca_gate * (ca_er - ca_cyt), 0.0)
                J_net = J_release - J_serca
                
                ca_er = jnp.maximum(0.0, ca_er - (0.83 / 0.17) * J_net * dt)
                ca_cicr = jnp.maximum(0.0, ca_cicr + J_net * dt)
                
                decay = jnp.exp(-dt / 200.0)
                effcai = effcai * decay + (cai_raw + ca_cicr - 70e-6) * 200.0 * (1.0 - decay)
                
                dep = jnp.where(effcai > theta_d, 1.0, 0.0)
                pot = jnp.where(effcai > theta_p, 1.0, 0.0)
                drho = (-rho*(1.0-rho)*(0.5-rho) + pot*gb['gamma_p']*(1.0-rho) - dep*gb['gamma_d']*rho) / 70000.0
                rho = jnp.clip(rho + dt * drho, 0.0, 1.0)
                
                return (Priming, ca_er, ca_cicr, is_above, event_peak, effcai, rho), None
                
            return scan_step
        return scan_factory

if __name__ == "__main__":
    DualHillCICRModel().run()