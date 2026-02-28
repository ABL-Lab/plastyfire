#!/usr/bin/env python3
"""
Primed CICR variant (sigmoid gate).
Rewritten entirely in JAX for massive parallel vectorization.
"""

import jax.numpy as jnp
from cicr_common import CICRModel

class PrimedCICRModel(CICRModel):
    DESCRIPTION = "Primed CICR (JAX Batched)"
    
    FIT_PARAMS = [
        ("gamma_d", 50.0, 200.0), ("gamma_p", 150.0, 300.0),
        ("a00", 1.0, 5.0), ("a01", 1.0, 5.0),
        ("a10", 1.0, 5.0), ("a11", 1.0, 5.0),
        ("a20", 1.0, 10.0), ("a21", 1.0, 5.0),
        ("a30", 1.0, 10.0), ("a31", 1.0, 5.0),
        ("tau_IP3", 100.0, 5000.0),
        ("IP3_threshold", 0.5, 20.0),
        ("a_evt", 0.0, 0.05), ("b_evt", 0.0, 0.05),
        ("a_er", 0.0, 0.05), ("b_er", 0.0, 0.05),
        ("g_serca", 0.01, 50.0),
        ("g_release", 0.001, 10.0),
        ("sat_cap", 0.001, 0.5),
        ("sigmoid_slope", 1.0, 20.0),
        ("K_act", 0.0001, 0.005),
    ]

    DEFAULT_PARAMS = {
        "gamma_d": 101.5, "gamma_p": 216.2,
        "a00": 1.002, "a01": 1.954, "a10": 1.159, "a11": 2.483,
        "a20": 1.127, "a21": 2.456, "a30": 5.236, "a31": 1.782,
        "tau_IP3": 500.0, "IP3_threshold": 3.0,
        "a_evt": 0.003, "b_evt": 0.003, "a_er": 0.005, "b_er": 0.005,
        "g_serca": 1.9565, "g_release": 0.5, "sat_cap": 0.05,
        "sigmoid_slope": 5.0, "K_act": 0.0005,
    }

    def unpack_params(self, x):
        gb = {'gamma_d': x[0], 'gamma_p': x[1], 'a00': x[2], 'a01': x[3],
              'a10': x[4], 'a11': x[5], 'a20': x[6], 'a21': x[7], 'a30': x[8], 'a31': x[9]}
        cicr = {'tau_IP3': x[10], 'IP3_threshold': x[11],
                'a_evt': x[12], 'b_evt': x[13], 'a_er': x[14], 'b_er': x[15],
                'g_serca': x[16], 'g_release': x[17], 'sat_cap': x[18],
                'sigmoid_slope': x[19], 'K_act': x[20]}
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
                arg = -cicr['sigmoid_slope'] * (Priming - cicr['IP3_threshold'])
                P_gate = jnp.where(arg > 500.0, 0.0, jnp.where(arg < -500.0, 1.0, 1.0 / (1.0 + jnp.exp(arg))))
                
                ca_cyt = cai_raw + ca_cicr
                ca_cicr_um = 1000.0 * ca_cicr
                J_serca_raw = cicr['g_serca'] * (ca_cicr_um**2) / (0.1**2 + ca_cicr_um**2)
                J_serca = jnp.where(ca_cicr > gamma_er, jnp.minimum(J_serca_raw, ca_cicr / dt), 0.0)
                J_release = jnp.where(ca_er > ca_cyt, cicr['g_release'] * P_gate * (ca_cyt / (ca_cyt + cicr['K_act'])) * (ca_er - ca_cyt), 0.0)
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
    PrimedCICRModel().run()