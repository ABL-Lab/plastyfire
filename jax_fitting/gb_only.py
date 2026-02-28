#!/usr/bin/env python3
"""
GB-Only Baseline Model (Continuous JAX).
Turns off CICR completely to fit only the 10 Graupner-Brunel parameters.
Use this to find a solid baseline before freezing GB and tuning CICR.
"""

import numpy as np
import jax.numpy as jnp
from cicr_common import CICRModel

def _apply_gb_only_debug(cai_1syn, t, cp, cq, is_apical, dp):
    """Plain-NumPy debug simulation. Returns zero for all CICR variables."""
    n = len(cai_1syn)
    cai_total = np.copy(cai_1syn)
    
    return {
        "cai_total": cai_total,
        "priming":   np.zeros(n), # No IP3
        "ca_er":     np.zeros(n), # No ER pool
        "ca_cicr":   np.zeros(n), # No CICR flux
    }

class GBOnlyModel(CICRModel):
    DESCRIPTION = "GB-Only Baseline (No CICR)"
        
    FIT_PARAMS = [
            # Only the 10 GB parameters are exposed to the optimizer
            ("gamma_d", 50.0, 500.0), ("gamma_p", 150.0, 500.0),
            ("a00", 0.01, 20.0), ("a01", 0.01, 20.0),
            ("a10", 0.01, 20.0), ("a11", 0.01, 20.0),
            ("a20", 0.01, 20.0), ("a21", 0.01, 20.0),
            ("a30", 0.01, 20.0), ("a31", 0.01, 20.0),
        ]
    
    DEFAULT_PARAMS = {
        "gamma_d": 181.4, "gamma_p": 209.9,
        "a00": 1.03, "a01": 1.94, "a10": 1.92, "a11": 3.56,
        "a20": 3.16, "a21": 2.69, "a30": 7.73, "a31": 2.74,
    }

    def unpack_params(self, x):
        # Unpack only the GB variables, return empty dict for CICR
        gb = {'gamma_d': x[0], 'gamma_p': x[1], 'a00': x[2], 'a01': x[3],
              'a10': x[4], 'a11': x[5], 'a20': x[6], 'a21': x[7], 'a30': x[8], 'a31': x[9]}
        return gb, {}

    def get_debug_sim_fn(self):
        return _apply_gb_only_debug

    def get_step_factory(self):
        def scan_factory(params, syn_params):
            gb, _ = params
            c_pre, c_post, is_apical = syn_params
            
            # GB Thresholds
            theta_d = jnp.where(is_apical, gb['a20']*c_pre + gb['a21']*c_post, gb['a00']*c_pre + gb['a01']*c_post)
            theta_p = jnp.where(is_apical, gb['a30']*c_pre + gb['a31']*c_post, gb['a10']*c_pre + gb['a11']*c_post)
            
            cai_rest = 70e-6 
            tau_effca = 278.3177658387
            
            def scan_step(carry, inputs):
                # We unpack carry, but ignore all CICR states (index 0, 1, 2, 3, 4)
                _, _, _, _, _, effcai, rho = carry
                cai_raw, dt = inputs
                
                # 1. Graupner-Brunel Plasticity Integrator (NO CICR CALCIUM ADDED)
                decay_eff = jnp.exp(-dt / tau_effca)
                
                # Notice `ca_cicr_new` is completely removed from this calculation
                effcai_new = effcai * decay_eff + (cai_raw - cai_rest) * tau_effca * (1.0 - decay_eff)
                
                pot = jnp.where(effcai_new > theta_p, 1.0, 0.0)
                dep = jnp.where(effcai_new > theta_d, 1.0, 0.0)
                
                # JAX handles cubic terms cleanly.
                drho = (-rho*(1.0-rho)*(0.5-rho) + pot*gb['gamma_p']*(1.0-rho) - dep*gb['gamma_d']*rho) / 70000.0
                rho_new  = jnp.where(dt > 0, jnp.clip(rho + dt * drho, 0.0, 1.0), rho)
                
                # Return zeros for the CICR state variables to keep shapes consistent
                return (0.0, 0.0, 0.0, 0.0, 0.0, effcai_new, rho_new), None
            return scan_step
        return scan_factory

if __name__ == "__main__":
    GBOnlyModel().run()