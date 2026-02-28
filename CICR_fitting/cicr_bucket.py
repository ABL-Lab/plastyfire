#!/usr/bin/env python3
"""
Theme Park Bucket (Integrate-and-Fire) CICR variant (Continuous JAX).
Models the ER as a tipping bucket: integrates calcium to a threshold, 
triggers a saturated discrete dump, and decays with a long time constant.
"""

import numpy as np
import jax.numpy as jnp
from cicr_common import CICRModel

def _apply_bucket_cicr_debug(cai_1syn, t, cp, cq, is_apical, dp):
    """Plain-NumPy Bucket CICR debug simulation for diagnostic plotting."""
    
    # Bucket Parameters
    v_fill    = float(dp["v_fill"])
    tau_leak  = float(dp["tau_leak"])
    threshold = float(dp["threshold"])
    amp_dump  = float(dp["amp_dump"])
    tau_dump  = float(dp["tau_dump"])
    
    cai_rest  = 70e-6
    
    n = len(cai_1syn)
    cai_total   = np.copy(cai_1syn)
    bucket_out  = np.zeros(n)
    ca_cicr_out = np.zeros(n)

    bucket_level = 0.0
    ca_cicr = 0.0

    for i in range(n - 1):
        dt = t[i + 1] - t[i]
        if dt <= 0.0:
            cai_total[i + 1]   = cai_1syn[i + 1] + ca_cicr
            bucket_out[i + 1]  = bucket_level
            ca_cicr_out[i + 1] = ca_cicr
            continue

        cai_raw = cai_1syn[i]

        # 1. Fill the leaky bucket
        cai_evoked = max(0.0, cai_raw - cai_rest)
        bucket_level = bucket_level * np.exp(-dt / tau_leak) + (v_fill * cai_evoked) * dt
        
        # 2. Tipping the bucket
        is_tipping = 1.0 if bucket_level >= threshold else 0.0
        if is_tipping:
            bucket_level = 0.0 # Bucket empties
            
        # 3. The Splash (STRICT NON-ADDITIVE SATURATION)
        decayed_ca = ca_cicr * np.exp(-dt / tau_dump)
        splash_ca = is_tipping * amp_dump
        ca_cicr = max(decayed_ca, splash_ca) # Cannot exceed amp_dump

        cai_total[i + 1]   = cai_1syn[i + 1] + ca_cicr
        bucket_out[i + 1]  = bucket_level
        ca_cicr_out[i + 1] = ca_cicr

    return {
        "cai_total": cai_total,
        "priming":   bucket_out,   # Plotted on Row 3
        "ca_er":     np.zeros(n),  # Unused in macro-model
        "ca_cicr":   ca_cicr_out,  # Plotted on Row 4
    }


class BucketCICRModel(CICRModel):
    DESCRIPTION = "Theme Park Bucket (Integrate-and-Fire) CICR"
        
    FIT_PARAMS = [
            # --- GB PARAMS LOCKED (Tuned perfectly for LTD) ---
            # ("gamma_d", 50.0, 300.0), ("gamma_p", 150.0, 400.0),
            # ("a00", 0.01, 20.0), ("a01", 0.01, 20.0),
            # ("a10", 0.01, 20.0), ("a11", 0.01, 20.0),
            # ("a20", 0.01, 20.0), ("a21", 0.01, 20.0), 
            # ("a30", 0.01, 20.0), ("a31", 0.01, 20.0), 
            
            # --- CICR PARAMS ACTIVE (Tuning to act as an LTP Booster) ---
            ("v_fill", 1.0, 5000.0),      
            ("tau_leak", 50.0, 2000.0),   
            ("threshold", 0.1, 10.0),     
            # Widened amp_dump slightly so it can cross the huge theta_p = 1.27 ceiling
            ("amp_dump", 0.0001, 0.02),  
            ("tau_dump", 10.0, 5000.0),  
        ]
    
    DEFAULT_PARAMS = {
        # Locked GB Params (from your perfect 0.7922 LTD fit)
        "gamma_d": 286.3765, "gamma_p": 292.4128,
        "a00": 0.3552, "a01": 4.1632, "a10": 9.8053, "a11": 6.7515,
        "a20": 9.5501, "a21": 17.6940, "a30": 17.6844, "a31": 5.6505,
        
        # Initial CICR Bucket Params
        "v_fill": 1000.0, "tau_leak": 200.0, "threshold": 1.0, 
        "amp_dump": 0.005, "tau_dump": 100.0
    }

    def unpack_params(self, x):
        # 1. Start with everything locked to default (Perfect LTD GB params)
        p = dict(self.DEFAULT_PARAMS)
        
        # 2. Overwrite only the parameters the optimizer is currently fitting (The Bucket)
        for name, val in zip(self.PARAM_NAMES, x):
            p[name] = val
            
        # 3. Route them to their specific dictionaries
        gb = {k: p[k] for k in ['gamma_d', 'gamma_p', 'a00', 'a01', 'a10', 
                                'a11', 'a20', 'a21', 'a30', 'a31']}
        cicr = {k: p[k] for k in p if k not in gb}
        
        return gb, cicr

    def get_debug_sim_fn(self):
        return _apply_bucket_cicr_debug

    def get_step_factory(self):
        def scan_factory(params, syn_params):
            gb, cicr = params
            c_pre, c_post, is_apical = syn_params
            
            # GB Thresholds - Forcing theta_p to be biologically higher than theta_d
            theta_d = jnp.where(is_apical, gb['a20']*c_pre + gb['a21']*c_post, gb['a00']*c_pre + gb['a01']*c_post)
            raw_theta_p = jnp.where(is_apical, gb['a30']*c_pre + gb['a31']*c_post, gb['a10']*c_pre + gb['a11']*c_post)
            
            # STRICT ENFORCEMENT: Potentiation requires more calcium than depression
            theta_p = jnp.maximum(raw_theta_p, theta_d + 0.1)
            
            cai_rest = 70e-6 
            
            def scan_step(carry, inputs):
                bucket_level, _dummy1, ca_cicr, _dummy2, _dummy3, effcai, rho = carry
                cai_raw, dt = inputs
                
                # 1. Fill the bucket 
                cai_evoked = jnp.maximum(0.0, cai_raw - cai_rest)
                bucket_new = bucket_level * jnp.exp(-dt / cicr['tau_leak']) + (cicr['v_fill'] * cai_evoked) * dt
                
                # 2. Tipping the bucket
                is_tipping = jnp.where(bucket_new >= cicr['threshold'], 1.0, 0.0)
                bucket_new = jnp.where(is_tipping == 1.0, 0.0, bucket_new)
                
                # 3. The Splash (STRICT NON-ADDITIVE SATURATION)
                decayed_ca = ca_cicr * jnp.exp(-dt / cicr['tau_dump'])
                splash_ca = is_tipping * cicr['amp_dump']
                ca_cicr_new = jnp.maximum(decayed_ca, splash_ca)

                # 4. Graupner-Brunel Plasticity Integrator
                decay_eff = jnp.exp(-dt / 200.0)
                effcai_new = effcai * decay_eff + (cai_raw + ca_cicr_new - cai_rest) * 200.0 * (1.0 - decay_eff)
                
                dep = jnp.where(effcai_new > theta_d, 1.0, 0.0)
                pot = jnp.where(effcai_new > theta_p, 1.0, 0.0)
                
                # Equation: tau * d_rho/dt = -rho(1-rho)(0.5-rho) + gamma_p(1-rho)Theta_p - gamma_d*rho*Theta_d
                drho = (-rho*(1.0-rho)*(0.5-rho) + pot*gb['gamma_p']*(1.0-rho) - dep*gb['gamma_d']*rho) / 70000.0
                rho_new = jnp.clip(rho + dt * drho, 0.0, 1.0)
                
                return (bucket_new, 0.0, ca_cicr_new, 0.0, 0.0, effcai_new, rho_new), None
            return scan_step
        return scan_factory

if __name__ == "__main__":
    BucketCICRModel().run()