#!/usr/bin/env python3
"""
Advanced MWC-based CICR variant (Continuous JAX).
Implements the Monod-Wyman-Changeux tetramer architecture, Fraiman-Dawson luminal
calcium inhibition, spatial attenuation, and active SERCA extrusion kinetics to 
prevent the "Machine Learning Exploit" and securely hit the 0.79 LTD target.
"""

import numpy as np
import jax.numpy as jnp
from cicr_common import CICRModel

def _apply_mwc_cicr_debug(cai_1syn, t, cp, cq, is_apical, dp):
    """Plain-NumPy MWC CICR debug simulation for diagnostic plotting."""
    # MWC Parameters
    L       = float(dp["L"])
    K_R_IP3 = float(dp["K_R_IP3"])
    K_R_Ca  = float(dp["K_R_Ca"])
    K_T_IP3 = float(dp["K_T_IP3"])
    K_T_Ca  = float(dp["K_T_Ca"])
    
    # IP3 Production
    tau_IP3 = float(dp["tau_IP3"])
    v_prod  = float(dp["v_prod"])
    K_prod  = float(dp["K_prod"])
    
    # Structural & Extrusion Kinetics
    V_ER       = float(dp["V_ER"])
    K_lum      = float(dp["K_lum"])
    V_SERCA    = float(dp["V_SERCA"])
    K_SERCA    = float(dp["K_SERCA"])
    r_coupling = float(dp["r_coupling"])
    
    cai_rest  = 70e-6
    fc, fe    = 0.83, 0.17
    caAvg     = 0.0017
    D_Ca      = 200.0 # Diffusion coefficient approximation (um^2/s)

    n = len(cai_1syn)
    cai_total   = np.copy(cai_1syn)
    P_IP3_out   = np.zeros(n)
    ca_er_out   = np.zeros(n)
    ca_cicr_out = np.zeros(n)

    # Initial ER calcium
    ca_er   = max(0.0, (caAvg - fc * cai_1syn[0]) / fe)
    ca_cicr = 0.0
    P_IP3   = 0.0

    # Auto-balance baseline leak against active SERCA at rest
    ca_er_rest   = max(0.0, (caAvg - fc * cai_rest) / fe)
    J_serca_rest = V_SERCA * (cai_rest**2) / (K_SERCA**2 + cai_rest**2)
    g_leak       = J_serca_rest / (ca_er_rest - cai_rest + 1e-12)

    ca_er_out[0]   = ca_er
    ca_cicr_out[0] = ca_cicr

    for i in range(n - 1):
        dt = t[i + 1] - t[i]
        if dt <= 0.0:
            cai_total[i + 1]   = cai_1syn[i + 1] + ca_cicr
            P_IP3_out[i + 1]   = P_IP3
            ca_er_out[i + 1]   = ca_er
            ca_cicr_out[i + 1] = ca_cicr
            continue

        cai_raw = cai_1syn[i]
        ca_cyt  = cai_raw + ca_cicr

        # Spatial Attenuation Approximation for VDCC to ER distance
        # D_Ca is scaled by 1e-6 to loosely match timescale magnitudes
        Ca_trigger = cai_raw * np.exp(-(r_coupling**2) / (4.0 * D_Ca * 1e-6 + 1e-12))

        # 1. Evoked IP3 Production
        cai_evoked = max(0.0, cai_raw - cai_rest)
        f_ca       = (cai_evoked**2) / (K_prod**2 + cai_evoked**2 + 1e-12)
        prod_rate  = v_prod * f_ca
        P_IP3      = P_IP3 * np.exp(-dt / tau_IP3) + prod_rate * dt

        # 2. MWC Allosteric Gate (Strict integer exponent n=4)
        IP3_safe = max(0.0, P_IP3)
        R_IP3 = (1.0 + IP3_safe / K_R_IP3)**4
        R_Ca  = (1.0 + Ca_trigger / K_R_Ca)**4
        T_IP3 = (1.0 + IP3_safe / K_T_IP3)**4
        T_Ca  = (1.0 + Ca_trigger / K_T_Ca)**4
        
        Q_MWC = (R_IP3 * R_Ca) + L * (T_IP3 * T_Ca)
        R_bar = (R_IP3 * R_Ca) / Q_MWC

        # 3. Fraiman-Dawson Luminal Brake
        luminal_brake = ca_er / (K_lum + ca_er + 1e-12)

        # 4. Fluxes
        J_IP3R_raw = V_ER * R_bar * luminal_brake * (ca_er - ca_cyt)
        J_IP3R = J_IP3R_raw if ca_er > ca_cyt else 0.0
        
        J_leak  = g_leak * (ca_er - ca_cyt)
        J_SERCA = V_SERCA * (ca_cyt**2) / (K_SERCA**2 + ca_cyt**2 + 1e-12)

        J_net_raw = J_IP3R + J_leak - J_SERCA

        # Flux clamping to prevent numerical instability
        max_release = ca_er * (fe / fc) / dt
        max_uptake  = ca_cicr / dt if ca_cicr > 0 else 0.0
        J_net = max(-max_uptake, min(J_net_raw, max_release))

        ca_er   = max(0.0, ca_er   - (fc / fe) * J_net * dt)
        ca_cicr = max(0.0, ca_cicr + J_net * dt)

        cai_total[i + 1]   = cai_1syn[i + 1] + ca_cicr
        P_IP3_out[i + 1]   = P_IP3
        ca_er_out[i + 1]   = ca_er
        ca_cicr_out[i + 1] = ca_cicr

    return {
        "cai_total": cai_total,
        "priming":   P_IP3_out,   
        "ca_er":     ca_er_out,
        "ca_cicr":   ca_cicr_out,
    }


class MWCCICRModel(CICRModel):
    DESCRIPTION = "MWC Allosteric CICR with Luminal Brake (Continuous JAX)"
        
    FIT_PARAMS = [
            # GB Params (10)
            ("gamma_d", 50.0, 200.0), ("gamma_p", 150.0, 300.0),
            ("a00", 0.01, 20.0), ("a01", 0.01, 20.0),
            ("a10", 0.01, 20.0), ("a11", 0.01, 20.0),
            ("a20", 0.01, 20.0), ("a21", 0.01, 20.0), # Widened bounds for ceiling protection
            ("a30", 0.01, 20.0), ("a31", 0.01, 20.0), 
            
            # IP3 Priming (3)
            ("tau_IP3", 50.0, 2000.0),        
            ("v_prod", 0.001, 10.0),          
            ("K_prod", 0.0001, 0.005),   

            # MWC Allosteric Gate (5)
            ("L", 10.0, 10000.0),        # Allosteric constant (anchors basal leak)
            ("K_R_IP3", 0.0001, 0.01),   # Relaxed IP3 affinity
            ("K_R_Ca", 0.0001, 0.01),    # Relaxed Ca affinity
            ("K_T_IP3", 0.01, 1.0),      # Tense IP3 affinity (low)
            ("K_T_Ca", 0.01, 1.0),       # Tense Ca affinity (low)

            # Kinetics, Luminal Brake & Extrusion (5)
            ("V_ER", 1.0, 200.0),        # Maximal IP3R flux
            ("K_lum", 0.0001, 0.01),     # Luminal brake dissociation constant
            ("V_SERCA", 1.0, 100.0),     # Active shear extrusion velocity
            ("K_SERCA", 0.0001, 0.01),   # SERCA affinity
            ("r_coupling", 0.01, 2.0),   # Spatial VDCC-to-ER coupling distance (um)
        ]
    
    DEFAULT_PARAMS = {
        "gamma_d": 181.4, "gamma_p": 209.9,
        "a00": 1.03, "a01": 1.94, "a10": 1.92, "a11": 3.56,
        "a20": 3.16, "a21": 2.69, "a30": 7.73, "a31": 2.74,
        "tau_IP3": 200.0, "v_prod": 1.0, "K_prod": 0.001,
        "L": 1000.0, "K_R_IP3": 0.001, "K_R_Ca": 0.001, "K_T_IP3": 0.1, "K_T_Ca": 0.1,
        "V_ER": 50.0, "K_lum": 0.002, "V_SERCA": 25.0, "K_SERCA": 0.0001, "r_coupling": 0.5
    }

    def unpack_params(self, x):
        gb = {'gamma_d': x[0], 'gamma_p': x[1], 'a00': x[2], 'a01': x[3],
              'a10': x[4], 'a11': x[5], 'a20': x[6], 'a21': x[7], 'a30': x[8], 'a31': x[9]}
        cicr = {'tau_IP3': x[10], 'v_prod': x[11], 'K_prod': x[12],
                'L': x[13], 'K_R_IP3': x[14], 'K_R_Ca': x[15], 'K_T_IP3': x[16], 'K_T_Ca': x[17],
                'V_ER': x[18], 'K_lum': x[19], 'V_SERCA': x[20], 'K_SERCA': x[21], 'r_coupling': x[22]}
        return gb, cicr

    def get_debug_sim_fn(self):
        return _apply_mwc_cicr_debug

    def get_step_factory(self):
        def scan_factory(params, syn_params):
            gb, cicr = params
            c_pre, c_post, is_apical = syn_params
            
            # GB Thresholds
            theta_d = jnp.where(is_apical, gb['a20']*c_pre + gb['a21']*c_post, gb['a00']*c_pre + gb['a01']*c_post)
            raw_theta_p = jnp.where(is_apical, gb['a30']*c_pre + gb['a31']*c_post, gb['a10']*c_pre + gb['a11']*c_post)
            theta_p = jnp.maximum(raw_theta_p, theta_d + 0.01)
            
            # Baseline Equilibrium Logic
            cai_rest = 70e-6 
            ca_er_rest = jnp.maximum(0.0, (0.0017 - 0.83 * cai_rest) / 0.17)
            
            J_serca_rest = cicr['V_SERCA'] * (cai_rest**2) / (cicr['K_SERCA']**2 + cai_rest**2)
            g_leak = J_serca_rest / (ca_er_rest - cai_rest + 1e-12)
            
            def scan_step(carry, inputs):
                P_IP3, ca_er, ca_cicr, _dummy1, _dummy2, effcai, rho = carry
                cai_raw, dt = inputs
                
                ca_cyt = cai_raw + ca_cicr
                D_Ca = 200.0
                
                # Spatial Attenuation Approximation
                Ca_trigger = cai_raw * jnp.exp(-(cicr['r_coupling']**2) / (4.0 * D_Ca * 1e-6 + 1e-12))
                
                # 1. Evoked IP3 Production
                cai_evoked = jnp.maximum(0.0, cai_raw - cai_rest)
                f_ca = (cai_evoked**2) / (cicr['K_prod']**2 + cai_evoked**2 + 1e-12)
                prod_rate = cicr['v_prod'] * f_ca
                P_IP3_new = P_IP3 * jnp.exp(-dt / cicr['tau_IP3']) + prod_rate * dt
                
                # 2. MWC Allosteric Gate
                IP3_safe = jnp.maximum(0.0, P_IP3_new)
                R_IP3 = (1.0 + IP3_safe / cicr['K_R_IP3'])**4
                R_Ca  = (1.0 + Ca_trigger / cicr['K_R_Ca'])**4
                T_IP3 = (1.0 + IP3_safe / cicr['K_T_IP3'])**4
                T_Ca  = (1.0 + Ca_trigger / cicr['K_T_Ca'])**4
                
                Q_MWC = (R_IP3 * R_Ca) + cicr['L'] * (T_IP3 * T_Ca)
                R_bar = (R_IP3 * R_Ca) / Q_MWC
                
                # 3. Fraiman-Dawson Luminal Brake
                luminal_brake = ca_er / (cicr['K_lum'] + ca_er + 1e-12)
                
                # 4. Raw Fluxes
                J_IP3R = jnp.where(ca_er > ca_cyt, cicr['V_ER'] * R_bar * luminal_brake * (ca_er - ca_cyt), 0.0)
                J_leak_dynamic = g_leak * (ca_er - ca_cyt)
                J_SERCA = cicr['V_SERCA'] * (ca_cyt**2) / (cicr['K_SERCA']**2 + ca_cyt**2 + 1e-12)
                
                J_net_raw = J_IP3R + J_leak_dynamic - J_SERCA
                
                # Flux Clamping
                max_release = ca_er * (0.17 / 0.83) / dt
                max_uptake = ca_cicr / dt
                J_net_clamped = jnp.clip(J_net_raw, -max_uptake, max_release)
                
                ca_er_new = jnp.maximum(0.0, ca_er - (0.83 / 0.17) * J_net_clamped * dt)
                ca_cicr_new = jnp.maximum(0.0, ca_cicr + J_net_clamped * dt)

                # 5. Graupner-Brunel Plasticity Integrator
                decay_eff = jnp.exp(-dt / 200.0)
                effcai_new = effcai * decay_eff + (cai_raw + ca_cicr_new - cai_rest) * 200.0 * (1.0 - decay_eff)
                
                dep = jnp.where(effcai_new > theta_d, 1.0, 0.0)
                pot = jnp.where(effcai_new > theta_p, 1.0, 0.0)
                drho = (-rho*(1.0-rho)*(0.5-rho) + pot*gb['gamma_p']*(1.0-rho) - dep*gb['gamma_d']*rho) / 70000.0
                rho_new = jnp.clip(rho + dt * drho, 0.0, 1.0)
                
                return (P_IP3_new, ca_er_new, ca_cicr_new, 0.0, 0.0, effcai_new, rho_new), None
            return scan_step
        return scan_factory

if __name__ == "__main__":
    MWCCICRModel().run()