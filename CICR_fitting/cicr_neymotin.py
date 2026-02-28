#!/usr/bin/env python3
"""
Neymotin-Li-Rinzel CICR variant (Continuous JAX).
Uses the biologically rigorous tetramer IP3R model (m^3 * n^3 * h^3).
Automatically balances ER leak against SERCA at rest and features 
strict mass-conservation bounds to prevent Euler integration overflow.
"""

import numpy as np
import jax.numpy as jnp
from cicr_common import CICRModel


def _apply_neymotin_cicr_debug(cai_1syn, t, cp, cq, is_apical, dp):
    """Plain-NumPy Neymotin CICR debug simulation for diagnostic plotting."""
    tau_IP3   = float(dp["tau_IP3"])
    v_prod    = float(dp["v_prod"])
    K_prod    = float(dp["K_prod"])
    k_IP3     = float(dp["k_IP3"])
    k_act     = float(dp["k_act"])
    k_inh     = float(dp["k_inh"])
    tau_h     = float(dp["tau_h"])
    g_release = float(dp["g_release"])
    g_serca   = float(dp["g_serca"])

    cai_rest  = 70e-6
    K_serca   = 0.0001   # 0.1 µM
    fc, fe    = 0.83, 0.17
    caAvg     = 0.0017

    n = len(cai_1syn)
    cai_total   = np.copy(cai_1syn)
    P_IP3_out   = np.zeros(n)
    ca_er_out   = np.zeros(n)
    ca_cicr_out = np.zeros(n)

    # Initial ER calcium
    ca_er   = max(0.0, (caAvg - fc * cai_1syn[0]) / fe)
    ca_cicr = 0.0

    # Auto-balance ER leak against SERCA at rest
    ca_er_rest   = max(0.0, (caAvg - fc * cai_rest) / fe)
    J_serca_rest = g_serca * (cai_rest**2) / (K_serca**2 + cai_rest**2)
    g_leak       = J_serca_rest / (ca_er_rest - cai_rest + 1e-12)

    P_IP3  = 0.0
    h_gate = 0.0  

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

        # 1. Evoked IP3 Production
        cai_evoked = max(0.0, cai_raw - cai_rest)
        f_ca       = (cai_evoked**2) / (K_prod**2 + cai_evoked**2 + 1e-12)
        prod_rate  = v_prod * f_ca
        P_IP3      = P_IP3 * np.exp(-dt / tau_IP3) + prod_rate * dt

        # 2. Fast Activation Gates (m, n)
        P_safe = max(0.0, P_IP3)
        m_inf  = P_safe / (P_safe + k_IP3 + 1e-12)
        n_inf  = ca_cyt / (ca_cyt + k_act + 1e-12)

        # 3. Slow Inactivation Gate (h)
        h_inf = k_inh / (k_inh + ca_cyt + 1e-12)
        if h_gate == 0.0:
            h_gate = h_inf
        decay_h = np.exp(-dt / tau_h)
        h_gate  = h_gate * decay_h + h_inf * (1.0 - decay_h)

        # 4. Open Probability (Tetramer: m³ · n³ · h³)
        P_open = (m_inf**3) * (n_inf**3) * (h_gate**3)

        # 5. Fluxes & Absolute Mass-Conservation Clamping
        J_IP3R_raw = g_release * P_open * (ca_er - ca_cyt) if ca_er > ca_cyt else 0.0
        
        max_release = ca_er * (fe / fc) / dt
        J_IP3R_clamped = min(J_IP3R_raw, max_release)

        J_SERCA_raw = g_serca * (ca_cyt**2) / (K_serca**2 + ca_cyt**2 + 1e-12)
        J_SERCA_active = max(0.0, J_SERCA_raw - J_serca_rest)
        
        max_uptake = ca_cicr / dt if ca_cicr > 0 else 0.0
        J_SERCA_clamped = min(J_SERCA_active, max_uptake)

        J_leak_dynamic = g_leak * (ca_er - ca_cyt)

        # Total unconstrained net flux
        J_net_total = J_IP3R_clamped + J_leak_dynamic - J_SERCA_clamped

        # STRICT BOUND: Never move more than the compartments actually hold
        J_net_safe = max(-max_uptake, min(J_net_total, max_release))

        ca_er   = max(0.0, ca_er   - (fc / fe) * J_net_safe * dt)
        ca_cicr = max(0.0, ca_cicr + J_net_safe * dt)

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


class NeymotinCICRModel(CICRModel):
    DESCRIPTION = "Neymotin IP3R CICR (Continuous JAX)"
        
    FIT_PARAMS = [
            # # GB Params (10)
            # ("gamma_d", 50.0, 200.0), ("gamma_p", 150.0, 300.0),
            # ("a00", 0.01, 5.0), ("a01", 0.01, 5.0),
            # ("a10", 0.01, 5.0), ("a11", 0.01, 5.0),
            # ("a20", 0.01, 5.0), ("a21", 0.01, 5.0),
            # ("a30", 0.01, 5.0), ("a31", 0.01, 5.0),
            
            # Neymotin CICR Params (9)
            ("tau_IP3", 50.0, 500.0),         
            ("v_prod", 0.001, 10.0),          
            ("K_prod", 0.0001, 0.005),        
            ("k_IP3", 0.0001, 0.005),         
            ("k_act", 0.0001, 0.005),         
            ("k_inh", 0.0001, 0.005),         
            ("tau_h", 500.0, 3000.0),         
            ("g_release", 0.1, 100.0),       
            ("g_serca", 5.0, 100.0),          
        ]
    
    DEFAULT_PARAMS = {
        "gamma_d": 152.5869, "gamma_p": 215.7528,
        "a00": 0.2136, "a01": 3.3134, "a10": 0.8019, "a11": 2.8953,
        "a20": 19.2455, "a21": 13.9981, "a30": 7.7053, "a31": 1.0006,
        "tau_IP3": 1120.3375, "v_prod": 1.5286, "K_prod": 0.0026,
        "k_IP3": 0.0042, "k_act": 0.0012, "k_inh": 0.0018, 
        "tau_h": 1100.5462, "g_release": 35.6268, "g_serca": 56.4306,
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

    def get_debug_sim_fn(self):
        return _apply_neymotin_cicr_debug

    def get_step_factory(self):
        def scan_factory(params, syn_params):
            gb, cicr = params
            c_pre, c_post, is_apical = syn_params
            
            theta_d = jnp.where(is_apical, gb['a20']*c_pre + gb['a21']*c_post, gb['a00']*c_pre + gb['a01']*c_post)
            raw_theta_p = jnp.where(is_apical, gb['a30']*c_pre + gb['a31']*c_post, gb['a10']*c_pre + gb['a11']*c_post)
            theta_p = jnp.maximum(raw_theta_p, theta_d + 0.01)
            
            cai_rest = 70e-6 
            ca_er_rest = jnp.maximum(0.0, (0.0017 - 0.83 * cai_rest) / 0.17)
            K_serca = 0.0001 
            
            J_serca_rest = cicr['g_serca'] * (cai_rest**2) / (K_serca**2 + cai_rest**2)
            g_leak = J_serca_rest / (ca_er_rest - cai_rest + 1e-12)
            
            def scan_step(carry, inputs):
                P_IP3, ca_er, ca_cicr, h_gate, _dummy, effcai, rho = carry
                cai_raw, dt = inputs
                
                ca_cyt = cai_raw + ca_cicr
                
                # 1. Evoked IP3 Production
                cai_evoked = jnp.maximum(0.0, cai_raw - cai_rest)
                f_ca = (cai_evoked**2) / (cicr['K_prod']**2 + cai_evoked**2 + 1e-12)
                prod_rate = cicr['v_prod'] * f_ca
                P_IP3_new = P_IP3 * jnp.exp(-dt / cicr['tau_IP3']) + prod_rate * dt
                
                # 2. Fast Activation Gates (m, n)
                P_safe = jnp.maximum(0.0, P_IP3_new)
                m_inf = P_safe / (P_safe + cicr['k_IP3'] + 1e-12)
                n_inf = ca_cyt / (ca_cyt + cicr['k_act'] + 1e-12)
                
                # 3. Slow Inactivation Gate (h)
                h_inf = cicr['k_inh'] / (cicr['k_inh'] + ca_cyt + 1e-12)
                h_gate_safe = jnp.where(h_gate == 0.0, h_inf, h_gate)
                
                decay_h = jnp.exp(-dt / cicr['tau_h'])
                h_new = h_gate_safe * decay_h + h_inf * (1.0 - decay_h)
                
                # 4. Open Probability 
                P_open = (m_inf**3) * (n_inf**3) * (h_new**3)
                
                # 5. Fluxes & Absolute Mass-Conservation Clamping
                J_IP3R_raw = jnp.where(ca_er > ca_cyt, cicr['g_release'] * P_open * (ca_er - ca_cyt), 0.0)
                
                max_release = ca_er * (0.17 / 0.83) / dt
                J_IP3R_clamped = jnp.minimum(J_IP3R_raw, max_release)
                
                J_SERCA_raw = cicr['g_serca'] * (ca_cyt**2) / (K_serca**2 + ca_cyt**2 + 1e-12)
                J_SERCA_active = jnp.maximum(0.0, J_SERCA_raw - J_serca_rest)
                
                max_uptake = ca_cicr / dt
                J_SERCA_clamped = jnp.minimum(J_SERCA_active, max_uptake)
                
                J_leak_dynamic = g_leak * (ca_er - ca_cyt)
                
                # Total unconstrained net flux
                J_net_total = J_IP3R_clamped + J_leak_dynamic - J_SERCA_clamped
                
                # STRICT BOUND: Never move more than the compartments actually hold
                J_net_safe = jnp.clip(J_net_total, -max_uptake, max_release)
                
                ca_er_new = jnp.maximum(0.0, ca_er - (0.83 / 0.17) * J_net_safe * dt)
                ca_cicr_new = jnp.maximum(0.0, ca_cicr + J_net_safe * dt)

                # 6. Graupner-Brunel Plasticity Integrator with Potentiation Veto
                decay_eff = jnp.exp(-dt / 200.0)
                effcai_new = effcai * decay_eff + (cai_raw + ca_cicr_new - cai_rest) * 200.0 * (1.0 - decay_eff)
                
                pot = jnp.where(effcai_new > theta_p, 1.0, 0.0)
                dep = jnp.where((effcai_new > theta_d) & (pot == 0.0), 1.0, 0.0)
                
                drho = (-rho*(1.0-rho)*(0.5-rho) + pot*gb['gamma_p']*(1.0-rho) - dep*gb['gamma_d']*rho) / 70000.0
                rho_new = jnp.clip(rho + dt * drho, 0.0, 1.0)
                
                return (P_IP3_new, ca_er_new, ca_cicr_new, h_new, 0.0, effcai_new, rho_new), None
            return scan_step
        return scan_factory

if __name__ == "__main__":
    NeymotinCICRModel().run()