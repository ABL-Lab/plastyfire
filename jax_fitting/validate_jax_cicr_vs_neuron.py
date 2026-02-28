#!/usr/bin/env python3
"""Validate JAX replay against NEURON traces for Primed-Junctional CICR.

Tests GB+CICR JAX implementation against NEURON's trace outputs.
"""

import os, sys, pickle, argparse
import numpy as np
import pandas as pd
from pathlib import Path
from functools import partial

import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

DHURUVA_PARAMS_V3 = {
    "gamma_d": 50.4690,
    "gamma_p": 161.2494,
    "a00": 0.5169,
    "a01": 4.6079,
    "a10": 1.8423,
    "a11": 4.9230,
    "a20": 0.1282,
    "a21": 0.4070,
    "a30": 0.3265,
    "a31": 4.9275,
    "tau_prime": 436.0976,
    "k_prime": 44310.3178,
    "K_prime_limit": 0.0003,
    "K_trig": 0.0013,
    "n_trig": 5.0744,
    "Vmax_CICR": 0.0022,
    "tau_CICR": 313.2375,
    "tau_eff": 100.7488,
}

BASE_DIR = "/project/rrg-emuller/dhuruva/plastyfitting"
CAI_REST = 70e-6


def load_pkl(path):
    with open(path, "rb") as f:
        return pickle.load(f)

def get_thresholds(c_pre, c_post, is_apical, a):
    if is_apical:
        td = a["a20"] * c_pre + a["a21"] * c_post
        tp = a["a30"] * c_pre + a["a31"] * c_post
    else:
        td = a["a00"] * c_pre + a["a01"] * c_post
        tp = a["a10"] * c_pre + a["a11"] * c_post
    return td, tp

@partial(jax.jit, static_argnums=(4, ))
def compute_cicr_jax(cai, t, td, tp, cicr_params, rho0):
    """JAX implementation of Primed-Junctional CICR."""
    
    cicr_params = dict(cicr_params)
    
    dt_trace = jnp.diff(t, prepend=t[0])
    dt_trace = jnp.where(dt_trace <= 0, 1e-6, dt_trace)

    def scan_step(carry, inputs):
        prime_state, ca_cicr, effcai, rho = carry
        cai_raw, dt = inputs

        cai_evoked = jnp.maximum(0.0, cai_raw - CAI_REST)

        # 1. Soft-Clipped ER Priming
        cai_clip = cicr_params['K_prime_limit'] * jnp.tanh(cai_evoked / (cicr_params['K_prime_limit'] + 1e-12))
        
        rate_load = cicr_params['k_prime'] * cai_clip
        inv_tau_p = 1.0 / cicr_params['tau_prime']
        tau_eff_p = 1.0 / (rate_load + inv_tau_p + 1e-12)
        P_inf = rate_load * tau_eff_p
        P_new = P_inf + (prime_state - P_inf) * jnp.exp(-dt / tau_eff_p)

        # 2. Junctional Trigger
        trigger = (cai_evoked**cicr_params['n_trig']) / (cicr_params['K_trig']**cicr_params['n_trig'] + cai_evoked**cicr_params['n_trig'] + 1e-12)
        
        # 3. Synergistic CICR Release (uses P_new instead of prime_state to be stable internally)
        J_rate = cicr_params['Vmax_CICR'] * P_new * trigger

        ca_cicr_decay = jnp.exp(-dt / cicr_params['tau_CICR'])
        ca_cicr_new = ca_cicr * ca_cicr_decay + J_rate * cicr_params['tau_CICR'] * (1.0 - ca_cicr_decay)

        # 4. Chindemi c* Leaky Integrator
        ca_total = cai_evoked + ca_cicr_new
        eff_decay = jnp.exp(-dt / cicr_params['tau_eff'])
        effcai_new = effcai * eff_decay + ca_total * cicr_params['tau_eff'] * (1.0 - eff_decay) 

        # 5. Graupner-Brunel Plasticity Core
        pot = jnp.where(effcai_new > tp, 1.0, 0.0)
        dep = jnp.where(effcai_new > td, 1.0, 0.0)

        drho = (-rho*(1.0-rho)*(0.5-rho)
                + pot*cicr_params['gamma_p']*(1.0-rho)
                - dep*cicr_params['gamma_d']*rho) / 70000.0
        
        rho_new = jnp.where(dt > 0, jnp.clip(rho + dt * drho, 0.0, 1.0), rho)

        return (P_new, ca_cicr_new, effcai_new, rho_new), effcai_new

    init_carry = (0.0, 0.0, 0.0, rho0)
    (P_final, ca_cicr_final, effcai_final, rho_final), effcai_history = jax.lax.scan(scan_step, init_carry, (cai, dt_trace))
    
    return jnp.max(effcai_history), rho_final
    
    return jnp.max(effcai_history), rho_final


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trace-files", nargs='+', default=[
        os.path.join(BASE_DIR, "trace_results/DHURUVA_PARAMS_V3/180164-197248/10Hz_-10ms/simulation_traces.pkl"),
        os.path.join(BASE_DIR, "trace_results/DHURUVA_PARAMS_V3/180164-197248/10Hz_5ms/simulation_traces.pkl"),
        os.path.join(BASE_DIR, "trace_results/DHURUVA_PARAMS_V3/180164-197248/10Hz_10ms/simulation_traces.pkl"),
    ])
    args = parser.parse_args()

    cicr_params_tuple = tuple(DHURUVA_PARAMS_V3.items())
    cicr_params = DHURUVA_PARAMS_V3

    for pkl_path in args.trace_files:
        path_obj = Path(pkl_path)
        if not path_obj.exists():
            print(f"ERROR: {pkl_path} not found")
            sys.exit(1)

    print(f"\n{'=' * 80}")
    print(f"  Validating Primed-Junctional CICR (DHURUVA_PARAMS_V3) step sizes")
    print(f"{'=' * 80}")
    
    dt_test_list = [1.0, 10.0, 100.0, 1000.0]
    results = {dt: {'rho_diffs': [], 'eff_ratios': []} for dt in dt_test_list}

    for pkl_idx, pkl_path in enumerate(args.trace_files):
        data = load_pkl(str(pkl_path))
        t_list = data["t"]
        t = jnp.asarray(t_list, dtype=jnp.float64)

        cai_raw = np.asarray(data["cai_CR"], dtype=np.float64)
        if cai_raw.ndim == 1: cai_raw = cai_raw.reshape(1, -1)
        elif cai_raw.shape[0] == len(t) and cai_raw.shape[1] != len(t): cai_raw = cai_raw.T
        
        rho_gb = np.asarray(data.get("rho_GB", np.zeros_like(cai_raw)), dtype=np.float64)
        if rho_gb.ndim == 1: rho_gb = rho_gb.reshape(1, -1)
        elif rho_gb.shape[0] == len(t) and rho_gb.shape[1] != len(t): rho_gb = rho_gb.T

        effcai_gb = np.asarray(data.get("effcai_GB", np.zeros_like(cai_raw)), dtype=np.float64)
        if effcai_gb.ndim == 1: effcai_gb = effcai_gb.reshape(1, -1)
        elif effcai_gb.shape[0] == len(t) and effcai_gb.shape[1] != len(t): effcai_gb = effcai_gb.T

        sp = data.get("synprop", {})
        cpre = np.asarray(sp.get("Cpre", []), dtype=np.float64)
        cpost = np.asarray(sp.get("Cpost", []), dtype=np.float64)
        locs = sp.get("loc", ["basal"] * cai_raw.shape[0])
        is_api = [loc == "apical" for loc in locs]
        n_syn = cai_raw.shape[0]

        for s in range(min(5, n_syn)):
            td, tp = get_thresholds(cpre[s], cpost[s], is_api[s], DHURUVA_PARAMS_V3)
            tp = max(tp, td + 0.01)

            rho0 = rho_gb[s, 0]
            rho_n = rho_gb[s, -1]
            eff_n = effcai_gb[s].max() if len(effcai_gb[s]) > 0 else 0.0
            
            for test_dt in dt_test_list:
                # Interpolate trace and time
                t_interp = np.arange(t[0], t[-1], test_dt)
                if len(t_interp) < 2:
                    continue # Too coarse
                cai_syn_interp = jnp.interp(t_interp, t, cai_raw[s])
                
                eff_j, rho_j = compute_cicr_jax(cai_syn_interp, t_interp, td, tp, frozenset(cicr_params.items()), rho0)
                
                eff_j = float(eff_j)
                rho_j = float(rho_j)

                ratio = eff_j / eff_n if eff_n > 1e-9 else 1.0
                diff = abs(rho_j - rho_n)
                
                if eff_n > 1e-3:
                    results[test_dt]['eff_ratios'].append(ratio)
                results[test_dt]['rho_diffs'].append(diff)

    print("\n  Summary of Interpolation Timesteps Tradeoff:")
    print(f"  {'-'*60}")
    print(f"  {'dt (ms)':<10} | {'eff_ratio (JAX/NEURON)':<25} | {'Mean |rho_diff|':<15}")
    print(f"  {'-'*60}")
    
    for test_dt in dt_test_list:
        if len(results[test_dt]['eff_ratios']) == 0:
            continue
        mean_eff = np.mean(results[test_dt]['eff_ratios'])
        std_eff = np.std(results[test_dt]['eff_ratios'])
        mean_rho = np.mean(results[test_dt]['rho_diffs'])
        
        print(f"  {test_dt:<10.3f} | {mean_eff:.4f} +/- {std_eff:.4f}         | {mean_rho:.4f}")

if __name__ == "__main__":
    main()
