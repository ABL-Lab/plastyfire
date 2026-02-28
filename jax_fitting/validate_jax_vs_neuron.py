#!/usr/bin/env python3
"""Validate JAX replay against NEURON traces (Chindemi params).

Tests GB-only with JAX to match NEURON's fine temporal resolution.
"""

import os, sys, pickle, argparse
import numpy as np
import pandas as pd
from pathlib import Path
from functools import partial

import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

CHINDEMI_PARAMS = {
    "gamma_d": 101.5, "gamma_p": 216.2,
    "a00": 1.002, "a01": 1.954, "a10": 1.159, "a11": 2.483,
    "a20": 1.127, "a21": 2.456, "a30": 5.236, "a31": 1.782,
}

BASE_DIR = "/project/rrg-emuller/dhuruva/plastyfitting"
TRACE_DIR = os.path.join(BASE_DIR, "trace_results/Chindemi_params")
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

@partial(jax.jit, static_argnums=(5, 6, 7))
def compute_gb_jax(cai, t, td, tp, rho0, tau, gd, gp):
    """GB plasticity matched to current JAX pipelines."""
    inv_tau = 1.0 / (70.0 * 1000.0) # tau_ind_GB is 70s
    rho_star = 0.5
    
    dt_trace = jnp.diff(t, prepend=t[0])
    dt_trace = jnp.where(dt_trace <= 0, 1e-6, dt_trace)

    def scan_step(carry, inputs):
        effcai, rho = carry
        cai_raw, dt = inputs
        
        decay_eff = jnp.exp(-dt / tau)
        
        # JAX uses cai_raw (which is current cai[i])
        effcai_new = effcai * decay_eff + (cai_raw - CAI_REST) * tau * (1.0 - decay_eff)
        
        pot = jnp.where(effcai_new > tp, 1.0, 0.0)
        dep = jnp.where(effcai_new > td, 1.0, 0.0)
        
        drho = (-rho*(1.0-rho)*(rho_star-rho) + pot*gp*(1.0-rho) - dep*gd*rho) * inv_tau
        rho_new  = jnp.where(dt > 0, jnp.clip(rho + dt * drho, 0.0, 1.0), rho)
        
        return (effcai_new, rho_new), effcai_new

    init_carry = (0.0, rho0)
    (effcai_final, rho_final), effcai_history = jax.lax.scan(scan_step, init_carry, (cai, dt_trace))
    
    return jnp.max(effcai_history), rho_final


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-pairs", type=int, default=3)
    parser.add_argument("--protocol", type=str, default="10Hz_10ms")
    parser.add_argument("--max-dt", type=float, default=1.0,
                        help="Interpolation step for CICR sub-stepping (ms)")
    args = parser.parse_args()

    trace_path = Path(TRACE_DIR)
    if not trace_path.exists():
        print(f"ERROR: {TRACE_DIR} not found")
        sys.exit(1)

    pair_dirs = sorted(d for d in trace_path.iterdir() if d.is_dir())[:args.max_pairs]
    gd, gp = CHINDEMI_PARAMS["gamma_d"], CHINDEMI_PARAMS["gamma_p"]

    configs = [
        ("GB-only tau=278.3 (JAX)", 278.3177658387),
    ]

    for label, tau in configs:
        print(f"\n{'=' * 60}")
        print(f"  {label}")
        print(f"{'=' * 60}")

        all_rj, all_rn, all_er = [], [], []

        for i_p, pair_dir in enumerate(pair_dirs):
            pkl_path = pair_dir / args.protocol / "simulation_traces.pkl"
            if not pkl_path.exists():
                continue

            data = load_pkl(str(pkl_path))
            t_list = data["t"]
            t = jnp.asarray(t_list, dtype=jnp.float64)

            cai = np.asarray(data["cai_CR"], dtype=np.float64)
            if cai.ndim == 1:
                cai = cai.reshape(1, -1)
            elif cai.shape[0] == len(t) and cai.shape[1] != len(t):
                cai = cai.T
            cai = jnp.asarray(cai)

            rho_gb = np.asarray(data.get("rho_GB", np.zeros_like(cai)), dtype=np.float64)
            if rho_gb.ndim == 1:
                rho_gb = rho_gb.reshape(1, -1)
            elif rho_gb.shape[0] == len(t) and rho_gb.shape[1] != len(t):
                rho_gb = rho_gb.T

            effcai_gb = np.asarray(data.get("effcai_GB", np.zeros_like(cai)), dtype=np.float64)
            if effcai_gb.ndim == 1:
                effcai_gb = effcai_gb.reshape(1, -1)
            elif effcai_gb.shape[0] == len(t) and effcai_gb.shape[1] != len(t):
                effcai_gb = effcai_gb.T

            sp = data.get("synprop", {})
            cpre = np.asarray(sp.get("Cpre", []), dtype=np.float64)
            cpost = np.asarray(sp.get("Cpost", []), dtype=np.float64)
            locs = sp.get("loc", ["basal"] * cai.shape[0])
            is_api = [loc == "apical" for loc in locs]
            n_syn = cai.shape[0]

            print(f"\n  {pair_dir.name} ({n_syn} syn)")
            for s in range(min(5, n_syn)):
                td, tp = get_thresholds(cpre[s], cpost[s], is_api[s], CHINDEMI_PARAMS)
                rho0 = rho_gb[s, 0]
                rho_n = rho_gb[s, -1]
                eff_n = effcai_gb[s].max()

                eff_j, rho_j = compute_gb_jax(
                    cai[s], t, td, tp, rho0, tau, gd, gp)
                
                # Convert JAX back to regular scalar
                eff_j = float(eff_j)
                rho_j = float(rho_j)

                ratio = eff_j / eff_n if eff_n > 1e-9 else 1.0
                all_rj.append(rho_j)
                all_rn.append(rho_n)
                all_er.append(ratio)

                tag = "A" if is_api[s] else "B"
                print(f"    s{s}({tag}) rho:{rho0:.1f}->N{rho_n:.2f}/J{rho_j:.2f}  eff:N{eff_n:.4f}/J{eff_j:.4f} r={ratio:.3f}")
            
            # Now compute EPSP ratio explicitly for the pair
            rho_j_config = []
            rho_n_config = []
            for s in range(n_syn):
                rho_j_config.append(1 if all_rj[len(all_rj) - n_syn + s] >= 0.5 else 0)
                rho_n_config.append(1 if all_rn[len(all_rn) - n_syn + s] >= 0.5 else 0)
            
            print(f"    EPSP active syn JAX: {rho_j_config} | NEURON: {rho_n_config}")
            
            # Predict EPSP Ratio if basis exists
            try:
                pre_gid, post_gid = pair_dir.name.split("-")
                basis_csv = Path(BASE_DIR) / "basis_results" / f"basis_{pre_gid}_{post_gid}.csv"
                if basis_csv.exists():
                    df = pd.read_csv(basis_csv)
                    df["config_list"] = df["config"].apply(lambda x: [int(i) for i in x.split(",")])
                    baseline_row = df[df["config_list"].apply(lambda x: sum(x) == 0)]
                    baseline_mean = baseline_row["mean"].values[0]
                    
                    singleton_means = []
                    for i in range(n_syn):
                        row = df[df["config_list"].apply(lambda x: sum(x) == 1 and x[i] == 1)]
                        singleton_means.append(row["mean"].values[0])
                        
                    def get_epsp_mean(config):
                        epsp = baseline_mean
                        for i, val in enumerate(config):
                            if val >= 0.5:
                                epsp += (singleton_means[i] - baseline_mean)
                        return epsp
                        
                    jax_epsp = get_epsp_mean(rho_j_config)
                    neuron_epsp = get_epsp_mean(rho_n_config)
                    jax_ratio = jax_epsp / baseline_mean
                    neuron_ratio = neuron_epsp / baseline_mean
                    
                    print(f"    Predicted EPSP Ratio JAX: {jax_ratio:.3f} | NEURON: {neuron_ratio:.3f}")
                else:
                    print(f"    [EPSP] Basis file not found at {basis_csv}")
            except Exception as e:
                print(f"    [EPSP] Exception: {e}")


        if all_rj: 
            j, n, r = np.array(all_rj), np.array(all_rn), np.array(all_er) 
            print(f"\n  Summary: |rho_diff|={np.mean(np.abs(j - n)):.4f}" 
                  f"  eff_ratio={np.nanmean(r):.3f}+/-{np.nanstd(r):.3f}")


if __name__ == "__main__":
    main()
