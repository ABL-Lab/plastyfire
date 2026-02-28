#!/usr/bin/env python3
"""
Deep diagnostic: runs the ACTUAL JAX scan_step and reports what really happens.

Unlike plot_cicr_diagnostic.py (which uses its own _compute_rho with standard GB),
this script traces through the real model step by step, revealing:
  - Per-synapse final rho and rho0 for each protocol
  - Fraction of synapses that cross the 0.5 boundary (pot ↔ depot)
  - The actual predicted EPSP ratios (matching the loss function)
  - Summary of what the optimizer sees

Usage:
    python diagnose_model.py --params '{"gamma_d": 69.2, ...}'
    python diagnose_model.py   # uses DEFAULT_PARAMS
"""

import os, sys, json, argparse, pickle
import numpy as np

import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from cicr_common import preload_all_data, collate_protocol_to_jax, EXPERIMENTAL_TARGETS
from cicr_gemini_lse import DynamicTauEligibilityModel


def diagnose(model, params_dict=None, protocols=None, max_pairs=None):
    if protocols is None:
        protocols = ["10Hz_10ms", "10Hz_-10ms"]
    dp = dict(model.DEFAULT_PARAMS)
    if params_dict is not None:
        dp.update(params_dict)

    x = np.array([dp[n] for n in model.PARAM_NAMES])
    params = model.unpack_params(jnp.array(x))

    print(f"\n{'='*70}")
    print(f"  DEEP DIAGNOSTIC — {model.DESCRIPTION}")
    print(f"{'='*70}")
    print(f"\n  Parameters:")
    for name in model.PARAM_NAMES:
        print(f"    {name:20s} = {dp[name]:.6f}")

    protocol_data = preload_all_data(max_pairs=max_pairs, protocols=protocols)

    step_factory = model.get_step_factory()
    init_fn = model.get_init_fn()

    # Build vectorized simulation (same pattern as setup_jax in cicr_common.py)
    def sim_synapse(cai_trace, t_trace, c_pre, c_post, is_apical, rho0, params):
        dt_trace = jnp.diff(t_trace, prepend=t_trace[0])
        dt_trace = jnp.where(dt_trace <= 0, 1e-6, dt_trace)
        init = init_fn(cai_trace[0], rho0)
        scan_fn = step_factory(params, (c_pre, c_post, is_apical))
        final, _ = jax.lax.scan(scan_fn, init, (cai_trace, dt_trace))
        return final[-1]  # rho_final

    # vmap over synapses (c_pre, c_post, is_apical, rho0, cai vary; t_trace, params shared)
    vmap_syn = jax.vmap(sim_synapse, in_axes=(0, None, 0, 0, 0, 0, None))

    # vmap over pairs (everything varies except params)
    def sim_pair(cai_p, t_p, cpre_p, cpost_p, isapi_p, rho0_p, valid_p, params):
        rho_final = vmap_syn(cai_p, t_p, cpre_p, cpost_p, isapi_p, rho0_p, params)
        return jnp.where(valid_p, rho_final, jnp.nan)

    vmap_pair = jax.vmap(sim_pair, in_axes=(0, 0, 0, 0, 0, 0, 0, None))

    @jax.jit
    def run_protocol(collated, params):
        return vmap_pair(
            collated['cai'], collated['t'],
            collated['c_pre'], collated['c_post'],
            collated['is_apical'], collated['rho0'],
            collated['valid'], params
        )

    for proto in protocols:
        pairs = protocol_data.get(proto, [])
        if not pairs:
            print(f"\n  [{proto}] No data loaded, skipping.")
            continue

        target = EXPERIMENTAL_TARGETS.get(proto, float('nan'))
        print(f"\n{'─'*70}")
        print(f"  Protocol: {proto}   Target: {target:.4f}   Pairs loaded: {len(pairs)}")
        print(f"{'─'*70}")

        # Collate all pairs into padded JAX arrays (uses existing infrastructure)
        collated = collate_protocol_to_jax(pairs)

        # Run ALL synapses across ALL pairs in one vectorized JIT call
        rho_final_all = np.array(run_protocol(collated, params))  # (n_pairs, max_syn)
        valid = np.array(collated['valid'])                        # (n_pairs, max_syn)
        rho0_all = np.array(collated['rho0'])
        sm_all = np.array(collated['singletons'])
        bmean_all = np.array(collated['baseline'])

        # Flatten valid synapses for statistics
        all_rho0 = rho0_all[valid]
        all_rho_final = rho_final_all[valid]

        # Compute per-pair EPSP ratios
        pair_ratios = []
        for pi in range(len(pairs)):
            v = valid[pi]
            rho0_p = rho0_all[pi, v]
            rho_f_p = rho_final_all[pi, v]
            sm_p = sm_all[pi, v]
            bmean = bmean_all[pi]

            contribs_before = np.where(rho0_p >= 0.5, sm_p - bmean, 0.0)
            contribs_after = np.where(rho_f_p >= 0.5, sm_p - bmean, 0.0)
            epsp_before = bmean + contribs_before.sum()
            epsp_after = bmean + contribs_after.sum()
            ratio = epsp_after / epsp_before if epsp_before > 0 else float('nan')
            pair_ratios.append(ratio)
        n_total = len(all_rho0)

        initially_pot = all_rho0 >= 0.5
        initially_dep = ~initially_pot
        finally_pot = all_rho_final >= 0.5

        n_pot_before = initially_pot.sum()
        n_pot_after = finally_pot.sum()
        n_depotentiated = (initially_pot & ~finally_pot).sum()
        n_potentiated = (initially_dep & finally_pot).sum()
        n_unchanged_pot = (initially_pot & finally_pot).sum()
        n_unchanged_dep = (initially_dep & ~finally_pot).sum()

        pair_ratios = np.array(pair_ratios)
        pred_ratio = np.nanmean(pair_ratios)

        print(f"\n  Synapse-level summary ({n_total} synapses across {len(pairs)} pairs):")
        print(f"    Initially potentiated (rho0>=0.5):  {n_pot_before:4d} / {n_total}")
        print(f"    Finally potentiated   (rho>=0.5):   {n_pot_after:4d} / {n_total}")
        print(f"    Depotentiated (pot→dep):            {n_depotentiated:4d}")
        print(f"    Potentiated   (dep→pot):            {n_potentiated:4d}")
        print(f"    Unchanged pot (pot→pot):            {n_unchanged_pot:4d}")
        print(f"    Unchanged dep (dep→dep):            {n_unchanged_dep:4d}")

        print(f"\n  rho distribution (initially potentiated synapses):")
        if n_pot_before > 0:
            rho_pot = all_rho_final[initially_pot]
            print(f"    min={rho_pot.min():.4f}  p25={np.percentile(rho_pot,25):.4f}  "
                  f"median={np.median(rho_pot):.4f}  p75={np.percentile(rho_pot,75):.4f}  "
                  f"max={rho_pot.max():.4f}")
            for thresh in [0.3, 0.4, 0.45, 0.5, 0.6, 0.8, 0.99]:
                frac = (rho_pot < thresh).sum() / len(rho_pot) * 100
                print(f"    % below {thresh:.2f}: {frac:6.1f}%")

        print(f"\n  rho distribution (initially depotentiated synapses):")
        if initially_dep.sum() > 0:
            rho_dep = all_rho_final[initially_dep]
            print(f"    min={rho_dep.min():.4f}  p25={np.percentile(rho_dep,25):.4f}  "
                  f"median={np.median(rho_dep):.4f}  p75={np.percentile(rho_dep,75):.4f}  "
                  f"max={rho_dep.max():.4f}")
            for thresh in [0.3, 0.4, 0.45, 0.5, 0.6, 0.8, 0.99]:
                frac = (rho_dep >= thresh).sum() / len(rho_dep) * 100
                print(f"    % above {thresh:.2f}: {frac:6.1f}%")

        print(f"\n  Pair-level EPSP ratios ({len(pair_ratios)} pairs):")
        print(f"    min={np.nanmin(pair_ratios):.4f}  mean={pred_ratio:.4f}  "
              f"max={np.nanmax(pair_ratios):.4f}  std={np.nanstd(pair_ratios):.4f}")
        print(f"    Target: {target:.4f}   Error: {pred_ratio - target:+.4f}")

    print(f"\n{'='*70}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--params", type=str,
                        default=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                             "best_params_primed-junctional_cicr_synergy_20260221_232049.json"),
                        help="JSON string or file path with parameter dict")
    parser.add_argument("--protocols", nargs="+", default=["2Hz_5ms", "5Hz_5ms", "10Hz_10ms", "10Hz_-10ms", "50Hz_10ms"])
    parser.add_argument("--max-pairs", type=int, default=None)
    parser.add_argument("--model", type=str, default="PrimedJunctionalCICRModel",
                        help="Name of the model class to diagnose (e.g. PrimedJunctionalCICRModel)")
    args = parser.parse_args()

    import importlib
    
    # Try multiple potential files where models are defined
    model_modules = [
        "cicr_primed_junctional",
        "cicr_gemini_lse",
        "cicr_gemini_clipped",
        "cicr_gemini_syn",
        "cicr_dual_pathway",
        "cicr_dual_hill",
        "cicr_mwc_advanced",
        "cicr_neymotin",
        "cicr_bucket"
    ]
    
    model_class = None
    for mod_name in model_modules:
        try:
            module = importlib.import_module(mod_name)
            if hasattr(module, args.model):
                model_class = getattr(module, args.model)
                break
        except ImportError:
            pass

    if model_class is None:
        print(f"Error: Could not find model class '{args.model}' in any known script.")
        print(f"Ensure the class name is correct (e.g. ClippedLatentSTAPCDModel)")
        sys.exit(1)
        
    model = model_class()

    params = None
    if args.params:
        if os.path.isfile(args.params):
            with open(args.params) as f:
                params = json.load(f)
        else:
            params = json.loads(args.params)

    diagnose(model, params_dict=params, protocols=args.protocols, max_pairs=args.max_pairs)