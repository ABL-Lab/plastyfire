#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import numpy as np
import logging
from run_basis_pair import run_simulation_trial, count_synapses

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(description="Verify Linearity Logic")
    parser.add_argument("--basis-csv", required=True, help="Path to existing basis results CSV")
    parser.add_argument("--csv", required=True, help="Path to index CSV (to find config path)")
    parser.add_argument("--sim-config-override", help="Override sim config path")
    args = parser.parse_args()

    # 1. Load Basis Results
    df_basis = pd.read_csv(args.basis_csv)
    pre_gid = df_basis['pre_gid'].iloc[0]
    post_gid = df_basis['post_gid'].iloc[0]
    logger.info(f"Verifying pair {pre_gid}->{post_gid}")
    
    # 2. Re-resolve Sim Config (Recycled logic from validate_single_pair)
    if args.sim_config_override:
        sim_config = args.sim_config_override
    else:
        # Look up in index CSV
        df_index = pd.read_csv(args.csv)
        # Find row
        row = df_index[(df_index['pregid'] == pre_gid) & (df_index['postgid'] == post_gid)]
        if row.empty:
            logger.error("Pair not found in index CSV")
            return
        sim_path = row['path'].values[0]
        # Resolve config
        from pathlib import Path
        sim_dir = Path(sim_path).parent
        # Try generic or specific
        sim_config = sim_dir / "simulation_config.json"
        if not sim_config.exists():
            # Try finding any config in tree logic (simplified)
            matches = list(sim_dir.parent.glob("**/simulation_config.json"))
            if matches:
                sim_config = matches[0]
            else:
                logger.error(f"Could not resolve config from {sim_path}")
                return
    
    logger.info(f"Using config: {sim_config}")
    
    # 3. Parse Basis Values
    # We expect config 0,0,0... (Baseline)
    # And 1,0,0... etc (Singletons)
    
    # Parse configs to lists
    df_basis['config_list'] = df_basis['config'].apply(lambda x: [int(i) for i in x.split(',')])
    n_synapses = len(df_basis['config_list'].iloc[0])
    
    # Identify Baseline
    baseline_row = df_basis[df_basis['config_list'].apply(lambda x: sum(x) == 0)]
    if baseline_row.empty:
        logger.error("Baseline config (0,0,...) not found in basis results")
        return
    baseline_amp = baseline_row['mean'].values[0]
    logger.info(f"Baseline Amplitude: {baseline_amp:.4f}")
    
    # Identify Singletons
    singleton_amps = []
    for i in range(n_synapses):
        # Create pattern for synapse i
        # Config where sum is 1 and index i is 1
        row = df_basis[df_basis['config_list'].apply(lambda x: sum(x) == 1 and x[i] == 1)]
        if row.empty:
            logger.error(f"Singleton config for synapse {i} not found")
            return
        singleton_amps.append(row['mean'].values[0])
    
    logger.info(f"Singleton Amplitudes: {[f'{x:.4f}' for x in singleton_amps]}")
    
    # 4. Construct Test Case
    # Let's pick a combination of 2 random synapses (if n > 1)
    if n_synapses < 2:
        logger.info("Not enough synapses to test linearity (need > 1). Testing Max config instead.")
        test_config = [1] * n_synapses
    else:
        # Pick first two
        test_config = [0] * n_synapses
        test_config[0] = 1
        test_config[1] = 1
    
    test_config_str = ','.join(map(str, test_config))
    logger.info(f"Testing Config: {test_config_str}")
    
    # 5. Predict
    # Prediction = Baseline + Sum(Single_i - Baseline) for i where config[i]==1
    prediction = baseline_amp
    for i in range(n_synapses):
        if test_config[i] == 1:
            prediction += (singleton_amps[i] - baseline_amp)
            
    logger.info(f"Predicted Amplitude: {prediction:.4f}")
    
    # 6. Run Actual Simulation
    # Run 1 trial (trial 0)
    args_sim = (str(sim_config), int(pre_gid), int(post_gid), test_config_str, "S1nonbarrel_neurons", 0)
    
    logger.info("Running simulation...")
    # This might fail if run_basis_pair depends on other things, but should return amplitude
    _, _, actual_amp = run_simulation_trial(args_sim)
    
    logger.info(f"Actual Amplitude: {actual_amp:.4f}")
    
    # 7. Compare
    diff = abs(actual_amp - prediction)
    pct_error = (diff / abs(actual_amp)) * 100 if actual_amp != 0 else 0.0
    
    logger.info(f"Difference: {diff:.4f}")
    logger.info(f"Percent Error: {pct_error:.2f}%")
    
    if pct_error < 5.0: # Arbitrary threshold
        logger.info("SUCCESS: Linearity holds (within 5% error)")
    else:
        logger.warning("WARNING: Linearity error > 5%")

if __name__ == "__main__":
    main()
