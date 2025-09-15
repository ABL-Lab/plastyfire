#!/usr/bin/env python3
"""
Comprehensive EPSP ratio comparison between actual simulations and expected values
"""

import os
import pickle
import numpy as np
import pandas as pd
from pathlib import Path
import glob
import logging

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def compute_epsp_safe(pkl_file, window=100):
    """Safely compute EPSP values with bounds checking"""
    with open(pkl_file, 'rb') as f:
        data = pickle.load(f)
    try:
        t, v, spikes = data['t'], data['v'], data['pre_spikes']
    except:
        t, v, spikes = data['t'], data['v'], data['prespikes']
    
    # Extract rho values
    if len(data['rho_GB']) > 100:
        data['rho_GB'] = np.transpose(data['rho_GB'])
    
    initial_rho = [0 if k[0] < 0.5 else 1 for k in data['rho_GB']]
    final_rho = [0 if k[-1] < 0.5 else 1 for k in data['rho_GB']]
    
    valid_spikes = spikes[spikes < (t[-1] - window)][:60]
    
    # Import get_epsp_vector function
    try:
        from plastyfire.ephysutils import get_epsp_vector
    except ImportError:
        logger.error("Could not import get_epsp_vector from plastyfire.ephysutils")
        return None, 0, len(valid_spikes), len(spikes), t, v, data, initial_rho, final_rho
    
    if len(t) > 237823:  # Long simulation - compute before/after ratio
        epsp_values_before = get_epsp_vector(t, v, spikes[:60], window)
        epsp_values_after = get_epsp_vector(t, v, spikes[-60:], window)
        epsp_ratio = np.mean(epsp_values_after) / np.mean(epsp_values_before)
        avg_epsp = epsp_ratio
        epsp_values = 0
    else:
        epsp_values = get_epsp_vector(t, v, valid_spikes, window)  
        avg_epsp = np.mean(epsp_values) if len(epsp_values) > 0 else 0.0

    return epsp_values, avg_epsp, len(valid_spikes), len(spikes), t, v, data, initial_rho, final_rho

def map_protocol_to_folder(protocol_id):
    """Map protocol_id to simulation folder structure"""
    mapping = {
        'mrk97_01': '2Hz_5ms',
        'mrk97_02': '5Hz_5ms', 
        'mrk97_07': '10Hz_10ms',
        'mrk97_08': '10Hz_-10ms',
        'sjh06_02': '50Hz_10ms'
    }
    return mapping.get(protocol_id, None)

def create_epsp_comparison_csv():
    """Create comprehensive EPSP ratio comparison CSV for both connection types and cell pairs"""
    
    # Base paths
    base_path = Path("/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire")
    simulation_df_path = base_path / "plastyfire" / "simulation_epsp_df.csv"
    paired_recordings_path = base_path / "biodata" / "paired_recordings.csv"
    refitting_base = base_path / "refitting_results" / "fitting" / "n1" / "seed19091997"
    
    # Load existing data
    sim_df = pd.read_csv(simulation_df_path)
    paired_df = pd.read_csv(paired_recordings_path)
    
    # Define connection types and cell pairs to process
    connection_configs = [
        {"connection_type": "L23PC_L5TTPC", "cell_pair": "15864-199162"},
        {"connection_type": "L5TTPC_L5TTPC", "cell_pair": "205559-199162"}
    ]
    
    # Results storage
    all_comparison_results = []
    
    # Process each connection type and cell pair
    for config in connection_configs:
        connection_type = config["connection_type"]
        cell_pair = config["cell_pair"]
        
        logger.info(f"Processing {connection_type} with cell pair {cell_pair}")
        
        # Process each unique combination for this connection type
        unique_combinations = sim_df[['protocol_id', 'param_hash']].drop_duplicates()
        
        for _, row in unique_combinations.iterrows():
            protocol_id = row['protocol_id']
            param_hash = row['param_hash']
            
            logger.info(f"Processing {protocol_id} with hash {param_hash} for {connection_type}")

            # Map protocol to folder
            folder_name = map_protocol_to_folder(protocol_id)
            if not folder_name:
                logger.warning(f"No folder mapping for protocol {protocol_id}")
                continue
                
            # Find simulation file
            sim_pattern = refitting_base / connection_type / "simulations" / cell_pair / folder_name / f"simulation_{param_hash}.pkl"
            
            if not sim_pattern.exists():
                logger.warning(f"Simulation file not found: {sim_pattern}")
                continue
            
            # Extract rho values from simulation
            try:
                _, actual_ratio, _, _, t, v, data, initial_rho, final_rho = compute_epsp_safe(str(sim_pattern))
                logger.info(f"Initial rho: {initial_rho}, Final rho: {final_rho}")
                
                # Extract spikes from simulation data
                try:
                    spikes = data['pre_spikes']
                except:
                    spikes = data['prespikes']
                
                # Compute actual EPSP values from simulation data
                # For long simulations, we need to compute before/after values separately
                if len(t) > 237823:  # Long simulation - compute before/after values
                    try:
                        from plastyfire.ephysutils import get_epsp_vector
                        valid_spikes = spikes[spikes < (t[-1] - 100)][:60]
                        actual_epsp_before = np.mean(get_epsp_vector(t, v, valid_spikes, 100)) if len(valid_spikes) > 0 else 0.0
                        actual_epsp_after = np.mean(get_epsp_vector(t, v, spikes[-60:], 100)) if len(spikes) >= 60 else 0.0
                    except ImportError:
                        logger.warning("Could not import get_epsp_vector, using placeholder values")
                        actual_epsp_before = 0.0
                        actual_epsp_after = 0.0
                else:
                    logger.error("Short simulation not supported")
                    continue
                    
            except Exception as e:
                logger.error(f"Error processing simulation file {sim_pattern}: {e}")
                continue
                
            # Find corresponding ephys_data files using rho values
            ephys_dir = refitting_base / "ephys_data"
            
            # Convert rho lists to strings for filename matching
            initial_rho_str = '_'.join(map(str, initial_rho))
            final_rho_str = '_'.join(map(str, final_rho))
            
            # Find ephys files - use cell pair in filename
            ephys_before_pattern = f"ephys_data_{cell_pair.replace('-', '_')}_{initial_rho_str}.pkl"
            ephys_after_pattern = f"ephys_data_{cell_pair.replace('-', '_')}_{final_rho_str}.pkl"
            
            ephys_before_path = ephys_dir / ephys_before_pattern
            ephys_after_path = ephys_dir / ephys_after_pattern
            
            if not ephys_before_path.exists():
                logger.warning(f"Ephys before file not found: {ephys_before_path}")
                continue
                
            if not ephys_after_path.exists():
                logger.warning(f"Ephys after file not found: {ephys_after_path}")
                continue
                
            # Compute expected EPSP ratios from ephys data
            try:
                _, expected_epsp_before, _, _, _, _, _ = compute_epsp_safe(str(ephys_before_path))[:7]
                _, expected_epsp_after, _, _, _, _, _ = compute_epsp_safe(str(ephys_after_path))[:7]
                
                if expected_epsp_before > 0:
                    expected_epsp_ratio = expected_epsp_after / expected_epsp_before
                else:
                    expected_epsp_ratio = 0
                    
                # Calculate percentage difference
                if expected_epsp_ratio > 0:
                    difference_percentage = ((actual_ratio - expected_epsp_ratio) / expected_epsp_ratio) * 100
                else:
                    difference_percentage = 0
                    
                # Store results
                all_comparison_results.append({
                    'connection_type': connection_type,
                    'cell_pair': cell_pair,
                    'protocol_id': protocol_id,
                    'param_hash': param_hash,
                    'actual_epsp_ratio': actual_ratio,
                    'expected_epsp_ratio': expected_epsp_ratio,
                    'difference_in_percentage': difference_percentage,
                    'initial_rho': str(initial_rho),
                    'final_rho': str(final_rho),
                    'expected_epsp_before': expected_epsp_before,
                    'expected_epsp_after': expected_epsp_after,
                    'actual_epsp_before': actual_epsp_before,
                    'actual_epsp_after': actual_epsp_after
                })
                
                logger.info(f"✓ {connection_type} {cell_pair} {protocol_id}: Actual={actual_ratio:.4f}, Expected={expected_epsp_ratio:.4f}, Diff={difference_percentage:.2f}%")
                
            except Exception as e:
                logger.error(f"Error computing expected EPSP ratio for {connection_type} {cell_pair} {protocol_id}: {e}")
                continue
    
    # Create DataFrame and save
    if all_comparison_results:
        results_df = pd.DataFrame(all_comparison_results)
        output_path = base_path / "epsp_ratio_comparison.csv"
        results_df.to_csv(output_path, index=False)
        logger.info(f"Results saved to {output_path}")
        
        # Print summary statistics
        print("\n=== EPSP Ratio Comparison Summary ===")
        print(f"Total comparisons: {len(results_df)}")
        print(f"Connection types processed: {results_df['connection_type'].unique()}")
        print(f"Cell pairs processed: {results_df['cell_pair'].unique()}")
        print(f"Mean absolute difference: {results_df['difference_in_percentage'].abs().mean():.2f}%")
        print(f"Max difference: {results_df['difference_in_percentage'].abs().max():.2f}%")
        print(f"Min difference: {results_df['difference_in_percentage'].abs().min():.2f}%")
        
        # Print summary by connection type
        for connection_type in results_df['connection_type'].unique():
            conn_df = results_df[results_df['connection_type'] == connection_type]
            print(f"\n{connection_type}:")
            print(f"  Comparisons: {len(conn_df)}")
            print(f"  Mean absolute difference: {conn_df['difference_in_percentage'].abs().mean():.2f}%")
        
        print(f"\nFirst few results:")
        print(results_df[['connection_type', 'cell_pair', 'protocol_id', 'actual_epsp_ratio', 'expected_epsp_ratio', 'difference_in_percentage', 'actual_epsp_before', 'actual_epsp_after']].head())
        
        return results_df
    else:
        logger.warning("No valid comparisons found")
        return None

if __name__ == "__main__":
    results = create_epsp_comparison_csv()
