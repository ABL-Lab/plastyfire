#!/usr/bin/env python3
import os
import sys
import logging
import argparse
import multiprocessing
import numpy as np
import pandas as pd
import traceback
import bluecellulab
from pathlib import Path
from libsonata import SpikeReader
from plastyfire.epg import ParamsGenerator
from bluepysnap import Simulation as BluePySnapSimulation
from conntility.io.synapse_report import get_presyn_mapping

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Constants and Mappings
SYNPROPS = ["Cpre", "Cpost", "loc", "Use0_TM", "Dep_TM", "Fac_TM", "Nrrp_TM", "gmax0_AMPA", "gmax_NMDA",
            "volume_CR", "synapseID", "theta_d_GB", "theta_p_GB"]
PARAM_MAP = {"Use_d_TM": "Use_d", "Use_p_TM": "Use_p", "Use0_TM": "Use",
             "Dep_TM": "Dep", "Fac_TM": "Fac", "Nrrp_TM": "Nrrp"}

def _set_local_params(synapse, fit_params, extra_params, c_pre=0., c_post=0.):
    """Sets synaptic parameters in bluecellulab"""
    for key, val in extra_params.items():
        if key in PARAM_MAP:
            setattr(synapse.hsynapse, PARAM_MAP[key], val)
        else:
            if key == "loc":
                continue
            setattr(synapse.hsynapse, key, val)
    # Note: Threshold setting logic removed/simplified as in original if fit_params is None
    # If necessary, copy the full logic from synapse_combinations_generator.py

def _map_syn_idx(sim_config, post_gid, syn_idx, edge_pop):
    return get_presyn_mapping(BluePySnapSimulation(sim_config).circuit, edge_pop,
                              pd.MultiIndex.from_tuples([(post_gid, syn_id) for syn_id in syn_idx]))

def calculate_epsp_amplitude(t, v, pre_spikes):
    """
    Calculate EPSP amplitude.
    Assumes positive deflection.
    Baseline: Average voltage in 10ms window BEFORE the first spike.
    Peak: Max voltage in window [spike_time, spike_time + 100ms] (or end of trace).
    If multiple spikes, this simple logic might need adjustment, but for basis generation 
    usually simple protocols are used. 
    """
    # Default to 0 if no spikes or empty trace
    if len(pre_spikes) == 0:
        return 0.0
    
    # Use the first spike for baseline alignment
    first_spike = pre_spikes[0]
    
    # Baseline window: [first_spike - 10, first_spike]
    
    # Ensure indices are valid
    baseline_indices = np.where((t >= first_spike - 10.0) & (t < first_spike))[0]
    
    if len(baseline_indices) > 0:
        baseline_v = np.mean(v[baseline_indices])
    else:
        # Fallback: first 10ms of simulation
        baseline_v = np.mean(v[:min(len(v), 100)])
        
    # Search for peak after the spike
    # We look for the global max in the trace after the spike
    # Taking the max over the whole trace after spike might cover multiple EPSPs if they exist.
    # For a basis, we usually want the response to the train or pulse.
    # The user manual calculation might be more complex, but (Max - Baseline) is a safe proxy for "Magnitude".
    
    # Look for max V after first spike
    post_spike_indices = np.where(t >= first_spike)[0]
    if len(post_spike_indices) == 0:
        return 0.0
        
    peak_v = np.max(v[post_spike_indices])
    
    amplitude = peak_v - baseline_v
    return amplitude

def run_simulation_trial(args):
    """
    Worker function to run a single simulation trial.
    Returns: (rho_config, trial, epsp_amplitude)
    """
    sim_config_path, pre_gid, post_gid, rho_config, node_pop, trial = args
    
    try:
        # Initialize bluecellulab
        # Note: Imports are inside to ensure fresh state if using multiprocessing 'spawn' or similar, 
        # though 'fork' is standard on Linux.
        # Original script re-imported bluecellulab here.
        import bluecellulab
        bluecellulab.set_verbose(2)
        
        # Load mechanisms - Hardcoded path from original script
        mechanisms_path = "/project/ctb-emuller/dhuruva/DEES_cell_packages/"
        bluecellulab.neuron.load_mechanisms(mechanisms_path)
        
        # Set global parameters
        h = bluecellulab.neuron.h
        h.cao_CR_GluSynapse = 2.0
        h.minis_single_vesicle_GluSynapse = 0.0
        h.init_depleted_GluSynapse = 0.0
        h.tau_effca_GB_GluSynapse = 278.3177658387
        h.gamma_d_GB_GluSynapse = 101.5387594661
        h.gamma_p_GB_GluSynapse = 216.1841700668
        
        # Random seed
        np.random.seed(trial)
        
        sim = bluecellulab.CircuitSimulation(sim_config_path, base_seed=trial)
        
        workdir = os.path.dirname(sim_config_path)
        pre_spikes = SpikeReader(os.path.join(workdir, "prespikes.h5"))[node_pop].get_dict()["timestamps"]
        
        sim.instantiate_gids([(node_pop, post_gid)], 
                           add_synapses=True, 
                           add_minis=False, 
                           add_pulse_stimuli=True,
                           intersect_pre_gids=[(node_pop, pre_gid)],
                           pre_spike_trains={(node_pop, pre_gid): pre_spikes})
        
        cell = sim.cells[(node_pop, post_gid)]
        
        # Uninsert SK_E2
        for sec in cell.somatic + cell.axonal:
            sec.uninsert("SK_E2")
            
        # Synapse Recording Setup (only if needed, here we simplify to check config)
        syn_idx = []
        for syn_id, synapse in cell.synapses.items():
            syn_idx.append(syn_id[1])

        # Parameter Generation and Setting
        bluepysnap_sim = BluePySnapSimulation(sim_config_path)
        EXTRA_RECIPE_PATH = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/biodata/recipe.csv"
        edge_pop = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"
        pgen = ParamsGenerator(bluepysnap_sim.circuit, node_pop, edge_pop, EXTRA_RECIPE_PATH)
        syn_extra_params = pgen.generate_params(pre_gid, post_gid)
        
        df = _map_syn_idx(sim_config_path, post_gid, syn_idx, edge_pop)
        
        for syn_id, synapse in cell.synapses.items():
            if syn_extra_params is not None:
                global_syn_id = df.loc[df["local_syn_idx"] == syn_id[1]].index[0]
                _set_local_params(synapse, None, syn_extra_params[global_syn_id])
                
        # Apply Rho Configuration
        rho_values = [int(x) for x in rho_config.split(",")]
        syn_idx_counter = 0
        for syn_id, synapse in cell.synapses.items():
            if syn_idx_counter < len(rho_values):
                rho_val = rho_values[syn_idx_counter]
                if rho_val >= 0.5:
                    synapse.hsynapse.rho0_GB = 1.0
                    synapse.hsynapse.Use = synapse.hsynapse.Use_p
                    synapse.hsynapse.gmax_AMPA = synapse.hsynapse.gmax_p_AMPA
                else:
                    synapse.hsynapse.rho0_GB = 0.0
                    synapse.hsynapse.Use = synapse.hsynapse.Use_d
                    synapse.hsynapse.gmax_AMPA = synapse.hsynapse.gmax_d_AMPA
            syn_idx_counter += 1
            
        # Disable plasticity
        for syn_id, synapse in cell.synapses.items():
            synapse.hsynapse.theta_d_GB = 1000000
            synapse.hsynapse.theta_p_GB = 1000000
            
        # Run Simulation
        h.cvode_active(1)
        sim.run(240000, cvode=True) # Check if this duration is required or if it can be shortened
        
        # Extract Data
        t = np.array(sim.get_time())
        v = np.array(sim.get_voltage_trace((node_pop, post_gid)))
        
        # Compute EPSP Amplitude
        epsp_amp = calculate_epsp_amplitude(t, v, pre_spikes)
        
        # Clean up
        del sim
        
        return (rho_config, trial, epsp_amp)
        
    except Exception as e:
        logger.error(f"Error in trial {trial}: {e}")
        traceback.print_exc()
        return (rho_config, trial, float('nan'))

def count_synapses(sim_config_path, pre_gid, post_gid, node_pop="S1nonbarrel_neurons"):
    try:
        import bluecellulab
        sim = bluecellulab.CircuitSimulation(sim_config_path)
        sim.instantiate_gids([(node_pop, post_gid)], 
                           add_synapses=True, 
                           add_minis=False, 
                           add_pulse_stimuli=False,
                           intersect_pre_gids=[(node_pop, pre_gid)])
        cell = sim.cells[(node_pop, post_gid)]
        count = len(cell.synapses)
        del sim
        return count
    except Exception as e:
        logger.error(f"Error counting synapses: {e}")
        return 0

def generate_basis_configs(synapse_count):
    configs = []
    # 1. Baseline: All 0s
    configs.append([0] * synapse_count)
    # 2. Single Potentiation
    for i in range(synapse_count):
        c = [0] * synapse_count
        c[i] = 1
        configs.append(c)
    # 3. Maximum: All 1s
    configs.append([1] * synapse_count)
    return [','.join(map(str, c)) for c in configs]

def main():
    parser = argparse.ArgumentParser(description="Run EPSP Basis Simulations")
    parser.add_argument("--pre-gid", type=int, required=True)
    parser.add_argument("--post-gid", type=int, required=True)
    parser.add_argument("--sim-config", type=str, required=True)
    parser.add_argument("--output-csv", type=str, required=True)
    parser.add_argument("--num-trials", type=int, default=10)
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args()
    
    # 1. Resolve Simulation Config
    sim_config = args.sim_config
    if not os.path.exists(sim_config):
        # User might have provided a directory or a non-existent file path
        # Try to find it in common locations relative to the seed/base dir
        base_path = Path(sim_config).parent
        potential_paths = [
            base_path / "L5TTPC_L5TTPC/simulations/single_cells/simulation_config.json",
            base_path / "simulations/single_cells/simulation_config.json",
            base_path / "simulation_config.json"
        ]
        found = False
        for p in potential_paths:
            if p.exists():
                logger.info(f"Config not found at {sim_config}, found valid config at {p}")
                sim_config = str(p)
                found = True
                break
        
        if not found:
             # Try finding ANY simulation_config.json in the directory tree
             logger.info(f"Searching for any simulation_config.json in {base_path}...")
             matches = list(base_path.glob("**/simulation_config.json"))
             
             # Filter matches to prioritize those that look like pair simulation directories
             # Heuristic: look for pre-post pair in path
             pair_str = f"{args.pre_gid}-{args.post_gid}"
             pair_matches = [p for p in matches if pair_str in str(p)]
             
             if pair_matches:
                 # If multiple, maybe pick one with "10Hz" or just the first
                 # Let's pick the first one
                 sim_config = str(pair_matches[0])
                 logger.info(f"Auto-resolved pair-specific config to: {sim_config}")
                 found = True
             elif matches:
                 # Fallback to any config found (e.g. single_cells)
                 # BUT single_cells might not have prespikes.h5
                 # So we prefer one that is NOT single_cells if possible?
                 # Actually, let's just take the first one and hope, but log warning
                 sim_config = str(matches[0])
                 logger.info(f"Auto-resolved config to (potentially generic): {sim_config}")
                 found = True

        if not found:
            logger.error(f"Could not find valid simulation_config.json starting from {args.sim_config}")
            return

        if not found:
            logger.error(f"Could not find valid simulation_config.json starting from {args.sim_config}")
            return

    # 2. Count Synapses
    n_synapses = count_synapses(sim_config, args.pre_gid, args.post_gid)
    if n_synapses == 0:
        logger.warning(f"No synapses for {args.pre_gid}->{args.post_gid} using config {sim_config}")
        return

    logger.info(f"Generating basis for {args.pre_gid}->{args.post_gid} ({n_synapses} synapses)")
    
    # 3. Generate Configs
    basis_configs = generate_basis_configs(n_synapses)
    
    # 4. Prepare Tasks
    tasks = []
    for config in basis_configs:
        for t in range(args.num_trials):
            tasks.append((sim_config, args.pre_gid, args.post_gid, config, "S1nonbarrel_neurons", t))
            
    # 5. Run Parallel
    n_workers = args.workers if args.workers else min(multiprocessing.cpu_count(), len(tasks))
    logger.info(f"Running {len(tasks)} tasks with {n_workers} workers")
    
    results = []
    with multiprocessing.Pool(n_workers) as pool:
        for res in pool.map(run_simulation_trial, tasks):
            results.append(res)
            
    # 5. Save to CSV
    # Aggregate: Mean EPSP per config
    if not results:
        logger.warning("No results to save.")
        return

    df_raw = pd.DataFrame(results, columns=["config", "trial", "epsp"])
    
    # Calculate stats
    stats = df_raw.groupby("config")["epsp"].agg(["mean", "std", "count"]).reset_index()
    stats["pre_gid"] = args.pre_gid
    stats["post_gid"] = args.post_gid
    
    # Ensure output directory exists
    Path(args.output_csv).parent.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Saving results to {args.output_csv}")
    stats.to_csv(args.output_csv, index=False)
    
if __name__ == "__main__":
    main()
