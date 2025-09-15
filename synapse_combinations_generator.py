#!/usr/bin/env python3

import os
import json
import pandas as pd
import itertools
import logging
import multiprocessing
import numpy as np
import pickle
from pathlib import Path
import bluecellulab
from libsonata import SpikeReader
import plastyfire.ephysutils as ephysutils
from conntility.io.synapse_report import get_presyn_mapping
import traceback
from bluepysnap import Simulation as BluePySnapSimulation

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# bluecellulab.set_verbose(2)
# bluecellulab.neuron.load_mechanisms("/project/ctb-emuller/dhuruva/DEES_cell_packages/")
DEBUG = False  # Default to False, will be overridden by modelfitter
SYNPROPS = ["Cpre", "Cpost", "loc", "Use0_TM", "Dep_TM", "Fac_TM", "Nrrp_TM", "gmax0_AMPA", "gmax_NMDA",
            "volume_CR", "synapseID", "theta_d_GB", "theta_p_GB"]
MOD_PROPS = ["gamma_d_GB", "gamma_p_GB"]  # 'tau_exp_GB'
SYNREC = ["rho_GB", "Use_GB", "gmax_AMPA", "cai_CR", "vsyn", "ica_NMDA", "ica_VDCC", "effcai_GB"]
# because of the constant ping-pong between GluSynapse.mod, SONATA parameter names
# and synapse helpers both in `neurodamus` and in bluecellulab.synapses.synapse_types/GluSynapse() mimicking it
# some variable names have to be patched (to match the current state of GluSynapse.mod)
PARAM_MAP = {"Use_d_TM": "Use_d", "Use_p_TM": "Use_p", "Use0_TM": "Use",
             "Dep_TM": "Dep", "Fac_TM": "Fac", "Nrrp_TM": "Nrrp"}

from plastyfire.epg import ParamsGenerator


def _set_local_params(synapse, fit_params, extra_params, c_pre=0., c_post=0.):
    """Sets synaptic parameters in bluecellulab"""
    for key, val in extra_params.items():  # update basic synapse parameters
        if key in PARAM_MAP:
            setattr(synapse.hsynapse, PARAM_MAP[key], val)
        else:
            if key  == "loc":
                continue
            setattr(synapse.hsynapse, key, val)
    if fit_params is not None:  # update thresholds
        if all(key in fit_params for key in ["a00", "a01"]) and extra_params["loc"] == "basal":
            # set basal depression threshold
            synapse.hsynapse.theta_d_GB = fit_params["a00"] * c_pre + fit_params["a01"] * c_post
        if all(key in fit_params for key in ["a10", "a11"]) and extra_params["loc"] == "basal":
            # set basal potentiation threshold
            synapse.hsynapse.theta_p_GB = fit_params["a10"] * c_pre + fit_params["a11"] * c_post
        if all(key in fit_params for key in ["a20", "a21"]) and extra_params["loc"] == "apical":
            # set apical depression threshold
            synapse.hsynapse.theta_d_GB = fit_params["a20"] * c_pre + fit_params["a21"] * c_post
        if all(key in fit_params for key in ["a30", "a31"]) and extra_params["loc"] == "apical":
            # set apical potentiation threshold
            synapse.hsynapse.theta_p_GB = fit_params["a30"] * c_pre + fit_params["a31"] * c_post


def _map_syn_idx(sim_config, post_gid, syn_idx, edge_pop):
    """Compared to `bluepysnap` which has one global synapse ID for all synapses in the circuit,
    `bluecellulab` (just as `neurodamus`) re-indexes synapses for each postsynaptic cell starting at 0
    (which ID is used for seeding the synapses). This helper gets the mapping between the two indexing versions
    It's unfortunate that this is here (and gets called so many times)... but I couldn't find a better way"""
    return get_presyn_mapping(BluePySnapSimulation(sim_config).circuit, edge_pop,
                              pd.MultiIndex.from_tuples([(post_gid, syn_id) for syn_id in syn_idx]))


def get_unique_cell_pairs(csv_path):
    df = pd.read_csv(csv_path)
    unique_pairs = df[['pregid', 'postgid']].drop_duplicates()
    return unique_pairs.values.tolist()

def dummy_test_function(args):
    """Dummy function to test multiprocessing"""
    sim_config_path, pre_gid, post_gid, rho_config, node_pop, trial = args
    print(f"DEBUG: Processing {pre_gid}->{post_gid} with config {rho_config}, trial {trial} on process {os.getpid()}")
    return rho_config, f"dummy_result_{rho_config}_trial_{trial}"

def _set_local_params(synapse, fit_params, extra_params, c_pre=0., c_post=0.):
    """Sets synaptic parameters in bluecellulab"""
    for key, val in extra_params.items():  # update basic synapse parameters
        if key in PARAM_MAP:
            setattr(synapse.hsynapse, PARAM_MAP[key], val)
        else:
            if key  == "loc":
                continue
            setattr(synapse.hsynapse, key, val)
    if fit_params is not None:  # update thresholds
        if all(key in fit_params for key in ["a00", "a01"]) and extra_params["loc"] == "basal":
            # set basal depression threshold
            synapse.hsynapse.theta_d_GB = fit_params["a00"] * c_pre + fit_params["a01"] * c_post
        if all(key in fit_params for key in ["a10", "a11"]) and extra_params["loc"] == "basal":
            # set basal potentiation threshold
            synapse.hsynapse.theta_p_GB = fit_params["a10"] * c_pre + fit_params["a11"] * c_post
        if all(key in fit_params for key in ["a20", "a21"]) and extra_params["loc"] == "apical":
            # set apical depression threshold
            synapse.hsynapse.theta_d_GB = fit_params["a20"] * c_pre + fit_params["a21"] * c_post
        if all(key in fit_params for key in ["a30", "a31"]) and extra_params["loc"] == "apical":
            # set apical potentiation threshold
            synapse.hsynapse.theta_p_GB = fit_params["a30"] * c_pre + fit_params["a31"] * c_post


def get_epsp_value(args):
    sim_config_path, pre_gid, post_gid, rho_config, node_pop, trial = args
    
    # Create filename using rho config and trial number (replace commas with underscores for valid filename)
    safe_rho_config = rho_config.replace(',', '_')
    filename = f"ephys_data_{pre_gid}_{post_gid}_{safe_rho_config}_trial_{trial}.pkl"
    filepath = os.path.join(os.path.dirname(sim_config_path), '..', '..', '..', '..', 'ephys_data', filename)
    
    # Check if file already exists
    if os.path.exists(filepath):
        logger.info(f"File already exists, skipping simulation: {filepath}")
        return rho_config, 0.0  # Return placeholder value since we're not computing EPSP anyway
    
    logger.info(f"Starting EPSP computation for {pre_gid}->{post_gid}, config {rho_config}, trial {trial} on process {os.getpid()}")
    try:
        import bluecellulab
        
        # Initialize fresh bluecellulab environment for this process
        logger.info(f"Initializing fresh bluecellulab environment on process {os.getpid()}")
        import importlib
        #importlib.reload(bluecellulab)
        bluecellulab.set_verbose(2)
        bluecellulab.neuron.load_mechanisms("/project/ctb-emuller/dhuruva/DEES_cell_packages/")
        bluecellulab.neuron.h.cao_CR_GluSynapse = 2.0
        bluecellulab.neuron.h.minis_single_vesicle_GluSynapse = 0.0
        bluecellulab.neuron.h.init_depleted_GluSynapse = 0.0
        bluecellulab.neuron.h.tau_effca_GB_GluSynapse = 278.3177658387
        bluecellulab.neuron.h.gamma_d_GB_GluSynapse = 101.5387594661
        bluecellulab.neuron.h.gamma_p_GB_GluSynapse = 216.1841700668
        
        # Set random seed for this trial to ensure different results
        np.random.seed(trial)
        logger.info(f"Creating CircuitSimulation for {pre_gid}->{post_gid}")
        sim = bluecellulab.CircuitSimulation(sim_config_path, base_seed=trial)
        
        workdir = os.path.dirname(sim_config_path)
        logger.info(f"Loading spike data from {workdir}/prespikes.h5")
        pre_spikes = SpikeReader(os.path.join(workdir, "prespikes.h5"))[node_pop].get_dict()["timestamps"]
        
        logger.info(f"Instantiating GIDs for {pre_gid}->{post_gid} on process {os.getpid()}")
        try:
            sim.instantiate_gids([(node_pop, post_gid)], 
                               add_synapses=True, 
                               add_minis=False, 
                               add_pulse_stimuli=True,
                               intersect_pre_gids=[(node_pop, pre_gid)],
                               pre_spike_trains={(node_pop, pre_gid): pre_spikes})
            logger.info(f"GID instantiation completed for {pre_gid}->{post_gid} on process {os.getpid()}")
        except Exception as e:
            logger.error(f"GID instantiation failed for {pre_gid}->{post_gid} on process {os.getpid()}: {e}")
            raise
        
        cell = sim.cells[(node_pop, post_gid)]
        for sec in cell.somatic + cell.axonal:
            sec.uninsert("SK_E2")
        logger.info(f"Found {len(cell.synapses)} synapses for {pre_gid}->{post_gid}")
        syn_rec_lst = ["rho_GB"]
        syn_rec, syn_idx = {key: [] for key in syn_rec_lst}, []
        for syn_id, synapse in cell.synapses.items():
            syn_idx.append(syn_id[1])
            if len(syn_rec_lst) != 0:
                for key, lst in syn_rec.items():  # set up recordings
                    recorder = bluecellulab.neuron.h.Vector()
                    recorder.record(getattr(synapse.hsynapse, "_ref_%s" % key))
                    lst.append(recorder)
                    
        bluepysnap_sim = BluePySnapSimulation(sim_config_path)
        EXTRA_RECIPE_PATH = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/biodata/recipe.csv"
        edge_pop = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"
        pgen = ParamsGenerator(bluepysnap_sim.circuit, node_pop, edge_pop, EXTRA_RECIPE_PATH)
        syn_extra_params = pgen.generate_params(pre_gid, post_gid)

        df = _map_syn_idx(sim_config_path, post_gid, syn_idx, edge_pop)
        syn_props = {key: [] for key in SYNPROPS}
        for syn_id, synapse in cell.synapses.items():
            # Reduced logging: removed verbose synapse debug
            # logger.debug("Configuring synapse %d", syn_id[1])
            if syn_extra_params is not None:  # configure local parameters
                global_syn_id = df.loc[df["local_syn_idx"] == syn_id[1]].index[0]
                _set_local_params(synapse, None, syn_extra_params[global_syn_id])
            for key, lst in syn_props.items():  # store synapse properties
                if key == "Cpre":
                    lst.append(0)
                elif key == "Cpost":
                    lst.append(0)
                elif key == "loc":
                    lst.append(syn_extra_params[global_syn_id]["loc"])
                else:
                    if key in PARAM_MAP:
                        lst.append(getattr(synapse.hsynapse, PARAM_MAP[key]))
                    else:
                        lst.append(getattr(synapse.hsynapse, key))

        rho_values = [int(x) for x in rho_config.split(",")]
        logger.info(f"Applying rho configuration: {rho_config}")
        syn_idx = 0
        logger.info("rho_GB before")
        initial_rho_values = []
        for syn_id, synapse in cell.synapses.items():
            print(synapse.hsynapse)
            #print(dir(synapse.hsynapse))
            print(f"printing Use: {synapse.hsynapse.Use}")
            print(f"printing Use_p: {synapse.hsynapse.Use_p}")
            print(f"printing Use_d: {synapse.hsynapse.Use_d}")
            print(f"printing gmax_AMPA: {synapse.hsynapse.gmax_AMPA}")
            print(f"printing gmax_p_AMPA: {synapse.hsynapse.gmax_p_AMPA}")
            print(f"printing gmax_d_AMPA: {synapse.hsynapse.gmax_d_AMPA}")
            print(f"printing rho0_GB: {synapse.hsynapse.rho0_GB}")
            print(f"printing rho_GB: {synapse.hsynapse.rho_GB}")
            initial_rho_values.append(synapse.hsynapse.rho0_GB)
        logger.info(f"initial_rho_values: {initial_rho_values}")
        for syn_id, synapse in cell.synapses.items():
            rho_val = rho_values[syn_idx]
            #logger.info(f"rho_val: {rho_val}")
            if rho_val >= 0.5:
                #logger.info(f"rho_val greater than 0.5: {rho_val}")
                synapse.hsynapse.rho0_GB = 1.0
                synapse.hsynapse.Use = synapse.hsynapse.Use_p
                #logger.info(f"synapse.hsynapse.Use: {synapse.hsynapse.Use}")
                #logger.info(f"synapse.hsynapse.Use_p: {synapse.hsynapse.Use_p}")
                synapse.hsynapse.gmax_AMPA = synapse.hsynapse.gmax_p_AMPA
            else:
                #logger.info(f"rho_val less than 0.5: {rho_val}")
                synapse.hsynapse.rho0_GB = 0.0
                synapse.hsynapse.Use = synapse.hsynapse.Use_d
                synapse.hsynapse.gmax_AMPA = synapse.hsynapse.gmax_d_AMPA
            syn_idx += 1
        
        # Log final rho values after configuration
        logger.info("rho_GB after configuration:")
        final_rho_values = []
        for syn_id, synapse in cell.synapses.items():
            print(f"Synapse {syn_id[1]}: rho0_GB = {synapse.hsynapse.rho0_GB}, rho_GB = {synapse.hsynapse.rho_GB}")
            final_rho_values.append(synapse.hsynapse.rho0_GB)
        logger.info(f"final_rho_values: {final_rho_values}")
        
        #disable plasticity
        for syn_id, synapse in cell.synapses.items():
            synapse.hsynapse.theta_d_GB = 1000000
            synapse.hsynapse.theta_p_GB = 1000000

        logger.info(f"Starting simulation for {pre_gid}->{post_gid}")
        bluecellulab.neuron.h.cvode_active(1)
        sim.run(240000, cvode=True)
        #sim.run(2, cvode=True)
        logger.info(f"Simulation completed for {pre_gid}->{post_gid}")
        
        logger.info(f"Extracting voltage trace for {pre_gid}->{post_gid}")
        t = np.array(sim.get_time())
        v = np.array(sim.get_voltage_trace((node_pop, post_gid)))
        
        # Save t, v, pre_spikes as pickle file with rho config and trial as filename
        data_dict = {
            't': t,
            'v': v,
            'rho_GB': np.array(syn_rec['rho_GB']),
            'pre_spikes': pre_spikes,
            'pre_gid': pre_gid,
            'post_gid': post_gid,
            'rho_config': rho_config,
            'trial': trial
        }
        
        
        # Ensure the directory exists
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        
        logger.info(f"Saving ephys data to {filepath}")
        with open(filepath, 'wb') as f:
            pickle.dump(data_dict, f)
        
        logger.info(f"Data saved successfully for {pre_gid}->{post_gid}, config {rho_config}, trial {trial}")
        avg_epsp_value = 0.0  # Placeholder value since we're not computing EPSP
        del sim
        return rho_config, avg_epsp_value

        
    except Exception as e:
        #print stack trace
        print(traceback.format_exc())
        logger.error(f"Error computing EPSP for {pre_gid}->{post_gid}, config {rho_config}, trial {trial}: {e}")
        return rho_config, 0.0

def count_synapses_between_cells(sim_config_path, pre_gid, post_gid, node_pop="S1nonbarrel_neurons"):
    try:
        sim = bluecellulab.CircuitSimulation(sim_config_path)
        sim.instantiate_gids([(node_pop, post_gid)], 
                           add_synapses=True, 
                           add_minis=False, 
                           add_pulse_stimuli=False,
                           intersect_pre_gids=[(node_pop, pre_gid)])
        
        cell = sim.cells[(node_pop, post_gid)]
        synapse_count = len(cell.synapses)
        
        logger.info(f"Cell pair {pre_gid}->{post_gid}: {synapse_count} synapses")
        del sim
        return synapse_count
        
    except Exception as e:
        logger.error(f"Error counting synapses for pair {pre_gid}->{post_gid}: {e}")
        return 0

def generate_binary_combinations(n_synapses):
    if n_synapses == 0:
        return []
    
    combinations = []
    for combination in itertools.product([0, 1], repeat=n_synapses):
        combinations.append(','.join(map(str, combination)))
    
    return combinations

def create_synapse_combinations_json(csv_path, sim_config_path, output_path, num_trials=10):
    logger.info("Extracting unique cell pairs from CSV...")
    cell_pairs = get_unique_cell_pairs(csv_path)
    logger.info(f"Found {len(cell_pairs)} unique cell pairs")
    logger.info(f"Running {num_trials} trials for each combination")
    
    result_dict = {}
    
    for i, (pre_gid, post_gid) in enumerate(cell_pairs):
        logger.info(f"Processing pair {i+1}/{len(cell_pairs)}: {pre_gid} -> {post_gid}")
        
        synapse_count = count_synapses_between_cells(sim_config_path, pre_gid, post_gid)
        
        if synapse_count > 0:
            combinations = generate_binary_combinations(synapse_count)
            cell_pair_key = f"pre_post_cell_{pre_gid}_{post_gid}"
            
            # Filter out combinations that already have existing files for all trials
            logger.info(f"Checking for existing files before processing {len(combinations)} combinations...")
            pending_combinations = []
            existing_results = {}
            
            for combo in combinations:
                safe_rho_config = combo.replace(',', '_')
                # Check if all trial files exist for this combination
                all_trials_exist = True
                for trial in range(num_trials):
                    filename = f"ephys_data_{pre_gid}_{post_gid}_{safe_rho_config}_trial_{trial}.pkl"
                    filepath = os.path.join(os.path.dirname(sim_config_path), '..', '..', '..', '..', 'ephys_data', filename)
                    if not os.path.exists(filepath):
                        all_trials_exist = False
                        break
                
                if all_trials_exist:
                    logger.info(f"Found all trial files for combination {combo}")
                    existing_results[combo] = 0.0  # Placeholder value
                else:
                    pending_combinations.append(combo)
            
            logger.info(f"Found {len(existing_results)} existing files, {len(pending_combinations)} combinations need processing")
            
            if pending_combinations:
                logger.info(f"Computing EPSP values for {len(pending_combinations)} combinations using multiprocessing...")
                
                # Log CPU information
                total_cpus = multiprocessing.cpu_count()
                # Use more processes for parallel execution - limit to reasonable number to avoid resource conflicts
                total_tasks = len(pending_combinations) * num_trials
                available_cpus = min(total_tasks, total_cpus - 2)  # Leave 2 cores free
                available_cpus = max(available_cpus, 1)  # At least 1 process
                logger.info(f"System has {total_cpus} CPU cores available")
                logger.info(f"Using {available_cpus} cores for {len(pending_combinations)} combinations × {num_trials} trials = {total_tasks} total tasks")
                
                # Create args for all combinations and all trials
                args_list = []
                for combo in pending_combinations:
                    for trial in range(num_trials):
                        args_list.append((sim_config_path, pre_gid, post_gid, combo, "S1nonbarrel_neurons", trial))
            
                # DEBUG: Use dummy function to test multiprocessing
                use_dummy = False  # Set to False to use real function
                
                logger.info(f"Created {len(args_list)} total tasks for multiprocessing")
                with multiprocessing.Pool(processes=available_cpus) as pool:
                    if use_dummy:
                        print(f"DEBUG: Using dummy function with {available_cpus} processes")
                        results = pool.map(dummy_test_function, args_list)
                    else:
                        logger.info(f"Starting multiprocessing with {available_cpus} processes and 66000s timeout")
                        logger.info(f"Task breakdown: {len(pending_combinations)} combinations × {num_trials} trials = {len(args_list)} total tasks")
                        try:
                            # Use map_async with timeout to avoid hanging
                            async_result = pool.map_async(get_epsp_value, args_list)
                            results = async_result.get(timeout=66000)  # 66000 second timeout
                            logger.info("Multiprocessing completed successfully")
                        except multiprocessing.TimeoutError:
                            logger.error("Multiprocessing timed out after 66000 seconds")
                            pool.terminate()
                            pool.join()
                            raise Exception("Multiprocessing timed out - possible resource conflict")
                
                # Combine existing results with newly computed results
                # Group results by combination (ignoring trial number for the summary)
                computed_results = {}
                for combo, epsp_val in results:
                    if combo not in computed_results:
                        computed_results[combo] = epsp_val
                all_results = {**existing_results, **computed_results}
            else:
                logger.info("All combinations already have existing files, skipping multiprocessing")
                all_results = existing_results
            
            result_dict[cell_pair_key] = all_results
            
            logger.info(f"Generated {len(combinations)} combinations for {synapse_count} synapses")
        else:
            logger.warning(f"No synapses found for pair {pre_gid} -> {post_gid}")
    
    logger.info(f"Writing results to {output_path}")
    with open(output_path, 'w') as f:
        json.dump(result_dict, f, indent=2)
    
    total_combinations = sum(len(combinations) for combinations in result_dict.values())
    logger.info(f"Complete! Generated {total_combinations} total combinations for {len(result_dict)} cell pairs")
    
    return result_dict

def run_custom_combinations(csv_path, sim_config_path, custom_configs, output_path, num_trials=10):
    """Run only the specified custom configurations"""
    logger.info("Running custom synapse combinations...")
    logger.info(f"Custom configs to run: {custom_configs}")
    logger.info(f"Running {num_trials} trials for each combination")
    
    cell_pairs = get_unique_cell_pairs(csv_path)
    logger.info(f"Found {len(cell_pairs)} unique cell pairs")
    
    result_dict = {}
    
    # Convert custom configs to string format
    custom_combinations = []
    for config in custom_configs:
        combo_str = ','.join([str(int(x)) for x in config])
        custom_combinations.append(combo_str)
    
    logger.info(f"Custom combinations: {custom_combinations}")
    
    for i, (pre_gid, post_gid) in enumerate(cell_pairs):
        logger.info(f"Processing pair {i+1}/{len(cell_pairs)}: {pre_gid} -> {post_gid}")
        
        synapse_count = count_synapses_between_cells(sim_config_path, pre_gid, post_gid)
        
        if synapse_count > 0:
            # Validate that custom configs match synapse count
            valid_combinations = []
            for combo in custom_combinations:
                if len(combo.split(',')) == synapse_count:
                    valid_combinations.append(combo)
                else:
                    logger.warning(f"Skipping config {combo} - length {len(combo.split(','))} doesn't match {synapse_count} synapses")
            
            if not valid_combinations:
                logger.warning(f"No valid custom combinations for {synapse_count} synapses")
                continue
                
            cell_pair_key = f"pre_post_cell_{pre_gid}_{post_gid}"
            
            logger.info(f"Running {len(valid_combinations)} custom combinations for {synapse_count} synapses")
            
            # Log CPU information
            total_cpus = multiprocessing.cpu_count()
            # Use more processes for parallel execution - limit to reasonable number to avoid resource conflicts
            available_cpus = min(len(valid_combinations) * num_trials, total_cpus - 2)  # Leave 2 cores free
            available_cpus = max(available_cpus, 1)  # At least 1 process
            logger.info(f"System has {total_cpus} CPU cores available")
            logger.info(f"Using {available_cpus} cores for {len(valid_combinations)} combinations × {num_trials} trials = {len(valid_combinations) * num_trials} total tasks")
            
            # Create args for all combinations and all trials
            args_list = []
            for combo in valid_combinations:
                for trial in range(num_trials):
                    args_list.append((sim_config_path, pre_gid, post_gid, combo, "S1nonbarrel_neurons", trial))
            
            logger.info(f"Created {len(args_list)} total tasks for multiprocessing")
            with multiprocessing.Pool(processes=available_cpus) as pool:
                logger.info(f"Starting multiprocessing with {available_cpus} processes")
                logger.info(f"Task breakdown: {len(valid_combinations)} combinations × {num_trials} trials = {len(args_list)} total tasks")
                try:
                    results = pool.map(get_epsp_value, args_list)
                    logger.info("Custom combinations completed successfully")
                except Exception as e:
                    logger.error(f"Error in multiprocessing: {e}")
                    raise
            
            # Group results by combination (ignoring trial number for the summary)
            computed_results = {}
            for combo, epsp_val in results:
                if combo not in computed_results:
                    computed_results[combo] = epsp_val
            result_dict[cell_pair_key] = computed_results
            
        else:
            logger.warning(f"No synapses found for pair {pre_gid} -> {post_gid}")
    
    logger.info(f"Writing results to {output_path}")
    with open(output_path, 'w') as f:
        json.dump(result_dict, f, indent=2)
    
    total_combinations = sum(len(combinations) for combinations in result_dict.values())
    logger.info(f"Complete! Generated {total_combinations} custom combinations for {len(result_dict)} cell pairs")
    
    return result_dict

def process_connection_type(base_dir, connection_type, cell_pair, custom_run_configs=None, num_trials=10):
    """Process a specific connection type and cell pair"""
    logger.info(f"Processing {connection_type} with cell pair {cell_pair}")
    
    # Construct CSV path based on connection type
    csv_path = base_dir / f"refitting_results/fitting/n1/seed19091997/index_{connection_type}.csv"
    
    # Construct simulation config candidates based on connection type and cell pair
    sim_config_candidates = [
        base_dir / f"refitting_results/fitting/n1/seed19091997/{connection_type}/simulations/{cell_pair}/2Hz_-10ms/simulation_config.json",
        base_dir / f"refitting_results/fitting/n1/seed19091997/{connection_type}/simulations/{cell_pair}/2Hz_5ms/simulation_config.json",
        base_dir / f"refitting_results/fitting/n1/seed19091997/{connection_type}/simulations/{cell_pair}/5Hz_5ms/simulation_config.json",
        base_dir / f"refitting_results/fitting/n1/seed19091997/{connection_type}/simulations/{cell_pair}/10Hz_10ms/simulation_config.json",
        base_dir / f"refitting_results/fitting/n1/seed19091997/{connection_type}/simulations/{cell_pair}/50Hz_10ms/simulation_config.json",
        base_dir / f"refitting_results/fitting/n1/seed19091997/simulation_config.json",
        base_dir / "configs/simulation_config.json",
        base_dir / "simulation_config.json"
    ]
    
    sim_config_path = None
    for candidate in sim_config_candidates:
        if candidate.exists():
            sim_config_path = candidate
            break
    
    if sim_config_path is None:
        logger.error(f"Could not find simulation_config.json for {connection_type} {cell_pair}. Please ensure it exists.")
        return None
    
    if not csv_path.exists():
        logger.error(f"Input CSV file not found: {csv_path}")
        return None
    
    # Create output path with connection type and cell pair
    safe_connection_type = connection_type.replace("_", "_")
    safe_cell_pair = cell_pair.replace("-", "_")
    
    if custom_run_configs is not None:
        output_path = base_dir / f"custom_synapse_combinations_{safe_connection_type}_{safe_cell_pair}.json"
        logger.info(f"CUSTOM MODE: Running only specified configurations for {connection_type}: {custom_run_configs}")
        logger.info(f"Output will be written to: {output_path}")
        return run_custom_combinations(str(csv_path), str(sim_config_path), custom_run_configs, str(output_path), num_trials)
    else:
        output_path = base_dir / f"synapse_combinations_{safe_connection_type}_{safe_cell_pair}.json"
        logger.info(f"FULL MODE: Running all possible combinations for {connection_type}")
        logger.info(f"Output will be written to: {output_path}")
        return create_synapse_combinations_json(str(csv_path), str(sim_config_path), str(output_path), num_trials)

def main():
    base_dir = Path(__file__).parent
    
    # Define connection types and cell pairs to process
    connection_configs = [
       {"connection_type": "L5TTPC_L5TTPC", "cell_pair": "205559-199162"},
        #{"connection_type": "L23PC_L5TTPC", "cell_pair": "15864-199162"}
    ]
    
    # Define your custom configurations here - SET TO None TO RUN ALL COMBINATIONS
    custom_run_configs = [
        [1, 1, 0, 1, 0, 0, 1, 1],
        [1, 1, 0, 1, 1, 0, 1, 1]
    ]
    # custom_run_configs = [
    #     [1, 1, 1, 0],
    #     [1, 1, 0, 0]
    # ]

    # custom_run_configs = None  # Uncomment this line to run ALL combinations instead
    
    all_results = {}
    
    # Process each connection type and cell pair
    for config in connection_configs:
        connection_type = config["connection_type"]
        cell_pair = config["cell_pair"]
        
        logger.info(f"Processing {connection_type} with cell pair {cell_pair}")
        
        # Check if data exists for this connection type
        csv_path = base_dir / f"refitting_results/fitting/n1/seed19091997/index_{connection_type}.csv"
        if not csv_path.exists():
            logger.warning(f"CSV file not found for {connection_type}: {csv_path}")
            logger.warning(f"Skipping {connection_type} - no index file available")
            continue
        
        # Check if simulation directory exists
        sim_dir = base_dir / f"refitting_results/fitting/n1/seed19091997/{connection_type}/simulations/{cell_pair}"
        if not sim_dir.exists():
            logger.warning(f"Simulation directory not found for {connection_type} {cell_pair}: {sim_dir}")
            logger.warning(f"Skipping {connection_type} {cell_pair} - no simulation data available")
            continue
        
        # Process this connection type
        result = process_connection_type(base_dir, connection_type, cell_pair, custom_run_configs, num_trials=10)
        if result is not None:
            all_results[f"{connection_type}_{cell_pair}"] = result
            logger.info(f"Successfully processed {connection_type} {cell_pair}")
        else:
            logger.error(f"Failed to process {connection_type} {cell_pair}")
    
    # Print summary
    logger.info(f"Processing complete! Successfully processed {len(all_results)} connection types:")
    for key in all_results.keys():
        logger.info(f"  - {key}")
    
    return all_results

if __name__ == "__main__":
    main()