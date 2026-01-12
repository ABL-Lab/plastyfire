#!/usr/bin/env python3

import os
import json
import pandas as pd
import logging
import numpy as np
import pickle
import argparse
import traceback
from pathlib import Path
import bluecellulab
from libsonata import SpikeReader
import plastyfire.ephysutils as ephysutils
from conntility.io.synapse_report import get_presyn_mapping
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


def get_epsp_value(sim_config_path, pre_gid, post_gid, rho_config, node_pop, trial, output_dir=None, no_cache=False, recipe_path=None, synapse_ids_str=None, fit_params=None):
    """Run simulation and save pkl file with ephys data"""
    
    # Create filename using rho config and trial number (replace commas with underscores for valid filename)
    safe_rho_config = rho_config.replace(',', '_')
    filename = f"ephys_data_{pre_gid}_{post_gid}_{safe_rho_config}.pkl"
    
    if output_dir:
        filepath = os.path.join(output_dir, filename)
    else:
        filepath = os.path.join(os.path.dirname(sim_config_path), '..', '..', '..', '..', 'ephys_data', filename)
    
    # Check if file already exists
    if os.path.exists(filepath):
        if no_cache:
            logger.info(f"File exists but no_cache is True, removing: {filepath}")
            try:
                os.remove(filepath)
            except OSError:
                pass
        else:
            logger.info(f"File already exists, skipping simulation: {filepath}")
            return rho_config, 0.0  # Return placeholder value since we're not computing EPSP anyway
    
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
    # bluecellulab.neuron.h.tau_effca_GB_GluSynapse = 278.3177658387
    # bluecellulab.neuron.h.gamma_d_GB_GluSynapse = 101.5387594661
    # bluecellulab.neuron.h.gamma_p_GB_GluSynapse = 216.1841700668
    
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
    
    # Use custom recipe path if provided, otherwise use default
    recipe_file = recipe_path if recipe_path else "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/biodata/recipe.csv"
    
    edge_pop = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"
    pgen = ParamsGenerator(bluepysnap_sim.circuit, node_pop, edge_pop, recipe_file)
    syn_extra_params = pgen.generate_params(pre_gid, post_gid)

    df = _map_syn_idx(sim_config_path, post_gid, syn_idx, edge_pop)
    syn_props = {key: [] for key in SYNPROPS}
    for syn_id, synapse in cell.synapses.items():
        # Reduced logging: removed verbose synapse debug
        # logger.debug("Configuring synapse %d", syn_id[1])
        if syn_extra_params is not None:  # configure local parameters
            global_syn_id = df.loc[df["local_syn_idx"] == syn_id[1]].index[0]
            _set_local_params(synapse, fit_params, syn_extra_params[global_syn_id])
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
    
    # Parse synapse IDs if provided
    target_synapse_ids = None
    if synapse_ids_str:
        try:
            target_synapse_ids = [int(x) for x in synapse_ids_str.split(",")]
            logger.info(f"Using synapse ID matching with {len(target_synapse_ids)} IDs")
        except ValueError:
            logger.error(f"Failed to parse synapse IDs: {synapse_ids_str}")
    
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
        current_syn_id = syn_id[1]
        
        # Determine which rho value to use
        rho_val = 0 # Default to 0 (depressed) if not found
        
        if target_synapse_ids:
            # Match by ID
            if current_syn_id in target_synapse_ids:
                idx = target_synapse_ids.index(current_syn_id)
                if idx < len(rho_values):
                    rho_val = rho_values[idx]
                else:
                    logger.warning(f"Synapse ID {current_syn_id} found at index {idx} but rho_values has length {len(rho_values)}")
            else:
                logger.warning(f"Synapse ID {current_syn_id} not found in target IDs. Defaulting to 0.")
        else:
            # Fallback to index-based matching (legacy behavior)
            if syn_idx < len(rho_values):
                rho_val = rho_values[syn_idx]
            else:
                logger.warning(f"Synapse index {syn_idx} out of range for rho_values (len={len(rho_values)}). Defaulting to 0.")
        
        if rho_val >= 0.5:
            synapse.hsynapse.rho0_GB = 1.0
            synapse.hsynapse.Use = synapse.hsynapse.Use_p
            synapse.hsynapse.gmax_AMPA = synapse.hsynapse.gmax_p_AMPA
        else:
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
    
    del sim
    return None

def main():
    parser = argparse.ArgumentParser(description='Run simulation and save pkl files')
    parser.add_argument('sim_config_path', help='Path to simulation config file')
    parser.add_argument('pre_gid', type=int, help='Presynaptic GID')
    parser.add_argument('post_gid', type=int, help='Postsynaptic GID')
    parser.add_argument('rho_config', help='Rho configuration (comma-separated 0s and 1s)')
    parser.add_argument('--node_pop', default='S1nonbarrel_neurons', help='Node population name')
    parser.add_argument('--trial', type=int, default=0, help='Trial number')
    parser.add_argument('--output-dir', default=None, help='Output directory for ephys data')
    parser.add_argument('--no-cache', action='store_true', help='Ignore existing files and force re-computation')
    parser.add_argument('--recipe-path', default=None, help='Path to custom recipe.csv file')
    parser.add_argument('--synapse-ids', default=None, help='Comma-separated list of synapse IDs to match rho values')
    
    # Add fit params arguments
    fit_param_names = [
        "gamma_d_GB_GluSynapse", "gamma_p_GB_GluSynapse",
        "a00", "a01", "a10", "a11", "a20", "a21", "a30", "a31"
    ]
    for param in fit_param_names:
        parser.add_argument(f"--{param}", type=float, default=None, help=f"Fit parameter {param}")
    
    args = parser.parse_args()
    
    # Collect fit params
    fit_params = {}
    for param in fit_param_names:
        val = getattr(args, param)
        if val is not None:
            fit_params[param] = val
    
    logger.info(f"Running simulation with args: {args}")
    logger.info(f"Fit params: {fit_params}")
    
    result = get_epsp_value(
        args.sim_config_path,
        args.pre_gid,
        args.post_gid,
        args.rho_config,
        args.node_pop,
        args.trial,
        args.output_dir,
        args.no_cache,
        args.recipe_path,
        args.synapse_ids,
        fit_params
    )
    
    logger.info(f"Simulation completed with result: {result}")


if __name__ == "__main__":
    main()