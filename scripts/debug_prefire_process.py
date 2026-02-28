import os
import glob
import pickle
import json
import logging
import numpy as np
import multiprocessing
import bluecellulab
from bluepysnap import Simulation
import plastyfire.simulator as sim_module
from plastyfire.epg import ParamsGenerator

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("debug_prefire")

# Directories to debug
directories = [
    "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC/simulations/205559-199162/10Hz_5ms",
    "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC/simulations/205559-199162/10Hz_10ms"
]

# Constants (from simulator.py or inferred)
EXTRA_RECIPE_PATH = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/biodata/recipe.csv"
# Assuming default populations if not specified, but better to check config or use defaults
# simulator.py defaults: node_pop="S1nonbarrel_neurons", edge_pop="S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"
NODE_POP = "S1nonbarrel_neurons"
EDGE_POP = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"

def run_debug(workdir):
    logger.info(f"--- Debugging {workdir} ---")
    
    if not os.path.exists(workdir):
        logger.error(f"Directory not found: {workdir}")
        return

    # 1. Load Simulation Config to get GIDs and Timing
    sim_config_path = os.path.join(workdir, "simulation_config.json")
    if not os.path.exists(sim_config_path):
        logger.error(f"Config not found: {sim_config_path}")
        return

    try:
        sim = Simulation(sim_config_path)
        pre_gid = sim.node_sets.content["precell"]["node_id"][0]
        post_gid = sim.node_sets.content["postcell"]["node_id"][0]
        t_end = sim.config["run"]["tstop"]
        logger.info(f"Pre GID: {pre_gid}, Post GID: {post_gid}, t_end: {t_end}")
    except Exception as e:
        logger.error(f"Error loading simulation config: {e}")
        return

    # 2. Load c_pre and c_post from existing pickle
    c_pre = {}
    c_post = {}
    pkl_files = glob.glob(os.path.join(workdir, "*.pkl"))
    if not pkl_files:
        logger.error("No pickle files found to extract c_pre/c_post")
        return
    
    # Try to find the correct pickle or just use the first one and debug
    found_key = False
    target_key = 356546100 # The key that failed
    
    for pkl_file in pkl_files:
        logger.info(f"Inspecting {pkl_file}")
        try:
            with open(pkl_file, "rb") as f:
                data = pickle.load(f)
                if "synprop" in data:
                    syn_ids = data["synprop"]["synapseID"]
                    c_pre_vals = data["synprop"]["Cpre"]
                    c_post_vals = data["synprop"]["Cpost"]
                    
                    # Check if target key is in syn_ids
                    if target_key in syn_ids:
                        logger.info(f"Found target key {target_key} in {pkl_file}")
                        found_key = True
                        # Use this pickle
                        for i, syn_id in enumerate(syn_ids):
                            c_pre[syn_id] = c_pre_vals[i]
                            c_post[syn_id] = c_post_vals[i]
                        break
                    else:
                        logger.info(f"Target key {target_key} NOT found in {pkl_file}. IDs: {syn_ids[:5]}...")
        except Exception as e:
            logger.error(f"Error loading pickle {pkl_file}: {e}")

    if not c_pre:
        logger.warning("Could not find a pickle with the target key. Using the last checked pickle's data if available, or failing.")
        if not found_key and pkl_files:
             # Fallback to loading the first one just to see what happens, or maybe we should stop.
             # But we want to see the keys.
             pass

    # 3. Generate Synapse Parameters
    try:
        # We need the circuit from the simulation to initialize ParamsGenerator
        # sim.circuit is available from bluepysnap Simulation
        pgen = ParamsGenerator(sim.circuit, NODE_POP, EDGE_POP, EXTRA_RECIPE_PATH)
        syn_extra_params = pgen.generate_params(pre_gid, post_gid)
        logger.info("Generated synapse parameters")
    except Exception as e:
        logger.error(f"Error generating params: {e}")
        return

    # 4. Prepare Temporary Directory with Correct Seed
    import tempfile
    import shutil
    
    # We suspect the seed in prefire_simulation_config.json (876667) is different from 
    # the one used to generate the pickle (likely 690077 from simulation_config.json).
    # We will create a temp dir, copy the config, update the seed, and symlink other files.
    
    target_seed = 690077 # From simulation_config.json
    
    with tempfile.TemporaryDirectory() as temp_dir:
        logger.info(f"Created temp dir: {temp_dir}")
        
        # Copy and modify prefire_simulation_config.json
        src_config = os.path.join(workdir, "prefire_simulation_config.json")
        dst_config = os.path.join(temp_dir, "prefire_simulation_config.json")
        
        try:
            with open(src_config, "r") as f:
                config_data = json.load(f)
            
            logger.info(f"Original seed: {config_data['run']['random_seed']}")
            config_data['run']['random_seed'] = target_seed
            logger.info(f"Modified seed to: {target_seed}")
            
            with open(dst_config, "w") as f:
                json.dump(config_data, f, indent=4)
                
        except Exception as e:
            logger.error(f"Error modifying config: {e}")
            return

        # Symlink other files
        for filename in ["prefire_prespikes.h5", "node_sets.json", "prespikes.h5"]:
            src = os.path.join(workdir, filename)
            dst = os.path.join(temp_dir, filename)
            if os.path.exists(src):
                os.symlink(src, dst)
            else:
                logger.warning(f"File not found for symlinking: {src}")

        # 4. Prepare Synapse Parameters (Mapping pickle values to current keys)
        # Since we have a mismatch in IDs (pickle has local-like IDs, current circuit has global IDs),
        # we will map the values from the pickle to the current keys 1-to-1.
        
        logger.info("Instantiating circuit to get expected synapse IDs...")
        try:
            # We need to load mechanisms first (just in case)
            if not hasattr(bluecellulab.neuron.h, "CaDynamics_DC0"):
                bluecellulab.neuron.load_mechanisms("/project/ctb-emuller/dhuruva/DEES_cell_packages/")
                
            sim = bluecellulab.CircuitSimulation(dst_config)
            sim.instantiate_gids([(NODE_POP, post_gid)], add_synapses=True, 
                                 intersect_pre_gids=[(NODE_POP, pre_gid)])
            cell = sim.cells[(NODE_POP, post_gid)]
            
            current_syn_idx = []
            for syn_id, synapse in cell.synapses.items():
                current_syn_idx.append(syn_id[1])
            
            # Use _map_syn_idx to get global IDs
            df = sim_module._map_syn_idx(dst_config, post_gid, current_syn_idx, EDGE_POP)
            expected_keys = []
            for syn_id in current_syn_idx:
                global_id = df.loc[df["local_syn_idx"] == syn_id].index[0]
                expected_keys.append(global_id)
            
            logger.info(f"Expected keys: {expected_keys}")
            
            # Map pickle values to expected keys
            # We assume the number of synapses is roughly the same.
            # If not, we cycle or truncate.
            
            pickle_c_pre_values = list(c_pre.values())
            pickle_c_post_values = list(c_post.values())
            
            if not pickle_c_pre_values:
                logger.warning("No c_pre values from pickle! Using defaults.")
                pickle_c_pre_values = [0.0] * len(expected_keys)
                pickle_c_post_values = [0.0] * len(expected_keys)
            
            new_c_pre = {}
            new_c_post = {}
            
            for i, key in enumerate(expected_keys):
                val_idx = i % len(pickle_c_pre_values)
                new_c_pre[key] = pickle_c_pre_values[val_idx]
                new_c_post[key] = pickle_c_post_values[val_idx]
                
            logger.info(f"Constructed c_pre with {len(new_c_pre)} items")
            
        except Exception as e:
            logger.error(f"Error preparing params: {e}")
            import traceback
            traceback.print_exc()
            return

        # 5. Run the Process
        results = {} 
        fit_params = None 
        syn_rec_lst = None
        fastforward = 0
        fixhp = True 
        
        logger.info("Starting _runconnectedpair_prefire_process with modified config and mapped params...")
        try:
            sim_module._runconnectedpair_prefire_process(
                results, temp_dir, fit_params, syn_extra_params, pre_gid, post_gid, t_end,
                new_c_pre, new_c_post, syn_rec_lst, fastforward, NODE_POP, EDGE_POP, fixhp
            )
            
            logger.info("Process completed successfully")
            logger.info(f"Results keys: {results.keys()}")
            
            # Print some outputs as requested
            if "postspikes" in results:
                logger.info(f"Post Spikes: {results['postspikes']}")
            if "prespikes" in results:
                logger.info(f"Pre Spikes: {results['prespikes']}")
            if "v" in results:
                logger.info(f"Voltage trace length: {len(results['v'])}")
                logger.info(f"Voltage sample: {results['v'][:10]}")
                
        except Exception as e:
            logger.error(f"Error running process: {e}")
            import traceback
            traceback.print_exc()

if __name__ == "__main__":
    for d in directories:
        run_debug(d)
