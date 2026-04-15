import os
import json
import logging
import numpy as np
import tempfile
import shutil
import bluecellulab
from plastyfire.epg_dhuruva import ParamsGenerator
import plastyfire.simulator as sim_module

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("validate_seed")

# Base configuration
BASE_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC/simulations/205559-199162/10Hz_10ms"
EXTRA_RECIPE_PATH = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/biodata/recipe.csv"
NODE_POP = "S1nonbarrel_neurons"
EDGE_POP = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"

def run_c_pre_finder_with_seed(seed, workdir):
    logger.info(f"--- Running with seed {seed} ---")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        # Copy config and modify seed
        src_config = os.path.join(workdir, "simulation_config.json")
        dst_config = os.path.join(temp_dir, "simulation_config.json")
        
        with open(src_config, "r") as f:
            config = json.load(f)
        
        config["run"]["random_seed"] = seed
        # Ensure output dir exists in temp
        config["output"]["output_dir"] = os.path.join(temp_dir, "out")
        os.makedirs(config["output"]["output_dir"], exist_ok=True)
        
        # Point node_sets_file to original (it's absolute path in config usually, but let's check)
        # In the file I read earlier: "node_sets_file": ".../node_sets.json" (absolute)
        # So we don't need to copy it if it's absolute. 
        # But wait, if it's relative, we need to copy.
        # The cat output showed absolute path.
        
        with open(dst_config, "w") as f:
            json.dump(config, f, indent=4)
            
        # We also need to copy/symlink node_sets.json if the config uses relative path, 
        # but it used absolute.
        # However, bluecellulab might look for other things.
        
        # Setup Params
        try:
            from bluepysnap import Simulation
            snap_sim = Simulation(dst_config)
            pre_gid = snap_sim.node_sets.content["precell"]["node_id"][0]
            post_gid = snap_sim.node_sets.content["postcell"]["node_id"][0]
            
            sim = bluecellulab.CircuitSimulation(dst_config)
            
            # Generate params
            pgen = ParamsGenerator(sim.circuit, NODE_POP, EDGE_POP, EXTRA_RECIPE_PATH)
            syn_extra_params = pgen.generate_params(pre_gid, post_gid)
            
            # Run c_pre_finder
            # c_pre_finder(sim_config, fit_params, syn_extra_params, pre_gid, post_gid, ...)
            # It uses multiprocessing.Pool.
            
            logger.info("Running c_pre_finder...")
            c_pre = sim_module.c_pre_finder(
                dst_config, 
                fit_params=None, 
                syn_extra_params=syn_extra_params, 
                pre_gid=pre_gid, 
                post_gid=post_gid,
                node_pop=NODE_POP, 
                edge_pop=EDGE_POP, 
                fixhp=True
            )
            
            return c_pre
            
        except Exception as e:
            logger.error(f"Error with seed {seed}: {e}")
            import traceback
            traceback.print_exc()
            return None

def compare_results(res1, res2):
    if res1 is None or res2 is None:
        logger.error("Cannot compare, one result is None")
        return

    logger.info(f"Result 1 keys (Seed A): {list(res1.keys())}")
    logger.info(f"Result 2 keys (Seed B): {list(res2.keys())}")
    
    vals1 = sorted(list(res1.values()))
    vals2 = sorted(list(res2.values()))
    
    logger.info(f"Values 1 (sorted): {vals1}")
    logger.info(f"Values 2 (sorted): {vals2}")
    
    if np.allclose(vals1, vals2):
        logger.info("SUCCESS: The SET of c_pre values is identical (ignoring keys).")
    else:
        logger.info("FAILURE: The set of c_pre values differs.")

if __name__ == "__main__":
    # Seed A: 690077 (from simulation_config.json)
    # Seed B: 876667 (from prefire_simulation_config.json)
    
    res_a = run_c_pre_finder_with_seed(690077, BASE_DIR)
    res_b = run_c_pre_finder_with_seed(876667, BASE_DIR)
    
    compare_results(res_a, res_b)
