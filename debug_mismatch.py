import os
import glob
import pickle
import json
import logging
import bluecellulab

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("debug_mismatch")

workdir = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC/simulations/205559-199162/10Hz_10ms"

def check():
    # 1. Load a pickle and check IDs
    pkl_files = glob.glob(os.path.join(workdir, "*.pkl"))
    if not pkl_files:
        logger.error("No pickles")
        return
    
    pkl_file = pkl_files[0]
    logger.info(f"Checking pickle: {pkl_file}")
    
    pickle_ids = set()
    with open(pkl_file, "rb") as f:
        data = pickle.load(f)
        if "synprop" in data:
            ids = data["synprop"]["synapseID"]
            pickle_ids.update(ids)
            logger.info(f"Pickle has {len(ids)} synapses. Sample: {ids[:5]}")
            if 356546100 in ids:
                logger.info("Target ID 356546100 FOUND in pickle")
            else:
                logger.info("Target ID 356546100 NOT found in pickle")
    
    # 2. Instantiate circuit with seed 690077
    sim_config_path = os.path.join(workdir, "simulation_config.json")
    with open(sim_config_path, "r") as f:
        config = json.load(f)
    
    seed = config["run"]["random_seed"]
    logger.info(f"Config seed: {seed}")
    
    # Use bluepysnap to get GIDs
    from bluepysnap import Simulation
    snap_sim = Simulation(sim_config_path)
    pre_gid = snap_sim.node_sets.content["precell"]["node_id"][0]
    post_gid = snap_sim.node_sets.content["postcell"]["node_id"][0]
    logger.info(f"Pre: {pre_gid}, Post: {post_gid}")
    
    sim = bluecellulab.CircuitSimulation(sim_config_path)
    node_pop = "S1nonbarrel_neurons"
    
    sim.instantiate_gids([(node_pop, post_gid)], add_synapses=True, 
                         intersect_pre_gids=[(node_pop, pre_gid)])
    
    cell = sim.cells[(node_pop, post_gid)]
    circuit_ids = set()
    for syn_id, synapse in cell.synapses.items():
        try:
            # Check what attributes are available
            # logger.info(dir(synapse.hsynapse))
            sid = int(synapse.hsynapse.synapseID)
            circuit_ids.add(sid)
        except:
            pass
            
    logger.info(f"Circuit has {len(circuit_ids)} synapses. Sample: {list(circuit_ids)[:5]}")
    
    if 356546100 in circuit_ids:
        logger.info("Target ID 356546100 FOUND in circuit")
    else:
        logger.info("Target ID 356546100 NOT found in circuit")
        
    # Intersection
    common = pickle_ids.intersection(circuit_ids)
    logger.info(f"Common IDs: {len(common)}")

if __name__ == "__main__":
    check()
