#!/usr/bin/env python3
"""
Extract raw calcium traces (cai_CR) from the pre- and post-spike simulation 
protocols across multiple simulation pairs using parallel processing.
"""

import os
import sys
import pickle
import argparse
import multiprocessing
import glob

# Prevent bluecellulab from making too much noise
os.environ["GLOG_minloglevel"] = "3"

def worker(workdir, recipe_path, outdir, tau_eff):
    # Import inside the worker to ensure NEURON states are kept cleanly isolated to this process
    sys.path.append("/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire")
    from bluepysnap import Simulation
    from plastyfire.epg import ParamsGenerator
    from plastyfire.simulator import _c_pre_finder_process, _c_post_finder_process

    sim_config = os.path.join(workdir, "simulation_config.json")
    if not os.path.exists(sim_config):
        return f"Skipped {workdir} (no sim_config)"

    try:
        sim = Simulation(sim_config)
        pre_gid = sim.node_sets.content["precell"]["node_id"][0]
        post_gid = sim.node_sets.content["postcell"]["node_id"][0]
        node_pop = "S1nonbarrel_neurons"
        edge_pop = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"

        pgen = ParamsGenerator(sim.circuit, node_pop, edge_pop, recipe_path)
        syn_extra_params = pgen.generate_params(pre_gid, post_gid)

        width = sim.config["inputs"]["pulse0"]["width"]
        amp = sim.config["inputs"]["pulse0"]["amp_start"]
        stimulus = {"nspikes": 1, "freq": 0.1, "width": width, "offset": 1000, "amp": amp}

        # We pass ref_params to ensure the SIMULATED c_pre/c_post (effcai_GB) 
        # computed by NEURON matches our analytical baseline assumptions (no CICR).
        ref_params = {
            "tau_effca_GB_GluSynapse": tau_eff,
            "Vmax_CICR_GluSynapse": 0.0,
        }
        res_pre = _c_pre_finder_process(sim_config, ref_params, syn_extra_params, pre_gid, post_gid, node_pop, edge_pop, True)
        res_post = _c_post_finder_process(sim_config, ref_params, syn_extra_params, pre_gid, post_gid, stimulus, node_pop, edge_pop, True)

        out_dict = {
            "pre": res_pre,
            "post": res_post,
            "pre_gid": pre_gid,
            "post_gid": post_gid,
        }

        pair_name = os.path.basename(os.path.dirname(workdir))
        if pair_name == "single_target":
            pair_name = "custom"
        output_file = os.path.join(outdir, f"{pair_name}_threshold_traces.pkl")
        with open(output_file, "wb") as f:
            pickle.dump(out_dict, f)
        
        return f"Success: {output_file}"
    except Exception as e:
        return f"Error in {workdir}: {e}"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract raw cai_CR traces across multiple directories.")
    parser.add_argument("rootdir", help="Root directory containing simulation subdirectories or a specific workdir.")
    parser.add_argument("--outdir", "-o", required=True, help="Central directory to save the extracted traces.")
    parser.add_argument("--recipe-path", "-r", default="/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/biodata/recipe.csv", help="Path to custom recipe.csv file")
    parser.add_argument("--cores", "-c", type=int, default=40, help="Number of CPU cores to use")
    parser.add_argument("--tau", type=float, default=278.318, help="Tau effca value (ms) to compute simulated Cpre baseline.")
    
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    # We only want to run this once per pair, not for every protocol.
    # The structure is typically: rootdir / PairDir / ProtocolDir / simulation_config.json
    # Examples:
    #   .../180164-197248/10Hz_10ms/simulation_config.json
    #   .../180164-197248/10Hz_-10ms/simulation_config.json
    
    pair_to_workdir = {}
    
    if os.path.exists(os.path.join(args.rootdir, "simulation_config.json")):
        # If user passed a specific protocol dir directly
        pair_to_workdir["single_target"] = args.rootdir
    else:
        # Find all directories containing simulation_config.json for the 10Hz_10ms protocol explicitly
        for path in glob.glob(os.path.join(args.rootdir, "**", "10Hz_10ms", "simulation_config.json"), recursive=True):
            workdir = os.path.dirname(path)
            # The pair directory is the parent of the workdir
            pair_dir = os.path.basename(os.path.dirname(workdir))
            
            # If we haven't seen this pair yet, save the first protocol dir we find 
            # to run the cpre/cpost extraction on.
            if pair_dir not in pair_to_workdir:
                pair_to_workdir[pair_dir] = workdir
            
    if not pair_to_workdir:
        print(f"No simulation_config.json found in {args.rootdir}")
        sys.exit(1)
        
    print(f"Found {len(pair_to_workdir)} unique simulation pairs. Starting pool with {args.cores} cores...")
    
    # maxtasksperchild=1 is absolutely critical! 
    # It ensures that NEURON global state doesn't build up/crash the process across iterations
    with multiprocessing.Pool(processes=args.cores, maxtasksperchild=1) as pool:
        results = pool.starmap(worker, [(w, args.recipe_path, args.outdir, args.tau) for w in pair_to_workdir.values()])
        
    for r in results:
        print(r)
