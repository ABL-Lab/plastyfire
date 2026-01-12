
import os
import sys
import numpy as np
import pickle
import plastyfire.simulator as sim

# Configuration
workdir = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC/simulations/180164-197248/10Hz_-10ms"
fastforward = 280000.0
recipe_orig = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/biodata/recipe.csv"
recipe_mod = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/biodata/recipe_mod.csv"

existing_orig_path = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/evaluation_results/Chindemi_params/sim_results/simulation_b7a94c1b3d8f_mrk97_08_180164_197248.pkl"
existing_mod_path = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/evaluation_results/Chindemi_params_mod/sim_results/simulation_b7a94c1b3d8f_mrk97_08_180164_197248.pkl"

fit_params = {
    "gamma_d_GB_GluSynapse": 101.5,
    "gamma_p_GB_GluSynapse": 216.2,
    "a00": 1.002,
    "a01": 1.954,
    "a10": 1.159,
    "a11": 2.483,
    "a20": 1.127,
    "a21": 2.456,
    "a30": 5.236,
    "a31": 1.782,
    "tau_effca_GB_GluSynapse": 278.3177658387
}

def load_pkl(path):
    print(f"Loading {path}...")
    with open(path, 'rb') as f:
        data = pickle.load(f)
    return data['v']

def run_sim(recipe_path, label):
    print(f"Running simulation with {label} recipe: {recipe_path}")
    results = sim.runconnectedpair(workdir, fit_params=fit_params, fastforward=fastforward, recipe_path=recipe_path)
    return results['v']

def compare_traces(v1, v2, name1, name2):
    # Truncate to min length
    min_len = min(len(v1), len(v2))
    v1 = np.array(v1[:min_len])
    v2 = np.array(v2[:min_len])
    
    diff = np.abs(v1 - v2)
    max_diff = np.max(diff)
    mean_diff = np.mean(diff)
    
    print(f"{name1} vs {name2}:")
    print(f"  Max Diff: {max_diff:.4f} mV")
    print(f"  Mean Diff: {mean_diff:.4f} mV")
    if max_diff < 1e-6:
        print("  -> IDENTICAL")
    else:
        print("  -> DIFFERENT")
    print("-" * 30)

print("--- Loading Existing Results ---")
v_exist_orig = load_pkl(existing_orig_path)
v_exist_mod = load_pkl(existing_mod_path)

print("\n--- Running Fresh Simulations ---")
v_fresh_orig = run_sim(recipe_orig, "ORIGINAL")
v_fresh_mod = run_sim(recipe_mod, "MODIFIED")

print("\n--- Comparisons ---")

# 1. Existing Orig vs Existing Mod (Confirming user's suspicion)
compare_traces(v_exist_orig, v_exist_mod, "Existing Original", "Existing Modified")

# 2. Fresh Orig vs Existing Orig (Verifying Original Baseline)
compare_traces(v_fresh_orig, v_exist_orig, "Fresh Original", "Existing Original")

# 3. Fresh Mod vs Existing Mod (Verifying Modified Baseline)
compare_traces(v_fresh_mod, v_exist_mod, "Fresh Modified", "Existing Modified")

# 4. Fresh Orig vs Existing Mod (Checking if Existing Mod is actually Original)
compare_traces(v_fresh_orig, v_exist_mod, "Fresh Original", "Existing Modified")

# 5. Fresh Mod vs Fresh Orig (Re-verifying the difference)
compare_traces(v_fresh_mod, v_fresh_orig, "Fresh Modified", "Fresh Original")
