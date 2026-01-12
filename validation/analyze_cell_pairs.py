import os
import pandas as pd
from bluepysnap import Circuit
from collections import Counter

# Configuration
SIM_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC/simulations"
CIRCUIT_CONFIG = "/project/ctb-emuller/datasets/SSCx_plasticity/O1_2023a_Ecker/circuit_config.json"
NODE_POP = "S1nonbarrel_neurons"

def analyze_pairs():
    # Load circuit
    print(f"Loading circuit from {CIRCUIT_CONFIG}...")
    circuit = Circuit(CIRCUIT_CONFIG)
    nodes = circuit.nodes[NODE_POP]

    # Get list of simulation directories
    if not os.path.exists(SIM_DIR):
        print(f"Error: Directory {SIM_DIR} does not exist.")
        return

    sim_dirs = [d for d in os.listdir(SIM_DIR) if os.path.isdir(os.path.join(SIM_DIR, d)) and '-' in d]
    print(f"Found {len(sim_dirs)} simulation directories.")

    pair_counts = Counter()
    
    print("Analyzing pairs...")
    for d in sim_dirs:
        try:
            pre_gid, post_gid = map(int, d.split('-'))
            
            # Get m-types
            # bluepysnap nodes.get returns a Series or DataFrame. 
            # We need to handle potential indexing differences (0-based vs 1-based if any, but usually GIDs are direct indices if using .get with ids)
            # However, bluepysnap GIDs are usually 0-indexed in the file but 1-indexed in SONATA? 
            # Let's assume standard bluepysnap usage: get(gid, property)
            
            pre_mtype = nodes.get(pre_gid, properties="mtype")
            post_mtype = nodes.get(post_gid, properties="mtype")
            
            pair_key = f"{pre_mtype},{post_mtype}"
            pair_counts[pair_key] += 1
            
        except ValueError:
            print(f"Skipping invalid directory name: {d}")
        except Exception as e:
            print(f"Error processing {d}: {e}")

    # Report results
    print("\n--- Pair Analysis Results ---")
    total_pairs = sum(pair_counts.values())
    print(f"Total valid pairs analyzed: {total_pairs}")
    
    for pair, count in pair_counts.most_common():
        pre, post = pair.split(',')
        proportion = (count / total_pairs) * 100
        print(f"{pre},{post} = {count} pairs ({proportion:.1f}%)")

if __name__ == "__main__":
    analyze_pairs()
