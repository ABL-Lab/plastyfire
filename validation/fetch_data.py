import os
import bluepysnap
import pandas as pd
import numpy as np
from tqdm import tqdm

# Configuration
circuit_path = "/home/dhuruva/projects/ctb-emuller/datasets/SSCx_plasticity/O1_2023a_Ecker/circuit_config.json"
# simulations_dir = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC/simulations"
# output_file = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/validation/synaptic_params.csv"

# L23PC_L5TTPC
simulations_dir = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/refitting_results/fitting/n100/seed19091997/L23PC_L5TTPC/simulations"
output_file = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/validation/synaptic_params_L23PC_L5TTPC.csv"

# Load circuit
print(f"Loading circuit from {circuit_path}")
circuit = bluepysnap.Circuit(circuit_path)
edge = circuit.edges['S1nonbarrel_neurons__S1nonbarrel_neurons__chemical']

# Get list of connections from simulation directories
print(f"Scanning {simulations_dir} for connections...")
connections = []
if os.path.exists(simulations_dir):
    for item in os.listdir(simulations_dir):
        if os.path.isdir(os.path.join(simulations_dir, item)):
            try:
                # Assuming format source-target
                parts = item.split('-')
                if len(parts) == 2:
                    source_id = int(parts[0])
                    target_id = int(parts[1])
                    connections.append((source_id, target_id))
            except ValueError:
                print(f"Skipping invalid directory name: {item}")

print(f"Found {len(connections)} connections.")

# Fetch parameters
properties = [
    'u_syn',
    'depression_time',
    'facilitation_time',
    'n_rrp_vesicles',
    'conductance',
    'volume_CR',
    '@source_node',
    '@target_node'
]

all_synapses = []

print("Fetching synaptic parameters using pair_edges...")
for source_id, target_id in tqdm(connections):
    try:
        # pair_edges returns a DataFrame if properties is a list
        syns = edge.pair_edges(source_id, target_id, properties=properties)
        
        if not syns.empty:
            all_synapses.append(syns)
    except Exception as e:
        print(f"Error fetching for {source_id}-{target_id}: {e}")

if all_synapses:
    df = pd.concat(all_synapses, ignore_index=True)
    
    # Calculate G_NMDAR
    df['g_nmdar'] = 0.55 * df['conductance']
    
    print(f"Total synapses fetched: {len(df)}")
    print(df.head())
    
    df.to_csv(output_file, index=False)
    print(f"Saved parameters to {output_file}")
else:
    print("No synapses found matching the desired pairs!")
