
import pickle
import numpy as np
import os

file1_path = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/evaluation_results/Chindemi_params/sim_results/simulation_b7a94c1b3d8f_mrk97_08_180164_197248.pkl"
file2_path = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/evaluation_results/Chindemi_params_mod/sim_results/simulation_b7a94c1b3d8f_mrk97_08_180164_197248.pkl"

print(f"Comparing:\n1: {file1_path}\n2: {file2_path}")

with open(file1_path, "rb") as f:
    data1 = pickle.load(f)

with open(file2_path, "rb") as f:
    data2 = pickle.load(f)

# Compare keys
keys1 = set(data1.keys())
keys2 = set(data2.keys())

if keys1 != keys2:
    print("Different keys!")
    print(f"In 1 only: {keys1 - keys2}")
    print(f"In 2 only: {keys2 - keys1}")
else:
    print("Keys match.")

# Compare specific fields
fields_to_compare = ['t', 'v', 'rho_GB']

for field in fields_to_compare:
    if field in data1 and field in data2:
        val1 = data1[field]
        val2 = data2[field]
        
        if isinstance(val1, np.ndarray) and isinstance(val2, np.ndarray):
            if np.array_equal(val1, val2):
                print(f"Field '{field}' is IDENTICAL.")
            else:
                print(f"Field '{field}' is DIFFERENT.")
                diff = np.abs(val1 - val2)
                print(f"  Max difference: {np.max(diff)}")
                print(f"  Mean difference: {np.mean(diff)}")
        else:
            if val1 == val2:
                print(f"Field '{field}' is IDENTICAL.")
            else:
                print(f"Field '{field}' is DIFFERENT.")
                print(f"  Value 1: {val1}")
                print(f"  Value 2: {val2}")
    else:
        print(f"Field '{field}' missing in one or both files.")

# Check if rho_GB is actually different in a meaningful way
if 'rho_GB' in data1 and 'rho_GB' in data2:
    rho1 = data1['rho_GB']
    rho2 = data2['rho_GB']
    # rho_GB is likely a list of arrays (one per synapse) or a 2D array
    # Let's inspect the structure
    print(f"rho_GB type: {type(rho1)}")
    if isinstance(rho1, list):
        print(f"rho_GB length: {len(rho1)}")
        for i in range(min(len(rho1), 5)): # Check first 5 synapses
            r1 = rho1[i]
            r2 = rho2[i]
            if np.array_equal(r1, r2):
                 print(f"  Synapse {i} rho trace: IDENTICAL")
            else:
                 print(f"  Synapse {i} rho trace: DIFFERENT")
                 print(f"    Max diff: {np.max(np.abs(r1 - r2))}")

