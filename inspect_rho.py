import pickle
import numpy as np
import sys

filepath = sys.argv[1]
print(f"Inspecting {filepath}")

with open(filepath, "rb") as f:
    data = pickle.load(f)

rho_gb = data["rho_GB"]
print(f"Type of rho_GB: {type(rho_gb)}")

if isinstance(rho_gb, list):
    print(f"Length of rho_GB list: {len(rho_gb)}")
    if len(rho_gb) > 0:
        print(f"Type of first element: {type(rho_gb[0])}")
        if hasattr(rho_gb[0], "shape"):
            print(f"Shape of first element: {rho_gb[0].shape}")
        elif isinstance(rho_gb[0], list):
             print(f"Length of first element: {len(rho_gb[0])}")

if isinstance(rho_gb, np.ndarray):
    print(f"Shape of rho_GB array: {rho_gb.shape}")

# Check synprop
if "synprop" in data:
    synprop = data["synprop"]
    if "synapseID" in synprop:
        print(f"Number of synapseIDs: {len(synprop['synapseID'])}")
