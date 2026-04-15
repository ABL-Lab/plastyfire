import pickle
import numpy as np
import os

pkl_path = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/trace_results/DHURUVA_PARAMS_V12/180164-197248/10Hz_10ms/simulation_traces.pkl"
try:
    with open(pkl_path, "rb") as f:
        d = pickle.load(f)

    rho = np.asarray(d["rho_GB"])
    if rho.ndim == 2 and rho.shape[0] > 100: rho = rho.T
    print("Final rho trace values:", rho[:, -1])
    print("Final rho binarized:   ", [1 if r[-1] >= 0.5 else 0 for r in rho])
except FileNotFoundError:
    print(f"File not found: {pkl_path}")
    print("Available directories:", os.listdir("trace_results/DHURUVA_PARAMS_V12/"))
