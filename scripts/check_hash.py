import hashlib
import json

# Hardcoded Chindemi_params from evaluate_best_solution.py
Chindemi_params = {
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
}

FIT_PARAM_NAMES = [
    "gamma_d_GB_GluSynapse",
    "gamma_p_GB_GluSynapse",
    "a00",
    "a01",
    "a10",
    "a11",
    "a20",
    "a21",
    "a30",
    "a31",
]

param_values = [Chindemi_params[name] for name in FIT_PARAM_NAMES]
# Match evaluate_best_solution.py: param_hash = hashlib.md5(str(param_values).encode()).hexdigest()[:12]
param_hash = hashlib.md5(str(param_values).encode()).hexdigest()[:12]
print(f"Calculated hash (str): {param_hash}")
