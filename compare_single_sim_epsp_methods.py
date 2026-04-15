#!/usr/bin/env python3
"""Compare basis-predicted and full-trace EPSP ratios for one simulation file.

This validator isolates the difference between the two STDP curve paths:
  - basis method: infer EPSP ratio from initial/final rho_GB states via basis CSVs
  - full method: measure EPSP ratio directly from the recorded voltage trace

Example:
  python compare_single_sim_epsp_methods.py       --sim-path /home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/full_trace_results/CHINDEMI_PARAMS/184033-189853/10Hz_10ms/simulation_traces.pkl       --params chindemi
"""

import argparse
import os
import pickle
import sys

import numpy as np

REPO_DIR = os.path.dirname(os.path.abspath(__file__))
NEW_FITTING_DIR = os.path.join(REPO_DIR, "new_fitting")

if NEW_FITTING_DIR not in sys.path:
    sys.path.insert(0, NEW_FITTING_DIR)
if REPO_DIR not in sys.path:
    sys.path.insert(0, REPO_DIR)

from effcai_to_epsp import fetch_epsp_ratio
from plastyfire.ephysutils import Experiment
from submit_l5ttpc_traces import (
    CHINDEMI_PARAMS,
    DHURUVA_PARAMS,
    DHURUVA_PARAMS_V2,
    DHURUVA_PARAMS_V3,
    DHURUVA_PARAMS_V4,
    DHURUVA_PARAMS_V5,
    DHURUVA_PARAMS_V6,
    DHURUVA_PARAMS_V7,
    DHURUVA_PARAMS_V8,
    DHURUVA_PARAMS_V9,
    DHURUVA_PARAMS_V10,
    DHURUVA_PARAMS_V11,
    DHURUVA_PARAMS_V12,
    DHURUVA_PARAMS_V13,
    DHURUVA_PARAMS_V14,
    DHURUVA_PARAMS_V15,
    DHURUVA_PARAMS_V16,
    DHURUVA_PARAMS_V17,
    DHURUVA_PARAMS_V18,
    DHURUVA_PARAMS_V19,
)

DEFAULT_BASIS_DIRS = [
    os.path.join(REPO_DIR, "basis_results"),
    os.path.join(REPO_DIR, "basis_results_old"),
]
C01_DURATION_MIN = 4.0
C02_DURATION_MIN = 4.0
TEST_PULSE_PERIOD_S = 4.0
N_EPSP = 60

PARAMS_MAPPING = {
    "v1": DHURUVA_PARAMS,
    "v2": DHURUVA_PARAMS_V2,
    "v3": DHURUVA_PARAMS_V3,
    "v4": DHURUVA_PARAMS_V4,
    "v5": DHURUVA_PARAMS_V5,
    "v6": DHURUVA_PARAMS_V6,
    "v7": DHURUVA_PARAMS_V7,
    "v8": DHURUVA_PARAMS_V8,
    "v9": DHURUVA_PARAMS_V9,
    "v10": DHURUVA_PARAMS_V10,
    "v11": DHURUVA_PARAMS_V11,
    "v12": DHURUVA_PARAMS_V12,
    "v13": DHURUVA_PARAMS_V13,
    "v14": DHURUVA_PARAMS_V14,
    "v15": DHURUVA_PARAMS_V15,
    "v16": DHURUVA_PARAMS_V16,
    "v17": DHURUVA_PARAMS_V17,
    "v18": DHURUVA_PARAMS_V18,
    "v19": DHURUVA_PARAMS_V19,
    "chindemi": CHINDEMI_PARAMS,
}


def resolve_basis_dir(basis_dir=None):
    candidates = [basis_dir] if basis_dir is not None else DEFAULT_BASIS_DIRS
    for candidate in candidates:
        if candidate and os.path.isdir(candidate):
            return candidate
    checked = ", ".join(candidates)
    raise FileNotFoundError(f"No basis directory found. Checked: {checked}")


def infer_protocol_from_path(sim_path):
    freq_dt = os.path.basename(os.path.dirname(sim_path))
    pair_name = os.path.basename(os.path.dirname(os.path.dirname(sim_path)))
    try:
        pre_gid_str, post_gid_str = pair_name.split("-", 1)
        pre_gid = int(pre_gid_str)
        post_gid = int(post_gid_str)
    except ValueError as exc:
        raise ValueError(f"Could not parse pair ids from '{pair_name}'") from exc

    try:
        freq_hz_str, dt_ms_str = freq_dt.split("_", 1)
        freq_hz = float(freq_hz_str.replace("Hz", ""))
        dt_ms = float(dt_ms_str.replace("ms", ""))
    except ValueError as exc:
        raise ValueError(f"Could not parse frequency/dt from '{freq_dt}'") from exc

    return pre_gid, post_gid, pair_name, freq_hz, dt_ms


def normalize_rho_trace(rho_gb_raw):
    if isinstance(rho_gb_raw, dict):
        rho_trace = np.array(list(rho_gb_raw.values()), dtype=float)
    else:
        rho_trace = np.asarray(rho_gb_raw, dtype=float)
        if rho_trace.ndim == 2 and rho_trace.shape[0] > rho_trace.shape[1]:
            rho_trace = rho_trace.T
        elif rho_trace.ndim == 1:
            rho_trace = rho_trace.reshape(1, -1)
    return rho_trace


def compute_rho_summary(rho_trace):
    initial = np.array([1 if rho[0] >= 0.5 else 0 for rho in rho_trace], dtype=int)
    final = np.array([1 if rho[-1] >= 0.5 else 0 for rho in rho_trace], dtype=int)
    rho_end = np.array([float(rho[-1]) for rho in rho_trace], dtype=float)
    return {
        "initial_rho": initial.tolist(),
        "final_rho": final.tolist(),
        "rho_end": rho_end.tolist(),
        "n_0to1": int(np.sum((initial == 0) & (final == 1))),
        "n_0to0": int(np.sum((initial == 0) & (final == 0))),
        "n_1to1": int(np.sum((initial == 1) & (final == 1))),
        "n_1to0": int(np.sum((initial == 1) & (final == 0))),
        "n_synapses": int(len(initial)),
    }


def compute_basis_result(sim_data, pre_gid, post_gid, freq_hz, dt_ms, params, basis_dir, ratio_method):
    if "rho_GB" not in sim_data:
        raise KeyError("Simulation data is missing 'rho_GB', required for basis comparison.")

    rho_trace = normalize_rho_trace(sim_data["rho_GB"])
    rho_summary = compute_rho_summary(rho_trace)
    basis_result = fetch_epsp_ratio(
        pre_gid,
        post_gid,
        freq_hz,
        dt_ms,
        rho_summary["initial_rho"],
        rho_summary["final_rho"],
        basis_dir=basis_dir,
        n_trials=0,
        params=params,
        ratio_method=ratio_method,
    )
    return basis_result, rho_summary


def compute_full_result(sim_data):
    required = ["t", "v", "prespikes"]
    missing = [key for key in required if key not in sim_data]
    if missing:
        raise KeyError("Simulation data is missing required keys: " + ", ".join(missing))

    experiment = Experiment(
        data=sim_data,
        c01duration=C01_DURATION_MIN,
        c02duration=C02_DURATION_MIN,
        period=TEST_PULSE_PERIOD_S,
    )
    epsp_before, epsp_after, epsp_ratio, epsp_before_std, epsp_after_std = experiment.compute_epsp_ratio(
        n=N_EPSP,
        full=True,
    )
    return {
        "epsp_before_mean": epsp_before,
        "epsp_before_std": epsp_before_std,
        "epsp_after_mean": epsp_after,
        "epsp_after_std": epsp_after_std,
        "ratio_mean": epsp_ratio,
        "n_epsp": N_EPSP,
    }


def print_report(sim_path, pair_name, freq_hz, dt_ms, basis_dir, ratio_method, basis_result, full_result, rho_summary):
    basis_ratio = float(basis_result["ratio_mean"])
    full_ratio = float(full_result["ratio_mean"])
    abs_diff = abs(full_ratio - basis_ratio)
    rel_diff = np.nan if basis_ratio == 0 else 100.0 * (full_ratio - basis_ratio) / basis_ratio

    print("=" * 76)
    print("Single-Simulation EPSP Validation: basis vs full")
    print("=" * 76)
    print(f"Simulation file : {sim_path}")
    print(f"Pair            : {pair_name}")
    print(f"Frequency       : {freq_hz:.3f} Hz")
    print(f"Δt              : {dt_ms:.3f} ms")
    print(f"Basis dir       : {basis_dir}")
    print(f"Basis method    : fetch_epsp_ratio(..., ratio_method='{ratio_method}')")
    print(f"Full method     : Experiment.compute_epsp_ratio(n={N_EPSP}, method='amplitude')")

    print("\nRho state summary")
    print("-" * 76)
    print(f"Synapses        : {rho_summary['n_synapses']}")
    print(f"0 -> 1          : {rho_summary['n_0to1']}")
    print(f"0 -> 0          : {rho_summary['n_0to0']}")
    print(f"1 -> 1          : {rho_summary['n_1to1']}")
    print(f"1 -> 0          : {rho_summary['n_1to0']}")
    print(f"Final rho vals  : {', '.join(f'{val:.4f}' for val in rho_summary['rho_end'])}")

    print("\nBasis extrapolation result")
    print("-" * 76)
    print(f"EPSP before     : {basis_result['epsp_before_mean']:.8f} ± {basis_result['epsp_before_std']:.8f}")
    print(f"EPSP after      : {basis_result['epsp_after_mean']:.8f} ± {basis_result['epsp_after_std']:.8f}")
    print(f"EPSP ratio      : {basis_ratio:.8f}")

    print("\nFull trace result")
    print("-" * 76)
    print(f"EPSP before     : {full_result['epsp_before_mean']:.8f} ± {full_result['epsp_before_std']:.8f}")
    print(f"EPSP after      : {full_result['epsp_after_mean']:.8f} ± {full_result['epsp_after_std']:.8f}")
    print(f"EPSP ratio      : {full_ratio:.8f}")

    print("\nDifference")
    print("-" * 76)
    print(f"Absolute diff   : {abs_diff:.8f}")
    if np.isnan(rel_diff):
        print("Relative diff   : undefined (basis ratio is zero)")
    else:
        print(f"Relative diff   : {rel_diff:.3f}%  (full - basis) / basis")

    print("\nInterpretation")
    print("-" * 76)
    print("Basis uses only binary initial/final rho states and the pair basis CSV.")
    print("Full uses the recorded voltage trace and the actual measured C01/C02 EPSPs.")
    print("If these differ, the gap comes from the method itself, not from batch aggregation.")
    print("=" * 76)


def main():
    parser = argparse.ArgumentParser(description="Compare basis and full EPSP ratios for one simulation_traces.pkl")
    parser.add_argument("--sim-path", required=True, help="Path to simulation_traces.pkl")
    parser.add_argument(
        "--params",
        choices=sorted(PARAMS_MAPPING.keys()),
        default="chindemi",
        help="Parameter set passed to fetch_epsp_ratio (default: chindemi)",
    )
    parser.add_argument(
        "--basis-dir",
        default=None,
        help="Directory containing basis_<pre>_<post>.csv files (default: auto-detect basis_results, then basis_results_old)",
    )
    parser.add_argument(
        "--ratio-method",
        choices=["simple", "delta_method"],
        default="delta_method",
        help="Ratio estimator for the basis method (default: delta_method)",
    )
    args = parser.parse_args()

    sim_path = os.path.abspath(args.sim_path)
    if not os.path.exists(sim_path):
        raise FileNotFoundError(f"Simulation file does not exist: {sim_path}")

    basis_dir = resolve_basis_dir(args.basis_dir)
    pre_gid, post_gid, pair_name, freq_hz, dt_ms = infer_protocol_from_path(sim_path)
    params = PARAMS_MAPPING[args.params]

    with open(sim_path, "rb") as handle:
        sim_data = pickle.load(handle)

    basis_result, rho_summary = compute_basis_result(
        sim_data,
        pre_gid,
        post_gid,
        freq_hz,
        dt_ms,
        params,
        basis_dir,
        args.ratio_method,
    )
    full_result = compute_full_result(sim_data)
    print_report(
        sim_path,
        pair_name,
        freq_hz,
        dt_ms,
        basis_dir,
        args.ratio_method,
        basis_result,
        full_result,
        rho_summary,
    )


if __name__ == "__main__":
    main()
