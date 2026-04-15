#!/usr/bin/env python3
"""Compare basis-based and full-trace EPSP ratio estimates for one simulation.

This script reproduces the same two paths used by:
  - plot_stdp_curves.py        -> basis extrapolation from rho_GB states
  - plot_stdp_curves_full.py   -> direct EPSP measurement from t/v/prespikes

Example:
  python validate_basis_vs_full_epsp.py \
      --sim-path /project/ctb-emuller/dhuruva/plastyfire/full_trace_results/CHINDEMI_PARAMS/180164-197248/10Hz_5ms/simulation_traces.pkl \
      --params chindemi
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


BASIS_DIR = os.path.join(REPO_DIR, "basis_results_old")
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


def infer_pair_and_dt(sim_path):
    freq_dt = os.path.basename(os.path.dirname(sim_path))
    pair_name = os.path.basename(os.path.dirname(os.path.dirname(sim_path)))
    try:
        pre_gid_str, post_gid_str = pair_name.split("-", 1)
        pre_gid = int(pre_gid_str)
        post_gid = int(post_gid_str)
    except ValueError as exc:
        raise ValueError(
            f"Could not parse pair ids from parent directory '{pair_name}'."
        ) from exc

    if not freq_dt.startswith("10Hz_") or not freq_dt.endswith("ms"):
        raise ValueError(
            f"Could not parse Δt from directory '{freq_dt}'. Expected format like 10Hz_5ms."
        )
    dt_ms = float(freq_dt.replace("10Hz_", "").replace("ms", ""))
    return pre_gid, post_gid, pair_name, dt_ms


def normalize_rho_trace(rho_gb_raw):
    if isinstance(rho_gb_raw, dict):
        rho_trace = np.array(list(rho_gb_raw.values()))
    else:
        rho_trace = np.asarray(rho_gb_raw)
        if rho_trace.ndim == 2 and rho_trace.shape[0] > rho_trace.shape[1]:
            rho_trace = rho_trace.T
        elif rho_trace.ndim == 1:
            rho_trace = rho_trace.reshape(1, -1)
    return rho_trace


def compute_rho_state_summary(rho_trace):
    initial = np.array([1 if rho[0] >= 0.5 else 0 for rho in rho_trace], dtype=int)
    final = np.array([1 if rho[-1] >= 0.5 else 0 for rho in rho_trace], dtype=int)
    return {
        "initial_rho": initial.tolist(),
        "final_rho": final.tolist(),
        "n_0to1": int(np.sum((initial == 0) & (final == 1))),
        "n_0to0": int(np.sum((initial == 0) & (final == 0))),
        "n_1to1": int(np.sum((initial == 1) & (final == 1))),
        "n_1to0": int(np.sum((initial == 1) & (final == 0))),
        "n_synapses": int(len(initial)),
    }


def compute_basis_result(sim_data, pre_gid, post_gid, dt_ms, params, basis_dir, ratio_method):
    if "rho_GB" not in sim_data:
        raise KeyError("Simulation data is missing 'rho_GB', required for the basis method.")

    rho_trace = normalize_rho_trace(sim_data["rho_GB"])
    rho_summary = compute_rho_state_summary(rho_trace)

    basis_result = fetch_epsp_ratio(
        pre_gid,
        post_gid,
        10.0,
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
    required_keys = ["t", "v", "prespikes"]
    missing = [key for key in required_keys if key not in sim_data]
    if missing:
        raise KeyError(
            "Simulation data is missing keys required for the full method: "
            + ", ".join(missing)
        )

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
        "n_prespikes": len(sim_data["prespikes"]),
    }


def print_report(sim_path, pair_name, dt_ms, basis_result, full_result, rho_summary, basis_dir, ratio_method):
    basis_ratio = basis_result["ratio_mean"]
    full_ratio = full_result["ratio_mean"]
    abs_diff = abs(full_ratio - basis_ratio)
    rel_diff_pct = np.nan
    if basis_ratio != 0:
        rel_diff_pct = 100.0 * (full_ratio - basis_ratio) / basis_ratio

    print("=" * 72)
    print("EPSP Ratio Validation: basis vs full")
    print("=" * 72)
    print(f"Simulation file : {sim_path}")
    print(f"Pair            : {pair_name}")
    print(f"Δt              : {dt_ms:.3f} ms")
    print(f"Basis dir       : {basis_dir}")
    print(f"Basis method    : fetch_epsp_ratio(..., ratio_method='{ratio_method}')")
    print(f"Full method     : Experiment.compute_epsp_ratio(n={N_EPSP}, method='amplitude')")

    print("\nRho state summary")
    print("-" * 72)
    print(f"Synapses        : {rho_summary['n_synapses']}")
    print(f"0 -> 1          : {rho_summary['n_0to1']}")
    print(f"0 -> 0          : {rho_summary['n_0to0']}")
    print(f"1 -> 1          : {rho_summary['n_1to1']}")
    print(f"1 -> 0          : {rho_summary['n_1to0']}")

    print("\nBasis extrapolation result")
    print("-" * 72)
    print(f"EPSP before     : {basis_result['epsp_before_mean']:.8f} ± {basis_result['epsp_before_std']:.8f}")
    print(f"EPSP after      : {basis_result['epsp_after_mean']:.8f} ± {basis_result['epsp_after_std']:.8f}")
    print(f"EPSP ratio      : {basis_ratio:.8f}")

    print("\nFull trace result")
    print("-" * 72)
    print(f"EPSP before     : {full_result['epsp_before_mean']:.8f} ± {full_result['epsp_before_std']:.8f}")
    print(f"EPSP after      : {full_result['epsp_after_mean']:.8f} ± {full_result['epsp_after_std']:.8f}")
    print(f"EPSP ratio      : {full_ratio:.8f}")
    print(f"Pre-spikes      : {full_result['n_prespikes']}")

    print("\nDifference")
    print("-" * 72)
    print(f"Absolute diff   : {abs_diff:.8f}")
    if np.isnan(rel_diff_pct):
        print("Relative diff   : undefined (basis ratio is zero)")
    else:
        print(f"Relative diff   : {rel_diff_pct:.3f}%  (full - basis) / basis")

    print("\nInterpretation")
    print("-" * 72)
    print("The basis method uses only the initial and final binary rho states plus the pair's basis file.")
    print("The full method measures actual EPSPs from the voltage trace over the last 60 pulses of C01 and C02.")
    print("They can differ because multiple voltage-trace effects are collapsed in the basis method:")
    print("  1. within-state amplitude variability across synapses,")
    print("  2. finite-trace noise and baseline drift,")
    print("  3. non-binary behavior during the run that is reduced to initial/final rho >= 0.5,")
    print("  4. bias correction in ratio_method='delta_method' rather than a direct after/before ratio.")
    print("=" * 72)


def main():
    parser = argparse.ArgumentParser(
        description="Compare basis-based and full-trace EPSP ratios for one simulation pickle."
    )
    parser.add_argument("--sim-path", required=True, help="Path to simulation_traces.pkl")
    parser.add_argument(
        "--params",
        choices=sorted(PARAMS_MAPPING.keys()),
        default="chindemi",
        help="Parameter set passed to fetch_epsp_ratio (default: chindemi)",
    )
    parser.add_argument(
        "--basis-dir",
        default=BASIS_DIR,
        help="Directory containing basis_<pre>_<post>.csv files",
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

    pre_gid, post_gid, pair_name, dt_ms = infer_pair_and_dt(sim_path)
    params = PARAMS_MAPPING[args.params]

    with open(sim_path, "rb") as handle:
        sim_data = pickle.load(handle)

    basis_result, rho_summary = compute_basis_result(
        sim_data,
        pre_gid,
        post_gid,
        dt_ms,
        params,
        args.basis_dir,
        args.ratio_method,
    )
    full_result = compute_full_result(sim_data)
    print_report(
        sim_path,
        pair_name,
        dt_ms,
        basis_result,
        full_result,
        rho_summary,
        args.basis_dir,
        args.ratio_method,
    )


if __name__ == "__main__":
    main()