#!/usr/bin/env python3
"""
Validation script: re-run Chindemi params through the FIXED pairrunner_edges_fit.py
for pair 181002-188173 at both ±10ms and compare resulting rho against BCL.

Expected after fix:
  - theta_d/theta_p computed by optimizer should match edges.h5 values
  - final_rho from optimizer should match BCL bluecellulab_results/rho.h5
"""

import os
import sys
import pickle
import subprocess
import hashlib
import h5py
import numpy as np

PLASTYFIRE_ROOT = os.path.dirname(os.path.abspath(__file__))
EDGES_H5 = "/project/ctb-emuller/dhuruva/plastyfire/data/dhuruva_modified_edges.h5"
EDGE_POP  = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"
NODE_POP  = "S1nonbarrel_neurons"

CHINDEMI_PARAMS = [101.5, 216.2, 1.002, 1.954, 1.159, 2.483, 1.127, 2.456, 5.236, 1.782]
CHINDEMI_HASH   = hashlib.md5(str(CHINDEMI_PARAMS).encode()).hexdigest()[:12]
PARAM_NAMES = [
    "gamma_d_GB_GluSynapse", "gamma_p_GB_GluSynapse",
    "a00", "a01", "a10", "a11", "a20", "a21", "a30", "a31",
]

BASE = os.path.join(PLASTYFIRE_ROOT, "refitting_results/fitting/n100/seed19091997/"
                    "L5TTPC_L5TTPC_STDP/simulations")

PAIRS_DTS = [
    ("181002-188173", "10Hz_-10ms"),
    ("181002-188173", "10Hz_10ms"),
    ("187212-192870", "10Hz_-10ms"),   # pair with big BCL LTD vs LTP difference
    ("187212-192870", "10Hz_10ms"),
]

def run_pairrunner(sim_dir, output_pkl):
    """Run pairrunner_edges_fit.py with Chindemi params, writing to output_pkl."""
    pairrunner = os.path.join(PLASTYFIRE_ROOT, "plastyfire", "pairrunner_edges_fit.py")
    cmd = [sys.executable, pairrunner]
    for name, val in zip(PARAM_NAMES, CHINDEMI_PARAMS):
        cmd += [f"--{name}={val}"]
    cmd += [f"--param_hash=VALIDATE_{CHINDEMI_HASH}"]

    # pairrunner writes to cwd, so run from sim_dir
    env = os.environ.copy()
    env["PYTHONPATH"] = PLASTYFIRE_ROOT + ":" + env.get("PYTHONPATH", "")
    print(f"\n  Running pairrunner in {sim_dir} ...")
    result = subprocess.run(cmd, cwd=sim_dir, env=env,
                            capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  STDERR:\n{result.stderr[-3000:]}")
        raise RuntimeError(f"pairrunner failed (exit {result.returncode})")
    print(f"  stdout: {result.stdout[-500:]}")
    return os.path.join(sim_dir, f"simulation_edges_VALIDATE_{CHINDEMI_HASH}.pkl")


def read_bcl_rho(sim_dir):
    h5 = os.path.join(sim_dir, "bluecellulab_results", "rho.h5")
    with h5py.File(h5, "r") as f:
        data = f[f"report/{NODE_POP}/data"][()]
    initial = np.array([1 if v >= 0.5 else 0 for v in data[0]], dtype=float)
    final   = data[-1]
    final_bin = np.array([1 if v >= 0.5 else 0 for v in final], dtype=float)
    return initial, final, final_bin


def read_edges_theta(sim_dir):
    """Read theta_d/theta_p from edges.h5 for the global_ids in this sim."""
    # Get global_ids from old optimizer pkl
    old_pkl = os.path.join(sim_dir, f"simulation_edges_{CHINDEMI_HASH}.pkl")
    with open(old_pkl, "rb") as f:
        d = pickle.load(f)
    gids = np.array(d["global_ids"], dtype=np.int64)
    with h5py.File(EDGES_H5, "r") as f:
        td = f[f"edges/{EDGE_POP}/0/theta_d"][gids]
        tp = f[f"edges/{EDGE_POP}/0/theta_p"][gids]
    return gids, td, tp


def main():
    print(f"Chindemi hash: {CHINDEMI_HASH}")
    print(f"Validation hash: VALIDATE_{CHINDEMI_HASH}")
    print("="*70)

    all_ok = True
    for pair, dt_dir in PAIRS_DTS:
        sim_dir = os.path.join(BASE, pair, dt_dir)
        print(f"\n{'='*70}")
        print(f"Pair: {pair}  dt: {dt_dir}")

        # 1. Read BCL reference
        bcl_init, bcl_final_raw, bcl_final_bin = read_bcl_rho(sim_dir)
        print(f"  BCL initial_rho: {bcl_init}")
        print(f"  BCL final_raw:   {np.round(bcl_final_raw, 3)}")
        print(f"  BCL final_bin:   {bcl_final_bin}  mean={bcl_final_bin.mean():.3f}")

        # 2. Read edges.h5 thresholds (ground truth)
        gids, td_h5, tp_h5 = read_edges_theta(sim_dir)
        print(f"  edges.h5 theta_d: {np.round(td_h5, 4)}")
        print(f"  edges.h5 theta_p: {np.round(tp_h5, 4)}")

        # 3. Run pairrunner with fixed code
        try:
            out_pkl = run_pairrunner(sim_dir, f"simulation_edges_VALIDATE_{CHINDEMI_HASH}.pkl")
        except RuntimeError as e:
            print(f"  ERROR: {e}")
            all_ok = False
            continue

        if not os.path.exists(out_pkl):
            print(f"  ERROR: output pkl not found: {out_pkl}")
            all_ok = False
            continue

        # 4. Load new optimizer result
        with open(out_pkl, "rb") as f:
            new = pickle.load(f)
        new_init  = np.array(new["initial_rho"])
        new_final = np.array(new["final_rho"])
        new_final_bin = np.array([1 if v >= 0.5 else 0 for v in new_final], dtype=float)

        print(f"  NEW initial_rho: {np.round(new_init, 3)}")
        print(f"  NEW final_rho:   {np.round(new_final, 3)}")
        print(f"  NEW final_bin:   {new_final_bin}  mean={new_final_bin.mean():.3f}")

        # 5. Compare
        bin_match = np.array_equal(bcl_final_bin, new_final_bin)
        init_match = np.allclose(bcl_init, new_init, atol=0.01)
        print(f"\n  initial_rho match (BCL vs new): {init_match}")
        print(f"  final_bin   match (BCL vs new): {bin_match}  {'✓ PASS' if bin_match else '✗ FAIL'}")

        # Compare raw final rho (continuous vs continuous, both should be close now)
        raw_diff = np.abs(bcl_final_raw - new_final)
        print(f"  final_rho raw diff: max={raw_diff.max():.4f}  mean={raw_diff.mean():.4f}")

        if not bin_match:
            all_ok = False
            print(f"  MISMATCH details:")
            for i, (b, n, br, nr) in enumerate(zip(bcl_final_bin, new_final_bin,
                                                    bcl_final_raw, new_final)):
                match = "✓" if b == n else "✗"
                print(f"    syn{i}: BCL_bin={b:.0f} (raw={br:.3f})  NEW_bin={n:.0f} (raw={nr:.3f})  {match}")

    print(f"\n{'='*70}")
    print(f"Overall: {'ALL PASS ✓' if all_ok else 'SOME FAILED ✗'}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
