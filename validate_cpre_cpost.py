#!/usr/bin/env python3
"""
Validate that the fixed _find_cpre_cpost (in simulator_edges.py) produces
c_pre/c_post values consistent with compute_thresholds.py.

Method:
  - compute_thresholds.py stored its results as theta_d/theta_p in edges.h5
  - We run _find_cpre_cpost via simulator_edges.compute_cpre_cpost_for_workdir
    on a few pairs and compare the resulting theta_d/theta_p against edges.h5

Validation:
  - theta_d(fixed) = a*c_pre + b*c_post  (per synapse, basal vs apical)
  - theta_p(fixed) same
  - Compare both against edges.h5 values; flag >20% deviation

Expected: all synapse theta values within 20% of edges.h5 (within ~10%)
"""

import os, sys, h5py, numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

EDGES_H5  = "/project/ctb-emuller/dhuruva/plastyfire/data/dhuruva_modified_edges.h5"
EDGE_POP  = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"
BASE      = "refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations"

# Chindemi a-params (same as edges.h5 was built with)
FIT_PARAMS = {
    "tau_effca_GB_GluSynapse": 278.3177658387,
    "gamma_d_GB_GluSynapse":   101.5,
    "gamma_p_GB_GluSynapse":   216.2,
    "a00": 1.002, "a01": 1.954,
    "a10": 1.159, "a11": 2.483,
    "a20": 1.127, "a21": 2.456,
    "a30": 5.236, "a31": 1.782,
}

# Test a few pairs (one protocol dir per pair — c_pre/c_post are protocol-independent)
TEST_PAIRS = [
    ("181002-188173", "10Hz_-10ms"),
    ("187212-192870", "10Hz_-10ms"),
    ("196734-181209", "10Hz_-10ms"),
    ("180164-197248", "10Hz_10ms"),   # has apical synapses — key test for c_post fix
]

import plastyfire.simulator_edges as sim_mod

print("Validating c_pre/c_post: fixed _find_cpre_cpost vs edges.h5")
print("Computes theta_d/theta_p from new c_pre/c_post and checks against stored theta.")
print("="*90)

a00, a01 = FIT_PARAMS["a00"], FIT_PARAMS["a01"]
a10, a11 = FIT_PARAMS["a10"], FIT_PARAMS["a11"]
a20, a21 = FIT_PARAMS["a20"], FIT_PARAMS["a21"]
a30, a31 = FIT_PARAMS["a30"], FIT_PARAMS["a31"]

all_ok = True
for pair, dt in TEST_PAIRS:
    workdir = os.path.join(BASE, pair, dt)
    pre_gid, post_gid = [int(x) for x in pair.split("-")]

    print(f"\nPair {pair}")

    # --- Run the fixed mini-sims ---
    result = sim_mod.compute_cpre_cpost_for_workdir(workdir, fit_params=FIT_PARAMS)
    if result is None:
        print("  ERROR: compute_cpre_cpost_for_workdir returned None")
        all_ok = False
        continue

    c_pre_map  = result["c_pre"]
    c_post_map = result["c_post"]
    gids = sorted(c_pre_map.keys())

    print(f"  c_pre  mean={np.mean(list(c_pre_map.values())):.4e}  "
          f"min={min(c_pre_map.values()):.4e}  max={max(c_pre_map.values()):.4e}")
    print(f"  c_post mean={np.mean(list(c_post_map.values())):.4e}  "
          f"min={min(c_post_map.values()):.4e}  max={max(c_post_map.values()):.4e}")

    # --- Read edges.h5 theta_d/theta_p and sec_type ---
    with h5py.File(EDGES_H5, "r") as f:
        td_h5    = {int(g): float(f[f"edges/{EDGE_POP}/0/theta_d"][g])                    for g in gids}
        tp_h5    = {int(g): float(f[f"edges/{EDGE_POP}/0/theta_p"][g])                    for g in gids}
        sec_types = {int(g): int(f[f"edges/{EDGE_POP}/0/afferent_section_type"][g])       for g in gids}

    hdr = (f"  {'gid':>10s}  {'sec':>4s}  "
           f"{'c_pre':>10s}  {'c_post':>10s}  "
           f"{'td_new':>8s}  {'td_h5':>8s}  {'td_ratio':>8s}  "
           f"{'tp_new':>8s}  {'tp_h5':>8s}  {'tp_ratio':>8s}  {'ok':>4s}")
    print(hdr)

    pair_ok = True
    for g in gids:
        cp  = c_pre_map[g]
        cq  = c_post_map[g]
        is_apical = (sec_types[g] == 3)
        sec_label = "api" if is_apical else "bas"

        if is_apical:
            td_new = a20 * cp + a21 * cq
            tp_new = a30 * cp + a31 * cq
        else:
            td_new = a00 * cp + a01 * cq
            tp_new = a10 * cp + a11 * cq

        td_ref = td_h5[g]
        tp_ref = tp_h5[g]

        td_ratio = td_new / td_ref if td_ref > 0 else float("nan")
        tp_ratio = tp_new / tp_ref if tp_ref > 0 else float("nan")

        ok = (0.8 < td_ratio < 1.2) and (0.8 < tp_ratio < 1.2)
        if not ok:
            pair_ok = False
            all_ok  = False

        print(f"  {g:>10d}  {sec_label:>4s}  "
              f"{cp:>10.4e}  {cq:>10.4e}  "
              f"{td_new:>8.4f}  {td_ref:>8.4f}  {td_ratio:>8.3f}  "
              f"{tp_new:>8.4f}  {tp_ref:>8.4f}  {tp_ratio:>8.3f}  "
              f"{'✓' if ok else '✗'}")

    print(f"  → {'PASS ✓' if pair_ok else 'FAIL ✗'}")

print(f"\n{'='*90}")
print(f"Overall: {'ALL PASS ✓' if all_ok else 'SOME FAILED ✗'}")
sys.exit(0 if all_ok else 1)
