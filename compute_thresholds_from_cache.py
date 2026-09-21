#!/usr/bin/env python3
"""
Compute theta_d/theta_p from an existing cpre/cpost cache and inject into edges.h5.

Avoids re-running any BCL simulations — uses cached c_pre/c_post scalars
(per synapse global index) and computes thresholds using given a-params.

Usage:
    python compute_thresholds_from_cache.py \\
        --cache cpre_cpost_cache/ion_channels_tau278.pkl \\
        --output-dir threshold_results_ic_gen4/ \\
        --a00 1.002 --a01 2.254638 \\
        --a10 1.209857 --a11 2.396290 \\
        --a20 1.127 --a21 2.500181 \\
        --a30 4.214614 --a31 2.239758 \\
        --edges-h5 data/dhuruva_modified_edges.h5

    # Then inject:
    python inject_thresholds.py --results-dir threshold_results_ic_gen4/ \\
        --edges-h5 data/dhuruva_modified_edges.h5
"""

import argparse
import logging
import os
import pickle
from pathlib import Path

import h5py
import numpy as np

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

EDGES_H5 = "/project/ctb-emuller/dhuruva/plastyfire/data/dhuruva_modified_edges.h5"
EDGE_POP = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"
POP_PATH = f"edges/{EDGE_POP}/0"

BASAL_SECTION_TYPE  = 2
APICAL_SECTION_TYPE = 3

# Gen-4 ion-channel optimizer best (defaults)
DEFAULT_A = {
    "a00": 1.002,     "a01": 2.254638,  # basal depression
    "a10": 1.209857,  "a11": 2.396290,  # basal potentiation
    "a20": 1.127,     "a21": 2.500181,  # apical depression
    "a30": 4.214614,  "a31": 2.239758,  # apical potentiation
}


def load_section_types(edges_h5, global_indices):
    """Load afferent_section_type for the given global synapse indices."""
    logger.info(f"Loading section types for {len(global_indices):,} synapses …")
    with h5py.File(edges_h5, "r") as f:
        sec_arr = f[f"{POP_PATH}/afferent_section_type"][global_indices]
    return sec_arr


def compute_thresholds_for_pair(pair_key, pair_data, edges_h5, a, output_dir):
    """Compute and save theta_d/theta_p .npz for one pair."""
    pre_gid, post_gid = pair_key
    out_file = output_dir / f"{pre_gid}_{post_gid}.npz"
    if out_file.exists():
        return f"skip:{pre_gid}_{post_gid}"

    c_pre_dict  = pair_data["c_pre"]   # {global_idx: float}
    c_post_dict = pair_data["c_post"]  # {global_idx: float}

    # Intersect keys
    global_ids = sorted(set(c_pre_dict.keys()) & set(c_post_dict.keys()))
    if not global_ids:
        return f"no_synapses:{pre_gid}_{post_gid}"

    global_arr = np.array(global_ids, dtype=np.int64)
    c_pre_arr  = np.array([c_pre_dict[g]  for g in global_ids], dtype=np.float32)
    c_post_arr = np.array([c_post_dict[g] for g in global_ids], dtype=np.float32)

    # Load section types in one read
    sec_arr = load_section_types(edges_h5, global_arr)

    basal  = sec_arr == BASAL_SECTION_TYPE
    apical = sec_arr == APICAL_SECTION_TYPE
    other  = ~basal & ~apical

    td = np.full(len(global_arr), -1.0, dtype=np.float32)
    tp = np.full(len(global_arr), -1.0, dtype=np.float32)

    td[basal]  = a["a00"] * c_pre_arr[basal]  + a["a01"] * c_post_arr[basal]
    tp[basal]  = a["a10"] * c_pre_arr[basal]  + a["a11"] * c_post_arr[basal]
    td[apical] = a["a20"] * c_pre_arr[apical] + a["a21"] * c_post_arr[apical]
    tp[apical] = a["a30"] * c_pre_arr[apical] + a["a31"] * c_post_arr[apical]

    np.savez_compressed(
        out_file,
        global_idx=global_arr,
        theta_d=td,
        theta_p=tp,
    )
    return (
        f"ok:{pre_gid}_{post_gid} n={len(global_arr)} "
        f"basal={basal.sum()} apical={apical.sum()} other={other.sum()} "
        f"td_mean={td[td>0].mean():.4f} tp_mean={tp[tp>0].mean():.4f}"
    )


def main():
    parser = argparse.ArgumentParser(description="Compute theta_d/theta_p from cpre/cpost cache")
    parser.add_argument("--cache",      required=True, help="Path to cpre/cpost cache pkl")
    parser.add_argument("--output-dir", required=True, help="Dir to write per-pair .npz files")
    parser.add_argument("--edges-h5",   default=EDGES_H5)
    parser.add_argument("--a00", type=float, default=DEFAULT_A["a00"])
    parser.add_argument("--a01", type=float, default=DEFAULT_A["a01"])
    parser.add_argument("--a10", type=float, default=DEFAULT_A["a10"])
    parser.add_argument("--a11", type=float, default=DEFAULT_A["a11"])
    parser.add_argument("--a20", type=float, default=DEFAULT_A["a20"])
    parser.add_argument("--a21", type=float, default=DEFAULT_A["a21"])
    parser.add_argument("--a30", type=float, default=DEFAULT_A["a30"])
    parser.add_argument("--a31", type=float, default=DEFAULT_A["a31"])
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    a = {k: getattr(args, k) for k in DEFAULT_A}
    logger.info("A-params: " + ", ".join(f"{k}={v}" for k, v in a.items()))

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Loading cache: {args.cache}")
    with open(args.cache, "rb") as f:
        cache = pickle.load(f)
    logger.info(f"  {len(cache)} pairs in cache")

    if args.dry_run:
        logger.info("Dry run — no files written.")
        for key in list(cache.keys())[:3]:
            pre, post = key
            n = len(cache[key]["c_pre"])
            logger.info(f"  {pre}_{post}: {n} synapses")
        return

    ok, skipped, failed = 0, 0, 0
    for pair_key, pair_data in cache.items():
        result = compute_thresholds_for_pair(pair_key, pair_data, args.edges_h5, a, output_dir)
        logger.info(result)
        if result.startswith("ok:"):
            ok += 1
        elif result.startswith("skip:"):
            skipped += 1
        else:
            failed += 1

    logger.info(f"\nDone: {ok} computed, {skipped} skipped, {failed} failed")
    logger.info(f"Output: {output_dir}/")
    logger.info(f"Next: python inject_thresholds.py --results-dir {output_dir} --edges-h5 {args.edges_h5}")


if __name__ == "__main__":
    main()
