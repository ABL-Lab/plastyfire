#!/usr/bin/env python3
"""
Create dhuruva_modified_edges.h5 from the original SONATA edges.h5.

For excitatory synapses of a specific set of (pre_gid, post_gid) pairs,
replaces the original (random Bernoulli) rho0_GB with a conductance-based
median split, matching epg_dhuruva._get_ltpltd_params:

    For each (pre_gid, post_gid) pair:
        rho0 = 1 if conductance >= median(connection conductances)  else 0
        Use_d_TM  = u_syn^(1/k_u)  if rho0==1  else u_syn
        Use_p_TM  = u_syn           if rho0==1  else u_syn^k_u
        gmax_d_AMPA = conductance/k_gsyn  if rho0==1  else conductance
        gmax_p_AMPA = conductance         if rho0==1  else k_gsyn*conductance

Uses the SONATA target_to_source index for fast targeted lookups instead of
scanning all 407M edges.

Also patches each pair's simulation_config.json files with a metadata.synapses
block recording edge_id, rho0, conductance, gmax_d_AMPA, gmax_p_AMPA per synapse.

Usage:
    python create_dhuruva_edges.py --pairs-file data/pairs_n100.txt
    python create_dhuruva_edges.py --pairs-file data/pairs_n100.txt --k-u 0.2 --k-gsyn 2
    python create_dhuruva_edges.py --pairs-file data/pairs_n100.txt --dry-run
"""

import os
import json
import argparse
import numpy as np
import h5py

# ── Circuit defaults ──────────────────────────────────────────────────────────
CIRCUIT_BASE   = "/project/ctb-emuller/datasets/SSCx_plasticity/O1_2023a_Ecker"
EDGES_FILE     = os.path.join(CIRCUIT_BASE,
                    "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical/edges.h5")
CIRCUIT_CONFIG = os.path.join(CIRCUIT_BASE, "circuit_config.json")

EDGE_POP = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"

DEFAULT_OUTPUT = ("/project/ctb-emuller/dhuruva/plastyfire/data/"
                  "dhuruva_modified_edges.h5")
DEFAULT_CIRCUIT_OUT = ("/project/ctb-emuller/dhuruva/plastyfire/data/"
                       "dhuruva_circuit_config.json")
DEFAULT_PAIRS_FILE = ("/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/"
                      "data/pairs_n100.txt")

DEFAULT_SIMS_DIR = ("/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/"
                    "refitting_results/fitting/n100/seed19091997/"
                    "L5TTPC_L5TTPC_STDP/simulations")

# Fields to store locally (patched); all others become ExternalLinks
MODIFIED_FIELDS = {"rho0_GB", "Use_d_TM", "Use_p_TM", "gmax_d_AMPA", "gmax_p_AMPA"}

# Excitatory synapse type IDs (≥ 100 in SONATA convention)
EXC_TYPE_MIN = 100


# ── Helper: compute plasticity params from conductance + rho0 ─────────────────
def _plasticity_params(u: np.ndarray, cond: np.ndarray,
                       rho0: np.ndarray, k_u: float, k_gsyn: float):
    """Vectorised version of epg_dhuruva._get_ltpltd_params."""
    pot   = rho0.astype(bool)
    use_d  = np.where(pot, np.power(u, 1.0 / k_u), u)
    use_p  = np.where(pot, u,                        np.power(u, k_u))
    gmax_d = np.where(pot, cond / k_gsyn,            cond)
    gmax_p = np.where(pot, cond,                     k_gsyn * cond)
    return use_d, use_p, gmax_d, gmax_p


# ── Targeted per-pair lookup using SONATA index ───────────────────────────────
def find_and_compute_pairs(
    edges_file: str, pairs: list[tuple[int, int]], k_u: float, k_gsyn: float
) -> tuple[np.ndarray, dict[str, np.ndarray], dict]:
    """
    For each (pre_gid, post_gid) pair, uses the SONATA target_to_source index
    to look up only the relevant synapses, then computes rho0/gmax/Use via
    conductance median split (matching epg_dhuruva logic).

    Returns:
        edge_idx       : flat array of all modified edge indices
        modifications  : dict of full-pair concatenated arrays for writing
        per_pair_data  : dict keyed by (pre_gid, post_gid) with per-synapse info
                         for patching simulation_config.json files
    """
    pop      = f"edges/{EDGE_POP}"
    pop0     = f"{pop}/0"
    idx_base = f"{pop}/indices/target_to_source"

    all_edge_idx  = []
    all_rho0      = []
    all_use_d     = []
    all_use_p     = []
    all_gmax_d    = []
    all_gmax_p    = []
    per_pair_data = {}

    with h5py.File(edges_file, "r") as h:
        n2r   = h[f"{idx_base}/node_id_to_ranges"][:]    # (N_nodes, 2) uint64
        r2e   = h[f"{idx_base}/range_to_edge_id"][:]     # (N_ranges, 2) uint64
        src_ds   = h[f"{pop}/source_node_id"]
        stype_ds = h[f"{pop0}/syn_type_id"]
        cond_ds  = h[f"{pop0}/conductance"]
        usyn_ds  = h[f"{pop0}/u_syn"]

        for pre_gid, post_gid in pairs:
            r_start, r_end = n2r[post_gid]
            # Collect all edge indices incoming to post_gid
            edge_ranges = r2e[r_start:r_end]  # [[e_start, e_end), ...]

            pair_idx   = []
            pair_cond  = []
            pair_usyn  = []

            for e_start, e_end in edge_ranges:
                e_start, e_end = int(e_start), int(e_end)
                src   = src_ds[e_start:e_end]
                stype = stype_ds[e_start:e_end]
                mask  = (src == pre_gid) & (stype >= EXC_TYPE_MIN)
                where = np.where(mask)[0]
                if len(where):
                    abs_idx = where + e_start
                    pair_idx.append(abs_idx)
                    pair_cond.append(cond_ds[e_start:e_end][where].astype(np.float64))
                    pair_usyn.append(usyn_ds[e_start:e_end][where].astype(np.float64))

            if not pair_idx:
                print(f"  WARNING: no excitatory synapses found for pair {pre_gid}→{post_gid}")
                continue

            idx   = np.concatenate(pair_idx)
            cond  = np.concatenate(pair_cond)
            usyn  = np.concatenate(pair_usyn)
            n_syn = len(idx)

            # Conductance median split (same as epg_dhuruva)
            rho0 = (cond >= np.median(cond)).astype(np.int64)
            use_d, use_p, gmax_d, gmax_p = _plasticity_params(usyn, cond, rho0, k_u, k_gsyn)

            print(f"  pair {pre_gid}→{post_gid}: {n_syn} synapses, "
                  f"rho0=1: {int(rho0.sum())} ({100*rho0.mean():.0f}%)")

            all_edge_idx.append(idx)
            all_rho0.append(rho0)
            all_use_d.append(use_d.astype(np.float32))
            all_use_p.append(use_p.astype(np.float32))
            all_gmax_d.append(gmax_d.astype(np.float32))
            all_gmax_p.append(gmax_p.astype(np.float32))

            per_pair_data[(pre_gid, post_gid)] = {
                "edge_ids":    idx.tolist(),
                "rho0":        rho0.tolist(),
                "conductance": cond.tolist(),
                "gmax_d_AMPA": gmax_d.astype(np.float32).tolist(),
                "gmax_p_AMPA": gmax_p.astype(np.float32).tolist(),
                "Use_d_TM":    use_d.astype(np.float32).tolist(),
                "Use_p_TM":    use_p.astype(np.float32).tolist(),
            }

    edge_idx = np.concatenate(all_edge_idx)
    modifications = {
        "rho0_GB":     np.concatenate(all_rho0).astype(np.int64),
        "Use_d_TM":    np.concatenate(all_use_d),
        "Use_p_TM":    np.concatenate(all_use_p),
        "gmax_d_AMPA": np.concatenate(all_gmax_d),
        "gmax_p_AMPA": np.concatenate(all_gmax_p),
    }
    n_pairs = len(all_edge_idx)
    print(f"\n  Total: {n_pairs} pairs, {len(edge_idx):,} synapses modified")
    return edge_idx, modifications, per_pair_data


# ── Write the modified edges file ─────────────────────────────────────────────
def write_modified_edges(
    edges_file: str,
    output_path: str,
    edge_idx: np.ndarray,
    modifications: dict[str, np.ndarray],
    chunk_size: int = 2_000_000,
):
    """
    Creates output_path as an HDF5 file where:
      - All datasets NOT in MODIFIED_FIELDS → HDF5 ExternalLinks to original
      - MODIFIED_FIELDS → full local arrays (original values + pair-targeted patch)
    """
    abs_orig    = os.path.abspath(edges_file)
    pop_path    = f"edges/{EDGE_POP}"
    group0_path = f"{pop_path}/0"

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    print(f"  Reading original arrays and applying patch ...")
    with h5py.File(edges_file, "r") as src_h:
        n_tot = src_h[f"{pop_path}/source_node_id"].shape[0]

        patched = {}
        for field in MODIFIED_FIELDS:
            ds    = src_h[f"{group0_path}/{field}"]
            arr   = np.empty(n_tot, dtype=ds.dtype)
            for start in range(0, n_tot, chunk_size):
                end = min(start + chunk_size, n_tot)
                arr[start:end] = ds[start:end]
                print(f"    {field}: {end//chunk_size}/{(n_tot+chunk_size-1)//chunk_size}", flush=True)
            arr[edge_idx] = modifications[field].astype(arr.dtype)
            patched[field] = arr

        print(f"  Writing {output_path} ...")
        with h5py.File(output_path, "w") as out_h:
            pop_grp = out_h.require_group(pop_path)

            for key in src_h[pop_path].keys():
                if key == "0":
                    continue
                pop_grp[key] = h5py.ExternalLink(abs_orig, f"{pop_path}/{key}")

            grp0 = out_h.require_group(group0_path)
            for field in src_h[group0_path].keys():
                if field not in MODIFIED_FIELDS:
                    grp0[field] = h5py.ExternalLink(abs_orig, f"{group0_path}/{field}")

            hdf_chunks = (min(chunk_size, n_tot),)
            for field, arr in patched.items():
                print(f"    Writing {field} ...")
                grp0.create_dataset(field, data=arr, chunks=hdf_chunks,
                                    compression="gzip", compression_opts=4, shuffle=True)


# ── Patch simulation_config.json files with rho0 metadata ────────────────────
def patch_simulation_configs(per_pair_data: dict, sims_dir: str, dry_run: bool = False):
    """
    For each (pre_gid, post_gid), finds all simulation_config.json files under
    sims_dir/{pre_gid}-{post_gid}/*/simulation_config.json and adds a
    metadata.synapses block with edge_id, rho0, conductance, gmax_d/p, Use_d/p.
    """
    sims_base = os.path.expanduser(sims_dir)
    n_patched = 0
    n_missing = 0

    for (pre_gid, post_gid), syn_info in per_pair_data.items():
        pair_dir = os.path.join(sims_base, f"{pre_gid}-{post_gid}")
        if not os.path.isdir(pair_dir):
            print(f"  WARNING: sim dir not found: {pair_dir}")
            n_missing += 1
            continue

        metadata = {
            "synapses": {
                "edge_ids":    syn_info["edge_ids"],
                "rho0":        syn_info["rho0"],
                "conductance": syn_info["conductance"],
                "gmax_d_AMPA": syn_info["gmax_d_AMPA"],
                "gmax_p_AMPA": syn_info["gmax_p_AMPA"],
                "Use_d_TM":    syn_info["Use_d_TM"],
                "Use_p_TM":    syn_info["Use_p_TM"],
            }
        }

        # Patch all dt subdirs for this pair
        for dt_dir in sorted(os.listdir(pair_dir)):
            cfg_path = os.path.join(pair_dir, dt_dir, "simulation_config.json")
            if not os.path.isfile(cfg_path):
                continue
            sidecar = os.path.join(pair_dir, dt_dir, "synapse_metadata.json")
            if dry_run:
                print(f"  [DRY RUN] Would write {sidecar}")
                n_patched += 1
                continue
            with open(sidecar, "w") as f:
                json.dump(metadata, f, indent=2)
            n_patched += 1

    status = "[DRY RUN] " if dry_run else ""
    print(f"  {status}Wrote synapse_metadata.json to {n_patched} sim dirs "
          f"({n_missing} pair dirs missing)")



def write_circuit_config(original_config: str, new_edges_path: str, output: str):
    with open(original_config) as f:
        cfg = json.load(f)
    for edge_entry in cfg["networks"]["edges"]:
        if EDGE_POP in edge_entry.get("populations", {}):
            edge_entry["edges_file"] = new_edges_path
            break
    os.makedirs(os.path.dirname(output), exist_ok=True)
    with open(output, "w") as f:
        json.dump(cfg, f, indent=2)
    print(f"  Circuit config written: {output}")


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Create dhuruva_modified_edges.h5 with conductance-based rho0 for specific pairs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--pairs-file", default=DEFAULT_PAIRS_FILE,
                        help="Two-column text file of (pre_gid post_gid) pairs to modify")
    parser.add_argument("--k-u",    type=float, default=0.2)
    parser.add_argument("--k-gsyn", type=float, default=2.0)
    parser.add_argument("--output", default=DEFAULT_OUTPUT,
                        help=f"Output edges file (default: {DEFAULT_OUTPUT})")
    parser.add_argument("--circuit-out", default=DEFAULT_CIRCUIT_OUT)
    parser.add_argument("--no-circuit-config", action="store_true")
    parser.add_argument("--sims-dir", default=DEFAULT_SIMS_DIR,
                        help="Base simulations directory containing {pre_gid}-{post_gid} subdirs")
    parser.add_argument("--no-patch-configs", action="store_true",
                        help="Skip patching simulation_config.json files with rho0 metadata")
    parser.add_argument("--dry-run", action="store_true",
                        help="Compute modifications but do not write any files")
    args = parser.parse_args()

    pairs = []
    with open(args.pairs_file) as f:
        for line in f:
            line = line.strip()
            if line:
                pre, post = map(int, line.split())
                pairs.append((pre, post))

    print("=" * 65)
    print("create_dhuruva_edges.py  (pair-targeted mode)")
    print(f"  pairs      = {len(pairs)} pairs from {args.pairs_file}")
    print(f"  k_u        = {args.k_u}")
    print(f"  k_gsyn     = {args.k_gsyn}")
    print(f"  output     = {args.output}")
    print("=" * 65)

    print("\n[1/2] Looking up synapses for each pair via SONATA index ...")
    edge_idx, modifications, per_pair_data = find_and_compute_pairs(EDGES_FILE, pairs, args.k_u, args.k_gsyn)

    print("\n[1b] Patching simulation_config.json files with rho0 metadata ...")
    if not args.no_patch_configs:
        patch_simulation_configs(per_pair_data, args.sims_dir, dry_run=args.dry_run)

    if args.dry_run:
        print("\n[DRY RUN] Skipping file writes.")
        return

    print(f"\n[2/2] Writing modified edges file ...")
    write_modified_edges(EDGES_FILE, args.output, edge_idx, modifications)
    print(f"  Done: {args.output}")

    if not args.no_circuit_config:
        write_circuit_config(CIRCUIT_CONFIG, args.output, args.circuit_out)

    print("\n" + "=" * 65)
    print(f"Modified {len(modifications['rho0_GB']):,} synapses across {len(pairs)} pairs.")
    print("=" * 65)


if __name__ == "__main__":
    main()
