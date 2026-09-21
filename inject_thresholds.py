#!/usr/bin/env python3
"""
Inject pre-computed theta_d / theta_p thresholds into dhuruva_modified_edges.h5.

Reads per-post_gid .npz files produced by compute_thresholds.py and writes the
theta_d and theta_p arrays into the HDF5 file in place (using chunked writes).

Usage:
    python inject_thresholds.py --results-dir /path/to/compute/results

    # Dry run (reports statistics without modifying the file)
    python inject_thresholds.py --results-dir /path/to/compute/results --dry-run

    # Point at a different edges file
    python inject_thresholds.py --results-dir /path/to/compute/results \\
        --edges-h5 /path/to/edges.h5
"""

import argparse
import logging
import os
from pathlib import Path

import h5py
import numpy as np

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

EDGES_H5 = "/project/ctb-emuller/dhuruva/plastyfire/data/dhuruva_modified_edges.h5"
EDGE_POP = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"
POP_PATH = f"edges/{EDGE_POP}/0"


def collect_results(results_dir):
    """Read all .npz files and return (global_idx, theta_d, theta_p) arrays."""
    npz_files = sorted(Path(results_dir).glob("*.npz"))
    logger.info(f"Found {len(npz_files):,} .npz files in {results_dir}")

    all_idx  = []
    all_td   = []
    all_tp   = []
    skipped  = 0

    for path in npz_files:
        try:
            d = np.load(path)
            all_idx.append(d["global_idx"])
            all_td.append(d["theta_d"])
            all_tp.append(d["theta_p"])
        except Exception as e:
            logger.warning(f"  Could not read {path.name}: {e}")
            skipped += 1

    if not all_idx:
        raise RuntimeError("No valid .npz files found — nothing to inject.")

    global_idx = np.concatenate(all_idx).astype(np.int64)
    theta_d    = np.concatenate(all_td).astype(np.float32)
    theta_p    = np.concatenate(all_tp).astype(np.float32)

    # Deduplicate (last write wins) in case of overlapping batches
    _, uniq = np.unique(global_idx, return_index=True)
    global_idx, theta_d, theta_p = global_idx[uniq], theta_d[uniq], theta_p[uniq]

    logger.info(f"  Collected {len(global_idx):,} synapse threshold entries "
                f"({skipped} files skipped)")
    logger.info(f"  theta_d  range: [{theta_d.min():.4f}, {theta_d.max():.4f}]  "
                f"n_neg1={(theta_d == -1.0).sum():,}")
    logger.info(f"  theta_p  range: [{theta_p.min():.4f}, {theta_p.max():.4f}]  "
                f"n_neg1={(theta_p == -1.0).sum():,}")
    return global_idx, theta_d, theta_p


LOCAL_SIZE_FIELD = "rho0_GB"  # local dataset we can use to get n_total in r+ mode


def _get_or_create_local_dataset(f, field_name, n_total):
    """
    Return a writable HDF5 dataset for field_name under POP_PATH.
    If the current entry is an ExternalLink (or inaccessible), delete it and
    create a local chunked dataset with fillvalue=-1.0 so only written chunks
    are actually stored on disk.
    """
    path = f"{POP_PATH}/{field_name}"
    grp  = f[f"edges/{EDGE_POP}/0"]
    lnk  = grp.get(field_name, getlink=True)

    chunk_size = 65_536

    if lnk is None:
        logger.info(f"  {field_name}: creating new local chunked dataset "
                    f"(fillvalue=-1, n={n_total:,})")
        ds = f.create_dataset(
            path,
            shape=(n_total,),
            dtype="float32",
            chunks=(chunk_size,),
            fillvalue=-1.0,
        )
    elif isinstance(lnk, h5py.ExternalLink):
        logger.info(f"  {field_name}: replacing ExternalLink → local chunked dataset "
                    f"(fillvalue=-1, n={n_total:,})")
        del f[path]
        ds = f.create_dataset(
            path,
            shape=(n_total,),
            dtype="float32",
            chunks=(chunk_size,),
            fillvalue=-1.0,
        )
    else:
        ds = f[path]
        logger.info(f"  {field_name}: existing local dataset, shape={ds.shape}")
    return ds


def compare_existing(edges_h5, global_idx, theta_d, theta_p, tol=1e-4):
    """
    Read current theta_d/theta_p from edges.h5 for the given indices and compare
    to the new values.  Logs a WARNING (and returns False) if any values differ
    by more than `tol`, so the caller can decide whether to abort.
    Returns True if everything matches within tolerance.
    """
    logger.info("Comparing new thresholds against values currently in edges.h5 …")
    with h5py.File(edges_h5, "r") as f:
        pop = f[POP_PATH]
        if "theta_d" not in pop or "theta_p" not in pop:
            logger.warning("  theta_d/theta_p not yet present in edges.h5 — skipping comparison.")
            return True

        lnk_d = pop.get("theta_d", getlink=True)
        lnk_p = pop.get("theta_p", getlink=True)
        if isinstance(lnk_d, h5py.ExternalLink) or isinstance(lnk_p, h5py.ExternalLink):
            logger.warning("  theta_d/theta_p are ExternalLinks — cannot read current values for comparison.")
            return True

        cur_td = f[f"{POP_PATH}/theta_d"][global_idx]
        cur_tp = f[f"{POP_PATH}/theta_p"][global_idx]

    delta_d = np.abs(cur_td - theta_d)
    delta_p = np.abs(cur_tp - theta_p)

    n_mismatch_d = int((delta_d > tol).sum())
    n_mismatch_p = int((delta_p > tol).sum())

    if n_mismatch_d == 0 and n_mismatch_p == 0:
        logger.info(f"  All {len(global_idx):,} theta_d and theta_p values match within tol={tol:.0e}. Safe to overwrite.")
        return True

    logger.warning(
        f"  THRESHOLD MISMATCH DETECTED (tol={tol:.0e}):\n"
        f"    theta_d: {n_mismatch_d:,}/{len(global_idx):,} entries differ "
        f"(max delta={delta_d.max():.6f}, mean delta={delta_d.mean():.6f})\n"
        f"    theta_p: {n_mismatch_p:,}/{len(global_idx):,} entries differ "
        f"(max delta={delta_p.max():.6f}, mean delta={delta_p.mean():.6f})"
    )

    # Show the worst offenders to help diagnose which params changed
    worst_d_idx = int(np.argmax(delta_d))
    worst_p_idx = int(np.argmax(delta_p))
    logger.warning(
        f"  Worst theta_d mismatch: global_idx={global_idx[worst_d_idx]} "
        f"old={cur_td[worst_d_idx]:.6f} → new={theta_d[worst_d_idx]:.6f}"
    )
    logger.warning(
        f"  Worst theta_p mismatch: global_idx={global_idx[worst_p_idx]} "
        f"old={cur_tp[worst_p_idx]:.6f} → new={theta_p[worst_p_idx]:.6f}"
    )
    return False


def inject(edges_h5, global_idx, theta_d, theta_p, dry_run=False, force=False):
    """Write theta_d / theta_p values into edges.h5, handling ExternalLinks."""
    if dry_run:
        logger.info("DRY RUN — no file changes.")
        return

    sort_order = np.argsort(global_idx)
    global_idx = global_idx[sort_order]
    theta_d    = theta_d[sort_order]
    theta_p    = theta_p[sort_order]

    match = compare_existing(edges_h5, global_idx, theta_d, theta_p)
    if not match and not force:
        logger.error(
            "Aborting inject — thresholds differ from what is currently in edges.h5.\n"
            "If this is intentional (e.g. you updated a-params), re-run with --force to overwrite."
        )
        raise SystemExit(1)
    elif not match and force:
        logger.warning("--force specified — overwriting mismatched thresholds.")

    logger.info(f"Opening {edges_h5} for in-place write …")
    with h5py.File(edges_h5, "r+") as f:
        # Use a local field to get n_total — avoids following ExternalLinks in r+ mode
        n_total = len(f[f"{POP_PATH}/{LOCAL_SIZE_FIELD}"])
        logger.info(f"  Total synapses: {n_total:,}")

        for field_name, values in [("theta_d", theta_d), ("theta_p", theta_p)]:
            ds = _get_or_create_local_dataset(f, field_name, n_total)
            # Write values one index at a time (few hundred entries — fast enough)
            for idx, val in zip(global_idx.tolist(), values.tolist()):
                ds[idx] = val
            logger.info(f"  {field_name}: wrote {len(global_idx):,} values")

    logger.info(f"Done. Updated {len(global_idx):,} theta_d/theta_p entries in {edges_h5}")


def verify(edges_h5, global_idx, sample_n=10):
    """Spot-check a few entries to confirm the write succeeded."""
    rng = np.random.default_rng(42)
    sample = np.sort(rng.choice(global_idx, min(sample_n, len(global_idx)), replace=False))
    with h5py.File(edges_h5, "r") as f:
        td_vals = f[f"{POP_PATH}/theta_d"][sample]
        tp_vals = f[f"{POP_PATH}/theta_p"][sample]
    logger.info("Spot-check (global_idx → theta_d, theta_p):")
    for i, idx in enumerate(sample):
        logger.info(f"  [{idx}] theta_d={td_vals[i]:.4f}  theta_p={tp_vals[i]:.4f}")


def main():
    parser = argparse.ArgumentParser(
        description="Inject computed theta_d/theta_p into dhuruva_modified_edges.h5")
    parser.add_argument("--results-dir", required=True,
                        help="Directory containing per-post_gid .npz files from compute_thresholds.py")
    parser.add_argument("--edges-h5", default=EDGES_H5,
                        help=f"Path to edges HDF5 file (default: {EDGES_H5})")
    parser.add_argument("--dry-run", action="store_true",
                        help="Read results and report statistics but do not modify the file")
    parser.add_argument("--force", action="store_true",
                        help="Overwrite even when new thresholds differ from the current values in edges.h5")
    args = parser.parse_args()

    global_idx, theta_d, theta_p = collect_results(args.results_dir)
    inject(args.edges_h5, global_idx, theta_d, theta_p, dry_run=args.dry_run, force=args.force)

    if not args.dry_run:
        verify(args.edges_h5, global_idx)


if __name__ == "__main__":
    main()
