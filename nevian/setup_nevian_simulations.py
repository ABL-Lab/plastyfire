#!/usr/bin/env python3
"""
Setup Nevian & Sakmann 2006 simulation folders.
================================================
Creates workdir folders for each Nevian protocol under:
  /home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/nevian/simulations/

Each workdir is structured exactly as the existing L5TTPC simulations so that
the existing `pairrunner.py` + `simulator.runconnectedpair()` can run them
without any modification.

Per-workdir files created:
  - node_sets.json            (pre/post cell GIDs)
  - simulation_config.json    (main pairing run)
  - prefire_simulation_config.json  (pre-firing calibration run)
  - prespikes.h5              (presynaptic spike train for main run)
  - prefire_prespikes.h5      (presynaptic spike train for prefire run)

Protocol timing:
  Adapted from Nevian & Sakmann 2006:
    - 60 pairings at 0.1 Hz  (T = 10,000 ms between pairings)
    - 10-min baseline at 0.1 Hz (C01: 60 pre-only spikes)
    - 10-min post-induction at 0.1 Hz (C02: 60 post-only spikes)
    - dt convention: positive = post follows pre (LTP), negative = post precedes pre (LTD)

Usage:
    python setup_nevian_simulations.py [--pairs-dir DIR] [--output-dir DIR]
                                       [--protocols all|control|<id1,id2,...>]
                                       [--dry-run]
"""

import argparse
import json
import os
import sys
import numpy as np
import h5py
from pathlib import Path

# Add plastyfire to path
PLASTYFIRE_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PLASTYFIRE_ROOT))

from nevian_protocols import ALL_PROTOCOLS, CONTROL_PROTOCOLS, PROTOCOLS_BY_ID, NevianProtocol

# --- Paths ---
CIRCUIT_CONFIG   = "/project/ctb-emuller/datasets/SSCx_plasticity/O1_2023a_Ecker/circuit_config.json"
NODE_POP         = "S1nonbarrel_neurons"
TARGET           = "hex_O1Excitatory"
DEFAULT_PAIRS_DIR = str(PLASTYFIRE_ROOT / "refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC/simulations")
DEFAULT_OUT_DIR  = str(Path(__file__).parent / "simulations")
RANDOM_SEED      = 19091997

# --- Timing constants (matching paper: 60 pairings at 0.1 Hz) ---
N_PAIRINGS       = 60          # number of pre-post pairing repetitions
PAIRING_PERIOD   = 10_000.0    # ms between pairings (0.1 Hz)
OFFSET           = 1_000.0     # ms: initial delay before first spike
C01_DURATION_MIN = 10.0        # minutes of baseline (pre-only stimulation)
C02_DURATION_MIN = 10.0        # minutes of post-induction monitoring (pre-only)
C01_PERIOD       = 10_000.0    # ms between baseline spikes (0.1 Hz)
C02_PERIOD       = 10_000.0    # ms between post-induction spikes (0.1 Hz)
PULSE_WIDTH      = 3.0         # ms current pulse width (from existing configs)


def ms2min(ms: float) -> float:
    return ms / 60_000.0


def save_spikes_h5(h5f_name: str, node_pop: str, spike_times: np.ndarray, gid: int):
    """Save spike times in SONATA/libsonata HDF5 format (same as simwriter.save_spikes)."""
    n = len(spike_times)
    with h5py.File(h5f_name, "w") as h5f:
        grp = h5f.require_group(f"spikes/{node_pop}")
        grp.create_dataset("timestamps", data=spike_times)
        grp["timestamps"].attrs["units"] = "ms"
        grp.create_dataset("node_ids", data=np.full(n, gid, dtype=int))


def compute_spike_train(protocol: NevianProtocol,
                        pre_gid: int,
                        amp: float,
                        spike_delay_ms: float = 1.5) -> dict:
    """
    Compute all spike times for one protocol and return a dict with everything
    needed to fill `simulation_config.json` and `prespikes.h5`.

    Burst structure (Nevian paper):
      - n_aps APs at freq_hz
      - timing defined by dt_ms (closest AP to EPSP onset)
      - For dt < 0 (pre before EPSP): last AP in burst is at dt ms before EPSP
        → first AP in burst is at (dt - (n_aps-1)*ISI) ms before EPSP
      - For dt > 0 (post after EPSP): first AP in burst is at dt ms after EPSP
        → last AP in burst is at (dt + (n_aps-1)*ISI) ms after EPSP

    Timeline:
      0 ... OFFSET ... [C01 baseline] ... [pairing period] ... [C02 post] ... t_stop
      Pre spikes during baseline and post-induction phases: every C01/C02_PERIOD ms
      Pre spikes during pairing: relative to each pairing period repeat
      Post spikes (current pulses): timed so that the closest AP is dt_ms from EPSP
    """
    n_before = int(C01_DURATION_MIN * 60_000 / C01_PERIOD)
    n_after  = int(C02_DURATION_MIN * 60_000 / C02_PERIOD)

    # Baseline pre-spikes: single pre-spike every C01_PERIOD, used to measure Cpre
    before_pre = np.array([OFFSET + i * C01_PERIOD for i in range(n_before)])

    # Pairing phase start time
    t_pairing_start = OFFSET + n_before * C01_PERIOD

    # ISI within burst (ms)
    isi_burst = 1000.0 / protocol.freq_hz if protocol.freq_hz > 0 else 0.0

    # Compute time of first AP in the burst relative to EPSP (at t=0 within each pairing period)
    # dt_ms = time of closest AP to EPSP
    # For dt > 0: first AP is at dt_ms (all APs come after EPSP)
    # For dt < 0: last AP is at dt_ms (all APs come before EPSP), first at dt_ms - (n_aps-1)*isi
    if protocol.n_aps == 0:
        # EPSP-only: no post APs, no pairing pre-burst either (just baseline pre-spikes)
        first_ap_offset = None
    elif protocol.n_aps == 1:
        first_ap_offset = protocol.dt_ms
    else:
        if protocol.dt_ms > 0:
            # APs follow EPSP: first AP at dt_ms, ISI between each
            first_ap_offset = protocol.dt_ms
        else:
            # APs precede EPSP: last AP at dt_ms (closest), so first AP earlier
            first_ap_offset = protocol.dt_ms - (protocol.n_aps - 1) * isi_burst

    # Pre-spikes during pairing: one pre-spike per pairing period, at t=0 within period
    # (EPSP is elicited by the pre-spike; post APs are timed relative to this)
    # Pre-spike arrives at t_pairing_start + j*PAIRING_PERIOD
    # Post burst offset = first_ap_offset (relative to EPSP)
    pairing_pre = np.array([
        t_pairing_start + j * PAIRING_PERIOD + spike_delay_ms
        for j in range(N_PAIRINGS)
    ])

    # Post-induction pre-spikes
    t_post_start = t_pairing_start + N_PAIRINGS * PAIRING_PERIOD
    after_pre = np.array([t_post_start + i * C02_PERIOD for i in range(n_after)])

    # Full pre-spike train
    all_pre_spikes = np.concatenate([before_pre, pairing_pre, after_pre])

    # t_stop: end of post-induction period
    t_stop = t_post_start + n_after * C02_PERIOD

    # ----------------------------------------------------------------
    # Postsynaptic current pulses (one per AP per pairing period)
    # The pulse "delay" in SONATA inputs = absolute time of pulse onset
    # We generate one pulse entry per AP position per pairing period
    # ----------------------------------------------------------------
    inputs = {}
    pulse_idx = 0
    if first_ap_offset is not None and protocol.n_aps > 0:
        for j in range(N_PAIRINGS):
            # EPSP time in this period
            t_epsp = t_pairing_start + j * PAIRING_PERIOD
            for k in range(protocol.n_aps):
                t_ap = t_epsp + first_ap_offset + k * isi_burst
                inputs[f"pulse{pulse_idx}"] = {
                    "input_type": "current_clamp",
                    "module": "pulse",
                    "node_set": TARGET,
                    "delay": t_ap,
                    "duration": PAIRING_PERIOD,   # duration keeps the stimulus active
                    "amp_start": amp,
                    "width": PULSE_WIDTH,
                    "frequency": 1000.0 / PAIRING_PERIOD,
                }
                pulse_idx += 1

    # ----------------------------------------------------------------
    # Prefire spike train (same structure but without C01/C02 bookends)
    # Used by prefire_simulation_config.json to calibrate
    # ----------------------------------------------------------------
    prefire_pre = np.array([
        OFFSET + j * PAIRING_PERIOD + spike_delay_ms
        for j in range(N_PAIRINGS)
    ])
    prefire_t_stop = OFFSET + N_PAIRINGS * PAIRING_PERIOD + 1000.0

    prefire_inputs = {}
    pulse_idx = 0
    if first_ap_offset is not None and protocol.n_aps > 0:
        for j in range(N_PAIRINGS):
            t_epsp = OFFSET + j * PAIRING_PERIOD
            for k in range(protocol.n_aps):
                t_ap = t_epsp + first_ap_offset + k * isi_burst
                prefire_inputs[f"pulse{pulse_idx}"] = {
                    "input_type": "current_clamp",
                    "module": "pulse",
                    "node_set": TARGET,
                    "delay": t_ap,
                    "duration": PAIRING_PERIOD,
                    "amp_start": amp,
                    "width": PULSE_WIDTH,
                    "frequency": 1000.0 / PAIRING_PERIOD,
                }
                pulse_idx += 1

    return {
        "all_pre_spikes": all_pre_spikes,
        "pairing_pre_spikes": pairing_pre,
        "prefire_pre_spikes": prefire_pre,
        "inputs": inputs,
        "prefire_inputs": prefire_inputs,
        "t_stop": t_stop,
        "prefire_t_stop": prefire_t_stop,
        "t_pairing_start": t_pairing_start,
        "n_before": n_before,
        "n_after": n_after,
    }


def write_workdir(workdir: str,
                  pre_gid: int,
                  post_gid: int,
                  protocol: NevianProtocol,
                  amp: float = 2.0,
                  spike_delay_ms: float = 1.5,
                  dry_run: bool = False) -> None:
    """
    Create the workdir and all required files for one (protocol, pair) simulation.

    Parameters
    ----------
    workdir     : absolute path of the simulation directory to create
    pre_gid     : presynaptic cell GID
    post_gid    : postsynaptic cell GID
    protocol    : NevianProtocol instance
    amp         : current injection amplitude (nA) — use value from single_cells cache if available
    spike_delay_ms : AP latency after current pulse onset (ms)
    dry_run     : if True, only print what would be created, don't write files
    """
    if not dry_run:
        os.makedirs(workdir, exist_ok=True)

    # --- node_sets.json ---
    node_sets = {
        "precell":  {"node_id": [int(pre_gid)],  "population": NODE_POP},
        "postcell": {"node_id": [int(post_gid)], "population": NODE_POP},
    }
    node_sets_path = os.path.join(workdir, "node_sets.json")

    # --- Spike trains and pulse inputs ---
    train = compute_spike_train(protocol, pre_gid, amp, spike_delay_ms)

    # --- simulation_config.json ---
    sim_config = {
        "run": {"dt": 0.025, "tstop": train["t_stop"], "random_seed": RANDOM_SEED},
        "network": CIRCUIT_CONFIG,
        "node_sets_file": node_sets_path,
        "node_set": "postcell",
        "output": {"output_dir": os.path.join(workdir, "out")},
        "inputs": train["inputs"],
        "connection_overrides": [
            {"name": "plasticity", "source": TARGET, "target": TARGET,
             "modoverride": "GluSynapse", "weight": 1.0}
        ],
    }

    # --- prefire_simulation_config.json ---
    prefire_sim_config = {
        "run": {"dt": 0.025, "tstop": train["prefire_t_stop"], "random_seed": RANDOM_SEED},
        "network": CIRCUIT_CONFIG,
        "node_sets_file": node_sets_path,
        "node_set": "postcell",
        "output": {"output_dir": os.path.join(workdir, "out")},
        "inputs": train["prefire_inputs"],
        "connection_overrides": [
            {"name": "plasticity", "source": TARGET, "target": TARGET,
             "modoverride": "GluSynapse", "weight": 1.0}
        ],
    }

    if dry_run:
        print(f"[DRY RUN] Would create: {workdir}/")
        print(f"  node_sets.json       pre={pre_gid}, post={post_gid}")
        print(f"  prespikes.h5         {len(train['all_pre_spikes'])} spikes, t_stop={train['t_stop']:.0f}ms")
        print(f"  prefire_prespikes.h5 {len(train['prefire_pre_spikes'])} spikes")
        print(f"  simulation_config.json  inputs: {len(train['inputs'])} pulses")
        print(f"  prefire_simulation_config.json  inputs: {len(train['prefire_inputs'])} pulses")
        return

    # Write files
    with open(node_sets_path, "w") as f:
        json.dump(node_sets, f, indent=4)

    with open(os.path.join(workdir, "simulation_config.json"), "w") as f:
        json.dump(sim_config, f, indent=4)

    with open(os.path.join(workdir, "prefire_simulation_config.json"), "w") as f:
        json.dump(prefire_sim_config, f, indent=4)

    prespikes_path = os.path.join(workdir, "prespikes.h5")
    save_spikes_h5(prespikes_path, NODE_POP, train["all_pre_spikes"], pre_gid)

    prefire_prespikes_path = os.path.join(workdir, "prefire_prespikes.h5")
    save_spikes_h5(prefire_prespikes_path, NODE_POP, train["prefire_pre_spikes"], pre_gid)

    print(f"  Created: {workdir}")


def read_amp_from_cache(single_cells_dir: str, post_gid: int,
                        freq_hz: float = 50.0) -> float:
    """
    Try to read the pre-calibrated current amplitude from the single_cells cache
    (same format as written by simwriter.check_electrical_constraint).
    Falls back to 2.0 nA if not found.
    """
    import pickle
    pkl_path = os.path.join(single_cells_dir, f"{post_gid}.pkl")
    if os.path.isfile(pkl_path):
        try:
            with open(pkl_path, "rb") as f:
                data = pickle.load(f)
            # data is a dict keyed by frequency
            # pick the closest available frequency
            avail = list(data.keys())
            if not avail:
                return 2.0
            closest = min(avail, key=lambda x: abs(x - freq_hz))
            return float(data[closest]["amp"])
        except Exception as e:
            print(f"  Warning: could not read amp cache for gid {post_gid}: {e}")
    return 2.0  # default fallback


def discover_pairs(pairs_dir: str):
    """
    Discover (pre_gid, post_gid) pairs from an existing simulations directory.
    Each pair folder is named '{pre_gid}-{post_gid}'.
    """
    pairs = []
    for entry in sorted(os.listdir(pairs_dir)):
        if entry == "single_cells":
            continue
        parts = entry.split("-")
        if len(parts) == 2:
            try:
                pre_gid  = int(parts[0])
                post_gid = int(parts[1])
                pairs.append((pre_gid, post_gid))
            except ValueError:
                pass
    return pairs


def main():
    parser = argparse.ArgumentParser(
        description="Create Nevian & Sakmann 2006 simulation folders",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Create all control protocol folders using existing L5TTPC pairs
  python setup_nevian_simulations.py

  # Only control protocols, dry run first
  python setup_nevian_simulations.py --protocols control --dry-run

  # Specific protocols only
  python setup_nevian_simulations.py --protocols LTP_3ap_50hz_dt+10ms,LTD_3ap_50hz_dt-10ms

  # Use first 5 pairs only (for testing)
  python setup_nevian_simulations.py --max-pairs 5 --dry-run

  # Custom pairs directory
  python setup_nevian_simulations.py --pairs-dir /path/to/L5TTPC/simulations
        """
    )
    parser.add_argument("--pairs-dir", default=DEFAULT_PAIRS_DIR,
                        help=f"Directory containing pre-post pair folders (default: {DEFAULT_PAIRS_DIR})")
    parser.add_argument("--output-dir", default=DEFAULT_OUT_DIR,
                        help=f"Output root directory (default: {DEFAULT_OUT_DIR})")
    parser.add_argument("--protocols", default="control",
                        help="Protocols to set up: 'all', 'control', or comma-separated protocol IDs")
    parser.add_argument("--max-pairs", type=int, default=None,
                        help="Limit to first N pairs (useful for testing)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print what would be created without writing any files")
    parser.add_argument("--default-amp", type=float, default=2.0,
                        help="Default current amplitude if not in cache (nA, default: 2.0)")
    parser.add_argument("--spike-delay", type=float, default=1.5,
                        help="AP latency after current pulse onset (ms, default: 1.5)")
    args = parser.parse_args()

    # Select protocols
    if args.protocols == "all":
        protocols = ALL_PROTOCOLS
    elif args.protocols == "control":
        protocols = CONTROL_PROTOCOLS
    else:
        ids = [p.strip() for p in args.protocols.split(",")]
        protocols = []
        for pid in ids:
            if pid not in PROTOCOLS_BY_ID:
                print(f"ERROR: Unknown protocol '{pid}'. Available: {list(PROTOCOLS_BY_ID.keys())}")
                sys.exit(1)
            protocols.append(PROTOCOLS_BY_ID[pid])

    # Discover pairs
    pairs = discover_pairs(args.pairs_dir)
    if not pairs:
        print(f"ERROR: No pairs found in {args.pairs_dir}")
        sys.exit(1)
    if args.max_pairs:
        pairs = pairs[:args.max_pairs]

    single_cells_dir = os.path.join(args.pairs_dir, "single_cells")

    print(f"Pairs directory : {args.pairs_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Protocols       : {len(protocols)} protocols")
    print(f"Pairs           : {len(pairs)} pairs")
    print(f"Dry run         : {args.dry_run}")
    print()

    if not args.dry_run:
        os.makedirs(args.output_dir, exist_ok=True)

    total = 0
    for pre_gid, post_gid in pairs:
        pair_name = f"{pre_gid}-{post_gid}"
        # Read amplitude from cache (uses frequency of 50Hz as default lookup key)
        amp = read_amp_from_cache(single_cells_dir, post_gid,
                                  freq_hz=50.0) if os.path.isdir(single_cells_dir) else args.default_amp

        for protocol in protocols:
            # workdir: output_dir / pair_name / protocol_id
            workdir = os.path.join(args.output_dir, pair_name, protocol.protocol_id)
            write_workdir(
                workdir=workdir,
                pre_gid=pre_gid,
                post_gid=post_gid,
                protocol=protocol,
                amp=amp,
                spike_delay_ms=args.spike_delay,
                dry_run=args.dry_run,
            )
            total += 1

    print(f"\n{'[DRY RUN] Would create' if args.dry_run else 'Created'} {total} simulation folders.")
    if not args.dry_run:
        print(f"\nNext step: run submit_nevian.py to submit to SLURM.")


if __name__ == "__main__":
    main()
