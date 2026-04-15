#!/usr/bin/env python3
"""
Compute theta_d / theta_p thresholds for specific (pre_gid, post_gid) pairs.

One submitit job per pair. Each job:
  - C_pre : fires pre_gid once, records peak effcai_GB at its synapses on post_gid
  - C_post: injects current into post_gid (1 bAP), records peak effcai_GB at same synapses
  - theta_d/theta_p computed per synapse (basal vs apical)

Results saved to OUTPUT_DIR/{pre_gid}_{post_gid}.npz with keys:
    global_idx (int64[N]), theta_d (float32[N]), theta_p (float32[N])

Usage:
    python compute_thresholds.py \\
        --pairs-file /path/to/pairs.txt \\
        --output-dir /path/to/results

    # Dry run
    python compute_thresholds.py --pairs-file pairs.txt --output-dir /tmp/out --dry-run
"""

import argparse
import logging
import os
import sys
from pathlib import Path

import h5py
import numpy as np

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

# ── Paths ────────────────────────────────────────────────────────────────────
PLASTYFIRE_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire"
EDGES_H5       = "/project/ctb-emuller/dhuruva/plastyfire/data/dhuruva_modified_edges.h5"
CIRCUIT_CONFIG = "/project/ctb-emuller/dhuruva/plastyfire/data/dhuruva_circuit_config.json"
NODES_H5       = "/project/ctb-emuller/datasets/SSCx_plasticity/O1_2023a_Ecker/S1nonbarrel_neurons/nodes.h5"
NODE_POP       = "S1nonbarrel_neurons"
EDGE_POP       = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"
POP_BASE       = f"edges/{EDGE_POP}"
POP_PATH       = f"edges/{EDGE_POP}/0"

ALL_L5_MTYPES = [
    "L5_TPC:A", "L5_TPC:B", "L5_TPC:C", "L5_UPC",
    "L5_BP", "L5_BTC", "L5_CHC", "L5_DBC", "L5_LBC",
    "L5_MC", "L5_NBC", "L5_NGC", "L5_SBC",
]
EXC_TYPE_MIN = 100

TAU_EFFCA = 278.3177658387

# Threshold linear coefficients (Chindemi defaults)
A_COEFF = {
    "a00": 1.002, "a01": 1.954,   # basal  depression
    "a10": 1.159, "a11": 2.483,   # basal  potentiation
    "a20": 1.127, "a21": 2.456,   # apical depression
    "a30": 5.236, "a31": 1.782,   # apical potentiation
}

BASAL_SECTION_TYPE  = 2
APICAL_SECTION_TYPE = 3

C_PRE_T0_MS     = 1000.0
C_PRE_TSTOP_MS  = 1500.0   # same as original c_pre_finder_process

C_POST_T0_MS    = 1000.0
C_POST_WIDTH_MS = 3.0
C_POST_TSTOP_MS = 3000.0

SINGLE_CELLS_DIRS = [
    "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations/single_cells",
]


def load_l5_node_ids(mtypes=None):
    if mtypes is None:
        mtypes = ALL_L5_MTYPES
    with h5py.File(NODES_H5, "r") as h:
        pop = f"nodes/{NODE_POP}/0"
        mtype_ids = h[f"{pop}/mtype"][:]
        library   = h[f"{pop}/@library/mtype"][:]
    library_str = np.array([s.decode() if isinstance(s, bytes) else s for s in library])
    target_ids  = np.where(np.isin(library_str, mtypes))[0].astype(np.uint32)
    return np.where(np.isin(mtype_ids, target_ids))[0].astype(np.int32)


def _get_spike_amp(post_gid):
    import pickle
    for d in SINGLE_CELLS_DIRS:
        pkl = os.path.join(d, f"{post_gid}.pkl")
        if os.path.exists(pkl):
            with open(pkl, "rb") as f:
                data = pickle.load(f)
            key = list(data.keys())[0]
            return float(data[key]["amp"])
    return 1.5   # safe default for L5 TPC


def _get_pair_edges(pre_gid, post_gid, l5_ids):
    """Return (global_idx, sec_type) arrays for EXC synapses from pre→post."""
    with h5py.File(EDGES_H5, "r") as f:
        n2r = f[f"{POP_BASE}/indices/target_to_source/node_id_to_ranges"]
        r2e = f[f"{POP_BASE}/indices/target_to_source/range_to_edge_id"]
        r0, r1 = int(n2r[post_gid][0]), int(n2r[post_gid][1])
        if r0 == r1:
            return None, None
        edge_slices = r2e[r0:r1]
        src_parts, stype_parts, sec_parts, gidx_parts = [], [], [], []
        for e_start, e_end in edge_slices:
            e0, e1 = int(e_start), int(e_end)
            src_parts.append(f[f"{POP_BASE}/source_node_id"][e0:e1].astype(np.int32))
            stype_parts.append(f[f"{POP_PATH}/syn_type_id"][e0:e1])
            sec_parts.append(f[f"{POP_PATH}/afferent_section_type"][e0:e1])
            gidx_parts.append(np.arange(e0, e1, dtype=np.int64))

    src    = np.concatenate(src_parts)
    stype  = np.concatenate(stype_parts)
    sec    = np.concatenate(sec_parts)
    gidx   = np.concatenate(gidx_parts)

    mask = (src == pre_gid) & (stype >= EXC_TYPE_MIN)
    if not mask.any():
        return None, None
    return gidx[mask], sec[mask]


def _worker(pre_gid, post_gid, output_dir):
    """Compute C_pre/C_post for one pair and write theta_d/theta_p."""
    os.environ["GLOG_minloglevel"] = "3"
    sys.path.insert(0, PLASTYFIRE_DIR)

    import tempfile, json
    import numpy as np
    import h5py
    import pandas as pd
    import bluecellulab
    from bluepysnap import Simulation
    from conntility.io.synapse_report import get_presyn_mapping

    output_file = os.path.join(output_dir, f"{pre_gid}_{post_gid}.npz")
    if os.path.exists(output_file):
        return f"skip:{pre_gid}_{post_gid}"

    # ── 1. Load L5 ids and find edges for this pair ───────────────────────────
    l5_ids = load_l5_node_ids()
    global_idx, sec_arr = _get_pair_edges(pre_gid, post_gid, l5_ids)
    if global_idx is None:
        return f"no_edges:{pre_gid}_{post_gid}"

    global_idx_set = set(global_idx.tolist())
    sec_lookup = dict(zip(global_idx.tolist(), sec_arr.tolist()))

    # ── 2. Build temp sim config ──────────────────────────────────────────────
    tmp_dir = Path(tempfile.mkdtemp())
    sim_cfg = {
        "run": {"dt": 0.025, "tstop": C_PRE_TSTOP_MS, "random_seed": 12345},
        "network": CIRCUIT_CONFIG,
        "node_set": "hex_O1",
        "output": {"output_dir": str(tmp_dir / "out")},
        "connection_overrides": [
            {"name": "plasticity",
             "source": "hex_O1", "target": "hex_O1",
             "modoverride": "GluSynapse", "weight": 1.0}
        ]
    }
    (tmp_dir / "out").mkdir()
    sim_config_path = str(tmp_dir / "sim_config.json")
    with open(sim_config_path, "w") as fh:
        json.dump(sim_cfg, fh)

    bluecellulab.set_verbose(0)
    bluecellulab.neuron.load_mechanisms("/project/ctb-emuller/dhuruva/DEES_cell_packages/")
    bluecellulab.neuron.h.cvode.atolscale("v", 0.1)

    intersect_list = [(NODE_POP, int(pre_gid))]

    # ── 3. Map local → global syn idx ────────────────────────────────────────
    sim_probe = bluecellulab.CircuitSimulation(sim_config_path)
    sim_probe.instantiate_gids(
        [(NODE_POP, int(post_gid))],
        add_synapses=True, add_minis=False,
        intersect_pre_gids=intersect_list,
    )
    cell_probe = sim_probe.cells[(NODE_POP, int(post_gid))]
    local_syn_idxs = [sid[1] for sid in cell_probe.synapses]
    if not local_syn_idxs:
        return f"no_synapses:{pre_gid}_{post_gid}"

    mi = pd.MultiIndex.from_tuples([(int(post_gid), s) for s in local_syn_idxs])
    presyn_df = get_presyn_mapping(Simulation(sim_config_path).circuit, EDGE_POP, mi)
    presyn_df = presyn_df[presyn_df.index.isin(global_idx_set)]
    if len(presyn_df) == 0:
        return f"no_presyn_match:{pre_gid}_{post_gid}"

    local2global = dict(zip(presyn_df["local_syn_idx"], presyn_df.index))
    del sim_probe, cell_probe

    # ── 4. C_pre ──────────────────────────────────────────────────────────────
    bluecellulab.neuron.h.tau_effca_GB_GluSynapse = TAU_EFFCA

    sim_cpre = bluecellulab.CircuitSimulation(sim_config_path)
    sim_cpre.instantiate_gids(
        [(NODE_POP, int(post_gid))],
        add_synapses=True, add_minis=False,
        pre_spike_trains={(NODE_POP, int(pre_gid)): [C_PRE_T0_MS]},
        intersect_pre_gids=intersect_list,
    )
    cell_cpre = sim_cpre.cells[(NODE_POP, int(post_gid))]

    for syn_id, syn in cell_cpre.synapses.items():
        h = syn.hsynapse
        h.rho0_GB    = 1
        h.Use_p      = 1
        h.Use        = 1
        h.gmax0_AMPA = h.gmax_p_AMPA
        h.theta_d_GB = -1
        h.theta_p_GB = -1

    recs_cpre = {sid[1]: bluecellulab.neuron.h.Vector() for sid in cell_cpre.synapses}
    for sid, syn in cell_cpre.synapses.items():
        recs_cpre[sid[1]].record(syn.hsynapse._ref_effcai_GB, 1.0)

    sim_cpre.run(C_PRE_TSTOP_MS, cvode=True)

    c_pre_map = {}
    for local_idx, rec in recs_cpre.items():
        if local_idx not in local2global:
            continue
        arr = np.array(rec.to_python())
        i0 = max(0, int(C_PRE_T0_MS))
        c_pre_map[local2global[local_idx]] = float(arr[i0:].max()) if len(arr) > i0 else 0.0

    del sim_cpre, cell_cpre, recs_cpre

    # ── 5. C_post ─────────────────────────────────────────────────────────────
    amp = _get_spike_amp(post_gid)

    sim_cpost = bluecellulab.CircuitSimulation(sim_config_path)
    sim_cpost.instantiate_gids(
        [(NODE_POP, int(post_gid))],
        add_synapses=True, add_minis=False,
        intersect_pre_gids=intersect_list,
    )
    cell_cpost = sim_cpost.cells[(NODE_POP, int(post_gid))]

    for sec in cell_cpost.somatic + cell_cpost.axonal:
        sec.uninsert("SK_E2")
    for syn_id, syn in cell_cpost.synapses.items():
        syn.hsynapse.theta_d_GB = -1
        syn.hsynapse.theta_p_GB = -1

    recs_cpost = {sid[1]: bluecellulab.neuron.h.Vector() for sid in cell_cpost.synapses}
    for sid, syn in cell_cpost.synapses.items():
        recs_cpost[sid[1]].record(syn.hsynapse._ref_effcai_GB, 1.0)

    stim = bluecellulab.neuron.h.IClamp(0.5, sec=cell_cpost.soma)
    stim.delay = C_POST_T0_MS
    stim.dur   = C_POST_WIDTH_MS
    stim.amp   = amp

    sim_cpost.run(C_POST_TSTOP_MS, cvode=True)

    c_post_map = {}
    for local_idx, rec in recs_cpost.items():
        if local_idx not in local2global:
            continue
        arr = np.array(rec.to_python())
        c_post_map[local2global[local_idx]] = float(arr.max()) if len(arr) > 0 else 0.0

    del sim_cpost, cell_cpost, recs_cpost

    # ── 6. Compute theta_d / theta_p ─────────────────────────────────────────
    out_global, out_td, out_tp = [], [], []
    for gidx_i in presyn_df.index:
        c_pre  = c_pre_map.get(gidx_i, 0.0)
        c_post = c_post_map.get(gidx_i, 0.0)
        sec_t  = sec_lookup.get(int(gidx_i), -1)

        if sec_t == BASAL_SECTION_TYPE:
            td = A_COEFF["a00"] * c_pre + A_COEFF["a01"] * c_post
            tp = A_COEFF["a10"] * c_pre + A_COEFF["a11"] * c_post
        elif sec_t == APICAL_SECTION_TYPE:
            td = A_COEFF["a20"] * c_pre + A_COEFF["a21"] * c_post
            tp = A_COEFF["a30"] * c_pre + A_COEFF["a31"] * c_post
        else:
            td, tp = -1.0, -1.0

        out_global.append(gidx_i)
        out_td.append(td)
        out_tp.append(tp)

    np.savez_compressed(
        output_file,
        global_idx=np.array(out_global, dtype=np.int64),
        theta_d=np.array(out_td, dtype=np.float32),
        theta_p=np.array(out_tp, dtype=np.float32),
    )
    return f"ok:{pre_gid}_{post_gid} n_syn={len(out_global)} Cpre_mean={np.mean(list(c_pre_map.values())):.4f} Cpost_mean={np.mean(list(c_post_map.values())):.4f}"


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Compute theta_d/theta_p for specific pairs")
    parser.add_argument("--pairs-file",
                        help="Text file with 'pre_gid post_gid' per line")
    parser.add_argument("--single-pair", nargs=2, type=int, metavar=("PRE_GID", "POST_GID"),
                        help="Run worker directly for a single pair (no submitit)")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--timeout-hours", type=float, default=1.0)
    parser.add_argument("--mem-gb", type=int, default=16)
    parser.add_argument("--slurm-account", default="ctb-emuller")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Direct single-pair mode: run _worker inline (used by individual sbatch jobs)
    if args.single_pair:
        pre_gid, post_gid = args.single_pair
        npz = output_dir / f"{pre_gid}_{post_gid}.npz"
        if npz.exists():
            logger.info(f"Already done: {pre_gid}_{post_gid}.npz — skipping")
            return
        logger.info(f"Running single pair: {pre_gid} -> {post_gid}")
        _worker(pre_gid, post_gid, str(output_dir))
        return

    if not args.pairs_file:
        parser.error("--pairs-file is required unless --single-pair is used")

    pairs = []
    with open(args.pairs_file) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) == 2:
                pairs.append((int(parts[0]), int(parts[1])))
    logger.info(f"Loaded {len(pairs)} pairs from {args.pairs_file}")

    # Skip already done
    pairs = [(pre, post) for pre, post in pairs
             if not os.path.exists(output_dir / f"{pre}_{post}.npz")]
    logger.info(f"  {len(pairs)} pairs still to process")

    if args.dry_run:
        logger.info(f"Dry run: would submit {len(pairs)} jobs. Exiting.")
        return

    import submitit
    log_folder = output_dir / "logs"
    log_folder.mkdir(exist_ok=True)

    executor = submitit.AutoExecutor(folder=str(log_folder))
    executor.update_parameters(
        timeout_min=int(args.timeout_hours * 60),
        cpus_per_task=1,
        mem_gb=args.mem_gb,
        slurm_account=args.slurm_account,
        slurm_setup=[f"source {PLASTYFIRE_DIR}/setupenv.sh"],
    )

    jobs = []
    with executor.batch():
        for pre_gid, post_gid in pairs:
            j = executor.submit(_worker, int(pre_gid), int(post_gid), str(output_dir))
            jobs.append(j)

    logger.info(f"Submitted {len(jobs)} jobs. Monitor: squeue -u $USER")
    logger.info(f"Logs: {log_folder}/")
    logger.info(f"When done, run: python inject_thresholds.py --results-dir {output_dir}")


if __name__ == "__main__":
    main()
