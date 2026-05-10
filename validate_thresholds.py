#!/usr/bin/env python3
"""
Validate stored theta_d / theta_p threshold values against a fresh bluecellulab
c_pre / c_post calibration run via plastyfire.simulator.

For a chosen pair:
  1. Load the stored theta_d / theta_p from data/thresholds_n100/<pre>_<post>.npz
  2. Run c_pre_finder  -> c_pre  per synapse
  3. Run c_post_finder -> c_post per synapse
  4. Read section type (basal / apical) from the edges HDF5
  5. Recompute theta_d / theta_p using the Chindemi A-coefficients
  6. Print a side-by-side comparison table
  7. Save a summary figure

Usage:
    python validate_thresholds.py                        # first pair found
    python validate_thresholds.py --pre 180164 --post 197248
    python validate_thresholds.py --npz data/thresholds_n100/180164_197248.npz
"""

import argparse
import glob
import json
import logging
import os
import pickle
import sys
import tempfile

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

# ── Paths ─────────────────────────────────────────────────────────────────────
PLASTYFIRE_DIR   = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire"
CIRCUIT_CONFIG   = "/project/ctb-emuller/dhuruva/plastyfire/data/dhuruva_circuit_config.json"
EDGES_H5         = "/project/ctb-emuller/dhuruva/plastyfire/data/dhuruva_modified_edges.h5"
THRESHOLDS_DIR   = os.path.join(PLASTYFIRE_DIR, "data", "thresholds_n100")
SINGLE_CELLS_DIR = os.path.join(
    PLASTYFIRE_DIR,
    "refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations/single_cells",
)

EDGE_POP  = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"
NODE_POP  = "S1nonbarrel_neurons"
POP_BASE  = f"edges/{EDGE_POP}"
POP_PATH  = f"edges/{EDGE_POP}/0"

BASAL_SECTION_TYPE  = 2
APICAL_SECTION_TYPE = 3
EXC_TYPE_MIN        = 100

# Chindemi A-coefficients (same as compute_thresholds.py)
A_COEFF = {
    "a00": 1.002, "a01": 1.954,   # basal  depression
    "a10": 1.159, "a11": 2.483,   # basal  potentiation
    "a20": 1.127, "a21": 2.456,   # apical depression
    "a30": 5.236, "a31": 1.782,   # apical potentiation
}
TAU_EFFCA = 278.3177658387


# ── Helpers ───────────────────────────────────────────────────────────────────

def _pick_npz(pre_gid=None, post_gid=None, npz_path=None):
    if npz_path:
        return npz_path
    if pre_gid and post_gid:
        p = os.path.join(THRESHOLDS_DIR, f"{pre_gid}_{post_gid}.npz")
        if not os.path.exists(p):
            raise FileNotFoundError(f"No .npz at {p}")
        return p
    # auto-pick first file
    files = sorted(glob.glob(os.path.join(THRESHOLDS_DIR, "*.npz")))
    if not files:
        raise FileNotFoundError(f"No .npz files in {THRESHOLDS_DIR}")
    return files[0]


def _load_npz(npz_path):
    d = np.load(npz_path)
    return (d["global_idx"].astype(np.int64),
            d["theta_d"].astype(np.float64),
            d["theta_p"].astype(np.float64))


def _get_section_types(pre_gid, post_gid):
    """Return dict {global_idx: section_type} for EXC synapses pre->post."""
    with h5py.File(EDGES_H5, "r") as f:
        n2r = f[f"{POP_BASE}/indices/target_to_source/node_id_to_ranges"]
        r2e = f[f"{POP_BASE}/indices/target_to_source/range_to_edge_id"]
        r0, r1 = int(n2r[post_gid][0]), int(n2r[post_gid][1])
        if r0 == r1:
            return {}
        edge_slices = r2e[r0:r1]
        src_parts, stype_parts, sec_parts, gidx_parts = [], [], [], []
        for e_start, e_end in edge_slices:
            e0, e1 = int(e_start), int(e_end)
            src_parts.append(f[f"{POP_BASE}/source_node_id"][e0:e1].astype(np.int32))
            stype_parts.append(f[f"{POP_PATH}/syn_type_id"][e0:e1])
            sec_parts.append(f[f"{POP_PATH}/afferent_section_type"][e0:e1])
            gidx_parts.append(np.arange(e0, e1, dtype=np.int64))

    src   = np.concatenate(src_parts)
    stype = np.concatenate(stype_parts)
    sec   = np.concatenate(sec_parts)
    gidx  = np.concatenate(gidx_parts)

    mask = (src == pre_gid) & (stype >= EXC_TYPE_MIN)
    return {int(gidx[i]): int(sec[i]) for i in np.where(mask)[0]}


def _get_spike_amp(post_gid):
    """Read calibrated spike amplitude from single_cells pkl."""
    pkl = os.path.join(SINGLE_CELLS_DIR, f"{post_gid}.pkl")
    if os.path.exists(pkl):
        with open(pkl, "rb") as f:
            d = pickle.load(f)
        key = list(d.keys())[0]
        return float(d[key]["amp"]), float(d[key]["width"])
    # fallback: read from prefire_simulation_config if available
    cfg_glob = glob.glob(
        os.path.join(PLASTYFIRE_DIR,
                     f"refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP"
                     f"/simulations/*-{post_gid}/*/prefire_simulation_config.json")
    )
    if cfg_glob:
        with open(cfg_glob[0]) as f:
            cfg = json.load(f)
        amp   = cfg["inputs"]["pulse0"]["amp_start"]
        width = cfg["inputs"]["pulse0"]["width"]
        return float(amp), float(width)
    logger.warning("No spike amp found for %d, using default 1.5 nA", post_gid)
    return 1.5, 3.0


def _make_sim_config(tmp_dir):
    """Write a minimal bluecellulab SONATA sim config to tmp_dir."""
    out_dir = os.path.join(tmp_dir, "out")
    os.makedirs(out_dir, exist_ok=True)
    cfg = {
        "run": {
            "dt": 0.025,
            "tstop": 3000.0,
            "random_seed": 12345,
        },
        "network": CIRCUIT_CONFIG,
        "node_set": "hex_O1",
        "output": {"output_dir": out_dir},
        "conditions": {
            "extracellular_calcium": 2.0,
            "v_init": -80.0,
            "mechanisms": {
                "GluSynapse": {
                    "cao_CR": 2.0,
                    "gamma_d_GB": 101.5387594661,
                    "gamma_p_GB": 216.1841700668,
                    "minis_single_vesicle": False,
                }
            },
        },
        "connection_overrides": [
            {
                "name": "plasticity",
                "source": "hex_O1",
                "target": "hex_O1",
                "modoverride": "GluSynapse",
                "weight": 1.0,
            }
        ],
    }
    path = os.path.join(tmp_dir, "sim_config.json")
    with open(path, "w") as fh:
        json.dump(cfg, fh, indent=2)
    return path


def _compute_thresholds(c_pre_map, c_post_map, sec_type_map):
    """Recompute theta_d / theta_p from c_pre, c_post, section types."""
    td, tp = {}, {}
    for gidx, c_pre in c_pre_map.items():
        c_post  = c_post_map.get(gidx, 0.0)
        sec_t   = sec_type_map.get(int(gidx), -1)
        if sec_t == BASAL_SECTION_TYPE:
            td[gidx] = A_COEFF["a00"] * c_pre + A_COEFF["a01"] * c_post
            tp[gidx] = A_COEFF["a10"] * c_pre + A_COEFF["a11"] * c_post
        elif sec_t == APICAL_SECTION_TYPE:
            td[gidx] = A_COEFF["a20"] * c_pre + A_COEFF["a21"] * c_post
            tp[gidx] = A_COEFF["a30"] * c_pre + A_COEFF["a31"] * c_post
        else:
            td[gidx] = -1.0
            tp[gidx] = -1.0
    return td, tp


def _print_table(global_idx, stored_td, stored_tp, c_pre_map, c_post_map,
                 recomp_td, recomp_tp, sec_type_map):
    sec_name = {BASAL_SECTION_TYPE: "basal", APICAL_SECTION_TYPE: "apical"}
    hdr = (f"{'global_idx':>12}  {'loc':>6}  "
           f"{'c_pre':>10}  {'c_post':>10}  "
           f"{'stored_td':>10}  {'recomp_td':>10}  {'Δtd':>8}  "
           f"{'stored_tp':>10}  {'recomp_tp':>10}  {'Δtp':>8}")
    print("\n" + "=" * len(hdr))
    print(hdr)
    print("-" * len(hdr))

    max_delta_td = max_delta_tp = 0.0
    for gidx, s_td, s_tp in zip(global_idx, stored_td, stored_tp):
        cp  = c_pre_map.get(int(gidx), float("nan"))
        cpo = c_post_map.get(int(gidx), float("nan"))
        r_td = recomp_td.get(int(gidx), float("nan"))
        r_tp = recomp_tp.get(int(gidx), float("nan"))
        d_td = abs(r_td - s_td)
        d_tp = abs(r_tp - s_tp)
        max_delta_td = max(max_delta_td, d_td)
        max_delta_tp = max(max_delta_tp, d_tp)
        loc  = sec_name.get(sec_type_map.get(int(gidx), -1), "?")
        print(f"{gidx:>12}  {loc:>6}  "
              f"{cp:>10.6f}  {cpo:>10.6f}  "
              f"{s_td:>10.6f}  {r_td:>10.6f}  {d_td:>8.2e}  "
              f"{s_tp:>10.6f}  {r_tp:>10.6f}  {d_tp:>8.2e}")

    print("=" * len(hdr))
    print(f"\nMax |Δtheta_d| = {max_delta_td:.2e}   Max |Δtheta_p| = {max_delta_tp:.2e}")
    if max_delta_td < 1e-4 and max_delta_tp < 1e-4:
        print("✓  Stored thresholds match recomputed values within 1e-4")
    else:
        print("✗  MISMATCH detected — stored and recomputed values differ")


def _plot(global_idx, stored_td, stored_tp, recomp_td, recomp_tp,
          c_pre_map, c_post_map, pre_gid, post_gid, out_path):
    gidxs   = [int(g) for g in global_idx]
    s_td    = stored_td
    s_tp    = stored_tp
    r_td_v  = np.array([recomp_td.get(g, np.nan) for g in gidxs])
    r_tp_v  = np.array([recomp_tp.get(g, np.nan) for g in gidxs])
    cp_v    = np.array([c_pre_map.get(g, np.nan) for g in gidxs])
    cpo_v   = np.array([c_post_map.get(g, np.nan) for g in gidxs])
    x       = np.arange(len(gidxs))
    labels  = [str(g) for g in gidxs]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle(f"Threshold validation — pair {pre_gid}→{post_gid}", fontsize=13)

    # c_pre / c_post
    ax = axes[0, 0]
    ax.bar(x - 0.2, cp_v,  0.4, label="c_pre",  color="steelblue")
    ax.bar(x + 0.2, cpo_v, 0.4, label="c_post", color="tomato")
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Peak effcai_GB"); ax.set_title("c_pre / c_post per synapse")
    ax.legend(fontsize=9); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    # theta_d comparison
    ax = axes[0, 1]
    ax.bar(x - 0.2, s_td,   0.4, label="stored",    color="#66b3e6")
    ax.bar(x + 0.2, r_td_v, 0.4, label="recomputed", color="#ffaa55")
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("theta_d"); ax.set_title("theta_d: stored vs recomputed")
    ax.legend(fontsize=9); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    # theta_p comparison
    ax = axes[1, 0]
    ax.bar(x - 0.2, s_tp,   0.4, label="stored",    color="#66b3e6")
    ax.bar(x + 0.2, r_tp_v, 0.4, label="recomputed", color="#ffaa55")
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("theta_p"); ax.set_title("theta_p: stored vs recomputed")
    ax.legend(fontsize=9); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    # residuals
    ax = axes[1, 1]
    ax.plot(x, s_td - r_td_v, "o-", label="Δtheta_d (stored−recomp)", color="steelblue")
    ax.plot(x, s_tp - r_tp_v, "s-", label="Δtheta_p (stored−recomp)", color="tomato")
    ax.axhline(0, color="k", linestyle="--", alpha=0.5)
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Residual"); ax.set_title("Residuals")
    ax.legend(fontsize=9); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"Saved figure → {out_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pre",  type=int, default=None, help="Pre-synaptic GID")
    parser.add_argument("--post", type=int, default=None, help="Post-synaptic GID")
    parser.add_argument("--npz",  default=None, help="Direct path to .npz threshold file")
    parser.add_argument("--out",  default=None, help="Output PNG (default: threshold_validation_<pre>_<post>.png)")
    args = parser.parse_args()

    sys.path.insert(0, PLASTYFIRE_DIR)
    import bluecellulab
    bluecellulab.set_verbose(0)
    bluecellulab.neuron.load_mechanisms("/project/ctb-emuller/dhuruva/DEES_cell_packages/")
    bluecellulab.neuron.h.cvode.atolscale("v", 0.1)

    from plastyfire.simulator import c_pre_finder, c_post_finder

    # 1. Resolve npz and pair GIDs
    npz_path = _pick_npz(args.pre, args.post, args.npz)
    basename = os.path.splitext(os.path.basename(npz_path))[0]  # e.g. "180164_197248"
    pre_gid, post_gid = [int(x) for x in basename.split("_")]
    logger.info("Validating pair %d → %d  (npz: %s)", pre_gid, post_gid, npz_path)

    # 2. Load stored thresholds
    global_idx, stored_td, stored_tp = _load_npz(npz_path)
    logger.info("Loaded %d synapse entries from npz", len(global_idx))

    # 3. Read section types from edges file
    sec_type_map = _get_section_types(pre_gid, post_gid)
    logger.info("Found %d EXC synapses in edges file", len(sec_type_map))

    # 4. Spike amplitude for c_post calibration
    amp, width = _get_spike_amp(post_gid)
    stimulus = {"nspikes": 1, "freq": 0.1, "width": width, "offset": 1000.0, "amp": amp}
    logger.info("Post-cell spike stimulus: amp=%.4f nA, width=%.1f ms", amp, width)

    # 5. Build minimal sim config in a temp directory
    with tempfile.TemporaryDirectory() as tmp_dir:
        sim_config = _make_sim_config(tmp_dir)
        logger.info("Sim config written to %s", sim_config)

        fit_params = {"tau_effca_GB_GluSynapse": TAU_EFFCA}

        # 6. Run c_pre_finder (syn_extra_params=None: same overrides as compute_thresholds.py)
        logger.info("Running c_pre_finder ...")
        c_pre_map = c_pre_finder(
            sim_config, fit_params, None,
            pre_gid, post_gid,
            node_pop=NODE_POP, edge_pop=EDGE_POP, fixhp=True,
        )
        logger.info("c_pre done: %d synapses  mean=%.6f",
                    len(c_pre_map), np.mean(list(c_pre_map.values())))

        # 7. Run c_post_finder
        logger.info("Running c_post_finder ...")
        c_post_map = c_post_finder(
            sim_config, fit_params, None,
            pre_gid, post_gid, stimulus,
            node_pop=NODE_POP, edge_pop=EDGE_POP, fixhp=True,
        )
        if c_post_map is None:
            logger.error("c_post_finder returned None — cell failed to spike. "
                         "Check spike amplitude.")
            sys.exit(1)
        logger.info("c_post done: %d synapses  mean=%.6f",
                    len(c_post_map), np.mean(list(c_post_map.values())))

    # 8. Recompute thresholds
    recomp_td, recomp_tp = _compute_thresholds(c_pre_map, c_post_map, sec_type_map)

    # 9. Print comparison table
    _print_table(global_idx, stored_td, stored_tp,
                 c_pre_map, c_post_map,
                 recomp_td, recomp_tp, sec_type_map)

    # 10. Plot
    out_path = args.out or f"threshold_validation_{pre_gid}_{post_gid}.png"
    _plot(global_idx, stored_td, stored_tp, recomp_td, recomp_tp,
          c_pre_map, c_post_map, pre_gid, post_gid, out_path)


if __name__ == "__main__":
    main()
