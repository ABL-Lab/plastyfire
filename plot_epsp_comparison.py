#!/usr/bin/env python3

import argparse
import hashlib
import logging
import os
import pickle
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "plastyfire"))

from plastyfire.config import OptConfig
from plastyfire.ephysutils import Experiment

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

CSVF_NAME = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/biodata/paired_recordings.csv"
CONFIGS_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/configs"
BESTSOL_PKL = "/lustre06/project/6077694/dhuruva/plastyfire/bestsol.pkl"
CHECKPOINT_PKL = (
    "/lustre06/project/6077694/dhuruva/plastyfire/plastyfire/checkpoint_full.pkl"
)

FIT_PARAM_NAMES = [
    "gamma_d_GB_GluSynapse",
    "gamma_p_GB_GluSynapse",
    "a00",
    "a01",
    "a10",
    "a11",
    "a20",
    "a21",
    "a30",
    "a31",
]

MIN2MS = 60 * 1000.0

PROTOCOLS = ["mrk97_01", "mrk97_02", "mrk97_07", "mrk97_08", "sjh06_02"]

PROTOCOL_ORDER = ["sjh06_02", "mrk97_08", "mrk97_07", "mrk97_02", "mrk97_01"]


def load_best_solution():
    if os.path.exists(BESTSOL_PKL):
        logger.info(f"Loading best solution from {BESTSOL_PKL}")
        with open(BESTSOL_PKL, "rb") as f:
            best_params = pickle.load(f)
        return best_params, None
    elif os.path.exists(CHECKPOINT_PKL):
        logger.warning(
            f"bestsol.pkl not found, loading from checkpoint: {CHECKPOINT_PKL}"
        )
        with open(CHECKPOINT_PKL, "rb") as f:
            checkpoint = pickle.load(f)
        if "halloffame" in checkpoint and len(checkpoint["halloffame"]) > 0:
            hof = checkpoint["halloffame"]
            best_ind = hof[0]
            best_params = dict(zip(FIT_PARAM_NAMES, best_ind))
            fitness_scores = best_ind.fitness.values if hasattr(best_ind, 'fitness') else None
            logger.info(f"Extracted best individual from checkpoint: {best_params}")
            if fitness_scores:
                logger.info(f"Fitness scores: {fitness_scores}")
            return best_params, fitness_scores
        else:
            raise ValueError("No hall of fame found in checkpoint")
    else:
        raise FileNotFoundError(f"Neither {BESTSOL_PKL} nor {CHECKPOINT_PKL} found")


def compute_epsp_ratios(sim_dict, param_hash, config):
    workdir = sim_dict["workdir"]

    pkl_file = os.path.join(workdir, f"simulation_{param_hash}.pkl")
    if not os.path.exists(pkl_file):
        return None, "file_not_found"

    try:
        with open(pkl_file, "rb") as f:
            raw_results = pickle.load(f)

        if "pre_spikes" in raw_results and "prespikes" not in raw_results:
            raw_results["prespikes"] = raw_results["pre_spikes"]

        exp_handler = Experiment(
            raw_results,
            c01duration=config.C01_duration,
            c02duration=config.C02_duration,
            period=sim_dict["period"],
        )

        nepsp = int(config.C01_duration * MIN2MS / config.T)
        epsp_ratio = exp_handler.compute_epsp_ratio(nepsp)
        return epsp_ratio, None
    except RuntimeError as e:
        if "Postsynaptic cell spiking" in str(e):
            return None, "postsynaptic_spiking"
        return None, "runtime_error"
    except Exception as e:
        logger.warning(f"Error computing EPSP ratio: {e}")
        return None, "other_error"


def get_protocol_label(protocol_id, invitro_db):
    row = invitro_db[invitro_db["protocol_id"] == protocol_id].iloc[0]
    pre_mtype = row["pre_mtype"]
    post_mtype = row["post_mtype"]

    mtype_map = {
        "L23PC": "L2/3",
        "L5TTPC": "L5",
        "L4PC": "L4",
        "L4SS": "L4",
    }

    pre_label = mtype_map.get(pre_mtype, pre_mtype)
    post_label = mtype_map.get(post_mtype, post_mtype)

    freq = int(row["frequency_train"])
    dt = int(row["dt_train"])

    return f"{pre_label} to {post_label}, {freq} Hz, {dt} ms"


def main():
    parser = argparse.ArgumentParser(description="Plot EPSP ratio comparison")
    parser.add_argument("--seed", type=int, default=1234, help="Random seed")
    parser.add_argument(
        "--sample-size", type=int, default=100, help="Number of pairs to sample"
    )
    args = parser.parse_args()

    best_params, fitness_scores = load_best_solution()
    logger.info(f"Best parameters: {best_params}")

    param_values = [best_params[name] for name in FIT_PARAM_NAMES]
    param_hash = hashlib.md5(str(param_values).encode()).hexdigest()[:12]
    logger.info(f"Parameter hash: {param_hash}")

    # Load in vitro database
    invitro_db = pd.read_csv(CSVF_NAME)
    invitro_db = invitro_db.loc[invitro_db["protocol_id"].isin(PROTOCOLS)]

    # Try to load results from cache first
    cache_key = hashlib.md5(str(param_values).encode()).hexdigest()
    cache_file = os.path.join("/lustre06/project/6077694/dhuruva/plastyfire/plastyfire/.cache", f"{cache_key}.pkl")

    if os.path.exists(cache_file):
        logger.info(f"Loading EPSP ratios from cache: {cache_file}")
        with open(cache_file, 'rb') as f:
            cache_data = pickle.load(f)

        if 'resdb' in cache_data and isinstance(cache_data['resdb'], pd.DataFrame):
            results_df = cache_data['resdb']
            logger.info(f"Loaded {len(results_df)} EPSP ratios from cache")
        else:
            logger.error("Cache file does not contain expected 'resdb' DataFrame")
            return
    else:
        logger.warning(f"Cache file not found: {cache_file}")
        logger.info("Falling back to computing EPSP ratios from simulation files...")

        all_sims = []
        for elem in invitro_db.itertuples():
            config = OptConfig(
                os.path.join(CONFIGS_DIR, f"{elem.pre_mtype}_{elem.post_mtype}.yaml")
            )

            fastforward = config.fastforward
            if fastforward is None:
                fastforward = config.C01_duration * MIN2MS + config.nreps * config.T
            fastforward = float(fastforward)

            index_file = os.path.join(
                os.path.split(os.path.split(config.out_dir)[0])[0],
                f"index_{elem.pre_mtype}_{elem.post_mtype}.csv",
            )

            if not os.path.exists(index_file):
                logger.error(f"Index file not found: {index_file}")
                continue

            sim_idx = pd.read_csv(index_file)
            sim_idx.set_index(["frequency", "dt"], inplace=True)
            sim_idx.sort_index(inplace=True)

            paths = sim_idx.loc[elem.frequency_train, elem.dt_train]["path"]
            if isinstance(paths, str):
                paths = pd.Series([paths])

            np.random.seed(args.seed)
            sample_size = min(args.sample_size, len(paths))
            paths = paths.sample(sample_size, random_state=args.seed)

            logger.info(f"Protocol {elem.protocol_id}: Found {len(paths)} simulation paths")

            for path in paths:
                all_sims.append(
                    {
                        "protocol_id": elem.protocol_id,
                        "fastforward": fastforward,
                        "workdir": "/".join(path.split("/")[:-1]),
                        "config": config,
                        "period": elem.period_sweep,
                    }
                )

        logger.info(f"Total simulations: {len(all_sims)}")

        epsp_results = []
        skip_stats = {
            "postsynaptic_spiking": 0,
            "file_not_found": 0,
            "runtime_error": 0,
            "other_error": 0,
        }

        for sim_dict in all_sims:
            epsp_ratio, skip_reason = compute_epsp_ratios(
                sim_dict, param_hash, sim_dict["config"]
            )
            if epsp_ratio is not None:
                epsp_results.append(
                    {
                        "protocol_id": sim_dict["protocol_id"],
                        "epsp_ratio": epsp_ratio,
                    }
                )
            elif skip_reason:
                skip_stats[skip_reason] = skip_stats.get(skip_reason, 0) + 1

        results_df = pd.DataFrame(epsp_results)
        logger.info(
            f"Computed {len(results_df)} EPSP ratios out of {len(all_sims)} simulations"
        )
        if sum(skip_stats.values()) > 0:
            logger.info(f"Skipped simulations: {skip_stats}")

    insilico_stats = results_df.groupby("protocol_id")["epsp_ratio"].agg(
        ["mean", "sem", "std", "count"]
    )
    insilico_stats = pd.DataFrame(insilico_stats)
    insilico_stats.columns = [
        "mean_epsp_ratio_sim",
        "sem_epsp_ratio_sim",
        "std_epsp_ratio_sim",
        "count",
    ]

    comparison = pd.merge(
        invitro_db[
            [
                "protocol_id",
                "mean_epsp_ratio",
                "sem_epsp_ratio",
                "frequency_train",
                "dt_train",
                "pre_mtype",
                "post_mtype",
            ]
        ],
        insilico_stats,
        on="protocol_id",
    )

    logger.info(f"\nComparison:\n{comparison}")

    fig, ax = plt.subplots(figsize=(6, 8))

    y_positions = {}
    y_pos = 0
    for protocol_id in PROTOCOL_ORDER:
        if protocol_id in comparison["protocol_id"].values:
            y_positions[protocol_id] = y_pos
            y_pos += 1

    y_labels = []
    y_insilico = []
    y_invitro = []
    x_insilico = []
    x_invitro = []
    err_insilico = []
    err_invitro = []

    for protocol_id in PROTOCOL_ORDER:
        if protocol_id not in y_positions:
            continue

        row = comparison[comparison["protocol_id"] == protocol_id].iloc[0]
        label = get_protocol_label(protocol_id, invitro_db)
        y_labels.append(label)

        y_pos = y_positions[protocol_id]

        y_insilico.append(y_pos)
        x_insilico.append(row["mean_epsp_ratio_sim"])
        err_insilico.append(row["sem_epsp_ratio_sim"])

        y_invitro.append(y_pos)
        x_invitro.append(row["mean_epsp_ratio"])
        err_invitro.append(row["sem_epsp_ratio"])

    ax.errorbar(
        x_insilico,
        y_insilico,
        xerr=err_insilico,
        fmt="o",
        color="#1f77b4",
        label="in silico",
        capsize=5,
        capthick=2,
        markersize=8,
        elinewidth=2,
    )

    ax.errorbar(
        x_invitro,
        y_invitro,
        xerr=err_invitro,
        fmt="o",
        color="#ff7f0e",
        label="in vitro",
        capsize=5,
        capthick=2,
        markersize=8,
        elinewidth=2,
    )

    ax.axvline(1.0, color="gray", linestyle="--", linewidth=1, alpha=0.7)

    ax.set_yticks(list(range(len(y_labels))))
    ax.set_yticklabels(y_labels)
    ax.set_xlabel("EPSP ratio", fontsize=12)
    ax.set_xlim(0.75, 1.25)
    ax.set_ylim(-0.5, len(y_labels) - 0.5)
    ax.legend(loc="best", fontsize=10)
    ax.grid(True, alpha=0.3, axis="x")

    ax.text(
        0.02,
        0.98,
        "C",
        transform=ax.transAxes,
        fontsize=16,
        fontweight="bold",
        va="top",
        ha="left",
    )

    if fitness_scores is not None:
        fitness_text = "Fitness scores:\n"
        # Map fitness scores to protocol labels using PROTOCOLS order (not PROTOCOL_ORDER)
        # fitness_scores correspond to PROTOCOLS order from the optimization
        for i, protocol_id in enumerate(PROTOCOLS):
            if i < len(fitness_scores) and protocol_id in comparison["protocol_id"].values:
                label = get_protocol_label(protocol_id, invitro_db)
                fitness_text += f"  {label}: {fitness_scores[i]:.6f}\n"

        ax.text(
            0.02,
            0.02,
            fitness_text.strip(),
            transform=ax.transAxes,
            fontsize=8,
            va="bottom",
            ha="left",
            bbox=dict(boxstyle="round,pad=0.5", facecolor="wheat", alpha=0.8),
        )

    plt.tight_layout()

    output_file = "epsp_ratio_comparison.png"
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    logger.info(f"Plot saved to {output_file}")
    plt.close()


if __name__ == "__main__":
    main()
