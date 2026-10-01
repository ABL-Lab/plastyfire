#!/bin/bash
#SBATCH --job-name=eval_l23l5
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/eval_l23l5_%j.out
# Transfer test: score an L5-L5 fit on paired L2/3->L5 records (L23_PATHWAYS.md). No refit, CPU/numba, no jax.
# Sizing basis: v22r_evalall 21941602 = 4054 records, 160 GiB MaxRSS, 16.5 min, 4 CPU at ~19% eff (~1 core).
# Pilot 22068719 (td2_s1, 334 records): MaxRSS 15.9 GB, 1:50, 75% CPU eff -> 20G, 0:15.
# usage: FIT=<results prefix> SPLIT=rise|median sbatch glusynapse_v2/run_eval_l23l5.sh
set -euo pipefail
source .venv/bin/activate
export NUMBA_NUM_THREADS=1
ANALYTICAL_BASIS_DIR=basis_results_edges_ebner_l23l5_delta python glusynapse_v2/eval_l23l5.py \
    --fit glusynapse_v2/results/${FIT}.json --split ${SPLIT} \
    --save glusynapse_v2/results/l23l5_${FIT#reduced_gpu_subset_}_${SPLIT}
