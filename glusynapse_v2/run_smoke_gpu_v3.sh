#!/bin/bash
#SBATCH --job-name=v3_smoke
#SBATCH --account=def-emuller
#SBATCH --gpus-per-node=nvidia_h100_80gb_hbm3_1g.10gb:1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/v3_smoke_%j.out
# gpu_v3 kernel vs jax_v2 scan on synthetic data (all modes). No trace files: jax + numba + CUDA contexts only (pilot size 4G).
# 22107147: PASS (max abs 4.4e-16, all 7 modes; jax + numba-cuda in one process ok). 1.27 GB MaxRSS, 0:47 -> 4G is the floor Slurm needs; keep.
set -euo pipefail
source glusynapse_v2/env_v3.sh
nvidia-smi -L
python -u glusynapse_v2/tests/smoke_gpu_v3.py
