#!/bin/bash
# seff 22108766 (gate smoke): 435 MB, 0:19 -> 1G, 0:15 (rates checks add a few kernel compiles)
# SMOKE_ARGS=--real (C1/C2 on A0g_s3, GPU vs CPU scan, real L5 data): sbatch --mem=33G --time=00:15:00. Basis: synthetic 22114130
# 594 MB, 0:27; the L5 GPU build = v4 1g fits 22-26.4 GB MaxRSS (S already holds -ica_VDCC: no extra trace), CPU scan of L5
# alone after it is freed <= 19.0 GB (22120612) -> 26.4 x 1.25 = 33G; ~2 min load + ~2 min CPU + compiles ~6 min x 1.5 -> 0:15.
#SBATCH --job-name=v4_smoke
#SBATCH --account=def-emuller
#SBATCH --gpus-per-node=nvidia_h100_80gb_hbm3_1g.10gb:1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/v4_smoke_%j.out
# v4 post-rule prototype: thetas, fast kernel vs CPU reference and vs gpu_v3, readout (synthetic, no traces). Sized as run_smoke_gpu_v3.sh (1.27 GB, 0:47 -> 4G floor).
set -euo pipefail
source glusynapse_v2/env_v3.sh
nvidia-smi -L
python -m py_compile glusynapse_v2/rho_redesign/{fit_v4,eval_v4,rho_v4,rho_v4d,gpu_v4_rho}.py && echo "py_compile ok"
python -u glusynapse_v2/rho_redesign/${SMOKE:-smoke_gpu_v4.py} ${SMOKE_ARGS:-}
