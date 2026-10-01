#!/bin/bash
#SBATCH --job-name=v2_jax_td4
#SBATCH --account=def-emuller
#SBATCH --gpus-per-node=nvidia_h100_80gb_hbm3_3g.40gb:1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/v2_jax_td4_%j.out
# t_drive 4 GPU==CPU check on the 2 pilot pairs (vca dirs). Sizing: pilot extract npz are 2/24 of the subset (GPU fit 24 GB),
# Measured 22073582: 3.2 GB MaxRSS, 17 s, 82% CPU -> 4G. PASS (GPU == CPU to 4e-16).
# With no_block lanes, 22088066: 3.75 GB, 51 s -> 5G. PASS.
set -euo pipefail
source .venv/bin/activate
python -u glusynapse_v2/tests/test_jax_tdrive4.py
