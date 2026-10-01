#!/bin/bash
#SBATCH --job-name=v2_explore
#SBATCH --account=def-emuller      # rrg-emuller has no GPU allocation
#SBATCH --gpus-per-node=nvidia_h100_80gb_hbm3_3g.40gb:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=35G          # preview data, measured 21 GiB for the preview fits
#SBATCH --time=00:45:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/v2_explore_%j.out
# Score hand-made variants on the preview (explore_v2.py). env: BASE (fit json), VARIANTS (json file), SHOW.
set -euo pipefail
source .venv/bin/activate
python -u glusynapse_v2/explore_v2.py --base $BASE --variants @$VARIANTS --show "${SHOW:-}"
