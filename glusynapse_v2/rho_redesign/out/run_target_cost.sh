#!/bin/bash
#SBATCH --job-name=tgt_cost
#SBATCH --account=def-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/tgt_cost_%j.out
# Sizing basis: plot job 22114017 (213 MB MaxRSS, 3 s); this reads npz headers only.
# MEASURED 22115351: 1:13 wall, 1.99 GB MaxRSS of 2G (at the limit), 19% CPU -> next run 3G, 0:15.
set -euo pipefail
source glusynapse_v2/env_v3.sh
python -u /lustre09/project/6070394/dhuruva/plastyfire/glusynapse_v2/rho_redesign/out/target_cost.py > /lustre09/project/6070394/dhuruva/plastyfire/glusynapse_v2/rho_redesign/out/target_cost.csv
