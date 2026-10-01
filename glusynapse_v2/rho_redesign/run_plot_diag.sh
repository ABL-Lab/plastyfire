#!/bin/bash
# seff 22108767: 247 MB, 0:07
#SBATCH --job-name=plot_diag
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/plot_diag_%j.out
# Figures/tables from diag_rho outputs (csv ~6k rows + traces npz). Pilot size (no comparable run): 1 CPU, 4G, 0:15.
set -euo pipefail
source .venv/bin/activate
python -u glusynapse_v2/rho_redesign/plot_diag.py --tag ${TAG:-td4_s2}
