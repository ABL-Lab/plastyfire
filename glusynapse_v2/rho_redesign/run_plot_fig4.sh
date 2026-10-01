#!/bin/bash
# Sizing basis: comparable plot job 22108767 (plot_diag): MaxRSS 247 MB, 0:07 elapsed.
# 247 MB + 25% = 309 MB -> floor 2G; time 0:15 (min); 1 CPU (single-threaded matplotlib).
#SBATCH --job-name=plot_fig4
#SBATCH --account=def-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/plot_fig4_%j.out
set -euo pipefail
source .venv/bin/activate
python -u glusynapse_v2/rho_redesign/plot_fig4_l5_fits.py
