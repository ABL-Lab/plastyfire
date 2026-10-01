#!/bin/bash
#SBATCH --job-name=replot_fig7
#SBATCH --account=rrg-emuller
# Redraw fig7 (spine VDCC Ca labels) from results/vgate_amp_A0g_s3*.csv. CSV read + matplotlib only: 1 CPU, 2G, 0:15 (no comparable job; pilot).
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/replot_fig7_%j.out
set -euo pipefail
source .venv/bin/activate
R=glusynapse_v2/rho_redesign
python -u $R/replot_fig7.py $R/results/vgate_amp_A0g_s3 $R/figs/fig7_vgate_amp.png
cp $R/figs/fig7_vgate_amp.png /lustre09/project/6070394/dhuruva/Dhuruva_DEES_Plasticity_Copy/figures/glusyn_v2/fig7_vgate_amp.png
