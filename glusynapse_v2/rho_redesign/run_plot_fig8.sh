#!/bin/bash
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --mem=2G
#SBATCH --job-name=plot_fig8
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/plot_fig8_%j.out
# No measurement yet (small pandas/matplotlib job); check seff after run.
source /project/rrg-emuller/dhuruva/plastyfire/.venv/bin/activate
python /project/rrg-emuller/dhuruva/plastyfire/glusynapse_v2/rho_redesign/plot_fig8_joint_fit.py
