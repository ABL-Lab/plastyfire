#!/bin/bash
# Geometric Letzkus split csv (letzkus_geo.py); pandas only. Measured analogue 22136361 (letzkus_loc): 50 s, 1.99 GB -> 1 CPU, 2500M, 0:15.
#SBATCH --job-name=letzkus_geo
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=2500M
#SBATCH --time=00:15:00
#SBATCH --output=/scratch/dhuruva/letzkus_loc/letzkus_geo_%j.out
ROOT=/lustre09/project/6070394/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
python -u glusynapse_v2/rho_redesign/letzkus_geo.py /scratch/dhuruva/letzkus_loc/ours_synapses.csv ebner/pair_geometry_L23PC_L5TTPC.csv ebner/pair_geometry_L23PC_L5TTPC_geo.csv
