#!/bin/bash
# L2/3 -> L5 synapse location table vs Letzkus 2006 (letzkus_loc.py). Single-threaded, read-only on the circuit.
# Pilot sizing (no comparable job; geometry part of ebner/run_pair_geometry.sh without NEURON): 1 CPU, 2G, 0:15.
#SBATCH --job-name=letzkus_loc
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=/scratch/dhuruva/letzkus_loc/letzkus_loc_%j.out
ROOT=/lustre09/project/6070394/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
python -u glusynapse_v2/rho_redesign/letzkus_loc.py /scratch/dhuruva/letzkus_loc
