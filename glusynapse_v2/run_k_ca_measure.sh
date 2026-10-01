#!/bin/bash
#SBATCH --job-name=k_ca_measure
#SBATCH --account=def-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/k_ca_measure_%j.out
# K_ca for ljp 0 and ljp 25 (sj07 post-only traces). Sizing: no measurement yet (slurmdb down): one pair's 7 synapses x 20M float64
# MEASURED 22096823: ljp0 K 3.1344e-8 (pooled) / 3.1365e-8 (per-synapse); ljp25 K 3.4715e-7 / 3.3239e-7; seff unavailable (slurmdb down).
# = ~1.1 GB plus the pickle -> 6G; 48 pickles loaded one at a time. Correct from seff.
export OPENBLAS_NUM_THREADS=1
source .venv/bin/activate
python -u glusynapse_v2/k_ca_measure.py
