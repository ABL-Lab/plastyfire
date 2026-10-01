#!/bin/bash
#SBATCH --job-name=v2_unit_td4
#SBATCH --account=def-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/v2_unit_td4_%j.out
# t_drive 4 unit tests + t_drive 0-3 regression, one name per call (full suite OOMs on the login node).
# Measured 22091011: MaxRSS 248 MB, 7 s, 86% CPU on 1 CPU with 4G -> 1G (+25%: 310 MB, floor 1G).
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
source .venv/bin/activate
for n in t_drive4 t_drive3 t_drive_ no_block; do python glusynapse_v2/tests/test_units.py $n; done
