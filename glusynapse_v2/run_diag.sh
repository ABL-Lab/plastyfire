#!/bin/bash
#SBATCH --job-name=v2_diag
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=01:00:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/v2_diag_%j.out
source .venv/bin/activate
python -u glusynapse_v2/diag_no.py
