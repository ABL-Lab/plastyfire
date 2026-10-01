#!/bin/bash
# Summary of the ljp25 validation (reads 8 small jsons). No measurement yet: 1G / 0:15 is the minimum; resize after seff. Measured 22090852: 2 s, 44 MB -> 256M.
#SBATCH --job-name=ljp25_sum
#SBATCH --account=def-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=256M
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/glusynapse_v2/cell_audit/ljp25_validation/summary_%j.txt
set -euo pipefail
source .venv/bin/activate
python -u glusynapse_v2/cell_audit/ljp25_validation/summarize.py ljp25
