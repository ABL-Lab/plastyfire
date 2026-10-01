#!/bin/bash
#SBATCH --job-name=v22r_evalall
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G         # measured MaxRSS 160 GiB (21941602, GPU-fit json), +25%
#SBATCH --time=04:00:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/v22r_evalall_%j.out
# Score a fit (env FIT = results basename, default reduced_subset) once on all 120 pairs (4054 records).
set -euo pipefail
source .venv/bin/activate
V2=glusynapse_v2
python -u $V2/eval_v2.py --fit $V2/results/${FIT:-reduced_subset}.json \
    --dirs $V2/extracted/ebner_delta-prefire,$V2/extracted/markram_delta-prefire-tr --save $V2/results/${FIT:-reduced_subset}_allpairs ${EXTRA:-}
