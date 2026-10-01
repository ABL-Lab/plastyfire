#!/bin/bash
#SBATCH --job-name=bd_screen
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=1600M   # measured MaxRSS 1.28 GB
#SBATCH --time=02:00:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/bd_screen_%A_%a.out
# delta-bd screen (screen_bd.py). env: CANDS (json), PAIRS (comma list) or PAIRFILE + array index, TAG, NOPHARM
set -euo pipefail
source .venv/bin/activate
if [ -n "${TASKFILE:-}" ]; then read PAIRS CANDS < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $TASKFILE); fi
python -u glusynapse_v2/emodel_bd/screen_bd.py --pairs $PAIRS --cands $CANDS ${NOPHARM:+--no-pharm} \
    --out glusynapse_v2/emodel_bd/results/${TAG}_${SLURM_ARRAY_TASK_ID:-0}.csv
