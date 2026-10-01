#!/bin/bash
#SBATCH --job-name=basal_nakv
#SBATCH --account=def-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=1600M   # diag_burst MaxRSS 0.4-1.3 GB
#SBATCH --time=01:30:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/basal_nakv_%A_%a.out
# one task = one pair x one variant: TASKFILE line = "PAIR VARIANT"
set -euo pipefail
source .venv/bin/activate
read PAIR VAR < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $TASKFILE)
mkdir -p glusynapse_v2/cell_audit/basal_nakv/out
python -u glusynapse_v2/cell_audit/basal_nakv/grid.py --pair $PAIR --variants $VAR ${PROTOS:+--protos $PROTOS} \
    --out glusynapse_v2/cell_audit/basal_nakv/out/${PAIR}_${VAR}.json
