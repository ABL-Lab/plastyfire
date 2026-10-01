#!/bin/bash
#SBATCH --job-name=nmda_vdcc
#SBATCH --account=def-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=1600M   # diag_burst MaxRSS 1.0-1.3 GB
#SBATCH --time=02:00:00
#SBATCH --array=0-7
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/nmda_vdcc_%A_%a.out
set -euo pipefail
source .venv/bin/activate
D=glusynapse_v2/cell_audit/nmda_vdcc
PAIR=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $D/pairs.txt)
python -u glusynapse_v2/diag_burst.py --pair $PAIR --variants og_basal,s5,g2 --protos epsp,1ap,1ap@+10,1ap@-10 \
    --var-cands $D/var_cands.json --syn-cands $D/syn_cands.json \
    --out glusynapse_v2/results/diag_burst_${PAIR}_cal1.json
