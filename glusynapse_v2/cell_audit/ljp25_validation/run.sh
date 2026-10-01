#!/bin/bash
# Sabatini/Chindemi Fig 2 validation of the Ba->Ca VDCC correction (ljp_VDCC 25 mV, gca_bar x0.25/0.31/0.40 of 0.0744),
# same 8 pairs, protocols and diag_burst setup as cell_audit/nmda_vdcc (NMDA_VDCC_CALIBRATION.md s2), plus the unchanged
# delta cell as reference. Sizing: 22066821 (3 variants x 4 protos) max 5:02, MaxRSS 0.79 GB -> 4 variants ~6.7 min;
# 1000M, 0:15. Measured 22090528: 5:57-6:29, MaxRSS 0.94-0.99 GB (at the limit) -> 1250M.
#SBATCH --job-name=ljp25_val
#SBATCH --account=def-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=1250M
#SBATCH --time=00:15:00
#SBATCH --array=0-7
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/ljp25_val_%A_%a.out
set -euo pipefail
source .venv/bin/activate
D=glusynapse_v2/cell_audit/ljp25_validation
PAIR=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $D/pairs.txt)
python -u glusynapse_v2/diag_burst.py --pair $PAIR --variants delta,ljp25_g025,ljp25_g031,ljp25_g040 --protos epsp,1ap,1ap@+10,1ap@-10 \
    --syn-cands $D/syn_cands.json --out $D/out/diag_burst_${PAIR}_ljp25.json
