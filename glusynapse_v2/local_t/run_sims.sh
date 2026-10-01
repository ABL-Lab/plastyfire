#!/bin/bash
#SBATCH --job-name=local_t
#SBATCH --account=def-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=1100M   # pilot 22072098: MaxRSS 0.73 GB (diag_burst 22066821: 0.43-0.79 GB) + 25%
#SBATCH --time=00:15:00   # pilot 22072098: 5:40 for 18 sims + 50%, next 15 min
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/local_t_%A_%a.out
# measured 22072313 (24 pairs): elapsed 5:30-11:25, MaxRSS 0.44-0.89 GB, CPU eff ~98%; total 2.4 CPU·h
# one task = one pair (line SLURM_ARRAY_TASK_ID+1 of PAIRFILE), variants og,s5, all protocols
set -euo pipefail
source .venv/bin/activate
PAIR=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" ${PAIRFILE:-glusynapse_v2/local_t/pairs.txt})
mkdir -p glusynapse_v2/local_t/out
python -u glusynapse_v2/local_t/sims.py --pair $PAIR --out glusynapse_v2/local_t/out/${PAIR}.npz
