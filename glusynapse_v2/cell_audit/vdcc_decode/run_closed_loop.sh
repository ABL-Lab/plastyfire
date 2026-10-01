#!/bin/bash
#SBATCH --job-name=vdcc_cl
#SBATCH --account=def-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=1100M       # local_t 22072313 (same sims, 2 variants x 9 protos): MaxRSS 0.44-0.89 GB + 25%
#SBATCH --time=00:15:00   # 22072313: 5:30-11:25 for 18 sims; here 4 sims -> <= 3 min, 15 min is the floor
#SBATCH --array=0-1
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/vdcc_cl_%A_%a.out
set -euo pipefail
source .venv/bin/activate
D=glusynapse_v2/cell_audit/vdcc_decode
PAIR=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $D/pairs_cl.txt)
python -u $D/closed_loop.py --pair $PAIR --variants s25 --protos ap1,burst,epsp,sj20@-10 --out $D/out/${PAIR}_s25.npz
