#!/bin/bash
#SBATCH --job-name=diag_ltd
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/diag_ltd_%j.out
# LTD_DIAG.md: per-synapse LTD mechanism + 34 fixed-parameter trigger / veto / d_min / theta_d variants of the v5c rule
# (diag_ltd.py), one pathway model per job. CPU/numba, single-threaded.
# usage: FIT=<json> MODEL=l5|l23|l23l23 sbatch --export=ALL,FIT=...,MODEL=... glusynapse_v2/rho_redesign/run_diag_ltd.sh
# MEASURED round 1 (19 variants, 3G): 22155436 l5 1:26, 22155437 l23 2:08, 22155438 l23l23 4:17, MaxRSS 2.99 GB each
# (at the 3G limit). Round 2 (34 variants + veto queue): 4G, 0:15 (4:17 x 34/19 x 1.5 = 11.5 min).
# MEASURED round 2 (4G): 22156062 l5 1:45 (CPU 79 %), 22156064 l23l23 3:58, MaxRSS 3.99 GB = at the 4G limit again
# (page cache counted; completed) -> next run --mem=5G.
set -euo pipefail
source .venv/bin/activate
export NUMBA_NUM_THREADS=1
R=glusynapse_v2/rho_redesign
FIT=${FIT:-$R/results/v5_S1_V5c_s6_39.json}
MODEL=${MODEL:-l5}
OUT=${OUT:-/scratch/dhuruva/ltd_diag}
python -m py_compile $R/diag_ltd.py
python -u $R/diag_ltd.py --fit $FIT --model $MODEL --out $OUT
