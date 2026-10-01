#!/bin/bash
#SBATCH --job-name=diag_l23v5
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/diag_l23v5_%j.out
# L23L23_DIAG.md section 7: per-record CPU replay of the v5c rule (gpu_v5_rho mode 2 / 3) at the params of a v5 fit json
# + 13 fixed-parameter variants, one pathway model per job. CPU/numba, single-threaded.
# usage: FIT=<json> MODEL=l23l23|l5|l23 sbatch --export=ALL,FIT=...,MODEL=... glusynapse_v2/rho_redesign/run_diag_l23l23_v5.sh
# Sizing basis: diag_l23l23 run 2 22135934/5/6 (same loader, 11 variants, l23l23 / l5 / l23): 3:49 / 1:59 / 2:39,
# MaxRSS ~2.09 GB at 2G (page-cache ceiling, completed). Here 13 lighter variants + a per-synapse Ca inversion
# (one trace at a time): same 2G, 0:15 (3:49 x ~1.3 x 1.5 = 7.5 min).
set -euo pipefail
source .venv/bin/activate
export NUMBA_NUM_THREADS=1
R=glusynapse_v2/rho_redesign
FIT=${FIT:-$R/results/v5_V5c_s7.json}
MODEL=${MODEL:-l23l23}
OUT=${OUT:-/scratch/dhuruva/l23l23_diag_v5}
python -m py_compile $R/diag_l23l23_v5.py
python -u $R/diag_l23l23_v5.py --fit $FIT --model $MODEL --out $OUT
