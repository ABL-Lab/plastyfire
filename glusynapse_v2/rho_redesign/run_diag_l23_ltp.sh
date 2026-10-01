#!/bin/bash
#SBATCH --job-name=diag_l23_ltp
#SBATCH --account=def-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/diag_l23_ltp_%j.out
# Round-2 diagnosis (ROUND2.md): L2/3->L5 +10 over-potentiation vs L5 +10 LTP under A0g_s3. CPU/numba, no jax.
# usage: sbatch glusynapse_v2/rho_redesign/run_diag_l23_ltp.sh
# Sizing basis: diag_distal 22115348 (223 L2/3 records + vdcc + shaft_cai): 15.74 GB MaxRSS, 8:12, 95% CPU.
# This job: 276 L2/3 records (110 1AP+10, 55 3AP+10, 111 S&H 50 Hz+10) + 120 L5 records (5 protos x 24) = 396
# (L5 batch freed before L2/3 loads). mem = 15.74 x 396/223 x 1.25 = 34.9 -> 35G; time = 8:12 x 396/223 x 1.5 = 21.8
# -> 0:30; 1 CPU (single-threaded numba).
# MEASURED: 22120179 15.86 GB MaxRSS, 2:28 elapsed, 66% CPU (35G = 45% used) -> next run --mem=20G --time=0:15.
set -euo pipefail
source .venv/bin/activate
export NUMBA_NUM_THREADS=1
V2=glusynapse_v2; R=$V2/rho_redesign
FIT=${FIT:-$R/results/v4_A0g_s3.json}
python -u $R/diag_l23_ltp.py --fit $FIT --save $R/results/l23ltp_A0g_s3 --fig $R/figs/fig6_l23_ltp.png
