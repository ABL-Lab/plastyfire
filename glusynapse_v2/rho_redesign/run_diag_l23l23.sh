#!/bin/bash
#SBATCH --job-name=diag_l23l23
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/diag_l23l23_%j.out
# L23L23_DIAG.md: per-record CPU replay of the C1Ajn_s5 rule (fit json v4_C1Ajz_pilot = C1Ajn_s5 params, 3-pathway
# score) + parameter-free rule variants, one pathway model per job. CPU/numba, no jax, single-threaded.
# usage: MODEL=l23l23|l5|l23 sbatch --export=ALL,MODEL=l23l23 glusynapse_v2/rho_redesign/run_diag_l23l23.sh
# Sizing basis: diag_dltd 22132795 (same per-record replay, one BatchV2 per npz, 48 records of ~14 MB): 910 MB MaxRSS,
# 0:12, 100% CPU. Memory is per record (largest npz here 22 MB vs 14 MB) + 9-variant state (KB) + per-synapse rows
# (<= 6.3k x 60) -> ~1.2 GB x 1.25 -> 1600M. Time: <= 0.37 s/record (diag_l23_ltp 22120179: 396 records 2:28) + replay
# ~15 ms/synapse/9 variants: l23l23 1527 rec 6229 syn ~11 min, l5 502 rec ~4 min, l23 504 rec ~5 min; x1.5 -> 0:30 for
# l23l23 (l5 / l23 submitted with --time=0:15).
# MEASURED 22135726/7/8 (l23l23 / l5 / l23, 9 variants; crashed in the chi2 print after the targets table): 4:12 / 2:00 /
# 2:29, MaxRSS 1.63 GB = at the 1600M limit -> 2G; 11 variants ~ +20% time -> 0:15 for all three. Fill seff here.
set -euo pipefail
source .venv/bin/activate
export NUMBA_NUM_THREADS=1
V2=glusynapse_v2; R=$V2/rho_redesign
FIT=${FIT:-$R/results/v4_C1Ajz_pilot.json}
MODEL=${MODEL:-l23l23}
OUT=${OUT:-/scratch/dhuruva/l23l23_diag}
python -m py_compile $R/diag_l23l23.py
python -u $R/diag_l23l23.py --fit $FIT --model $MODEL --out $OUT
