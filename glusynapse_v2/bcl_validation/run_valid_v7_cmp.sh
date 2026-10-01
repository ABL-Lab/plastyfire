#!/bin/bash
# Compare step of run_valid_v7.sh (submit afterok / afterany, same FIT / ECB_REF / TAG / STAGE env): compare_prefire_v7.py
# (offline = CPU port of the gpu_v7_rho kernel on the split1 records and bases) -> /scratch/dhuruva/bcl_valid_v7_${TAG}/
# cmp_${STAGE}_{records,targets}.csv and .png/.pdf. shard_*.jsonl in that dir are merged first.
# Sizing basis: V5c compare pilot 22144695 MEASURED 0:20, 0.93 GB -> 1 CPU, 2G (pilot: 0.93 x 1.25 = 1.2, rounded),
#   0:15. Full: V5c merge + compare 22147508 MEASURED 2:00, 20.2 GB (1244 records; peak = the largest pathway, records
#   of one pathway are loaded together); W2 largest pathway L2/3->L2/3 10.8 GB of npz vs 7.3 GB (L5, V5c) -> 20.2 x 1.48
#   = 29.9 GB x 1.25 -> sbatch --mem=38G; 2:00 x 2.4 (2978 vs 1244 records) x 1.5 = 7:12 -> --time=0:15.
#SBATCH --job-name=valid_v7_cmp
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/valid_v7_cmp_%j.out
set -euo pipefail
: "${FIT:?FIT=<results json>}" "${ECB_REF:?ECB_REF=0|2}" "${TAG:?TAG=<name>}"
export ECB_REF
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1
S1=/scratch/dhuruva/split1
export BCL_EMODEL=delta-split1
export ANALYTICAL_BASIS_DIR=$S1/basis_l5l5 L23_BASIS_DIR=$S1/basis_l23l5 L23L23_BASIS_DIR=$S1/basis_l23l23 CEXP_DIR=$S1/cexp
D=/scratch/dhuruva/bcl_valid_v7_${TAG}; STAGE=${STAGE:-pilot}
if ls $D/shard_*.jsonl >/dev/null 2>&1; then cat $D/shard_*.jsonl | sort -u > $D/merged.jsonl; RES=$D/merged.jsonl
else RES=$D/live_${STAGE}.jsonl; fi
python -u glusynapse_v2/bcl_validation/compare_prefire_v7.py --fit $FIT --results $RES --save $D/cmp_${STAGE}
