#!/bin/bash
# Compare step of run_valid_v8.sh (submit afterok / afterany, same FIT / ECB_REF / TAG / STAGE env): compare_prefire_v8.py
# (offline = CPU port of the gpu_v8_rho G1 kernel on the split1 records and bases; gate_src 0 fits = compare_prefire_v7)
# -> /scratch/dhuruva/bcl_valid_v8_${TAG}/cmp_${STAGE}_{records,targets}.csv and .png/.pdf. shard_*.jsonl are merged first.
# Sizing basis: V5c compare pilot 22144695 MEASURED 0:20, 0.93 GB; v8 also keeps shaft_cai (one more float trace per
#   record, as large as vdcc) -> 0.93 x 1.25 = 1.2 -> 2G (pilot), 1 CPU, 0:15. Full: v7 compare 22174649 timed out at
#   0:15 with 15.7 GB (resubmitted 25G 1:00); v8 adds the shaft trace -> measure the v7 full compare first, then size
#   with sbatch --mem / --time overrides.
#SBATCH --job-name=valid_v8_cmp
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/valid_v8_cmp_%j.out
set -euo pipefail
: "${FIT:?FIT=<results json>}" "${ECB_REF:?ECB_REF=0|2}" "${TAG:?TAG=<name>}"
export ECB_REF
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1
S1=/scratch/dhuruva/split1
export BCL_EMODEL=delta-split1
export ANALYTICAL_BASIS_DIR=$S1/basis_l5l5 L23_BASIS_DIR=$S1/basis_l23l5 L23L23_BASIS_DIR=$S1/basis_l23l23 CEXP_DIR=$S1/cexp
D=/scratch/dhuruva/bcl_valid_v8_${TAG}; STAGE=${STAGE:-pilot}
if ls $D/shard_*.jsonl >/dev/null 2>&1; then cat $D/shard_*.jsonl | sort -u > $D/merged.jsonl; RES=$D/merged.jsonl
else RES=$D/live_${STAGE}.jsonl; fi
python -u glusynapse_v2/bcl_validation/compare_prefire_v8.py --fit $FIT --results $RES --save $D/cmp_${STAGE}
