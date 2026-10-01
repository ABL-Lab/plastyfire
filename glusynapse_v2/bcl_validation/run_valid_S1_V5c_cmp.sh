#!/bin/bash
# Compare step of run_valid_S1_V5c.sh (submit afterany): compare_prefire_v5.py (offline = CPU port of the gpu_v5_rho
# kernel on the split1 extracted records and bases) -> results/S1_V5c_${STAGE}_{targets,records}.csv and .png/.pdf.
# Sizing basis: V5c_s7 compare pilot 22138637 MEASURED 15 s, MaxRSS 0.94 GB -> 1 CPU, 2G, 0:15.
# Full: V5c_s7 compare used 18-22G, 0:30; resize from the pilot / split1 record count.
#SBATCH --job-name=valid_S1_V5c_cmp
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/valid_S1_V5c_cmp_%j.out
set -euo pipefail
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1
S1=/scratch/dhuruva/split1
export BCL_EMODEL=delta-split1
export ANALYTICAL_BASIS_DIR=$S1/basis_l5l5
export L23_BASIS_DIR=$S1/basis_l23l5
D=/scratch/dhuruva/bcl_valid_S1_V5c; if ls $D/shard_*.jsonl >/dev/null 2>&1; then cat $D/shard_*.jsonl | sort -u > $D/merged.jsonl; RES=$D/merged.jsonl; else RES=$D/live.jsonl; fi   # shards: merge (exact duplicate lines = shared base)
python -u glusynapse_v2/bcl_validation/compare_prefire_v5.py --fit glusynapse_v2/rho_redesign/results/v5_S1_V5c_s6_39.json \
    --results $RES \
    --save glusynapse_v2/bcl_validation/results/S1_V5c_${STAGE:-pilot}
