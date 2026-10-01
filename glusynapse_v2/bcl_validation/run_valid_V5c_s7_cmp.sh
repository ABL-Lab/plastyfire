#!/bin/bash
# Compare step of run_valid_V5c_s7.sh (submit afterok): compare_prefire_v5.py (offline = CPU port of the gpu_v5_rho
# kernel) -> results/V5c_s7_${STAGE}_{targets,records}.csv and V5c_s7_${STAGE}.png/.pdf; the log holds the markdown
# verdict table for BCL_V5c_s7.md.
# Sizing: v4 compare pilot 22135323 (38 records) MEASURED 17 s, MaxRSS 0.95 GB -> 1 CPU, 2G, 0:15 (numba JIT of the
#   v5 port adds seconds). Full (1295 records): v4 plan 18G, 0:30; resize from this pilot's seff.
# MEASURED pilot 22138637: 15 s, MaxRSS 0.94 GB. Full: 1 CPU, 18G (v4 plan), 0:30. Raised to 22G while pending (v4 full compare 22135838 used 17.99 GB of 18G).
#SBATCH --job-name=valid_V5c_s7_cmp
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/valid_V5c_s7_cmp_%j.out
set -euo pipefail
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1
python -u glusynapse_v2/bcl_validation/compare_prefire_v5.py --fit glusynapse_v2/rho_redesign/results/v5_V5c_s7.json \
    --results /scratch/dhuruva/bcl_valid_V5c_s7/live.jsonl \
    --save glusynapse_v2/bcl_validation/results/V5c_s7_${STAGE:-pilot}
