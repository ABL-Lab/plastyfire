#!/bin/bash
# Compare step of run_valid_C1Ajn_s5.sh (submit afterok): compare_prefire_v4.py -> results/C1Ajn_s5_${STAGE}_{targets,records}.csv
# and C1Ajn_s5_${STAGE}.png/.pdf; the log holds the markdown verdict table for BCL_C1Ajn_s5.md.
# Sizing (pilot, 36 records): no standalone measurement; equiv 22131696 ran compare_v4 (5 records) inside its 4.66 GB
#   peak; eval_v4 on all 111 L2/3 pairs: 15.9 GB for 5 fits (22108820). Pilot 3G, 0:15. Full: resize from the pilot seff.
# MEASURED pilot 22135323 (38 records): 17 s, MaxRSS 0.95 GB. Full (1295 records, ~34x): BatchV2 ~12-20 MB/record (pilot;
# MEASURED full 22135838 (1295 records): 4:36, MaxRSS 17.99 GB of 18G (at the limit) -> 22G, 0:15 next time.
#   eval_v4 111 L2/3 pairs 15.9 GB) -> L5 716 records 9-14 GB peak, +25% -> 18G; ~10 min -> 0:30.
#   sbatch --mem=18G --time=00:30:00 --export=ALL,STAGE=full run_valid_C1Ajn_s5_cmp.sh
#SBATCH --job-name=valid_C1Ajn_cmp
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=00:15:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/valid_C1Ajn_cmp_%j.out
set -euo pipefail
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1
python -u glusynapse_v2/bcl_validation/compare_prefire_v4.py --fit glusynapse_v2/rho_redesign/results/v4_C1Ajn_s5.json \
    --results /scratch/dhuruva/bcl_valid_C1Ajn_s5/live.jsonl \
    --save glusynapse_v2/bcl_validation/results/C1Ajn_s5_${STAGE:-pilot}
