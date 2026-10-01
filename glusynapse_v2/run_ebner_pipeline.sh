#!/bin/bash
#SBATCH --job-name=v2_ebner_A
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=64
#SBATCH --mem=500G
#SBATCH --time=24:00:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/v2_ebner_A_%j.out
# Job A, after ebner T5 (delta-prefire) and Markram delta-prefire-tr:
#   [1] windowed extraction (shaft_cai + vdcc)  [2] data checks  [3] cooker (post-only) fit on all protocols (T11)
# Job B (run_ebner_v2fit.sh, afterok A): v2.2 fit on all protocols (T17). Each step stops the chain on failure.
set -euo pipefail
source .venv/bin/activate
PF=$PWD; V2=$PF/glusynapse_v2
SIMS=$PF/refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
EB=$V2/extracted/ebner_delta-prefire; MK=$V2/extracted/markram_delta-prefire-tr

echo "=== [0] unit tests $(date)"; python $V2/tests/test_units.py
echo "=== [1] extract $(date)"
python -u $V2/extract_v2.py --param-hash delta-prefire --index-csv $SIMS/../../index_Ebner2019_L5TTPC_L5TTPC.csv \
    --window --workers 30 --skip-existing --out $EB   # ~10 GB per worker (40M-sample traces)
python -u $V2/extract_v2.py --param-hash delta-prefire-tr --workers 60 --skip-existing --out $MK
echo "=== [2] validate $(date)"
python -u $V2/validate_v2.py --dirs $EB,$MK --a $PF/fit_results/delta-cooker.json --pairs 10
echo "=== [3] cooker fit, all protocols $(date)"
python -u $V2/fit_v2.py --dirs $EB,$MK --model post --rho full --workers 60 --maxiter 100 \
    --x0 $PF/fit_results/delta-cooker.json --save $V2/results/fit_post_all
echo "=== done $(date)"
