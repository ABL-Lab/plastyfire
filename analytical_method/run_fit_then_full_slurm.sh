#!/bin/bash
#SBATCH --job-name=fit_then_full
#SBATCH --account=rrg-emuller
#SBATCH --time=03:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=60
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/fit_%x_%j.out
#SBATCH --error=/project/rrg-emuller/dhuruva/plastyfire/logs/fit_%x_%j.err
set -euo pipefail
# Analytical fit, then submit the full-protocol BCL validation of its best params
# the moment the fit finishes (run_full_stdp_slurm.sh, no basis).
#
#   sbatch --job-name=myfit --export=ALL,NAME=myfit,FIT_FLAGS="--cpre-only --fit-gamma" \
#       analytical_method/run_fit_then_full_slurm.sh
#
# Common flags (--de --workers 60 --maxiter 600 --stride 1 --fit-dt=-10,5,10
# --anchor-tails) are always passed; FIT_FLAGS adds the method-specific ones.
: "${NAME:?set NAME}"
FIT_FLAGS="${FIT_FLAGS:-}"

PF=/project/rrg-emuller/dhuruva/plastyfire
cd $PF/analytical_method
source $PF/.venv/bin/activate
export PYTHONPATH=$PF
export ANALYTICAL_EXTRACTED=$PF/analytical_method/extracted_sabrina_n120_d0p25
export ANALYTICAL_BASIS_DIR=$PF/basis_results_edges_sabrina_n120
JSON=$PF/fit_results/${NAME}.json
mkdir -p $PF/fit_results

echo "=== fit $NAME: $FIT_FLAGS ($(date)) ==="
python -u fit.py --de --workers 60 --maxiter 600 --stride 1 \
    --fit-dt=-10,5,10 --anchor-tails $FIT_FLAGS --save-json $JSON

echo "=== submit full-protocol validation ($(date)) ==="
sbatch --job-name=full_$NAME --export=ALL,NAME=$NAME,PARAMS_JSON=$JSON \
    $PF/run_full_stdp_slurm.sh
