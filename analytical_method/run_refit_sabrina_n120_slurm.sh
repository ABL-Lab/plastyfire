#!/bin/bash
#SBATCH --job-name=refit_sabrina_n120
#SBATCH --account=rrg-emuller
#SBATCH --time=06:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=60
#SBATCH --output=/home/dhuruva/ws/plastyfire/logs/refit_sabrina_n120_%j.out
#SBATCH --error=/home/dhuruva/ws/plastyfire/logs/refit_sabrina_n120_%j.err
set -euo pipefail

PF=/home/dhuruva/ws/plastyfire
cd $PF/analytical_method
source $PF/.venv/bin/activate
export PYTHONPATH=$PF

# n120 inputs; kept separate from the n100 data in extracted/
export ANALYTICAL_EXTRACTED=$PF/analytical_method/extracted_sabrina_n120_d0p25
export ANALYTICAL_BASIS_DIR=$PF/basis_results_edges_sabrina_n120

echo "=== [1] extract calcium traces, 0.25 ms ($(date)) ==="
python -u extract.py --workers 60 --decim 0.25 --param-hash b8c7ff3ecf0a \
    --sims $PF/refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations \
    --cache $PF/cpre_cpost_cache/sabrina_n120_defit2.pkl \
    --out $ANALYTICAL_EXTRACTED
ls $ANALYTICAL_EXTRACTED | wc -l; du -sh $ANALYTICAL_EXTRACTED

echo "=== [2] validate offline model vs BCL DE fit #2 run ($(date)) ==="
python -u validate.py --ref defit2_n120

echo "=== [3] refit: cpre-only, DE, stride 1 ($(date)) ==="
python -u fit.py --cpre-only --fit-gamma --de \
    --workers 60 --maxiter 600 --stride 1 --fit-dt=-10,5,10 --anchor-tails

echo "=== done ($(date)) ==="
