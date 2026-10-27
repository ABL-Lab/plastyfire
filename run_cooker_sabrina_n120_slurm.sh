#!/bin/bash
#SBATCH --job-name=cooker_sabrina_n120
#SBATCH --account=rrg-emuller
#SBATCH --time=02:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=60
#SBATCH --output=/home/dhuruva/ws/plastyfire/logs/cooker_sabrina_n120_%j.out
#SBATCH --error=/home/dhuruva/ws/plastyfire/logs/cooker_sabrina_n120_%j.err
set -euo pipefail

cd /home/dhuruva/ws/plastyfire
source /home/dhuruva/ws/plastyfire/.venv/bin/activate

SIMS=refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations

echo "=== BCL prefire STDP, cooker ($(date)) ==="
python -u run_de_fit2_pool.py --cooker --workers 60 --results-dir $SIMS \
    --cache cpre_cpost_cache/sabrina_n120_defit2.pkl

echo "=== plot ($(date)) ==="
python -u plot_stdp_ic_best.py --results-dir $SIMS --param-hash cooker \
    --basis-dir basis_results_edges_sabrina_n120 \
    --label "cooker, n120 (BCL prefire)" --out stdp_cooker_sabrina_n120.png
