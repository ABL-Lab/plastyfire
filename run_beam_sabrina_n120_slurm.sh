#!/bin/bash
#SBATCH --job-name=beam_sabrina_n120
#SBATCH --account=rrg-emuller
#SBATCH --time=02:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=60
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/beam_sabrina_n120_%j.out
#SBATCH --error=/project/rrg-emuller/dhuruva/plastyfire/logs/beam_sabrina_n120_%j.err
set -euo pipefail

cd /project/rrg-emuller/dhuruva/plastyfire
source /project/rrg-emuller/dhuruva/plastyfire/.venv/bin/activate

SIMS=refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations

echo "=== BCL prefire STDP, beam ($(date)) ==="
python -u run_de_fit2_pool.py --beam --workers 60 --results-dir $SIMS \
    --cache cpre_cpost_cache/sabrina_n120_defit2.pkl

echo "=== plot ($(date)) ==="
python -u plot_stdp_ic_best.py --results-dir $SIMS --param-hash beam \
    --basis-dir basis_results_edges_sabrina_n120 \
    --label "beam, n120 (BCL prefire)" --out stdp_beam_sabrina_n120.png
