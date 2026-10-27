#!/bin/bash
#SBATCH --job-name=defit2_sabrina_n120
#SBATCH --account=rrg-emuller
#SBATCH --time=04:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=60
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/defit2_sabrina_n120_%j.out
#SBATCH --error=/project/rrg-emuller/dhuruva/plastyfire/logs/defit2_sabrina_n120_%j.err
set -euo pipefail

cd /project/rrg-emuller/dhuruva/plastyfire
source /project/rrg-emuller/dhuruva/plastyfire/.venv/bin/activate

CIRCUIT=data/dhuruva_modified_ion_channels_circuit_config.json
EDGES=data/dhuruva_modified_edges.h5
CACHE=cpre_cpost_cache/sabrina_n120_defit2.pkl
THR=threshold_results_defit2_sabrina_n120
SIMS=refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
mkdir -p cpre_cpost_cache

echo "=== [4a] c_pre/c_post cache ($(date)) ==="
python -u precompute_cpre_cpost.py --params defit2 \
    --results-dir refitting_results \
    --circuit-config $CIRCUIT --edges-h5 $EDGES \
    --output $CACHE --workers 60

echo "=== [4b] thresholds + inject ($(date)) ==="
python -u compute_thresholds_from_cache.py --cache $CACHE --output-dir $THR \
    --a00 1.003498 --a01 2.902478 --a10 1.644558 --a11 2.764812 \
    --a20 1.003498 --a21 2.902478 --a30 1.644558 --a31 2.764812 \
    --edges-h5 $EDGES
python -u inject_thresholds.py --results-dir $THR --edges-h5 $EDGES --force

echo "=== [5] BCL prefire STDP, DE fit #2 ($(date)) ==="
python -u run_de_fit2_pool.py --workers 60 --results-dir $SIMS --cache $CACHE

echo "=== done ($(date)) ==="
