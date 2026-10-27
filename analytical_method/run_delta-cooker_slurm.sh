#!/bin/bash
#SBATCH --job-name=delta-cooker
#SBATCH --account=rrg-emuller
#SBATCH --time=04:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=60
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/delta-cooker_%j.out
#SBATCH --error=/project/rrg-emuller/dhuruva/plastyfire/logs/delta-cooker_%j.err
set -euo pipefail
# delta-cooker: cooker refit on the delta emodel (SSCx-AAD-delta-emodels), end to end:
#   c_pre/c_post cache -> BCL prefire (traces) -> extract -> [wait for basis]
#   -> cooker fit -> submit full-protocol validation of the fit.
# The basis runs separately (one sbatch per pair, basis_results_edges_sabrina_n120_delta).
PF=/project/rrg-emuller/dhuruva/plastyfire
cd $PF
source $PF/.venv/bin/activate
export PYTHONPATH=$PF

NAME=delta-cooker
CIRCUIT=$PF/data/dhuruva_delta_circuit_config.json
CACHE=$PF/cpre_cpost_cache/sabrina_n120_delta.pkl
BASIS=$PF/basis_results_edges_sabrina_n120_delta
EXTRACTED=$PF/analytical_method/extracted_sabrina_n120_delta_d0p25
SIMS=$PF/refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
PREFIRE_HASH=delta-prefire    # cooker (gamma-fit) params, only used to generate the calcium traces
N_PAIRS=$(ls -d $SIMS/*-* | wc -l)

echo "=== [1] c_pre/c_post cache, delta ($(date)) ==="
python -u precompute_cpre_cpost.py --params defit2 --results-dir refitting_results \
    --circuit-config $CIRCUIT --edges-h5 data/dhuruva_modified_edges.h5 \
    --output $CACHE --workers 60

echo "=== [2] BCL prefire, delta ($(date)) ==="
python -u run_de_fit2_pool.py --cooker --param-hash $PREFIRE_HASH --workers 60 \
    --results-dir $SIMS --cache $CACHE --circuit-config $CIRCUIT

echo "=== [3] extract calcium traces, 0.25 ms ($(date)) ==="
cd $PF/analytical_method
python -u extract.py --workers 60 --decim 0.25 --param-hash $PREFIRE_HASH \
    --sims $SIMS --cache $CACHE --out $EXTRACTED

echo "=== [4] wait for basis ($(date)) ==="
while [ $(ls $BASIS/basis_*.csv 2>/dev/null | wc -l) -lt $N_PAIRS ] \
      && squeue -u $USER -h -o '%j' | grep -q '^basis_'; do sleep 30; done
echo "basis csvs: $(ls $BASIS/basis_*.csv | wc -l) / $N_PAIRS"

echo "=== [5] cooker fit ($(date)) ==="
export ANALYTICAL_EXTRACTED=$EXTRACTED ANALYTICAL_BASIS_DIR=$BASIS
mkdir -p $PF/fit_results
python -u fit.py --de --workers 60 --maxiter 600 --stride 1 \
    --fit-dt=-10,5,10 --anchor-tails --save-json $PF/fit_results/$NAME.json

echo "=== [6] submit full-protocol validation ($(date)) ==="
sbatch --job-name=full_$NAME \
    --export=ALL,NAME=$NAME,PARAMS_JSON=$PF/fit_results/$NAME.json,CIRCUIT=$CIRCUIT,CACHE=$CACHE,BASIS=$BASIS \
    $PF/run_full_stdp_slurm.sh
echo "=== done ($(date)) ==="
