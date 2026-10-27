#!/bin/bash
# T24: c_pre/c_post cache (defit2 mini-sims, og-delta) for the 120 L23PC->L5TTPC pairs, own pkl.
# Sizing: v2_cache_sv 22005793, 120 pairs x 16 workers: 6 min, MaxRSS 13.3 GB -> 17G.
#SBATCH --job-name=cache_l23l5
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=16
#SBATCH --mem=17G
#SBATCH --time=01:00:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/cache_l23l5_%j.out
set -euo pipefail
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
python -u precompute_cpre_cpost.py --params defit2 --results-dir refitting_results \
    --sims-dir refitting_results/fitting/n120/seed20262009/Ebner2019_L23PC_L5TTPC/simulations \
    --circuit-config data/dhuruva_delta_circuit_config.json --edges-h5 data/dhuruva_modified_edges.h5 \
    --output cpre_cpost_cache/ebner_l23l5_delta.pkl --workers 16
