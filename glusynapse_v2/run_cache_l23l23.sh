#!/bin/bash
# MEASURED 22128321: 2:30, MaxRSS 10.2 GB, CPU eff 72% on 16 -> next 13G, 0:15.
# c_pre/c_post cache, L2/3 PC -> L2/3 PC (Zilberter) on the new l23l23 edges (og-delta, <=120 pairs, 16 workers).
# Sizing from cache_l23l5_rs 22085683 (120 pairs, 16 workers): 5:11, MaxRSS 19.8 GiB, CPU eff 80% -> 25G (+25%), 0:15.
# Same pair count; L2/3 post cells are smaller than L5 TTPCs, so no upward scaling.
#SBATCH --job-name=cache_l23l23
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=16
#SBATCH --mem=25G
#SBATCH --time=00:15:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/cache_l23l23_%j.out
set -euo pipefail
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
python -u precompute_cpre_cpost.py --params defit2 --results-dir refitting_results \
    --sims-dir refitting_results/fitting/n120/seed20262009/Zilberter2009_L23PC_L23PC/simulations \
    --circuit-config data/dhuruva_delta_l23l23_circuit_config.json --edges-h5 data/dhuruva_modified_edges_l23l23.h5 \
    --output cpre_cpost_cache/zilberter_l23l23_delta_rs.pkl --workers 16
