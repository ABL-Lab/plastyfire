#!/bin/bash
# T29: c_pre/c_post cache on the re-split L23->L5 edges (og-delta, 120 pairs, 16 workers).
# T29: rho0 re-split L23->L5 (new edges). Sizing: cache 22064518 (same 120 pairs, 16 workers): 4:41, MaxRSS 17.0 GiB -> 22G, 15 min.
# MEASURED 22085683: 5:11, MaxRSS 19.8 GiB (90% of 22G), CPU eff 80%.
#SBATCH --job-name=cache_l23l5_rs
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=16
#SBATCH --mem=22G
#SBATCH --time=00:15:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/cache_l23l5_rs_%j.out
set -euo pipefail
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
python -u precompute_cpre_cpost.py --params defit2 --results-dir refitting_results \
    --sims-dir refitting_results/fitting/n120/seed20262009/Ebner2019_L23PC_L5TTPC/simulations \
    --circuit-config data/dhuruva_delta_l23l5_circuit_config.json --edges-h5 data/dhuruva_modified_edges_l23l5.h5 \
    --output cpre_cpost_cache/ebner_l23l5_delta_rs.pkl --workers 16
