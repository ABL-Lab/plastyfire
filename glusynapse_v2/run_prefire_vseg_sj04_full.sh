#!/bin/bash
# T31: og-delta prefire (hash delta-prefire-vseg, v_seg) of sjostrom04_dltd_step250ms on the 24 subset24 pairs, 12 workers
# (2 rounds). Basis: pilot 22127911 (2 workers, 6:33, MaxRSS 9.99 GB at the 10G limit -> 6 GB/worker, 261-390 s/sim):
# 12 x 6 = 72 GB + 25% -> 90G. 2 rounds x ~6.5 min = 13 min + 50% -> 0:30. --skip-existing skips the 2 pilot pairs.
#SBATCH --job-name=vseg_sj04_full
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=12
# MEASURED 22128324 (24 pairs, 12 workers): 11:41, 56.2 GB MaxRSS (62%), 72% CPU -> next 70G, 0:30.
#SBATCH --mem=90G
#SBATCH --time=00:30:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/vseg_sj04_full_%j.out
set -u
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
S=refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
PAIRS=$(cat glusynapse_v2/subset24_pairs.txt)
python -u run_de_fit2_pool.py --cooker --param-hash delta-prefire-vseg --results-dir $S \
    --cache cpre_cpost_cache/sabrina_n120_delta.pkl --circuit-config data/dhuruva_delta_circuit_config.json \
    --trace-vars cai_CR,shaft_cai,ica_VDCC,v_seg --skip-existing --pairs $PAIRS --workers 12 --protocols sjostrom04_dltd_step250ms
