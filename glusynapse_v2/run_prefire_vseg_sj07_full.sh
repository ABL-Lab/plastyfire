#!/bin/bash
# T30 full: og-delta prefire (hash delta-prefire-vseg, v_seg) of the 4 Sjostrom 2007 step protocols on all 24
# subset24 pairs = 94 workdirs (pair 189436-193575 has no _pair / _pre_only: the pre cell fails the 10-12 AP search);
# the 8 pilot sims are skipped (--skip-existing) -> 86 sims.
# Sizing from pilot 22089134 (8 workers: 10:43, MaxRSS 31.66 GB = 3.96 GB/worker, sims mean 378 s, max 643 s):
# 86 x 378 s = 9.0 CPU-h; 16 workers -> 34 min + 11 min tail = ~45 min, + 50% -> 1:15.
# Memory 16 x 3.96 = 63 GB + 25% -> 79G.
# MEASURED 22089771: 86/86 ok, 21:44, MaxRSS 78.99 GB (100% of 79G, no OOM kill but at the limit: 4.9 GB/worker),
#   CPU eff 87%, 5.0 CPU-h. A rerun needs 16 x 4.9 x 1.25 = 98G, time 22 min + 50% -> 0:45.
#SBATCH --job-name=vseg_sj07
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=16
#SBATCH --mem=79G
#SBATCH --time=01:15:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/vseg_sj07_%j.out
set -u
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
S=refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
PAIRS=$(cat glusynapse_v2/subset24_pairs.txt)
PROTOS=sjostrom07_step200ms_pair,sjostrom07_step200ms_pre_only,sjostrom07_step200ms_post_only,sjostrom07_step200ms_1.2nA_pair
echo "pairs $PAIRS"; echo "protocols $PROTOS"
python -u run_de_fit2_pool.py --cooker --param-hash delta-prefire-vseg --results-dir $S \
    --cache cpre_cpost_cache/sabrina_n120_delta.pkl --circuit-config data/dhuruva_delta_circuit_config.json \
    --trace-vars cai_CR,shaft_cai,ica_VDCC,v_seg --skip-existing --pairs $PAIRS --workers 16 --protocols $PROTOS
