#!/bin/bash
# T30 pilot: og-delta prefire (hash delta-prefire-vseg, v_seg recorded) of the 4 Sjostrom 2007 step protocols
# on 2 subset24 pairs = 8 sims, 8 workers (one wave).
# Sizing from the vseg pilot 22073044 (26 workers, 17:25, MaxRSS 101.8 GB = 3.9 GB/worker, 502 s biological
# time per sim): 8 x 3.9 = 31 GB + 25% -> 39G. Step sims are 302 s biological, but have ~300 pre + ~300 post
# APs; the time is taken as <= the 17:25 pilot, + 50% -> 0:30.
# MEASURED 22089134: 8/8 ok, 10:43, MaxRSS 31.66 GB (81% of 39G, 3.96 GB/worker), CPU eff 53% (one wave, sims 146-643 s,
#   mean 378 s). Guardrail passed: 10 APs per step in every step (12-13 at 1.2 nA).
#SBATCH --job-name=vseg_sj07_pilot
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=8
#SBATCH --mem=39G
#SBATCH --time=00:30:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/vseg_sj07_pilot_%j.out
set -u
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
S=refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
PAIRS=${PAIRS:-180351-198084,181455-195199}
PROTOS=sjostrom07_step200ms_pair,sjostrom07_step200ms_pre_only,sjostrom07_step200ms_post_only,sjostrom07_step200ms_1.2nA_pair
echo "pairs $PAIRS"; echo "protocols $PROTOS"
python -u run_de_fit2_pool.py --cooker --param-hash delta-prefire-vseg --results-dir $S \
    --cache cpre_cpost_cache/sabrina_n120_delta.pkl --circuit-config data/dhuruva_delta_circuit_config.json \
    --trace-vars cai_CR,shaft_cai,ica_VDCC,v_seg --skip-existing --pairs $PAIRS --workers 8 --protocols $PROTOS
