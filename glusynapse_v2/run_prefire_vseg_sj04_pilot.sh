#!/bin/bash
# T31 pilot: og-delta prefire (hash delta-prefire-vseg, v_seg recorded) of sjostrom04_dltd_step250ms on 2 subset24 pairs = 2 sims,
# 2 workers. Basis: sj07 pilot 22089134 (8 workers, 10:43, MaxRSS 31.66 GB = 3.96 GB/worker, sims 146-643 s, 302 s bio) and vseg
# pilot 22073044 (502 s bio per sim, 17:25 total, 3.9 GB/worker); this sim is 502 s bio (50 x 10 s) but with ~0 post APs.
# Memory 2 x 3.96 = 7.9 GB + 25% -> 10G. Time <= 17:25 + 50% = 26 min -> 0:30. Record seff here after the run.
# MEASURED 22127911: 6:33, MaxRSS 9.99 GB on 2 workers (~5-6 GB/worker, at the 10G limit), 261-390 s per pair, both ok, 0 APs in steps.
#SBATCH --job-name=vseg_sj04_pilot
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --time=00:30:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/vseg_sj04_pilot_%j.out
set -u
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
S=refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
PAIRS=${PAIRS:-180351-198084,181455-195199}
echo "pairs $PAIRS"
python -u run_de_fit2_pool.py --cooker --param-hash delta-prefire-vseg --results-dir $S \
    --cache cpre_cpost_cache/sabrina_n120_delta.pkl --circuit-config data/dhuruva_delta_circuit_config.json \
    --trace-vars cai_CR,shaft_cai,ica_VDCC,v_seg --skip-existing --pairs $PAIRS --workers 2 --protocols sjostrom04_dltd_step250ms
