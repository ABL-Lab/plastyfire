#!/bin/bash
# delta-ljp25 part C: Sjostrom 2007 step protocols (4), 24 pairs, hash delta-ljp25-prefire-vseg, v_seg; pilot sims skipped.
# Sizing: vseg sj07 22089771 (4.9 GB/worker, 21:44 on 16, 79 GB at limit) -> 16 x 4.9 x 1.25 = 98G; pilot sj07 mean 324 s/sim -> ~22 min, +50% -> 0:45.
# MEASURED 22093948: 86/86 ok, 20:41 (limit 0:45). MaxRSS/seff unavailable (slurmdb down): fill in.
#SBATCH --job-name=ljp25_sj07
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=16
#SBATCH --mem=98G
#SBATCH --time=00:45:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/ljp25_sj07_%j.out
set -u
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
GLOB=$(python -c "import json; print(json.dumps(json.load(open('glusynapse_v2/spine/delta_ljp25.json'))['globals']))")
S=refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
PAIRS=$(cat glusynapse_v2/subset24_pairs.txt)
PROTOS=sjostrom07_step200ms_pair,sjostrom07_step200ms_pre_only,sjostrom07_step200ms_post_only,sjostrom07_step200ms_1.2nA_pair
echo "pairs $PAIRS"; echo "protocols $PROTOS"
python -u run_de_fit2_pool.py --cooker --param-hash delta-ljp25-prefire-vseg --glusyn-globals "$GLOB" --results-dir $S \
    --cache cpre_cpost_cache/sabrina_n120_delta-ljp25.pkl --circuit-config data/dhuruva_delta_circuit_config.json \
    --trace-vars cai_CR,shaft_cai,ica_VDCC,v_seg --skip-existing --pairs $PAIRS --workers 16 --protocols $PROTOS
