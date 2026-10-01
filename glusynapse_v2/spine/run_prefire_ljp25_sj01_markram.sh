#!/bin/bash
# delta-ljp25 part B: Sjostrom 2001 (10) + Markram (7), 24 pairs, hash delta-ljp25-prefire-vseg, v_seg; 136 sims per array task.
# Sizing: vseg 22073523 (MaxRSS 103-110 GiB on 16 -> 140G; 18-25 min) and pilot 22091529 (sj01 mean 245 s, Markram 89 s): 80x245+56x89 = 24600 s / 16 = 26 min + tail, +50% -> 0:45.
# MEASURED 22093947: tasks 119/119/136 ok, 0 failed, 20:53 / 22:48 / 27:45 (limit 0:45). MaxRSS/seff unavailable (slurmdb down): fill in.
#SBATCH --job-name=ljp25_sj01mk
#SBATCH --account=rrg-emuller
#SBATCH --array=0-2
#SBATCH --cpus-per-task=16
#SBATCH --mem=140G
#SBATCH --time=00:45:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/ljp25_sj01mk_%A_%a.out
set -u
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
GLOB=$(python -c "import json; print(json.dumps(json.load(open('glusynapse_v2/spine/delta_ljp25.json'))['globals']))")
S=refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
PAIRS=$(python -c "
p=open('glusynapse_v2/subset24_pairs.txt').read().strip().split(',')
print(','.join(p[${SLURM_ARRAY_TASK_ID:-0}::3]))")
PROTOS=10Hz_-50ms,10Hz_-30ms,10Hz_-10ms,10Hz_5ms,10Hz_10ms,10Hz_30ms,10Hz_50ms,sjostrom_0.1hz_dt+10ms,sjostrom_0.1hz_dt-10ms,sjostrom_10hz_dt+10ms,sjostrom_10hz_dt-10ms,sjostrom_20hz_dt+10ms,sjostrom_20hz_dt-10ms,sjostrom_40hz_dt+10ms,sjostrom_40hz_dt-10ms,sjostrom_50hz_dt+10ms,sjostrom_50hz_dt-10ms
echo "pairs $PAIRS"; echo "protocols $PROTOS"
python -u run_de_fit2_pool.py --cooker --param-hash delta-ljp25-prefire-vseg --glusyn-globals "$GLOB" --results-dir $S \
    --cache cpre_cpost_cache/sabrina_n120_delta-ljp25.pkl --circuit-config data/dhuruva_delta_circuit_config.json \
    --trace-vars cai_CR,shaft_cai,ica_VDCC,v_seg --skip-existing --pairs $PAIRS --workers 16 --protocols $PROTOS
