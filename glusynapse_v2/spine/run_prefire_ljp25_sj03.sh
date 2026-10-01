#!/bin/bash
# delta-ljp25 (spine/delta_ljp25.json: ljp_VDCC 25, gca_bar 0.02306) part A: Sjostrom 2003 (9) + 2 r50, 24 pairs, hash delta-ljp25-prefire-vseg, v_seg; 88 sims per array task.
# Sizing: vseg 22073522 (MaxRSS 131-164 GiB on 32 workers -> 195G) and pilot 22091529 (64 sims ok, 23:41; sj03 mean 663 s/sim): 88 x 663 / 32 = 30 min + tail = ~35 min, +50% -> 1:00.
# MEASURED 22093946: tasks 77/77/88 ok, 0 failed, 19:57 / 24:06 / 27:50 (limit 1:00, could be 0:45). MaxRSS/seff unavailable (slurmdb down): fill in.
#SBATCH --job-name=ljp25_sj03
#SBATCH --account=rrg-emuller
#SBATCH --array=0-2
#SBATCH --cpus-per-task=32
#SBATCH --mem=195G
#SBATCH --time=01:00:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/ljp25_sj03_%A_%a.out
set -u
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
GLOB=$(python -c "import json; print(json.dumps(json.load(open('glusynapse_v2/spine/delta_ljp25.json'))['globals']))")
S=refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
PAIRS=$(python -c "
p=open('glusynapse_v2/subset24_pairs.txt').read().strip().split(',')
print(','.join(p[${SLURM_ARRAY_TASK_ID:-0}::3]))")
PROTOS=$(python -c "
import pandas as pd
o=[]
for f in ('index_Sjostrom2003_L5TTPC_L5TTPC.csv','index_Sjostrom2003b_L5TTPC_L5TTPC.csv'):
    o+=list(dict.fromkeys(pd.read_csv('refitting_results/fitting/n120/seed20262009/'+f).protocol_id))
print(','.join(x for x in o if x not in ('sjostrom_burst5x20hz_dt-120ms','sjostrom_burst5x20hz_dt-200ms')))")
echo "pairs $PAIRS"; echo "protocols $PROTOS"
python -u run_de_fit2_pool.py --cooker --param-hash delta-ljp25-prefire-vseg --glusyn-globals "$GLOB" --results-dir $S \
    --cache cpre_cpost_cache/sabrina_n120_delta-ljp25.pkl --circuit-config data/dhuruva_delta_circuit_config.json \
    --trace-vars cai_CR,shaft_cai,ica_VDCC,v_seg --skip-existing --pairs $PAIRS --workers 32 --protocols $PROTOS
