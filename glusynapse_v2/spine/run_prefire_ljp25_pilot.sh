#!/bin/bash
# MEASURED 22091529: 64/64 ok, 23:41 on 32 CPUs, globals ljp 25 / gca_bar 0.02306 confirmed in pairrunner fit_params; MaxRSS/seff not available (slurmdb down), fill in.
# delta-ljp25 pilot: 2 pairs x all protocols of sj03 (9+2 r50), sj01 (10), Markram (7), sj07 (4) = up to 2 x 32 (sj07 pre_only may be absent) sims,
# hash delta-ljp25-prefire-vseg, v_seg recorded. Sizing: 32 workers x 4.9 GB (sj07 22089771) x1.25 for v_seg-heavy sj03 (4.8-5.5 GB/worker in vseg parts) -> 32 x 5.5 = 176 + 25% = 220G.
# Time: ~583 s/sim sj03, 136-490 s sj01/Markram, 378 s sj07; per pair ~ 11x583+17x~190+4x378 = ~11000 s; 2 pairs / 32 CPUs -> ~12-20 min -> 0:45.
#SBATCH --job-name=ljp25_pilot
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=32
#SBATCH --mem=220G
#SBATCH --time=00:45:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/ljp25_pilot_%j.out
set -u
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
S=refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
GLOB=$(python -c "import json; print(json.dumps(json.load(open('glusynapse_v2/spine/delta_ljp25.json'))['globals']))")
PAIRS=${PAIRS:-180351-198084,181455-195199}
PROTOS=$(python -c "
import pandas as pd
o=[]
for f in ('index_Sjostrom2003_L5TTPC_L5TTPC.csv','index_Sjostrom2003b_L5TTPC_L5TTPC.csv'):
    o+=list(dict.fromkeys(pd.read_csv('refitting_results/fitting/n120/seed20262009/'+f).protocol_id))
print(','.join(x for x in o if x not in ('sjostrom_burst5x20hz_dt-120ms','sjostrom_burst5x20hz_dt-200ms')))")
PROTOS=$PROTOS,10Hz_-50ms,10Hz_-30ms,10Hz_-10ms,10Hz_5ms,10Hz_10ms,10Hz_30ms,10Hz_50ms,sjostrom_0.1hz_dt+10ms,sjostrom_0.1hz_dt-10ms,sjostrom_10hz_dt+10ms,sjostrom_10hz_dt-10ms,sjostrom_20hz_dt+10ms,sjostrom_20hz_dt-10ms,sjostrom_40hz_dt+10ms,sjostrom_40hz_dt-10ms,sjostrom_50hz_dt+10ms,sjostrom_50hz_dt-10ms,sjostrom07_step200ms_pair,sjostrom07_step200ms_pre_only,sjostrom07_step200ms_post_only,sjostrom07_step200ms_1.2nA_pair
echo "pairs $PAIRS"; echo "globals $GLOB"; echo "protocols $PROTOS"
python -u run_de_fit2_pool.py --cooker --param-hash delta-ljp25-prefire-vseg --glusyn-globals "$GLOB" --results-dir $S \
    --cache cpre_cpost_cache/sabrina_n120_delta-ljp25.pkl --circuit-config data/dhuruva_delta_circuit_config.json \
    --trace-vars cai_CR,shaft_cai,ica_VDCC,v_seg --skip-existing --pairs $PAIRS --workers ${WORKERS:-32} --protocols $PROTOS
