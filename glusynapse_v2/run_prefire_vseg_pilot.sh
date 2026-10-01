#!/bin/bash
# t_drive 3 pilot: og-delta re-prefire with v_seg recorded, new hash delta-prefire-vseg (nothing overwritten).
# 2 pairs x (11 Sjostrom 2003 + 2 r50 burst) = 26 sims, 26 workers (one wave, 1 sim per CPU).
# Sizing from sj03_prefire 22057545: 32 workers, 143-154 GB MaxRSS (4.8 GB/worker), 19-26 min for 88 sims (~10.5 min/sim).
# v_seg adds a 4th trace var: x1.33 -> 26 x 4.8 x 1.33 = 166 GB, +25% = 200G. Time: ~11-12 min + 50% -> 0:30.
#SBATCH --job-name=vseg_pilot
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=26
#SBATCH --mem=200G
#SBATCH --time=00:30:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/vseg_pilot_%j.out
set -u
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
S=refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
PAIRS=${PAIRS:-180351-198084,181455-195199}
PROTOS=$(python -c "
import pandas as pd
o=[]
for f in ('index_Sjostrom2003_L5TTPC_L5TTPC.csv','index_Sjostrom2003b_L5TTPC_L5TTPC.csv'):
    o+=list(dict.fromkeys(pd.read_csv('refitting_results/fitting/n120/seed20262009/'+f).protocol_id))
print(','.join(o))")
echo "pairs $PAIRS"; echo "protocols $PROTOS"
python -u run_de_fit2_pool.py --cooker --param-hash delta-prefire-vseg --results-dir $S \
    --cache cpre_cpost_cache/sabrina_n120_delta.pkl --circuit-config data/dhuruva_delta_circuit_config.json \
    --trace-vars cai_CR,shaft_cai,ica_VDCC,v_seg --skip-existing --pairs $PAIRS --workers ${WORKERS:-26} --protocols $PROTOS
