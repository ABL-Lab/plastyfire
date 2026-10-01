#!/bin/bash
# NOT SUBMITTED (est. ~48 CPU h > 30, orchestrator logs first). t_drive 3 full re-prefire, part A: Sjostrom 2003 (9 protocols; the nreps-15 bursts are
# superseded by r50) + 2 r50 bursts, all 24 subset24 pairs, hash delta-prefire-vseg. 3 array tasks x 8 pairs x 11 protocols = 88 sims.
# MEASURED 22073522: MaxRSS 131-164 GiB of 195G (84%), 20-27 min, ~80% CPU on 32; 242 new + 22 pilot sims. Keep 195G, 0:45.
# Sizing: pilot 22073044 (v_seg, 26 sims): mean 583 s/sim for these 11 protocols -> 88 x 583 / 32 = 27 min + tail = ~30 min;
# MaxRSS 3.9 GB/worker (pilot), 4.8 GB/worker in 22057545 (154 GB / 32) -> 32 x 4.8 = 154 GB + 25% = 195G. Time 30 min + 50% = 0:45.
#SBATCH --job-name=vseg_sj03
#SBATCH --account=rrg-emuller
#SBATCH --array=0-2
#SBATCH --cpus-per-task=32
#SBATCH --mem=195G
#SBATCH --time=00:45:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/vseg_sj03_%A_%a.out
set -u
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
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
python -u run_de_fit2_pool.py --cooker --param-hash delta-prefire-vseg --results-dir $S \
    --cache cpre_cpost_cache/sabrina_n120_delta.pkl --circuit-config data/dhuruva_delta_circuit_config.json \
    --trace-vars cai_CR,shaft_cai,ica_VDCC,v_seg --skip-existing --pairs $PAIRS --workers 32 --protocols $PROTOS
