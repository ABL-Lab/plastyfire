#!/bin/bash
# T25 step 6: og-delta prefire of the 11 Sjostrom 2003 protocols (configs/Sjostrom2003_L5TTPC_L5TTPC.yaml) on
# subset24, same hash/cache/trace vars as the Ebner delta-prefire set, so they extract next to it.
# Array task i runs pairs i, i+3, ... (8 pairs x 11 protocols each). Memory: sv_prefire 22005844 peaked at
# 346 GB with 46 workers on 500-1000 s protocols; 32 workers here -> ~240 GB + 25%.
#SBATCH --job-name=sj03_prefire
#SBATCH --account=rrg-emuller
#SBATCH --array=0-2
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=04:00:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/sj03_prefire_%A_%a.out
set -u
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
IDX=refitting_results/fitting/n120/seed20262009/index_Sjostrom2003_L5TTPC_L5TTPC.csv
S=refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
PAIRS=$(python -c "
p=open('glusynapse_v2/subset24_pairs.txt').read().strip().split(',')
print(','.join(p[${SLURM_ARRAY_TASK_ID:-0}::3]))")
PROTOS=$(python -c "
import pandas as pd; d=pd.read_csv('$IDX'); print(','.join(dict.fromkeys(d.protocol_id)))")
echo "pairs $PAIRS"; echo "protocols $PROTOS"
python -u run_de_fit2_pool.py --cooker --param-hash delta-prefire --results-dir $S \
    --cache cpre_cpost_cache/sabrina_n120_delta.pkl --circuit-config data/dhuruva_delta_circuit_config.json \
    --trace-vars cai_CR,shaft_cai,ica_VDCC --skip-existing --pairs $PAIRS --workers 32 --protocols $PROTOS
