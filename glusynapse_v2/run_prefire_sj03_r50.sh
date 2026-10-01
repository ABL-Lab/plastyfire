#!/bin/bash
# T25 r50: og-delta prefire of the 2 burst r50 protocols (configs/Sjostrom2003b_L5TTPC_L5TTPC.yaml) on subset24; copy of run_prefire_sj03.sh.
# Array task i runs pairs i, i+3, ... (8 pairs x 2 protocols). Memory: 16 concurrent tasks x ~7.5 GB (sv_prefire 22005844: 346 GB / 46 workers) = 120 GB + ~15%.
#SBATCH --job-name=sj03r50_prefire
#SBATCH --account=rrg-emuller
#SBATCH --array=0-2
#SBATCH --cpus-per-task=16
#SBATCH --mem=76G
#SBATCH --time=00:30:00
# seff 22068340: 12-15 min wall, 61 GB MaxRSS, 62-68% CPU eff on 16 -> 76G / 0:30.
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/sj03r50_prefire_%A_%a.out
set -u
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
IDX=refitting_results/fitting/n120/seed20262009/index_Sjostrom2003b_L5TTPC_L5TTPC.csv
S=refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
PAIRS=$(python -c "
p=open('glusynapse_v2/subset24_pairs.txt').read().strip().split(',')
print(','.join(p[${SLURM_ARRAY_TASK_ID:-0}::3]))")
PROTOS=$(python -c "
import pandas as pd; d=pd.read_csv('$IDX'); print(','.join(dict.fromkeys(d.protocol_id)))")
echo "pairs $PAIRS"; echo "protocols $PROTOS"
python -u run_de_fit2_pool.py --cooker --param-hash delta-prefire --results-dir $S \
    --cache cpre_cpost_cache/sabrina_n120_delta.pkl --circuit-config data/dhuruva_delta_circuit_config.json \
    --trace-vars cai_CR,shaft_cai,ica_VDCC --skip-existing --pairs $PAIRS --workers 16 --protocols $PROTOS
