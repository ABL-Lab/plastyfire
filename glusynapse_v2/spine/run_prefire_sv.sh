#!/bin/bash
# delta-sv prefire (spine-VDCC variant, glusynapse_v2/spine/delta_sv.json) on subset24: the Ebner 2019 protocols
# (as ebner/run_ebner_prefire.sh) then the Markram set (as delta-prefire-tr), new hashes, so no delta output is touched.
# Array task i runs pairs i, i+3, ... (8 pairs each). Memory from ebner_prefire 21902422 (265 GB MaxRSS, 46 workers).
#SBATCH --job-name=sv_prefire
#SBATCH --account=rrg-emuller
#SBATCH --array=0-2
#SBATCH --cpus-per-task=48
#SBATCH --mem=330G
#SBATCH --time=06:00:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/sv_prefire_%A_%a.out
set -u
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
IDX=refitting_results/fitting/n120/seed20262009/index_Ebner2019_L5TTPC_L5TTPC.csv
S=refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
PAIRS=$(python -c "
p=open('glusynapse_v2/subset24_pairs.txt').read().strip().split(',')
print(','.join(p[${SLURM_ARRAY_TASK_ID:-0}::3]))")
PROTOS=$(python -c "
import pandas as pd; d=pd.read_csv('$IDX'); print(','.join(dict.fromkeys(d.protocol_id)))")
OTHER=$(echo $PROTOS | tr , '\n' | grep -v '^letzkus' | paste -sd,)
LETZ=$(echo $PROTOS | tr , '\n' | grep '^letzkus' | paste -sd,)
GLOB=$(python -c "import json; print(json.dumps(json.load(open('glusynapse_v2/spine/delta_sv.json'))['globals']))")
COMMON="--cooker --cache cpre_cpost_cache/sabrina_n120_delta-sv.pkl --circuit-config data/dhuruva_delta_circuit_config.json \
 --trace-vars cai_CR,shaft_cai,ica_VDCC --skip-existing --pairs $PAIRS"
echo "pairs $PAIRS"; echo "globals $GLOB"
python -u run_de_fit2_pool.py $COMMON --param-hash delta-sv-prefire --glusyn-globals "$GLOB" --results-dir $S \
    --workers 40 --protocols $OTHER &
python -u run_de_fit2_pool.py $COMMON --param-hash delta-sv-prefire --glusyn-globals "$GLOB" --results-dir $S \
    --workers 6 --protocols $LETZ --allow-late-spikes &
wait
python -u run_de_fit2_pool.py $COMMON --param-hash delta-sv-prefire-tr --glusyn-globals "$GLOB" --results-dir $S --workers 46
