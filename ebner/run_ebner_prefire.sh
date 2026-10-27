#!/bin/bash
# T5: BCL prefire traces for the 28 Ebner 2019 protocols (configs/Ebner2019_L5TTPC_L5TTPC.yaml).
# Array task i runs every 8th pair (pairs i, i+8, ...), so the 8 tasks each get 15 pairs.
# MODE=other: 26 non-Letzkus protocols; MODE=letz: the 2 Letzkus 200 Hz protocols with
# --allow-late-spikes; MODE=both (default): both pools at once. Memory, not CPU, limits the
# worker count: 3 trace vars over 600-1000 s cost ~10-50 GB per task (19 synapses worst case).
#SBATCH --job-name=ebner_prefire
#SBATCH --account=rrg-emuller
#SBATCH --array=0-7
#SBATCH --cpus-per-task=48
#SBATCH --mem=740G
#SBATCH --time=12:00:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/ebner_prefire_%A_%a.out
set -u
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
IDX=refitting_results/fitting/n120/seed20262009/index_Ebner2019_L5TTPC_L5TTPC.csv
S=refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
PAIRS=$(python -c "
import pandas as pd; d=pd.read_csv('$IDX'); p=sorted({'%d-%d'%(a,b) for a,b in zip(d.pregid,d.postgid)})
print(','.join(p[${SLURM_ARRAY_TASK_ID:-0}::8]))")
PROTOS=$(python -c "
import pandas as pd; d=pd.read_csv('$IDX'); print(','.join(dict.fromkeys(d.protocol_id)))")
OTHER=$(echo $PROTOS | tr , '\n' | grep -v '^letzkus' | paste -sd,)
LETZ=$(echo $PROTOS | tr , '\n' | grep '^letzkus' | paste -sd,)
COMMON="--cooker --param-hash delta-prefire --results-dir $S --cache cpre_cpost_cache/sabrina_n120_delta.pkl \
 --circuit-config data/dhuruva_delta_circuit_config.json --trace-vars cai_CR,shaft_cai,ica_VDCC --skip-existing --pairs $PAIRS"
MODE=${MODE:-both}
[[ $MODE != letz ]] && python -u run_de_fit2_pool.py $COMMON --workers ${WORKERS_OTHER:-40} --protocols $OTHER &
[[ $MODE != other ]] && python -u run_de_fit2_pool.py $COMMON --workers ${WORKERS_LETZ:-6} --protocols $LETZ --allow-late-spikes &
wait
