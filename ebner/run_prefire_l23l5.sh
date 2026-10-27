#!/bin/bash
# T24: og-delta BCL prefire (hash delta-prefire, trace vars as sj03) of the 5 Ebner L23PC->L5TTPC protocols on
# the 120 pairs (~465 workdirs, 100-150 s sims). Memory: sj03_prefire 22057545 (32 workers) MaxRSS 154 GB -> 195G.
#SBATCH --job-name=prefire_l23l5
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=32
#SBATCH --mem=195G
#SBATCH --time=05:00:00
#SBATCH --output=/home/dhuruva/ws/plastyfire/logs/prefire_l23l5_%j.out
set -u
ROOT=/home/dhuruva/ws/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
R=refitting_results/fitting/n120/seed20262009
IDX=$R/index_Ebner2019_L23PC_L5TTPC.csv
PROTOS=$(python -c "
import pandas as pd; d=pd.read_csv('$IDX'); print(','.join(dict.fromkeys(d.protocol_id)))")
PAIRS=$(python -c "
import pandas as pd; d=pd.read_csv('$IDX'); print(','.join(dict.fromkeys(d.pair)))")
echo "protocols $PROTOS"
python -u run_de_fit2_pool.py --cooker --param-hash delta-prefire --results-dir $R/Ebner2019_L23PC_L5TTPC/simulations \
    --cache cpre_cpost_cache/ebner_l23l5_delta.pkl --circuit-config data/dhuruva_delta_circuit_config.json \
    --trace-vars cai_CR,shaft_cai,ica_VDCC --skip-existing --pairs $PAIRS --workers 32 --protocols $PROTOS
