#!/bin/bash
# T29 pilot: v_seg memory/time for the re-split L23->L5 prefire (hash delta-prefire-vseg-rs), 2 pairs x all protocols
# (~8-10 workdirs), 8 workers. Guess for the pilot only: 8 x 6.1 GB/worker (22064522: 195 GiB / 32, no v_seg) x 1.25 = 61 -> 64G.
# MEASURED 22087763: 5:59, 8 workers, 8 wd (2 letzkus 3ap guardrail fails), MaxRSS 31.3 GiB (3.9 GB/worker), CPU eff 55% (tail-bound).
# The pilot workdirs are kept (--skip-existing in the full run).
#SBATCH --job-name=prefire_rs_pilot
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=/home/dhuruva/ws/plastyfire/logs/prefire_l23l5_rs_pilot_%j.out
set -u
ROOT=/home/dhuruva/ws/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
R=refitting_results/fitting/n120/seed20262009
IDX=$R/index_Ebner2019_L23PC_L5TTPC.csv
PROTOS=$(python -c "
import pandas as pd; d=pd.read_csv('$IDX'); print(','.join(dict.fromkeys(d.protocol_id)))")
PAIRS=$(python -c "
import pandas as pd; d=pd.read_csv('$IDX'); print(','.join(list(dict.fromkeys(d.pair))[:2]))")
python -u run_de_fit2_pool.py --cooker --param-hash delta-prefire-vseg-rs --results-dir $R/Ebner2019_L23PC_L5TTPC/simulations \
    --cache cpre_cpost_cache/ebner_l23l5_delta_rs.pkl --circuit-config data/dhuruva_delta_l23l5_circuit_config.json \
    --trace-vars cai_CR,shaft_cai,ica_VDCC,v_seg --skip-existing --pairs $PAIRS --workers 8 --protocols $PROTOS
