#!/bin/bash
# T29: og-delta BCL prefire (hash delta-prefire-vseg-rs, trace vars incl. v_seg) of the 5 Ebner L23PC->L5TTPC protocols
# on the 120 pairs, re-split edges (new circuit config + cache). ~466 workdirs, pilot 22087763 (v_seg, 8 wd) ran in.
# Sizing: no-v_seg 22064522 (same workdirs, 32 workers): 41 min, MaxRSS 195 GiB (6.1 GB/worker, workers are reused, so
# higher than the 3.9 GB/worker of the 1-task-per-worker pilot). v_seg pilot: ~+17% CPU-s/wd vs 169 s -> ~48-56 min
# MEASURED 22089155: 37:18, MaxRSS 189.7 GiB of 224G, CPU eff 93%; 349 pkls (343+6 pilot), 117 guardrail fails. Next: 01:00.
# -> time 01:45 (56 min + 50%); mem 224G (195 GiB measured + ~15%, v_seg trace adds little). ~56 CPU-h allocated.
# Expect ~117 guardrail failures (dropped APs in Letzkus 3ap bursts), no --allow-late-spikes.
#SBATCH --job-name=prefire_l23l5_rs
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=32
#SBATCH --mem=224G
#SBATCH --time=01:45:00
#SBATCH --output=/home/dhuruva/ws/plastyfire/logs/prefire_l23l5_rs_%j.out
set -u
ROOT=/home/dhuruva/ws/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
R=refitting_results/fitting/n120/seed20262009
IDX=$R/index_Ebner2019_L23PC_L5TTPC.csv
PROTOS=$(python -c "
import pandas as pd; d=pd.read_csv('$IDX'); print(','.join(dict.fromkeys(d.protocol_id)))")
PAIRS=$(python -c "
import pandas as pd; d=pd.read_csv('$IDX'); print(','.join(dict.fromkeys(d.pair)))")
python -u run_de_fit2_pool.py --cooker --param-hash delta-prefire-vseg-rs --results-dir $R/Ebner2019_L23PC_L5TTPC/simulations \
    --cache cpre_cpost_cache/ebner_l23l5_delta_rs.pkl --circuit-config data/dhuruva_delta_l23l5_circuit_config.json \
    --trace-vars cai_CR,shaft_cai,ica_VDCC,v_seg --skip-existing --pairs $PAIRS --workers 32 --protocols $PROTOS
