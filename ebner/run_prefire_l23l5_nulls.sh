#!/bin/bash
# Full prefire of the Letzkus nulls (letzkus_nopost 120 wd, letzkus_3ap_200hz_dt-500ms 115 wd) under hash
# delta-prefire-vseg-rs-nulls (new; --skip-existing keeps the pilot pkls). 1ap offset left out (decision).
# Sizing: pilot 22123319: ~220 s/wd, 3.9 GB/worker with 1 task per worker. Full run reuses workers, which raises MaxRSS:
# 22089155 (32 workers, ~14.6 wd/worker) measured 6.1 GB/worker. Here ~14.7 wd/worker at 16 workers -> use 6.1 GB/worker:
# 16 x 6.1 = 98 GB + 25% = 125G. Time: 235 x 220 s / 16 = 54 min, + 50% = 81 min -> 1:30. CPU-h allocated 16 x 1.5 = 24 (<30).
# MEASURED 22123818: 39:18 wall, CPU eff 91%, MaxRSS 124.99 GB of 125G (AT THE LIMIT) -> next time 150G with 16 workers; 177 ok, 55 failed (all 3ap_dt-500ms spike-count guardrail: 201/300 spikes or 400-500 overfiring; nopost 0 failures).
#SBATCH --job-name=prefire_l23l5_nulls
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=16
#SBATCH --mem=125G
#SBATCH --time=01:30:00
#SBATCH --output=/home/dhuruva/ws/plastyfire/logs/prefire_l23l5_nulls_%j.out
set -u
ROOT=/home/dhuruva/ws/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
R=refitting_results/fitting/n120/seed20262009
IDX=$R/index_Ebner2019_L23PC_L5TTPC_nulls.csv
PAIRS=$(awk -F, 'NR>1 && ($7=="letzkus_nopost"||$7=="letzkus_3ap_200hz_dt-500ms"){print $1}' $IDX | awk '!s[$0]++' | paste -sd,)
python -u run_de_fit2_pool.py --cooker --param-hash delta-prefire-vseg-rs-nulls --results-dir $R/Ebner2019_L23PC_L5TTPC/simulations \
    --cache cpre_cpost_cache/ebner_l23l5_delta_rs.pkl --circuit-config data/dhuruva_delta_l23l5_circuit_config.json \
    --trace-vars cai_CR,shaft_cai,ica_VDCC,v_seg --skip-existing --pairs $PAIRS --workers 16 \
    --protocols letzkus_nopost,letzkus_3ap_200hz_dt-500ms
