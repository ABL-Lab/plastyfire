#!/bin/bash
# Pilot prefire of the Letzkus null controls (letzkus_nopost, letzkus_3ap_200hz_dt-500ms; 1ap offset dropped by decision)
# on 2 pairs that have both workdirs. New hash delta-prefire-vseg-rs-nulls, same cache/circuit config/trace vars as
# run_prefire_l23l5_rs.sh (22089155), so no existing pkl is touched. 4 workdirs, 4 workers (1 task per worker).
# Sizing: 22087763 pilot (8 wd, 1 task/worker): 3.9 GB/worker; ~170 s/wd (+17% v_seg) -> 4 x 3.9 = 15.6 GB + 25% = 20G;
# wall ~4 min incl. setup -> 0:15 (measured + 50%, rounded up). Record seff here after the run.
# MEASURED 22123319: 3:41 wall, 4 workers, MaxRSS 15.74 GB (3.9 GB/worker), CPU eff 61%, ~220 s/wd; 3 ok, 1 guardrail drop (10024-203348 3ap: 500 post spikes vs 300, calibration over-firing, per-pair drop accepted).
#SBATCH --job-name=prefire_nulls_pilot
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=00:15:00
#SBATCH --output=/home/dhuruva/ws/plastyfire/logs/prefire_nulls_pilot_%j.out
set -eu
ROOT=/home/dhuruva/ws/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
R=refitting_results/fitting/n120/seed20262009
IDX=$R/index_Ebner2019_L23PC_L5TTPC_nulls.csv
PAIRS=$(awk -F, 'NR>1 && ($7=="letzkus_nopost"||$7=="letzkus_3ap_200hz_dt-500ms"){c[$1]++} END{for(p in c) if(c[p]==2) print p}' $IDX | sort | head -2 | paste -sd,)
echo "pilot pairs: $PAIRS"
python -u run_de_fit2_pool.py --cooker --param-hash delta-prefire-vseg-rs-nulls --results-dir $R/Ebner2019_L23PC_L5TTPC/simulations \
    --cache cpre_cpost_cache/ebner_l23l5_delta_rs.pkl --circuit-config data/dhuruva_delta_l23l5_circuit_config.json \
    --trace-vars cai_CR,shaft_cai,ica_VDCC,v_seg --skip-existing --pairs $PAIRS --workers 4 \
    --protocols letzkus_nopost,letzkus_3ap_200hz_dt-500ms
ls $R/Ebner2019_L23PC_L5TTPC/simulations/*/*/*delta-prefire-vseg-rs-nulls* 2>/dev/null | wc -l
