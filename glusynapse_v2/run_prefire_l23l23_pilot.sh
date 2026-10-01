#!/bin/bash
# MEASURED 22128323: 26/26 ok, 0 guardrail fails, 5:00, MaxRSS 32.8 GB (4.1 GB/worker), CPU eff 90%; mean 85 s/wd (38-149 s).
# PILOT: og-delta BCL prefire (hash delta-prefire-vseg-rs, trace vars incl. v_seg) of the 13 Zilberter protocols on the
# first 2 pairs of the index (<= 26 workdirs), l23l23 edges + cache. Same hash as the full run, so --skip-existing reuses it.
# Sizing from the L2/3->L5 v_seg pilot 22087763 (8 workers, 8 wd, 100 s bio each): 5:59, MaxRSS 31.3 GiB (3.9 GB/worker),
# and the Sj07 pilot 22089134 (8 wd, 300 s bio): 10:43, 31.66 GB (3.96 GB/worker) -> per-worker memory does not grow with
# protocol duration: 8 workers, 40G (31.7 + 25%). Time: 40 x 5 s = 200 s bio per wd ~ 8.5 min; 26 wd / 8 workers = 4 rounds
# -> ~34 min + 50% -> 1:00. Measure here: per-wd CPU-s, GB/worker, guardrail fails per protocol -> size the full run.
#SBATCH --job-name=prefire_l23l23_pilot
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=01:00:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/prefire_l23l23_pilot_%j.out
set -u
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
R=refitting_results/fitting/n120/seed20262009
IDX=$R/index_Zilberter2009_L23PC_L23PC.csv
PROTOS=$(tail -n +2 $IDX | awk -F, '{print $7}' | awk '!s[$0]++' | paste -sd,)
PAIRS=$(tail -n +2 $IDX | cut -d, -f1 | awk '!s[$0]++' | head -2 | paste -sd,)
echo "pairs $PAIRS"; echo "protocols $PROTOS"
python -u run_de_fit2_pool.py --cooker --param-hash delta-prefire-vseg-rs --results-dir $R/Zilberter2009_L23PC_L23PC/simulations \
    --cache cpre_cpost_cache/zilberter_l23l23_delta_rs.pkl --circuit-config data/dhuruva_delta_l23l23_circuit_config.json \
    --trace-vars cai_CR,shaft_cai,ica_VDCC,v_seg --skip-existing --pairs $PAIRS --workers 8 --protocols $PROTOS
