#!/bin/bash
# FULL og-delta BCL prefire (hash delta-prefire-vseg-rs, v_seg traces) of the 13 Zilberter protocols, 120 pairs = 1555 workdirs
# (26 from pilot 22128323 skipped). l23l23 edges + cache.
# Sizing from pilot 22128323 (8 workers, 26 wd): 26/26 ok, 0 guardrail fails, 5:00, MaxRSS 32.8 GB = 4.1 GB/worker,
# CPU eff 90%, mean 85 s/wd (38-149 s). Workers are reused in a full run: 22089155 vs pilot 22087763 went 3.9 -> 5.9 GB/worker
# (x1.5) -> 32 workers x 4.1 x 1.5 = 197 GB + 25% -> 248G. Time 1529 wd x 85 s / 32 = 68 min + 50% -> 1:45 (~37 CPU-h).
#SBATCH --job-name=prefire_l23l23
#SBATCH --account=rrg-emuller
# MEASURED 22129908 (1529 wd, 32 workers): 1:13:57, 248.0 GB MaxRSS = AT the 248G limit, CPU 97%; 1501 ok, 28 guardrail fails (post-spike count 199/200 or 80/40, 6 pairs). Next full run 310G, 1:45.
#SBATCH --cpus-per-task=32
#SBATCH --mem=248G
#SBATCH --time=01:45:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/prefire_l23l23_%j.out
set -u
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
R=refitting_results/fitting/n120/seed20262009
IDX=$R/index_Zilberter2009_L23PC_L23PC.csv
PROTOS=$(tail -n +2 $IDX | awk -F, '{print $7}' | awk '!s[$0]++' | paste -sd,)
PAIRS=$(tail -n +2 $IDX | cut -d, -f1 | awk '!s[$0]++' | paste -sd,)
python -u run_de_fit2_pool.py --cooker --param-hash delta-prefire-vseg-rs --results-dir $R/Zilberter2009_L23PC_L23PC/simulations \
    --cache cpre_cpost_cache/zilberter_l23l23_delta_rs.pkl --circuit-config data/dhuruva_delta_l23l23_circuit_config.json \
    --trace-vars cai_CR,shaft_cai,ica_VDCC,v_seg --skip-existing --pairs $PAIRS --workers 32 --protocols $PROTOS
