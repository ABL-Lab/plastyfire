#!/bin/bash
# T27: full-protocol BCL (C01 + pairing + C02, soma-measured EPSP ratio, no basis) with delta-cooker params on the
# Ebner L23PC->L5TTPC pairs (og-delta emodel). Hash delta-cooker-full. 349 workdirs (the 117 that failed the
# delta-prefire post-spike guardrail are excluded, ebner/prefire_failed_l23l5.csv).
# Fastforward per protocol (snap rho after induction, before the first C02 pulse): sjostrom 390 s (induction ends
# 381.1 s, C02 starts 391 s), letzkus 340.9 s (induction ends 340.01 s, C02 starts 341 s). The default 280 s of the
# Sabrina runs would snap rho in the middle of the Ebner induction.
# Sizing basis (full_cooker 21723104, Sabrina, 805 wd, 60 workers: 1:55:39, MaxRSS 209.8 GB = 3.5 GB/worker,
# ~484 CPU-s/wd at tstop 520 s): Ebner tstop 630 s (sjostrom, 119 wd) / 580 s (letzkus, 230 wd) -> ~54 CPU-h,
# 60 workers -> ~58 min, x1.5 = 1:30; mem 60 x 3.5 GB x 1.25 = 265G.
# MEASURED (job 22071891): 348 ok + 1 fail, 55:08 elapsed, batch MaxRSS 236.6 GB (3.9 GB/worker; 265G request = 89%), TotalCPU 2d0h08 = 48.1 CPU-h,
# CPU efficiency 87%. ~400-850 s per workdir. Next time: mem 60 x 3.9 x 1.25 = 295G is not needed for the same set; time 0:55 x1.5 -> 1:30 was right.
#SBATCH --job-name=full_l23l5
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=60
#SBATCH --mem=265G
#SBATCH --time=01:30:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/full_l23l5_%j.out
set -u
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
R=refitting_results/fitting/n120/seed20262009
IDX=$R/index_Ebner2019_L23PC_L5TTPC.csv
PROTOS=$(python -c "
import pandas as pd; d=pd.read_csv('$IDX'); print(','.join(dict.fromkeys(d.protocol_id)))")
PAIRS=$(python -c "
import pandas as pd; d=pd.read_csv('$IDX'); print(','.join(dict.fromkeys(d.pair)))")
python -u run_de_fit2_pool.py --params-json fit_results/delta-cooker.json --param-hash delta-cooker-full --full \
    --results-dir $R/Ebner2019_L23PC_L5TTPC/simulations --cache cpre_cpost_cache/ebner_l23l5_delta.pkl \
    --circuit-config data/dhuruva_delta_circuit_config.json --skip-existing --workers ${WORKERS:-60} \
    --pairs $PAIRS --protocols $PROTOS --exclude-csv ebner/prefire_failed_l23l5.csv \
    --fastforward-map '{"sjostrom_50hz_dt+10ms": 390000, "letzkus_1ap_dt+10ms": 340900, "letzkus_3ap_200hz_dt+10ms": 340900, "letzkus_3ap_200hz_dt-10ms": 340900}' \
    ${DRY:-}
