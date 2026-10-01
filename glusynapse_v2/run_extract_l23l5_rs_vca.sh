#!/bin/bash
#SBATCH --job-name=vca_extract_l23
#SBATCH --account=rrg-emuller
# L2/3->L5 re-split prefire (ebner job 22089155, hash delta-prefire-vseg-rs, 349 pkls, 117 guardrail fails skipped as before)
# -> extracted/ebner_l23l5_delta-prefire-vseg-rs (with vev and cev*). Cache ebner_l23l5_delta_rs.pkl, edges via PLASTYFIRE_EDGES.
# Sizing: 22090550: ebner L5 records (94 files) in ~3.3 min on 4 workers with the slow ca_events, MaxRSS <= 73.8 GB on 4 workers.
#   ~232 files -> ~8 min slow-equivalent -> 0:15 (--skip-existing, resumable); 95G (73.8 + 25%).
#SBATCH --cpus-per-task=4
#SBATCH --mem=95G
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/vca_extract_l23l5_%j.out
set -euo pipefail
source .venv/bin/activate
PF=$PWD; V2=$PF/glusynapse_v2
export PLASTYFIRE_EDGES=$PF/data/dhuruva_modified_edges_l23l5.h5
SIMS=$PF/refitting_results/fitting/n120/seed20262009/Ebner2019_L23PC_L5TTPC/simulations
python -u $V2/extract_v2.py --param-hash delta-prefire-vseg-rs --sims $SIMS --cache $PF/cpre_cpost_cache/ebner_l23l5_delta_rs.pkl \
    --index-csv $SIMS/../../index_Ebner2019_L23PC_L5TTPC.csv --window --workers 4 --skip-existing \
    --out $V2/extracted/ebner_l23l5_delta-prefire-vseg-rs
echo "l23l5 $(ls $V2/extracted/ebner_l23l5_delta-prefire-vseg-rs | wc -l)"
