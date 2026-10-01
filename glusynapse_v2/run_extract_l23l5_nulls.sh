#!/bin/bash
#SBATCH --job-name=vca_extract_l23_nulls
#SBATCH --account=rrg-emuller
# Letzkus null controls prefire (job 22123818 + pilot 22123319, hash delta-prefire-vseg-rs-nulls, 177+ pkls) -> NEW npz files
# (names <pair>__<protocol>.npz, so the existing 349 files in ebner_l23l5_delta-prefire-vseg-rs are not touched; same
# extractor, so same fields incl. vev/cev). Cache and edges as run_extract_l23l5_rs_vca.sh.
# Sizing: 22090550: 94 files in ~3.3 min on 4 workers, MaxRSS <= 73.8 GB (peak scales with workers, not records).
#   ~180 files -> ~6.3 min + 50% = ~9.5 min -> 0:15 (--skip-existing, resumable); 95G (73.8 + 25%). 4 CPU = worker count.
#   MEASURED 22126528: 1:51 elapsed, 78.86 GB MaxRSS (121 written, 59 hit inode quota) -> next 99G, 0:15.
#   MEASURED 22128395 (remaining 59 files, --skip-existing): 0:52, 30.1 GB MaxRSS. A full rerun keeps 99G (peak 78.9 GB).
#SBATCH --cpus-per-task=4
#SBATCH --mem=95G
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/vca_extract_l23l5_nulls_%j.out
set -euo pipefail
source .venv/bin/activate
PF=$PWD; V2=$PF/glusynapse_v2
export PLASTYFIRE_EDGES=$PF/data/dhuruva_modified_edges_l23l5.h5
R=$PF/refitting_results/fitting/n120/seed20262009
SIMS=$R/Ebner2019_L23PC_L5TTPC/simulations
OUT=$V2/extracted/ebner_l23l5_delta-prefire-vseg-rs
python -u $V2/extract_v2.py --param-hash delta-prefire-vseg-rs-nulls --sims $SIMS --cache $PF/cpre_cpost_cache/ebner_l23l5_delta_rs.pkl \
    --index-csv $R/index_Ebner2019_L23PC_L5TTPC_nulls.csv --protocols letzkus_nopost,letzkus_3ap_200hz_dt-500ms \
    --window --workers 4 --skip-existing --out $OUT
echo "nulls npz: $(ls $OUT | grep -c -E '__letzkus_(nopost|3ap_200hz_dt-500ms)\.npz')"
