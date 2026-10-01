#!/bin/bash
# L2/3->L2/3 (Zilberter) prefire (hash delta-prefire-vseg-rs, 1555 wd) -> extracted/zilberter_l23l23_delta-prefire-vseg-rs
# (vev, cev* fields as the L2/3->L5 -vseg-rs dir). Cache zilberter_l23l23_delta_rs.pkl, edges via PLASTYFIRE_EDGES.
# Sizing from 22126528 (same extractor, 4 workers): 121 files in 1:51, MaxRSS 78.9 GB (peak scales with workers, not
# records) -> 99G; 1555 files x 0.92 s = ~24 min + 50% -> 0:45. 4 CPU = worker count. --skip-existing, resumable.
#SBATCH --job-name=vseg_extract_l23l23
#SBATCH --account=rrg-emuller
# MEASURED 22129909 (1527 files, 4 workers): 21:09, 98.99 GB MaxRSS = AT the 99G limit, CPU 83%. Next 124G, 0:45.
#SBATCH --cpus-per-task=4
#SBATCH --mem=99G
#SBATCH --time=00:45:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/vseg_extract_l23l23_%j.out
set -euo pipefail
source .venv/bin/activate
PF=$PWD; V2=$PF/glusynapse_v2
export PLASTYFIRE_EDGES=$PF/data/dhuruva_modified_edges_l23l23.h5
R=$PF/refitting_results/fitting/n120/seed20262009
SIMS=$R/Zilberter2009_L23PC_L23PC/simulations
OUT=$V2/extracted/zilberter_l23l23_delta-prefire-vseg-rs
python -u $V2/extract_v2.py --param-hash delta-prefire-vseg-rs --sims $SIMS --cache $PF/cpre_cpost_cache/zilberter_l23l23_delta_rs.pkl \
    --index-csv $R/index_Zilberter2009_L23PC_L23PC.csv --window --workers 4 --skip-existing --out $OUT
echo "l23l23 npz: $(ls $OUT | wc -l)"
