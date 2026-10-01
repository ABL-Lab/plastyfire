#!/bin/bash
#SBATCH --job-name=vca_extract_sj04
#SBATCH --account=rrg-emuller
# Sjostrom 2004 dLTD prefire (hash delta-prefire-vseg, 24 pairs) -> extracted/sj04_delta-prefire-vca. Clone of run_extract_sj07_vca.sh.
# Sizing: 22126528 (4 workers): 121 files in 1:51, MaxRSS 78.86 GB; peak scales with workers, not files -> 4 workers, 99G, 0:15.
#SBATCH --cpus-per-task=4
# MEASURED 22128325 (24 files, 4 workers): 1:09, 72.4 GB MaxRSS -> next 91G, 0:15.
#SBATCH --mem=99G
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/vca_extract_sj04_%j.out
set -euo pipefail
source .venv/bin/activate
PF=$PWD; V2=$PF/glusynapse_v2
SIMS=$PF/refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
PAIRS=$(cat $V2/subset24_pairs.txt)
python -u $V2/extract_v2.py --param-hash delta-prefire-vseg --index-csv $SIMS/../../index_Sjostrom2004_L5TTPC_L5TTPC.csv \
    --protocols sjostrom04_dltd_step250ms --window --workers 4 --skip-existing --out $V2/extracted/sj04_delta-prefire-vca --pairs $PAIRS
echo "sj04 $(ls $V2/extracted/sj04_delta-prefire-vca | wc -l)"
