#!/bin/bash
# delta-sv extraction (after run_prefire_sv.sh): subset24, same settings as run_ebner_pipeline.sh step [1], sv cache.
# Memory: Ebner windowed extraction ~10 GB per worker (run_ebner_pipeline.sh), 16 workers.
#SBATCH --job-name=sv_extract
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=03:00:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/sv_extract_%j.out
set -euo pipefail
source .venv/bin/activate
PF=$PWD; V2=$PF/glusynapse_v2
SIMS=$PF/refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
C=$PF/cpre_cpost_cache/sabrina_n120_delta-sv.pkl; PAIRS=$(cat $V2/subset24_pairs.txt)
python -u $V2/extract_v2.py --param-hash delta-sv-prefire --index-csv $SIMS/../../index_Ebner2019_L5TTPC_L5TTPC.csv \
    --window --workers 16 --skip-existing --cache $C --pairs $PAIRS --out $V2/extracted/ebner_delta-sv
python -u $V2/extract_v2.py --param-hash delta-sv-prefire-tr --workers 16 --skip-existing --cache $C --pairs $PAIRS \
    --out $V2/extracted/markram_delta-sv
echo "=== done $(date)"
