#!/bin/bash
# MEASURED 22091527: 120/120 pairs, 0 synapses at floor, finished ~11:49 (started ~11:41, about 8 min); MaxRSS/seff pending (slurmdb down).
#SBATCH --job-name=v2_cache_ljp25
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=16
#SBATCH --mem=16500M
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/v2_cache_ljp25_%j.out
# c_pre/c_post cache for delta-ljp25 (glusynapse_v2/spine/delta_ljp25.json), all 120 pairs like the delta-sv cache.
# Sizing: delta-sv 22005793 (120 pairs, 16 workers, 5:50, MaxRSS 13.3 G of 13 G) -> 16.5 G, 0:15.
set -euo pipefail
source .venv/bin/activate
GLOB=$(python -c "import json; print(json.dumps(json.load(open('glusynapse_v2/spine/delta_ljp25.json'))['globals']))")
echo "globals $GLOB"
python -u precompute_cpre_cpost.py --params defit2 --results-dir refitting_results \
    --circuit-config data/dhuruva_delta_circuit_config.json --edges-h5 data/dhuruva_modified_edges.h5 \
    --output cpre_cpost_cache/sabrina_n120_delta-ljp25.pkl --workers 16 --glusyn-globals "$GLOB"
