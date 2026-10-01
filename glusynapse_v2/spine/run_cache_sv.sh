#!/bin/bash
#SBATCH --job-name=v2_cache_sv
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/v2_cache_sv_%j.out
# c_pre/c_post cache for the spine-VDCC variant (delta emodel unchanged; GluSynapse GLOBALs GLOB).
# env: PAIRS, OUT, GLOB ('' = mod defaults, for the reproduction check), WORKERS.
set -euo pipefail
source .venv/bin/activate
python -u precompute_cpre_cpost.py --params defit2 --results-dir refitting_results \
    --circuit-config data/dhuruva_delta_circuit_config.json --edges-h5 data/dhuruva_modified_edges.h5 \
    --output $OUT --workers ${WORKERS:-2} --pairs $PAIRS ${GLOB:+--glusyn-globals "$GLOB"}
