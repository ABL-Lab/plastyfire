#!/bin/bash
# T25 step 6: windowed extraction of the Sjostrom 2003 prefire traces (run_prefire_sj03.sh) into an
# own npz dir extracted/sj03_delta-prefire (so fits running meanwhile never see a partial set).
# ~10 GB per worker on 40M-sample traces (run_ebner_pipeline.sh), 16 workers -> 160 G + 25%.
#SBATCH --job-name=sj03_extract
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=03:00:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/sj03_extract_%j.out
set -euo pipefail
source .venv/bin/activate
PF=$PWD; V2=$PF/glusynapse_v2
SIMS=$PF/refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
python -u $V2/extract_v2.py --param-hash delta-prefire --index-csv $SIMS/../../index_Sjostrom2003_L5TTPC_L5TTPC.csv \
    --window --workers 16 --skip-existing --out $V2/extracted/sj03_delta-prefire --pairs $(cat $V2/subset24_pairs.txt)
ls $V2/extracted/sj03_delta-prefire | wc -l
