#!/bin/bash
#SBATCH --job-name=v2_extract_markram
#SBATCH --account=rrg-emuller
#SBATCH --time=02:00:00
#SBATCH --mem=180G
#SBATCH --cpus-per-task=30
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/v2_extract_markram_%j.out
set -euo pipefail
PF=/project/rrg-emuller/dhuruva/plastyfire
cd $PF; source .venv/bin/activate; export PYTHONPATH=$PF
python -u glusynapse_v2/extract_v2.py --param-hash delta-cooker --workers 30 \
    --out $PF/glusynapse_v2/extracted/markram_delta-cooker
