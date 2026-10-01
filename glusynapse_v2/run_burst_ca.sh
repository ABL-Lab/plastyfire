#!/bin/bash
#SBATCH --job-name=v2_burstca
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G           # measured 3.3 GiB with 8 workers (21960073)
#SBATCH --time=01:00:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/v2_burstca_%j.out
# Burst Ca supralinearity at a pair's synapses (glusynapse_v2/burst_ca_test.py). env: PAIR, EMODEL, MODS
set -euo pipefail
source .venv/bin/activate
python -u glusynapse_v2/burst_ca_test.py --pair $PAIR --emodel ${EMODEL:-delta} --mods "${MODS:-}" \
    --workers 8 --out glusynapse_v2/results/burst_ca
