#!/bin/bash
#SBATCH --job-name=v2_nrntests
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/v2_nrntests_%j.out
# NEURON-level gates: mod behaviour (single synapse) and offline model vs mod.
set -euo pipefail
source .venv/bin/activate

python -u glusynapse_v2/tests/test_v2_single_synapse.py
python -u glusynapse_v2/tests/test_offline_gate.py
