#!/bin/bash
#SBATCH --job-name=v2_jaxtest
#SBATCH --account=def-emuller      # rrg-emuller has no GPU allocation
#SBATCH --gpus-per-node=h100:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G          # preview data 26-35 G + GPU inputs ~8 G
#SBATCH --time=01:00:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/v2_jaxtest_%j.out
set -euo pipefail
source .venv/bin/activate
nvidia-smi --query-gpu=name,memory.total --format=csv
python -u glusynapse_v2/tests/test_jax_gpu.py
