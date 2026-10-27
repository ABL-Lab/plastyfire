#!/bin/bash
#SBATCH --job-name=basis_sabrina_n120
#SBATCH --account=rrg-emuller
#SBATCH --time=06:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=60
#SBATCH --output=/home/dhuruva/ws/plastyfire/logs/basis_sabrina_n120_%j.out
#SBATCH --error=/home/dhuruva/ws/plastyfire/logs/basis_sabrina_n120_%j.err

cd /home/dhuruva/ws/plastyfire
source /home/dhuruva/ws/plastyfire/.venv/bin/activate

# 5 pairs at a time x 12 trial processes each = 60 cores
python -u submit_basis_edges.py --execution-mode cpu --workers 5 \
    --stdp-name Sabrina_L5TTPC_L5TTPC_STDP \
    --output-dir basis_results_edges_sabrina_n120 \
    --circuit-config data/dhuruva_modified_ion_channels_circuit_config.json \
    --skip-existing
