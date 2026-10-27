#!/bin/bash
#SBATCH --job-name=create_dhuruva_edges
#SBATCH --account=rrg-emuller
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/create_dhuruva_edges_%j.out
#SBATCH --error=/project/rrg-emuller/dhuruva/plastyfire/logs/create_dhuruva_edges_%j.err

cd /project/rrg-emuller/dhuruva/plastyfire
source /project/rrg-emuller/dhuruva/plastyfire/.venv/bin/activate

python -u create_dhuruva_edges.py --k-u 0.2 --k-gsyn 2.0 --pot-fraction 0.5 --no-patch-configs
