#!/bin/bash
# Synapse path distances + unitary EPSP, 10-90% rise time, all_protocols flag per L23PC->L5TTPC pair
# (ebner/pair_geometry.py). 22046061 (10 spikes): 19 min, MaxRSS 17.9 GB -> 23G; 30 spikes now.
#SBATCH --job-name=pair_geometry
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=16
#SBATCH --mem=23G
#SBATCH --time=02:00:00
#SBATCH --output=/home/dhuruva/ws/plastyfire/logs/pair_geometry_%j.out
ROOT=/home/dhuruva/ws/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
R=refitting_results/fitting/n120/seed20262009
python -u ebner/pair_geometry.py $R/index_Ebner2019_L23PC_L5TTPC.csv ebner/pair_geometry_L23PC_L5TTPC.csv \
 --epsp $R/Ebner2019_L23PC_L5TTPC/simulations 16
