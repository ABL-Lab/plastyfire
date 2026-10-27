#!/bin/bash
# T24: find_pairs + workdirs for configs/Ebner2019_L23PC_L5TTPC.yaml. The first attempt (22033992, 24G, 3 h)
# was OOM-killed at 24G (MaxRSS 25.2 GB, 30 pool workers) and then hung on the broken pool until the time
# limit; single_cells pkls of 94 cells are cached and reused. 22043887 (64G) ok, MaxRSS 36.4 GB.
# Reruns with --pairs-from reuse the cached stimuli: 22064506 took 1:46, MaxRSS 4.3 GB -> 6G.
#SBATCH --job-name=simwriter_l23l5
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=01:00:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/simwriter_l23l5_%j.out
#SBATCH --error=/project/rrg-emuller/dhuruva/plastyfire/logs/simwriter_l23l5_%j.err
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT/plastyfire && source $ROOT/.venv/bin/activate && export PYTHONPATH=$ROOT
python -u simwriter.py --config $ROOT/configs/Ebner2019_L23PC_L5TTPC.yaml \
    --pairs-from $ROOT/refitting_results/fitting/n120/seed20262009/index_Ebner2019_L23PC_L5TTPC.csv
