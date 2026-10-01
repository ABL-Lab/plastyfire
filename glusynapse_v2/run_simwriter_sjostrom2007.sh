#!/bin/bash
# T30: workdirs for configs/Sjostrom2007_L5TTPC_L5TTPC.yaml (4 step protocols, 24 subset24 pairs, prefire_only).
# The step amplitudes are new keys: 24 post + 24 pre gids are searched (~9 single-cell sims each, ~5 s per sim
# measured on the login node) -> 432 sims / 8 workers = ~4.5 min.
# Memory: 22043887 measured 36.4 GB MaxRSS with 30 search workers (~1.2 GB/worker); 22068339 (cached) 0.46 GB.
# 8 workers x 1.2 + 0.5 = 10.1 GB + 25% -> 13G. Time 4.5 min + 50% -> 0:15.
# MEASURED 22086838: 5:21, MaxRSS 1.20 GB (9% of 13G), CPU eff 80% on 8. Rerun with cached stimuli: 1 CPU / 1G / 0:15.
#SBATCH --job-name=simwriter_sj07
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=8
#SBATCH --mem=13G
#SBATCH --time=00:15:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/simwriter_sj07_%j.out
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT/plastyfire && source $ROOT/.venv/bin/activate && export PYTHONPATH=$ROOT
python -u simwriter.py --config $ROOT/configs/Sjostrom2007_L5TTPC_L5TTPC.yaml --workers 8
