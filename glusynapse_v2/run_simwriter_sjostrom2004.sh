#!/bin/bash
# T31: workdirs for configs/Sjostrom2004_L5TTPC_L5TTPC.yaml (1 protocol sjostrom04_dltd_step250ms, 24 subset24 pairs, prefire_only).
# New step key (250 ms, below 0.02 nA): 24 post gids searched, ~9-10 single-cell sims each (~5 s) = ~240 sims / 8 workers = ~3 min.
# Sizing basis 22086838 (sj07, 48 gids x 9 sims, 8 workers): 5:21, MaxRSS 1.20 GB, CPU eff 80% -> half the sims: 0:15 (3 min + 50%),
# 1.2 GB + 25% = 1.5 -> 2G, 8 CPU = worker count. Record seff here after the run.
# MEASURED 22127910: 1:39, MaxRSS 1.43 GB, CPU eff 81% on 8; 0 "no valid stimulus". Next: 2G, 0:15.
#SBATCH --job-name=simwriter_sj04
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=8
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/simwriter_sj04_%j.out
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT/plastyfire && source $ROOT/.venv/bin/activate && export PYTHONPATH=$ROOT
python -u simwriter.py --config $ROOT/configs/Sjostrom2004_L5TTPC_L5TTPC.yaml --workers 8
