#!/bin/bash
# Writes sjostrom04_dltd_step250ms_pre+100c (step_below_nA 0.10; existing ids unchanged, cached keys). Needs the new below0.1nA
# MEASURED 22136156: 1:47, 1.08 GB
# step calibration for the pairs' post cells. Basis: 22127910 (1 new key): 1:39, 1.43 GB, 81% CPU on 8 -> 2G, 0:15, 8 workers.
#SBATCH --job-name=simwriter_sj04_c
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=8
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/simwriter_sj04_c_%j.out
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT/plastyfire && source $ROOT/.venv/bin/activate && export PYTHONPATH=$ROOT
python -u simwriter.py --config $ROOT/configs/Sjostrom2004_L5TTPC_L5TTPC.yaml --workers 8
