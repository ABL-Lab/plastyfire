#!/bin/bash
#SBATCH --job-name=v2_readout
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=12
#SBATCH --mem=32g
#SBATCH --time=01:30:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/v2_readout_%j.out
source .venv/bin/activate
python glusynapse_v2/check_readout.py --pre-gid 180351 --post-gid 198084 --workers 12
python glusynapse_v2/check_readout.py --pre-gid 181116 --post-gid 188518 --workers 12
