#!/bin/bash
# T22 pilot: pipette groups (EPSP calibration) + simwriter for 3 L2/3 post cells of
# configs/Ebner2019_L23PC_L23PC_basal.yaml. Groups are cached in single_cells/pipette_<gid>.json,
# so the full run reuses them.
#SBATCH --job-name=nevian_pilot
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=/home/dhuruva/ws/plastyfire/logs/nevian_pilot_%j.out
ROOT=/home/dhuruva/ws/plastyfire
cd $ROOT/plastyfire && source $ROOT/.venv/bin/activate && export PYTHONPATH=$ROOT
python -u simwriter.py --config $ROOT/configs/Ebner2019_L23PC_L23PC_basal.yaml --pilot ${PILOT:-3}
