#!/bin/bash
# Workdirs for the 3 Letzkus null controls (configs/Ebner2019_L23PC_L5TTPC_nulls.yaml) on the existing 120 pairs.
# Same sims folder, index_label _nulls (own index csv + yaml), only the 3 new protocol ids are written: existing workdirs untouched.
# Index csv is written at seed20262009/index_Ebner2019_L23PC_L5TTPC_nulls.csv (not inside the label folder).
# MEASURED 22123227: 20 s, MaxRSS 315 MB, 6% CPU on 4 CPUs -> next time 1 CPU, 1G, 0:15. (Sizing basis was 22064506: 1:46, 4.3 GB.)
#SBATCH --job-name=simwriter_l23l5_nulls
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:15:00
#SBATCH --output=/home/dhuruva/ws/plastyfire/logs/simwriter_l23l5_nulls_%j.out
#SBATCH --error=/home/dhuruva/ws/plastyfire/logs/simwriter_l23l5_nulls_%j.err
ROOT=/home/dhuruva/ws/plastyfire
cd $ROOT/plastyfire && source $ROOT/.venv/bin/activate && export PYTHONPATH=$ROOT
python -u simwriter.py --config $ROOT/configs/Ebner2019_L23PC_L5TTPC_nulls.yaml
