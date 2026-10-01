#!/bin/bash
# t_drive 3 pilot: windowed extraction (with vev) of the delta-prefire-vseg pilot, new dirs, 2 pairs x 13 protocols.
# Sizing: 22057611 (16 workers, 210 GB = 13 GB/worker) and 22068341 (6 workers, 77 GB); v_seg trace adds ~x1.35:
# 4 workers x 17.5 GB = 70 GB + 25% = 90G. Time: 22068341 took 1 min for 48 npz; 26 npz + crossing detection -> 0:15.
#SBATCH --job-name=vseg_extract
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=4
#SBATCH --mem=90G
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/vseg_extract_%j.out
set -euo pipefail
source .venv/bin/activate
PF=$PWD; V2=$PF/glusynapse_v2
SIMS=$PF/refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
PAIRS=${PAIRS:-180351-198084,181455-195199}; W=${WORKERS:-4}
python -u $V2/extract_v2.py --param-hash delta-prefire-vseg --index-csv $SIMS/../../index_Sjostrom2003_L5TTPC_L5TTPC.csv \
    --window --workers $W --skip-existing --out $V2/extracted/sj03_delta-prefire-vseg --pairs $PAIRS
python -u $V2/extract_v2.py --param-hash delta-prefire-vseg --index-csv $SIMS/../../index_Sjostrom2003b_L5TTPC_L5TTPC.csv \
    --window --workers $W --skip-existing --out $V2/extracted/sj03r50_delta-prefire-vseg --pairs $PAIRS
ls $V2/extracted/sj03_delta-prefire-vseg $V2/extracted/sj03r50_delta-prefire-vseg | wc -l
