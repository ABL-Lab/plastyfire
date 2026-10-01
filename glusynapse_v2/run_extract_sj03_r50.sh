#!/bin/bash
# T25 r50: windowed extraction of the r50 burst prefire traces into extracted/sj03r50_delta-prefire (new dir; copy of run_extract_sj03.sh).
# Old run: 200 GiB MaxRSS for 264 npz; this is 48 npz -> 6 workers x ~12 GB = 80G.
#SBATCH --job-name=sj03r50_extract
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=6
#SBATCH --mem=80G
#SBATCH --time=00:15:00
# seff 22068341: 1:04 wall, 73.3 GB of 80G, 84% CPU eff on 6 -> keep 80G, 0:15.
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/sj03r50_extract_%j.out
set -euo pipefail
source .venv/bin/activate
PF=$PWD; V2=$PF/glusynapse_v2
SIMS=$PF/refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
python -u $V2/extract_v2.py --param-hash delta-prefire --index-csv $SIMS/../../index_Sjostrom2003b_L5TTPC_L5TTPC.csv \
    --window --workers 6 --skip-existing --out $V2/extracted/sj03r50_delta-prefire --pairs $(cat $V2/subset24_pairs.txt)
ls $V2/extracted/sj03r50_delta-prefire | wc -l
