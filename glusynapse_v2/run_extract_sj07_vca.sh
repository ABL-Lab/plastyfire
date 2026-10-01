#!/bin/bash
#SBATCH --job-name=vca_extract_sj07
#SBATCH --account=rrg-emuller
# t_drive 4: Sjostrom 2007 (hash delta-prefire-vseg, ebner job 22089771, 94 workdirs, 24 pairs) -> extracted/sj07_delta-prefire-vca.
# Sizing: 22090550 (sj03 220 + r50 48 + ebner 94 files): MaxRSS 73.8 GB on 4 workers (the per-worker memory follows the trace
#   length, unknown for sj07) -> 4 workers, 95G (73.8 + 25% = 92). Time: 268 files took ~26 min with the slow ca_events (4x slower now
#   blocked); 94 files ~ 9 min slow-equivalent -> 0:20. --skip-existing: rerun completes it if short.
#SBATCH --cpus-per-task=4
#SBATCH --mem=95G
#SBATCH --time=00:20:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/vca_extract_sj07_%j.out
set -euo pipefail
source .venv/bin/activate
PF=$PWD; V2=$PF/glusynapse_v2
SIMS=$PF/refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
PAIRS=$(cat $V2/subset24_pairs.txt)
python -u $V2/extract_v2.py --param-hash delta-prefire-vseg --index-csv $SIMS/../../index_Sjostrom2007_L5TTPC_L5TTPC.csv \
    --window --workers 4 --skip-existing --out $V2/extracted/sj07_delta-prefire-vca --pairs $PAIRS
echo "sj07 $(ls $V2/extracted/sj07_delta-prefire-vca | wc -l)"
