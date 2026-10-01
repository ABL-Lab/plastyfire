#!/bin/bash
#SBATCH --job-name=vseg_extract_full
#SBATCH --account=rrg-emuller
# MEASURED 22074074: MaxRSS 115 GB (at the limit!), 14:44, 93% CPU on 4 -> next run 145G, 0:30.
#SBATCH --cpus-per-task=4
#SBATCH --mem=145G
#SBATCH --time=00:30:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/vseg_extract_full_%j.out
# t_drive 3: extract the 24-pair v_seg re-prefire (22073522 sj03+r50, 22073523 sj01+Markram) into new *-vseg dirs.
# Sizing: pilot 22073365 (4 workers, 26 sj03/r50 records): 90 GB MaxRSS (memory is per worker, not per record), 1:48, 79% CPU
#   -> 4 workers, 115G. Full = ~672 records, ~26x the pilot -> ~45 min; +50% -> 1:15. Correct from seff.
set -euo pipefail
source .venv/bin/activate
PF=$PWD; V2=$PF/glusynapse_v2
SIMS=$PF/refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
PAIRS=$(cat $V2/subset24_pairs.txt); W=4
python -u $V2/extract_v2.py --param-hash delta-prefire-vseg --index-csv $SIMS/../../index_Sjostrom2003_L5TTPC_L5TTPC.csv \
    --window --workers $W --skip-existing --out $V2/extracted/sj03_delta-prefire-vseg --pairs $PAIRS
python -u $V2/extract_v2.py --param-hash delta-prefire-vseg --index-csv $SIMS/../../index_Sjostrom2003b_L5TTPC_L5TTPC.csv \
    --window --workers $W --skip-existing --out $V2/extracted/sj03r50_delta-prefire-vseg --pairs $PAIRS
python -u $V2/extract_v2.py --param-hash delta-prefire-vseg --index-csv $SIMS/../../index_Ebner2019_L5TTPC_L5TTPC.csv \
    --window --workers $W --skip-existing --out $V2/extracted/ebner_delta-prefire-vseg --pairs $PAIRS
python -u $V2/extract_v2.py --param-hash delta-prefire-vseg --workers $W --skip-existing \
    --out $V2/extracted/markram_delta-prefire-vseg --pairs $PAIRS
for d in sj03 sj03r50 ebner markram; do echo "$d $(ls $V2/extracted/${d}_delta-prefire-vseg | wc -l)"; done
