#!/bin/bash
#SBATCH --job-name=vca_extract_full
#SBATCH --account=rrg-emuller
# 22090550 (first attempt, all 4 dirs): TIMEOUT at 0:30 after sj03 220 + sj03r50 48 + ebner 94 files, MaxRSS 73.8 GB (of 145G), 94% CPU on 4.
#   The slow part was ca_events (3 K variants x full 20M-sample scans: ~4x slower than the vseg extract); now block-skipping
#   (model_v2.ca_events, identical result, test_t_drive4_events_fast). This rerun (--skip-existing) does the rest: ebner 146 + markram 168 files.
#   Memory kept at 145G: markram was the 115 GB peak in 22074074 and is not yet measured here. Original sizing basis 22074074: 115 GB, 14:44, 4 workers.
#   3 x n_syn x n_events floats and one float64 copy of one synapse's ica trace (160 MB): no change to memory.
#SBATCH --cpus-per-task=4
#SBATCH --mem=145G
#SBATCH --time=00:30:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/vca_extract_full_%j.out
# t_drive 4: re-extract the same delta-prefire-vseg outputs (vev AND cev/cev_lo/cev_hi) into new *-vca dirs.
# Sizing: pilot 22073365 (4 workers, 26 sj03/r50 records): 90 GB MaxRSS (memory is per worker, not per record), 1:48, 79% CPU
#   -> 4 workers, 115G. Full = ~672 records, ~26x the pilot -> ~45 min; +50% -> 1:15. Correct from seff.
set -euo pipefail
source .venv/bin/activate
PF=$PWD; V2=$PF/glusynapse_v2
SIMS=$PF/refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
PAIRS=$(cat $V2/subset24_pairs.txt); W=4
python -u $V2/extract_v2.py --param-hash delta-prefire-vseg --index-csv $SIMS/../../index_Sjostrom2003_L5TTPC_L5TTPC.csv \
    --window --workers $W --skip-existing --out $V2/extracted/sj03_delta-prefire-vca --pairs $PAIRS
python -u $V2/extract_v2.py --param-hash delta-prefire-vseg --index-csv $SIMS/../../index_Sjostrom2003b_L5TTPC_L5TTPC.csv \
    --window --workers $W --skip-existing --out $V2/extracted/sj03r50_delta-prefire-vca --pairs $PAIRS
python -u $V2/extract_v2.py --param-hash delta-prefire-vseg --index-csv $SIMS/../../index_Ebner2019_L5TTPC_L5TTPC.csv \
    --window --workers $W --skip-existing --out $V2/extracted/ebner_delta-prefire-vca --pairs $PAIRS
python -u $V2/extract_v2.py --param-hash delta-prefire-vseg --workers $W --skip-existing \
    --out $V2/extracted/markram_delta-prefire-vca --pairs $PAIRS
for d in sj03 sj03r50 ebner markram; do echo "$d $(ls $V2/extracted/${d}_delta-prefire-vca | wc -l)"; done
