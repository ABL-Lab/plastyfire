#!/bin/bash
#SBATCH --job-name=vca_extract_ljp25
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/vca_extract_ljp25_%j.out
# t_drive 4, ljp 25: the delta-ljp25-prefire-vseg prefire -> extracted/{sj03,sj03r50,ebner,markram,sj07}_delta-ljp25-prefire-vca
# MEASURED 22107521 (2026-09-30): MaxRSS 24.9 GB of 145G, 50:51 wall, 63% CPU on 4 workers -> next run 32G, 1:30.
# (same code path as run_extract_vca_full.sh; vev, cev, cev_lo, cev_hi). K = 0.2% of the median per-bAP -ica_VDCC peak on the ljp-25 sj07 post-only
# traces: pooled median 1.7358e-4 nA -> 3.4715e-7 nA (22096823, k_ca_measure.py; ljp 0 by the same measure 3.134e-8 vs 2.469e-8 from local_t used).
# Sizing basis (sacct 2026-09-30: 22074074 MaxRSS 120.6 GB, 22093100 82.6 GB / 9 min, 22090550 77.4 GB timed out at 30 min): 22090550 sj03+r50+ebner MaxRSS 73.8 GB on 4 workers; 22074074 markram-inclusive
#   peak 115 GB -> 145G as in run_extract_vca_full.sh. Time: ebner 146 + markram 168 files in <= 15.5 min (22093100 file times);
#   770 files total -> <= 38 min, +50% -> 1:00. --skip-existing: resumable.
set -euo pipefail
source .venv/bin/activate
PF=$PWD; V2=$PF/glusynapse_v2; H=delta-ljp25-prefire-vseg
SIMS=$PF/refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
PAIRS=$(cat $V2/subset24_pairs.txt); W=4; C=$PF/cpre_cpost_cache/sabrina_n120_delta-ljp25.pkl
X="python -u $V2/extract_v2.py --param-hash $H --cache $C --workers $W --skip-existing --k-ca 3.4715e-7 --pairs $PAIRS"
$X --index-csv $SIMS/../../index_Sjostrom2003_L5TTPC_L5TTPC.csv --window --out $V2/extracted/sj03_delta-ljp25-prefire-vca
$X --index-csv $SIMS/../../index_Sjostrom2003b_L5TTPC_L5TTPC.csv --window --out $V2/extracted/sj03r50_delta-ljp25-prefire-vca
$X --index-csv $SIMS/../../index_Ebner2019_L5TTPC_L5TTPC.csv --window --out $V2/extracted/ebner_delta-ljp25-prefire-vca
$X --out $V2/extracted/markram_delta-ljp25-prefire-vca
$X --index-csv $SIMS/../../index_Sjostrom2007_L5TTPC_L5TTPC.csv --window --out $V2/extracted/sj07_delta-ljp25-prefire-vca
for d in sj03 sj03r50 ebner markram sj07; do echo "$d $(ls $V2/extracted/${d}_delta-ljp25-prefire-vca | wc -l)"; done
