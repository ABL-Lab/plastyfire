#!/bin/bash
#SBATCH --job-name=v2_prevfit
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=32
#SBATCH --mem=180G
#SBATCH --time=08:00:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/v2_prevfit_%j.out
# Preview of T11/T12/T17 on the 8 preview pairs (26 Ebner protocols) + the same 8 pairs' Markram records.
set -euo pipefail
source .venv/bin/activate
V2=glusynapse_v2
PAIRS=180351-198084,181015-184976,182339-200396,182381-186941,182900-195902,183240-189246,183769-207759,186261-208490
D=$V2/extracted/ebner_preview,$V2/extracted/markram_delta-cooker
echo "=== post-only (cooker) refit $(date)"
python -u $V2/fit_v2.py --dirs $D --pairs $PAIRS --model post --workers 32 --maxiter 60 \
    --x0 fit_results/delta-cooker.json --save $V2/results/preview_fit_post
echo "=== v2 fit $(date)"
python -u $V2/fit_v2.py --dirs $D --pairs $PAIRS --model v2 --workers 32 --maxiter 80 \
    --free-filters theta_T,tau_T,theta_NO,tau_NO,tau_Z --filter-steps 5 \
    --x0 fit_results/delta-cooker.json --save $V2/results/preview_fit_v2
echo "=== done $(date)"
