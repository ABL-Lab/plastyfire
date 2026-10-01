#!/bin/bash
#SBATCH --job-name=v22r
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=48
#SBATCH --mem=250G
#SBATCH --time=12:00:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/v22r_%x_%j.out
# Reduced v2.2 (10 fitted: 4 a + A_mglu, A_NO + shared VDCC threshold, tau_NO, tau_Z, theta_Z; tau_T fixed 3 ms).
# MODE=preview: the 8 preview pairs (compare with the 12-parameter v2.2, chi2 84).
# MODE=subset:  24 random pairs with all 35 protocols (subset24_pairs.txt, seed 0).
set -euo pipefail
source .venv/bin/activate
V2=glusynapse_v2
python $V2/tests/test_units.py
if [ "$MODE" = preview ]; then
  D=$V2/extracted/ebner_preview,$V2/extracted/markram_delta-cooker
  PAIRS=180351-198084,181015-184976,182339-200396,182381-186941,182900-195902,183240-189246,183769-207759,186261-208490
  echo "=== eval_v2 check: 12-parameter v2.2 preview must give chi2 84.16 $(date)"
  python -u $V2/eval_v2.py --fit $V2/results/preview_fit_v22.json --dirs $D --pairs $PAIRS --save $V2/results/eval_check_preview_v22
else
  D=$V2/extracted/ebner_delta-prefire,$V2/extracted/markram_delta-prefire-tr
  PAIRS=$(cat $V2/subset24_pairs.txt)
fi
echo "=== reduced fit ($MODE) $(date)"
python -u $V2/fit_v2.py --dirs $D --pairs $PAIRS --model v2 --workers 48 --maxiter 80 \
    --filters '{"pre_drive": 1, "i_scale": 1e-5, "tau_T": 3.0}' --tie '{"theta_NOi": "theta_Ti"}' \
    --free-filters theta_Ti,tau_NO,tau_Z,theta_Z --filter-steps 5 \
    --x0 fit_results/delta-cooker.json --save $V2/results/reduced_$MODE
echo "=== done $(date)"
