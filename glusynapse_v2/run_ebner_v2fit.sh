#!/bin/bash
#SBATCH --job-name=v2_ebner_B
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=64
#SBATCH --mem=750G
#SBATCH --time=24:00:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/v2_ebner_B_%j.out
# Job B (tsk T17): v2.2 fit on all protocols; pre pathways driven by spine VDCC influx (pre_drive 1),
# NO gated by the pre-rate threshold theta_Z. Needs job A's extraction (run_ebner_pipeline.sh).
set -euo pipefail
source .venv/bin/activate
PF=$PWD; V2=$PF/glusynapse_v2
EB=$V2/extracted/ebner_delta-prefire; MK=$V2/extracted/markram_delta-prefire-tr
echo "=== [4] v2.2 fit, all protocols $(date)"
python -u $V2/fit_v2.py --dirs $EB,$MK --model v2 --rho full --workers 60 --maxiter 120 \
    --filters '{"pre_drive": 1, "i_scale": 1e-5}' \
    --free-filters theta_Ti,tau_T,theta_NOi,tau_NO,tau_Z,theta_Z --filter-steps 5 \
    --x0 $PF/fit_results/delta-cooker.json --save $V2/results/fit_v2_all
echo "=== done $(date)"
