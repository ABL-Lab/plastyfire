#!/bin/bash
#SBATCH --job-name=fit_cooker
#SBATCH --account=rrg-emuller
#SBATCH --time=03:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=60
#SBATCH --output=/home/dhuruva/ws/plastyfire/logs/fit_cooker_%j.out
#SBATCH --error=/home/dhuruva/ws/plastyfire/logs/fit_cooker_%j.err
set -euo pipefail

# "cooker": 4-param fit, theta_d = a00*c_pre + a01*c_post, theta_p = a10*c_pre + a11*c_post
# (apical tied to basal); gammas fixed at constants.py defaults (no --fit-gamma).
PF=/home/dhuruva/ws/plastyfire
cd $PF/analytical_method
source $PF/.venv/bin/activate
export PYTHONPATH=$PF
export ANALYTICAL_EXTRACTED=$PF/analytical_method/extracted_sabrina_n120_d0p25
export ANALYTICAL_BASIS_DIR=$PF/basis_results_edges_sabrina_n120

echo "=== cooker fit ($(date)) ==="
python -u fit.py --de --workers 60 --maxiter 600 --stride 1 \
    --fit-dt=-10,5,10 --anchor-tails
echo "=== done ($(date)) ==="
