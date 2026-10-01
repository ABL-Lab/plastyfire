#!/bin/bash
#SBATCH --job-name=eval_v4
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=21G
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/eval_v4_%j.out
# L2/3->L5 transfer (no refit) of v3/v4 fits with the v4 post rule, CPU/numba (eval_v4.py). usage: FITS="<json> <json> ..." sbatch ...
# Sizing: 22108820 (5 fits): 15.9 GB MaxRSS, 8:51 (~1.8 min per fit) -> 21G, 0:15 for up to 5 fits.
set -euo pipefail
source .venv/bin/activate
export NUMBA_NUM_THREADS=1
for F in $FITS; do
  B=$(basename $F .json)
  echo "=== $B"
  ANALYTICAL_BASIS_DIR=basis_results_edges_ebner_l23l5_delta_rs python -u glusynapse_v2/rho_redesign/eval_v4.py --l23 --fit $F \
      --dirs glusynapse_v2/extracted/ebner_l23l5_delta-prefire-vseg-rs --save glusynapse_v2/rho_redesign/results/l23_$B || echo "FAILED $B"
done
