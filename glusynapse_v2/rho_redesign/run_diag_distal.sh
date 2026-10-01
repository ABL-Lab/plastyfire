#!/bin/bash
#SBATCH --job-name=diag_distal
#SBATCH --account=def-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/diag_distal_%j.out
# Distal Letzkus diagnosis (DISTAL_LTP.md step 1), CPU/numba, no jax. usage: sbatch glusynapse_v2/rho_redesign/run_diag_distal.sh
# Sizing: from diag_rho l23 22108334 (same -rs dir, all 349 L2/3 records + vdcc): 14.99 GB MaxRSS (at its 15G limit),
# 1:29 elapsed, 80% CPU -> 15 x 1.25 = 19G; 1:29 x 1.5 -> 0:15. This job loads only the 3 Letzkus protocols (+ shaft_cai),
# MEASURED: 22115348 15.74 GB MaxRSS, 8:12 elapsed, 95% CPU (request 19G = 83%) -> next run --mem=20G --time=0:15 (8:12 x 1.5 = 12:18).
set -euo pipefail
source .venv/bin/activate
export NUMBA_NUM_THREADS=1
V2=glusynapse_v2; R=$V2/rho_redesign
FIT=${FIT:-$R/results/v4_A0_joint_s4.json}
ANALYTICAL_BASIS_DIR=basis_results_edges_ebner_l23l5_delta_rs python -u $R/diag_distal.py --fit $FIT \
  --dirs $V2/extracted/ebner_l23l5_delta-prefire-vseg-rs --save $R/results/distal_A0_joint_s4 --fig $R/figs/fig5_distal.png
