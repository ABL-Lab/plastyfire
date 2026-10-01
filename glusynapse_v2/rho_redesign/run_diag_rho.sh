#!/bin/bash
#SBATCH --job-name=diag_rho
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=17G
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/diag_rho_%j.out
# rho-rule diagnosis (RHO_REDESIGN.md step 1), CPU/numba, no jax. usage:
#   WHICH=l5  sbatch glusynapse_v2/rho_redesign/run_diag_rho.sh             (pl5sj07 -vca dirs, 502 records)
#   WHICH=l23 sbatch --mem=15G glusynapse_v2/rho_redesign/run_diag_rho.sh   (L2/3->L5 -rs dir, 349 records)
# Sizing: L5 from fit 22107773 (same 502 records + vdcc, 21.6 GB MaxRSS) -> 27G; L2/3 from eval 22108122 (349 records,
# 11.8 GB, 2:05) -> 15G; time 0:15.
# MEASURED: l5 22108333 13.2 GB, 8:05 (CPU ~1 core) -> next time --mem=17G; l23 22108334 15.0 GB (at the 15G limit), 1:29 -> 19G.
set -euo pipefail
source .venv/bin/activate
export NUMBA_NUM_THREADS=1
V2=glusynapse_v2; X=$V2/extracted; FIT=${FIT:-$V2/results/reduced_gpu_subset_pl5sj07_td4_s2.json}
TAG=${TAG:-td4_s2}
mkdir -p $V2/rho_redesign/out
if [ "$WHICH" = l5 ]; then
  python -u $V2/rho_redesign/diag_rho.py --fit $FIT --groups paired_l5,sjostrom07 \
    --dirs $X/ebner_delta-prefire-vca,$X/markram_delta-prefire-vca,$X/sj03_delta-prefire-vca,$X/sj03r50_delta-prefire-vca,$X/sj07_delta-prefire-vca \
    --save $V2/rho_redesign/out/l5_$TAG
else
  ANALYTICAL_BASIS_DIR=basis_results_edges_ebner_l23l5_delta_rs python -u $V2/rho_redesign/diag_rho.py --fit $FIT --l23 \
    --groups paired_l23l5 --dirs $X/ebner_l23l5_delta-prefire-vseg-rs --save $V2/rho_redesign/out/l23_$TAG
fi
