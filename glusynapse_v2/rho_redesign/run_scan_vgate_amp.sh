#!/bin/bash
#SBATCH --job-name=scan_vgate_amp
#SBATCH --account=def-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=33G
#SBATCH --time=00:45:00
# MEASURED 22120612 (12 theta_V x C1,C2 one pass): 19.0 GB MaxRSS, 3:27 wall, 82% CPU -> next 24G, 0:15.
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/scan_vgate_amp_%j.out
# ROUND2.md C1 / C2 theta_V scan (12 values x 2 candidates), no refit, A0g_s3. CPU/numba, no jax.
# usage: sbatch glusynapse_v2/rho_redesign/run_scan_vgate_amp.sh
# Sizing basis:
#  mem: the L5 batch (all 5 L5 dirs, 24 pairs, effcai + vdcc, 8.0 GB npz) is the peak; the same loader on the same dirs
#       in the v4 1g fits measured 26.4 GB MaxRSS (DECISIONS) -> 26.4 x 1.25 = 33G. The L2/3 batch (5.0 GB npz, no
#       shaft_cai) loads after the L5 one is freed; diag 22120179 had 15.9 GB for 3.8 GB of L2/3 npz + shaft_cai.
#  time: all 24 configs run in one numba pass per synapse (not 24 passes), so cost ~ load + one pass. 22120179 took
#       2:28 for 4.4 GB npz -> 13.0 GB here x 3 = 7.3 min, x 2 for the 24-config inner loop = ~15 min, x 1.5 = 22
#       -> 0:45 (allowing for the unmeasured inner loop). 1 CPU (single-threaded numba, no pool).
set -euo pipefail
source .venv/bin/activate
export NUMBA_NUM_THREADS=1
V2=glusynapse_v2; R=$V2/rho_redesign
FIT=${FIT:-$R/results/v4_A0g_s3.json}
python -u $R/scan_vgate_amp.py --fit $FIT --save $R/results/vgate_amp_A0g_s3 --fig $R/figs/fig7_vgate_amp.png
