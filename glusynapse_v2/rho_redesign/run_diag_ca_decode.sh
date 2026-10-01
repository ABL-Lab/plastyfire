#!/bin/bash
#SBATCH --job-name=diag_ca_decode
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=21G
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/diag_ca_decode_%j.out
# CA_DECODE.md: single-pool spine-Ca features (dCa/dt, fast-sensor grid, high-pass, sharpness) vs vd_int, A0g_s3.
# CPU/numba, no jax. usage: sbatch glusynapse_v2/rho_redesign/run_diag_ca_decode.sh
# Sizing basis: diag_l23_ltp 22120179 (same loads + run_path, 396 records): 15.86 GB MaxRSS, 2:28 elapsed, 66% CPU.
# Extra here: per-record cacr (<= ~10 MB transient) + per-synapse feature table (~2.5k rows x ~110 cols, < 10 MB)
# + numba JIT of two small kernels (~0.3 GB). mem = (15.86 + ~0.5) x 1.25 = 20.5 -> 21G.
# Time: 2:28 + decode pass (~22 numba kernel calls per pairing window, ~1e9 exp total, est. 1-3 min) ~ 5 min x 1.5
# -> 0:15. 1 CPU (single-threaded numba, NUMBA_NUM_THREADS=1).
# MEASURED 22133762: 7:28, 10.14 GB MaxRSS, 75% CPU -> next 13G, 0:15.
set -euo pipefail
source .venv/bin/activate
export NUMBA_NUM_THREADS=1
V2=glusynapse_v2; R=$V2/rho_redesign
FIT=${FIT:-$R/results/v4_A0g_s3.json}
python -m py_compile $R/diag_ca_decode.py
python -u $R/diag_ca_decode.py --fit $FIT --save $R/results/ca_decode_A0g_s3 --fig $R/figs/fig9_ca_decode.png \
    --md $R/CA_DECODE.md
