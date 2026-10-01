#!/bin/bash
#SBATCH --job-name=v3_eq
#SBATCH --account=def-emuller
#SBATCH --gpus-per-node=nvidia_h100_80gb_hbm3_2g.20gb:1
#SBATCH --cpus-per-task=1
#SBATCH --mem=45G
#SBATCH --time=00:30:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/v3_eq_%j.out
# gpu_v3 == jax_v2 (CHUNK 4, as the running fits) on the pl5sj07 setup, plus v2/v3 timings. TDRIVE=4 (default) or 3.
# Sizing: fit pilots 22098228/9 (B + v2 chunks) measured 35.9 GB MaxRSS, 3 min; v3 adds < 1 GB after its build -> 45G.
# Time: 3 min load + ~6 v2 calls at P~100 (1-2 min each on 2g) -> 0:30.
# 22105399 (failed at the first v3 kernel compile, before the v2 build): 24.0 GB MaxRSS, 1:57, 59% CPU. Loading takes 74 s.
# PASS 22107148 (td4, chi2 rel 1.2e-15): 40.1 GB MaxRSS, 8:13, 91% CPU; 22107522 (td3, 4.7e-16): 37.7 GB, 7:59 -> 45G / 0:15 is enough next time.
set -euo pipefail
source glusynapse_v2/env_v3.sh
nvidia-smi -L
python -u glusynapse_v2/tests/test_jax_v3.py
