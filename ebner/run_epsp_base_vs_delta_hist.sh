#!/bin/bash
# Plot ebner/epsp_base_vs_delta_hist.py (4 small csv). Sized from epsp_hist 22091030 (24 s, 2.1 GB) -> 3G, 0:15.
#SBATCH --job-name=epsp_bvd
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=00:15:00
#SBATCH --output=/home/dhuruva/ws/plastyfire/logs/epsp_bvd_%j.out
set -euo pipefail
ROOT=/home/dhuruva/ws/plastyfire
cd $ROOT && source .venv/bin/activate
python -u ebner/epsp_base_vs_delta_hist.py
