#!/bin/bash
# T29 step 7: EPSP histograms (csv reads, PairBasis, one h5 read of 798 rho0 values). 1 CPU, 22091030: 24 s, 2.1 GB -> 3G, 0:15; fix from seff.
#SBATCH --job-name=epsp_hist
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=00:15:00
#SBATCH --output=/home/dhuruva/ws/plastyfire/logs/epsp_hist_%j.out
set -euo pipefail
ROOT=/home/dhuruva/ws/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
python -u ebner/epsp_hist_l23l5_vs_l5l5.py
