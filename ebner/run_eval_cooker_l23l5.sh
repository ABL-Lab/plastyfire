#!/bin/bash
# Offline delta-cooker prediction on the L23PC->L5TTPC protocols (ebner/eval_cooker_l23l5.py).
# Sizing: v2_eval_all 160 GiB for ~3200 npz -> ~466 npz here ~23 GB + 25% -> 30G.
# Single-threaded (no precompute pool) -> 1 CPU; check seff after the run.
#SBATCH --job-name=eval_cooker_l23l5
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --time=01:00:00
#SBATCH --output=/home/dhuruva/ws/plastyfire/logs/eval_cooker_l23l5_%j.out
set -euo pipefail
ROOT=/home/dhuruva/ws/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
ls basis_results_edges_ebner_l23l5_delta/basis_*.csv | wc -l
python -u ebner/eval_cooker_l23l5.py
