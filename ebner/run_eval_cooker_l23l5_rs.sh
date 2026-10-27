#!/bin/bash
# Offline delta-cooker prediction on the L23PC->L5TTPC protocols (ebner/eval_cooker_l23l5_rs.py).
# Re-split (T29) inputs, ebner/eval_cooker_l23l5_rs.py. Sized from 22069136 (19.3 GiB, 1 min) -> 24G, 0:15.
#SBATCH --job-name=eval_cooker_l23l5_rs
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:15:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/eval_cooker_l23l5_rs_%j.out
set -euo pipefail
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
ls basis_results_edges_ebner_l23l5_delta_rs/basis_*.csv | wc -l
python -u ebner/eval_cooker_l23l5_rs.py
