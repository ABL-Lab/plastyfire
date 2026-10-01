#!/bin/bash
# Offline attribution of the live-vs-offline dpre difference of equiv 22131696 (diag_dpre_v4.py): 3 records, each
# loads one ~1.4 GB simulation_traces.pkl (sequentially) for the full-resolution -ica_VDCC.
# Sizing: no comparable job; one 1.4 GB pickle + one record's BatchV2 arrays (<0.5 GB) -> ~2.5 GB peak, +25% and
#   headroom for the unpickled dict -> 4G. Single-threaded: 1 CPU. ~3 x (pickle load ~1 min + features) -> 0:15.
#SBATCH --job-name=diag_dpre_v4
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/diag_dpre_v4_%j.out
set -euo pipefail
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1
python -u glusynapse_v2/bcl_validation/diag_dpre_v4.py --fit glusynapse_v2/rho_redesign/results/v4_C1Ajd_s5.json \
    --live /scratch/dhuruva/bcl_validation/equiv_C1Ajd_s5_22131696.jsonl
