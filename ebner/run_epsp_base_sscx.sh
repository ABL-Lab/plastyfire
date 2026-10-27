#!/bin/bash
# Unitary EPSP per pair, base SSCx (original emodels + edges) vs delta, ebner/epsp_base_sscx.py.
# Sized from pair_geometry 22064512 (same readout, 120 pairs, 30 spikes, 16 workers, ran within 23G / <2h).
# Task 0: L23->L5 base, 1: L5->L5 base. 22093162: task 2 (L5->L5 delta) done -> ebner/figs/epsp_delta_l5l5.csv; tasks 0-1 failed (ctb path, batch path).
#SBATCH --job-name=epsp_base
#SBATCH --account=rrg-emuller
#SBATCH --array=0-1
#SBATCH --cpus-per-task=16
#SBATCH --mem=23G
#SBATCH --time=01:30:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/epsp_base_%A_%a.out
set -euo pipefail
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
R=refitting_results/fitting/n120/seed20262009
BASE=$ROOT/data/sscx_base_circuit_config.json   # dataset config with $BASE_DIR moved ctb-emuller -> rrg-emuller (no ctb access)
DELTA=$ROOT/data/dhuruva_delta_circuit_config.json
IDX=($R/index_Ebner2019_L23PC_L5TTPC.csv $R/index_Ebner2019_L5TTPC_L5TTPC.csv $R/index_Ebner2019_L5TTPC_L5TTPC.csv)
CIRC=($BASE $BASE $DELTA)
OUT=(epsp_base_l23l5 epsp_base_l5l5 epsp_delta_l5l5)
i=$SLURM_ARRAY_TASK_ID
python -u ebner/epsp_base_sscx.py ${IDX[$i]} ${CIRC[$i]} ebner/figs/${OUT[$i]}.csv 16
