#!/bin/bash
# T6: effcai extraction of the Ebner 2019 prefire traces for the analytical fit.
# Window mode: 0.25 ms from 50 ms before to 1000 ms after every pre/post spike, 5 ms elsewhere
# (max |drho| vs flat 0.25 ms = 1.4e-7 on the smoke pair). Each worker unpickles a whole
# simulation_traces.pkl (all 3 trace vars, up to ~9 GB for 1000 s x 19 synapses), hence the memory.
#SBATCH --job-name=ebner_extract
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=24
#SBATCH --mem=500G
#SBATCH --time=12:00:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/ebner_extract_%j.out
set -eu
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
N=refitting_results/fitting/n120/seed20262009
python -u analytical_method/extract.py --workers ${WORKERS:-24} --window --decim 0.25 --skip-existing \
  --param-hash delta-prefire --sims $N/Sabrina_L5TTPC_L5TTPC_STDP/simulations \
  --index-csv $N/index_Ebner2019_L5TTPC_L5TTPC.csv --cache cpre_cpost_cache/sabrina_n120_delta.pkl \
  --out analytical_method/extracted_ebner2019_delta_win0p25
