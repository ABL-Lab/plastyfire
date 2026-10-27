#!/bin/bash
# T22 pilot: end-to-end BCL prefire on one pipette group (2 protocols), live c_pre/c_post (cache miss).
#SBATCH --job-name=nevian_prefire_pilot
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/home/dhuruva/ws/plastyfire/logs/nevian_prefire_pilot_%j.out
ROOT=/home/dhuruva/ws/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
python -u run_de_fit2_pool.py --cooker --param-hash delta-prefire \
 --results-dir refitting_results/fitting/n120/seed20262009/Ebner2019_L23PC_L23PC_basal/simulations \
 --cache cpre_cpost_cache/sabrina_n120_delta.pkl --circuit-config data/dhuruva_delta_circuit_config.json \
 --trace-vars cai_CR,shaft_cai,ica_VDCC --workers 2 --pairs ${PAIRS:-pip0-18559} \
 --protocols ${PROTOS:-nevian_3ap_50hz_dt+10ms,nevian_epsp_only}
