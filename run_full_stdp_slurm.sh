#!/bin/bash
#SBATCH --job-name=full_stdp
#SBATCH --account=rrg-emuller
#SBATCH --time=12:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=60
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/full_stdp_%x_%j.out
#SBATCH --error=/project/rrg-emuller/dhuruva/plastyfire/logs/full_stdp_%x_%j.err
set -euo pipefail
# Full-protocol STDP (C01 + induction + C02, fastforward at 280 s), EPSP ratio
# measured from the soma voltage -- no basis. Then curate the curve.
#
#   sbatch --job-name=full_fixed_beam --export=ALL,NAME=fixed_beam,PRESET=--fixed-beam run_full_stdp_slurm.sh
#   sbatch --job-name=full_x --export=ALL,NAME=x,PARAMS_JSON=/abs/path/best.json run_full_stdp_slurm.sh
: "${NAME:?set NAME}"

cd /project/rrg-emuller/dhuruva/plastyfire
source /project/rrg-emuller/dhuruva/plastyfire/.venv/bin/activate

SIMS=refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
CACHE=${CACHE:-cpre_cpost_cache/sabrina_n120_defit2.pkl}
BASIS=${BASIS:-basis_results_edges_sabrina_n120}
# CIRCUIT picks the emodel (e.g. data/dhuruva_delta_circuit_config.json)
CIRC_ARG=${CIRCUIT:+--circuit-config $CIRCUIT}
if [ -n "${PARAMS_JSON:-}" ]; then
    SRC="--params-json $PARAMS_JSON --param-hash $NAME"
else
    : "${PRESET:?set PRESET (e.g. --fixed-beam) or PARAMS_JSON}"
    SRC="$PRESET --param-hash $NAME"
fi

echo "=== full-protocol BCL STDP: $NAME ($(date)) ==="
python -u run_de_fit2_pool.py $SRC --full --workers 60 --results-dir $SIMS --cache $CACHE $CIRC_ARG

echo "=== curate ($(date)) ==="
python -u curate_stdp_full.py --results-dir $SIMS --param-hash $NAME \
    --basis-dir $BASIS \
    --label "$NAME, n120 (full protocol, BCL)" --out stdp_full_${NAME}_sabrina_n120.png
echo "=== done ($(date)) ==="
