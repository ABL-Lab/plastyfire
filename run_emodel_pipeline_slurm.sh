#!/bin/bash
#SBATCH --job-name=emodel_pipeline
#SBATCH --account=rrg-emuller
#SBATCH --time=10:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=60
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/pipeline_%x_%j.out
#SBATCH --error=/project/rrg-emuller/dhuruva/plastyfire/logs/pipeline_%x_%j.err
set -euo pipefail
# emodel -> fitted STDP params -> BCL full-protocol check -> neurodamus validation.
#
#   sbatch --job-name=delta-cooker \
#       --export=ALL,EMODEL_DIR=/project/rrg-emuller/dhuruva/DEES_cell_packages/SSCx-AAD-delta-emodels,EMODEL_TAG=delta,NAME=delta-cooker \
#       run_emodel_pipeline_slurm.sh
#
#   EMODEL_DIR  dir with the cADpyr_*.hoc templates (biophysical_neuron_models_dir)
#   EMODEL_TAG  short emodel label; keys everything that depends only on the emodel
#               (circuit config, c_pre/c_post cache, basis, prefire traces, extraction)
#   NAME        fit label (default <EMODEL_TAG>-cooker); keys the fit JSON, BCL pkls, ND out_<NAME>/
#   FIT_MAXITER DE generations (default 100)
#   FIT_FLAGS   extra fit.py flags (default none = cooker: 4 a-params, apical tied, gammas fixed;
#               e.g. "--cpre-only --fit-gamma" for beam)
#
# Every stage skips work whose output already exists, so rerunning resumes.
# Neurodamus reads thetas from the shared dhuruva_modified_edges.h5, so run one pipeline at a time
# through stage 11 (a second one would overwrite the injected thetas mid-validation).
# Stages: [1] circuit config  [2] c_pre/c_post cache  [3] basis (one sbatch per pair, async)
#   [4] BCL prefire traces  [5] extract  [6] wait basis  [7] analytical fit
#   [8] inject thetas + submit neurodamus array (async)  [9] BCL prefire with the fit (ND reference)
#   [10] BCL full protocol + curate  [11] wait neurodamus, compare ND vs BCL
: "${EMODEL_DIR:?set EMODEL_DIR}"
: "${EMODEL_TAG:?set EMODEL_TAG}"
NAME=${NAME:-${EMODEL_TAG}-cooker}
FIT_FLAGS=${FIT_FLAGS:-}
FIT_MAXITER=${FIT_MAXITER:-100}   # DE generations; delta-cooker plateaued by gen ~80

PF=/project/rrg-emuller/dhuruva/plastyfire
cd $PF
source $PF/.venv/bin/activate
export PYTHONPATH=$PF

SIMS=$PF/refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
EDGES=$PF/data/dhuruva_modified_edges.h5
CIRCUIT=$PF/data/dhuruva_${EMODEL_TAG}_circuit_config.json
CACHE=$PF/cpre_cpost_cache/sabrina_n120_${EMODEL_TAG}.pkl
BASIS=$PF/basis_results_edges_sabrina_n120_${EMODEL_TAG}
EXTRACTED=$PF/analytical_method/extracted_sabrina_n120_${EMODEL_TAG}_d0p25
PREFIRE_HASH=${EMODEL_TAG}-prefire   # cooker (gamma-fit) params, only used to generate calcium traces
JSON=$PF/fit_results/${NAME}.json
THR=$PF/threshold_results_${NAME}_sabrina_n120
N_PAIRS=$(ls -d $SIMS/*-* | wc -l)
stage() { echo; echo "=== [$1] $2 ($(date)) ==="; }
n_csv() { ls $BASIS/basis_*.csv 2>/dev/null | wc -l; }
echo "emodel $EMODEL_DIR ($EMODEL_TAG) | fit $NAME | flags '${FIT_FLAGS}' | $N_PAIRS pairs"

stage 1 "circuit config"
if [ ! -f $CIRCUIT ]; then
    python - "$EMODEL_DIR" "$CIRCUIT" <<'EOF'
import json, sys
c = json.load(open("data/dhuruva_modified_ion_channels_circuit_config.json"))
c["networks"]["nodes"][0]["populations"]["S1nonbarrel_neurons"]["biophysical_neuron_models_dir"] = sys.argv[1]
json.dump(c, open(sys.argv[2], "w"), indent=2)
EOF
fi
grep biophysical_neuron_models_dir $CIRCUIT | tail -1
grep -q "\"$EMODEL_DIR\"" $CIRCUIT || { echo "circuit config points elsewhere"; exit 1; }

stage 2 "c_pre/c_post cache"
[ -f $CACHE ] || python -u precompute_cpre_cpost.py --params defit2 --results-dir refitting_results \
    --circuit-config $CIRCUIT --edges-h5 $EDGES --output $CACHE --workers 60

stage 3 "basis (async, one job per pair)"
if [ $(n_csv) -lt $N_PAIRS ]; then
    python -u submit_basis_edges.py --execution-mode slurm --stdp-name Sabrina_L5TTPC_L5TTPC_STDP \
        --output-dir $BASIS --circuit-config $CIRCUIT --walltime 01:00:00 --skip-existing \
        2>&1 | grep -v -i 'memory request' || true
fi
echo "basis csvs now: $(n_csv) / $N_PAIRS"

stage 4 "BCL prefire traces"
python -u run_de_fit2_pool.py --cooker --param-hash $PREFIRE_HASH --workers 60 --skip-existing \
    --results-dir $SIMS --cache $CACHE --circuit-config $CIRCUIT

stage 5 "extract calcium traces, 0.25 ms"
if [ -z "$(ls $EXTRACTED 2>/dev/null)" ]; then
    (cd analytical_method && python -u extract.py --workers 60 --decim 0.25 \
        --param-hash $PREFIRE_HASH --sims $SIMS --cache $CACHE --out $EXTRACTED)
fi
echo "extracted: $(ls $EXTRACTED | wc -l) files"

stage 6 "wait for basis"
while [ $(n_csv) -lt $N_PAIRS ] && squeue -u $USER -h -o '%j' | grep -q '^basis_'; do sleep 30; done
echo "basis csvs: $(n_csv) / $N_PAIRS"

stage 7 "analytical fit -> $JSON"
export ANALYTICAL_EXTRACTED=$EXTRACTED ANALYTICAL_BASIS_DIR=$BASIS
mkdir -p $PF/fit_results
[ -f $JSON ] || (cd analytical_method && python -u fit.py --de --workers 60 --maxiter $FIT_MAXITER --stride 1 \
    --fit-dt=-10,5,10 --anchor-tails $FIT_FLAGS --save-json $JSON)
python -c "import json; print(json.load(open('$JSON'))['params'])"

stage 8 "inject thetas + submit neurodamus"
A=$(python -c "import json; p=json.load(open('$JSON'))['params']; print(' '.join(f'--{k} {p[k]}' for k in ('a00','a01','a10','a11','a20','a21','a30','a31')))")
python -u compute_thresholds_from_cache.py --cache $CACHE --output-dir $THR $A --edges-h5 $EDGES | tail -2
python -u inject_thresholds.py --results-dir $THR --edges-h5 $EDGES --force | grep -E 'Updated|wrote|Done' || true
ND_JOB=$(python -u submit_neurodamus_fit.py --tag $NAME --params-json $JSON --cache $CACHE \
    --circuit-config $CIRCUIT --results-dir $SIMS --skip-existing | tee /dev/stderr | sed -n 's/^ARRAY_JOB_ID=//p')
echo "neurodamus array: ${ND_JOB:-none}"

stage 9 "BCL prefire with the fit (neurodamus reference)"
python -u run_de_fit2_pool.py --params-json $JSON --param-hash $NAME --workers 60 --skip-existing \
    --results-dir $SIMS --cache $CACHE --circuit-config $CIRCUIT

stage 10 "BCL full protocol + curate"
# a standalone run_full_stdp_slurm.sh job for this NAME may already be going -- let it finish
while squeue -u $USER -h -o '%j' | grep -qx "full_$NAME"; do sleep 60; done
python -u run_de_fit2_pool.py --params-json $JSON --param-hash $NAME --full --workers 60 --skip-existing \
    --results-dir $SIMS --cache $CACHE --circuit-config $CIRCUIT
python -u curate_stdp_full.py --results-dir $SIMS --param-hash $NAME --basis-dir $BASIS \
    --label "$NAME, n120 (full protocol, BCL)" --out $PF/stdp_full_${NAME}_sabrina_n120.png

stage 11 "neurodamus vs BCL"
if [ -n "${ND_JOB:-}" ]; then
    while squeue -j $ND_JOB -h 2>/dev/null | grep -q .; do sleep 30; done
    sacct -j $ND_JOB -X -n -o State | sort | uniq -c
fi
python -u compare_nd_bcl_prefire.py --tag $NAME --results-dir $SIMS --basis-dir $BASIS
echo; echo "=== done ($(date)) ==="
echo "fit      : $JSON"
echo "BCL full : $PF/stdp_full_${NAME}_sabrina_n120.png"
echo "ND vs BCL: $PF/compare_nd_bcl_${NAME}.png"
