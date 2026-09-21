#!/bin/bash
# Submit pairrunner_edges_fit.py (with a-params) for the 30 pairs sampled by the
# optimizer (seed=1234, sample_size=30) × 2 protocols (+10ms, -10ms).
#
# These are the exact same 30 pairs the optimizer sees each evaluation.
# Output: simulation_edges_ff6d3f9f4dd2.pkl (a-params hash, different from no-a-params)
# Used as the fourth column in plot_three_way_comparison.py.
#
# Usage:  bash submit_chindemi_fit_30pairs.sh [--dry-run]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="$SCRIPT_DIR/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations"
PAIRRUNNER_FIT="$SCRIPT_DIR/plastyfire/pairrunner_edges_fit.py"
CACHE="$SCRIPT_DIR/cpre_cpost_cache/chindemi_tau278.pkl"

# 30 pairs sampled by optimizer with seed=1234, sample_size=30
PAIRS=(
    "180329-189548" "180687-187024" "180723-207048" "180834-205125"
    "181002-188173" "182369-181607" "182967-198780" "183796-186882"
    "184709-196154" "185337-193957" "186117-200624" "186198-204134"
    "187212-192870" "190353-181954" "191341-181892" "194532-187415"
    "196734-181209" "197652-203951" "199599-194835" "200586-209048"
    "201258-187071" "201281-199420" "202571-195398" "203038-207594"
    "204120-188861" "204883-200512" "205300-191467" "208904-200676"
    "208966-196199" "208976-188476"
)
# mrk97_07 = +10ms, mrk97_08 = -10ms
PROTOCOLS=("10Hz_10ms" "10Hz_-10ms")

DRY_RUN=0
[[ "$1" == "--dry-run" ]] && DRY_RUN=1

ENV_BLOCK='
unset PATH PYTHONPATH LD_LIBRARY_PATH CMAKE_PREFIX_PATH
unset LIBRARY_PATH JPKGDIR NEURODAMUS_DIR HOC_LIBRARY_PATH
hash -r
module --force purge
module load StdEnv/2023 scipy-stack/2024a gcc/12.3 openmpi/4.1.5 hdf5-mpi/1.14.2 cmake/3.31.0 mpi4py/4.0.3 python/3.11.5
export JPKGDIR=/project/def-emuller/opt/jupyterhub-pkgs/build-07.04.2026-py311
export PYTHONPATH=/cvmfs/soft.computecanada.ca/easybuild/python/site-packages
export PYTHONPATH=$PYTHONPATH:/cvmfs/soft.computecanada.ca/custom/python/site-packages
export PYTHONPATH=$PYTHONPATH:$JPKGDIR/lib/python3.11/site-packages/
export PYTHONPATH=$PYTHONPATH:$JPKGDIR/lib/python
export PYTHONPATH=$PYTHONPATH:$JPKGDIR/lib/python3.11/site-packages/neurodamus/core/hoc
export PYTHONPATH=$PYTHONPATH:$HOME/.local/lib/python3.11/site-packages
export HOC_LIBRARY_PATH=$JPKGDIR/lib/python3.11/site-packages/neurodamus/data/hoc:$JPKGDIR/share/neurodamus_neocortex/hoc
export LD_LIBRARY_PATH=$JPKGDIR/lib
export PATH=$PATH:/opt/software/slurm/bin/:$JPKGDIR/bin
export PYTHONPATH=$PYTHONPATH:/lustre06/project/6077694/dhuruva/plastyfire
'

# Chindemi params passed to pairrunner_edges_fit.py
# Let pairrunner auto-compute hash from sorted(fit_params.items()) → ff6d3f9f4dd2
FIT_ARGS="--gamma_d_GB_GluSynapse=101.5 --gamma_p_GB_GluSynapse=216.2 \
    --a00=1.002 --a01=1.954 --a10=1.159 --a11=2.483 \
    --a20=1.127 --a21=2.456 --a30=5.236 --a31=1.782 \
    --tau_effca_GB_GluSynapse=278.3177658387 \
    --cpre-cpost-cache=$CACHE"

submitted=0
for pair in "${PAIRS[@]}"; do
    for proto in "${PROTOCOLS[@]}"; do
        workdir="$RESULTS_DIR/$pair/$proto"
        [[ ! -f "$workdir/prefire_simulation_config.json" ]] && { echo "SKIP: $workdir"; continue; }

        job_name="chinfit_${proto}_${pair}"
        log_file="chinfit_${proto}_${pair}.log"

        if [[ $DRY_RUN -eq 1 ]]; then
            echo "[dry-run] $pair / $proto"
        else
            sbatch \
                --job-name="$job_name" \
                --account=ctb-emuller \
                --cpus-per-task=2 \
                --ntasks=1 \
                --mem=4g \
                --time=01:30:00 \
                --chdir="$workdir" \
                --output="$log_file" \
                --wrap="$ENV_BLOCK
python $PAIRRUNNER_FIT $FIT_ARGS"
            ((submitted++))
        fi
    done
done

if [[ $DRY_RUN -eq 0 ]]; then
    echo "Submitted $submitted jobs (expected hash: ff6d3f9f4dd2)"
else
    echo "[dry-run] Would submit $((${#PAIRS[@]} * ${#PROTOCOLS[@]})) jobs"
fi
