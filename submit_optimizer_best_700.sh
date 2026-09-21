#!/bin/bash
# Submit BCL prefire simulations (pairrunner_edges_fit.py) with optimizer best params
# for all 100 pairs × 7 protocols = 700 jobs.
#
# Optimizer best (hof[0] from checkpoint_edges.pkl, gen 21):
#   gamma_d=99.18727050237078  gamma_p=207.72948705417807
#   a-params: a00=1.0003856191061014  a01=1.2884800071747866
#             a10=1.1552690604404767  a11=2.4739773455167264
#             a20=1.0152594489565097  a21=2.4049220533677254
#             a30=5.613745269577143   a31=1.4748047556236923
#   tau=278.3177658387
#   → hash 23cd68d13dc3
#
# Output per workdir: simulation_edges_23cd68d13dc3.pkl
# Used by plot_stdp_optimizer_best.py to produce STDP curve vs in vitro.
#
# Usage:  bash submit_optimizer_best_700.sh [--dry-run]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="$SCRIPT_DIR/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations"
PAIRRUNNER_FIT="$SCRIPT_DIR/plastyfire/pairrunner_edges_fit.py"
CACHE="$SCRIPT_DIR/cpre_cpost_cache/chindemi_tau278.pkl"
PARAM_HASH="23cd68d13dc3"

PROTOCOLS=("10Hz_10ms" "10Hz_-10ms" "10Hz_5ms" "10Hz_30ms" "10Hz_-30ms" "10Hz_50ms" "10Hz_-50ms")

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

FIT_ARGS="--gamma_d_GB_GluSynapse=99.18727050237078 --gamma_p_GB_GluSynapse=207.72948705417807 \
    --a00=1.0003856191061014 --a01=1.2884800071747866 \
    --a10=1.1552690604404767 --a11=2.4739773455167264 \
    --a20=1.0152594489565097 --a21=2.4049220533677254 \
    --a30=5.613745269577143  --a31=1.4748047556236923 \
    --tau_effca_GB_GluSynapse=278.3177658387 \
    --param_hash=$PARAM_HASH \
    --cpre-cpost-cache=$CACHE"

submitted=0
skipped=0

for pair_dir in "$RESULTS_DIR"/*/; do
    pair=$(basename "$pair_dir")
    for proto in "${PROTOCOLS[@]}"; do
        workdir="$pair_dir$proto"
        [[ ! -f "$workdir/prefire_simulation_config.json" ]] && { ((skipped++)); continue; }

        job_name="optb_${proto}_${pair}"
        log_file="optbest_${proto}_${pair}.log"

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
    echo "Submitted $submitted jobs  (skipped $skipped missing workdirs)  hash: $PARAM_HASH"
else
    total=$(( ${#PROTOCOLS[@]} * 100 - skipped ))
    echo "[dry-run] Would submit ~$total jobs  (skipped $skipped missing workdirs)"
fi
