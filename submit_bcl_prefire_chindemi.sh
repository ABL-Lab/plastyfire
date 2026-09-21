#!/bin/bash
# Submit BCL prefire sims for Chindemi params (no a-params, no fastforward).
# Runs pairrunner_edges.py --params chindemi for all 100 pairs x 2 protocols.
# Results go to bluecellulab_results_prefire/ inside each workdir.
# This gives the "BCL prefire" leg of the 3-way comparison:
#   BCL (full 520s, explicit test pulses)  vs
#   BCL prefire (42s induction, basis EPSP) vs
#   Optimizer (same prefire but fitted params, basis EPSP)
#
# Usage:
#   ./submit_bcl_prefire_chindemi.sh [--dry-run] [--dependency=afterok:JOBID]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="$SCRIPT_DIR/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations"
PAIRRUNNER="$SCRIPT_DIR/plastyfire/pairrunner_edges.py"
# No --cpre-cpost-cache: chindemi has no a-params (theta_d/theta_p from edges.h5)
# No --fastforward: prefire tstop (42s) IS the induction; no fast-forward needed
BCL_SUBDIR="${BCL_SUBDIR_OVERRIDE:-bluecellulab_results_prefire}"
PROTOCOLS=("10Hz_10ms" "10Hz_-10ms" "10Hz_5ms" "10Hz_30ms" "10Hz_-30ms" "10Hz_50ms" "10Hz_-50ms")

DRY_RUN=0
DEPENDENCY=""

for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=1 ;;
        --dependency=*) DEPENDENCY="${arg#--dependency=}" ;;
    esac
done

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

submitted=0

for pair_dir in "$RESULTS_DIR"/*/; do
    pair=$(basename "$pair_dir")
    for proto in "${PROTOCOLS[@]}"; do
        workdir="$pair_dir$proto"
        [[ ! -f "$workdir/prefire_simulation_config.json" ]] && continue

        job_name="prefire_${proto}_${pair}"
        log_file="prefire_${proto}_${pair}.log"

        RUN_CMD="python $PAIRRUNNER --params chindemi --bcl-subdir $BCL_SUBDIR"

        SBATCH_ARGS=(
            --job-name="$job_name"
            --account=ctb-emuller
            --cpus-per-task=2
            --ntasks=1
            --mem=4g
            --time=01:30:00
            --chdir="$workdir"
            --output="$log_file"
        )
        [[ -n "$DEPENDENCY" ]] && SBATCH_ARGS+=(--dependency="$DEPENDENCY")

        if [[ $DRY_RUN -eq 1 ]]; then
            echo "[dry-run] $pair / $proto  ->  $workdir"
        else
            sbatch "${SBATCH_ARGS[@]}" --wrap="$ENV_BLOCK
$RUN_CMD"
            ((submitted++))
        fi
    done
done

if [[ $DRY_RUN -eq 0 ]]; then
    echo "Submitted $submitted jobs."
else
    echo "[dry-run] Would submit $submitted jobs."
fi
