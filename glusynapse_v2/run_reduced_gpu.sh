#!/bin/bash
#SBATCH --job-name=v22r_gpu
#SBATCH --account=def-emuller      # rrg-emuller has no GPU allocation
#SBATCH --gpus-per-node=h100:1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/v22r_gpu_%j.out
# Reduced v2.2 fit, as run_reduced.sh but with the GPU objective (jax_v2) and continuous filters.
# Pass --mem on the sbatch line: preview 40G (measured 27G), subset 85G (66 GiB), all 285G (226 GiB).
# MODE: preview | subset | all (all 120 pairs, streamed; CHUNK = GB per chunk, default 8).
# Variants via env: FILTERS (fixed filters json), TIE (json, '{}' = none), FREE (free filters), TAG (results suffix),
# TGROUPS (not GROUPS: bash builtin = numeric group ids; targets.py groups; T25 paired-only: paired_l5).
set -euo pipefail
source .venv/bin/activate
V2=glusynapse_v2
if [ "$MODE" = preview ]; then
  D=$V2/extracted/ebner_preview,$V2/extracted/markram_delta-cooker
  PAIRS=180351-198084,181015-184976,182339-200396,182381-186941,182900-195902,183240-189246,183769-207759,186261-208490
else
  D=$V2/extracted/ebner_delta-prefire,$V2/extracted/markram_delta-prefire-tr
  PAIRS=$(cat $V2/subset24_pairs.txt)
fi
D=${DIRS:-$D}                             # other extracted data (e.g. the delta-sv spine variant)
EXTRA=()
if [ "$MODE" = all ]; then                # all 120 pairs: traces exceed GPU memory, stream them in chunks
  PAIRS=; EXTRA=(--chunk-gb ${CHUNK:-8})
elif [ -n "${CHUNK:-}" ]; then            # stream a subset too (MIG slice)
  EXTRA=(--chunk-gb $CHUNK)
fi
if [ -n "${SEEDFITS:-}" ]; then          # fit jsons put into the initial DE population
  EXTRA+=(--seed-fits $SEEDFITS)
fi
if [ -n "${RESUME:-}" ]; then            # <save>_ckpt.npz of an earlier (timed-out) run
  EXTRA+=(--resume $RESUME)
fi
FILTERS=${FILTERS:-'{"pre_drive": 1, "i_scale": 1e-5, "tau_T": 3.0}'}
TIE=${TIE:-'{"theta_NOi": "theta_Ti"}'}
FREE=${FREE:-theta_Ti,tau_NO,tau_Z,theta_Z}
TAG=${TAG:-}
echo "filters $FILTERS tie $TIE free $FREE"
nvidia-smi --query-gpu=name,memory.total --format=csv
echo "=== reduced fit GPU ($MODE) $(date)"
python -u $V2/fit_v2.py --backend jax --groups ${TGROUPS:-markram,nevian,ebner} --dirs $D ${PAIRS:+--pairs $PAIRS} "${EXTRA[@]}" --conditions ${CONDS:-control,mglu_block,post_nmdar,nmdar_block} --model v2 --maxiter ${MAXITER:-80} --seed ${SEED:-1} \
    --filters "$FILTERS" --tie "$TIE" --free-filters $FREE \
    --x0 fit_results/delta-cooker.json --save $V2/results/reduced_gpu_$MODE$TAG
nvidia-smi --query-gpu=memory.used --format=csv || true
echo "=== done $(date)"
