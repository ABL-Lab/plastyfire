#!/bin/bash
# HOW TO RUN (v3 objective, FIT_SPEED.md). From /project/rrg-emuller/dhuruva/plastyfire:
#   env MODE=subset TGROUPS=paired_l5,sjostrom07 CONDS=control,mglu_block,post_nmdar,nmdar_block,no_block MAXITER=150 SEED=1 \
#     TAG=_pl5sj07_td4_s1 DIRS=<comma list of extracted dirs> FREE=theta_Te,tau_T,theta_Tg,dpre_min,theta_NOi,tau_NO,tau_Z,theta_Z \
#     FILTERS='{"pre_drive": 1, "i_scale": 1e-5, "t_drive": 4, "tau_E1": 100.0}' SEEDFITS=<fit json> [RESUME=<v2 or v3 _ckpt.npz>] \
#     sbatch glusynapse_v2/run_reduced_gpu_v3.sh
# Defaults below: 1g.10gb, 1 CPU, 33G, 0:30 for 150 generations (seff 22107773-6: 21.6-26.3 GB MaxRSS, 4.2-5.8 s/gen, ~2 min load).
# Env: glusynapse_v2/env_v3.sh (own venv .venv_v3 with numba-cuda; plastyfire/.venv stays for v2). Rule changes: edit both
# jax_v2.py and gpu_v3.py, then pass run_smoke_gpu_v3.sh and run_test_jax_v3.sh (TDRIVE=4 and 3).
# MEASURED 22110251/2 (ljp25, 150 gens): 12 min, 26.0 GB MaxRSS, 96% CPU -> 33G 0:30 stays right.
#SBATCH --job-name=v3_fit
#SBATCH --account=def-emuller
#SBATCH --gpus-per-node=nvidia_h100_80gb_hbm3_1g.10gb:1
#SBATCH --cpus-per-task=1
#SBATCH --mem=33G
#SBATCH --time=00:30:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/v3_fit_%j.out
# Drop-in for run_reduced_gpu.sh with the v3 objective (fit_v3.py + gpu_v3.py; FIT_SPEED.md). Same env vars:
# MODE (preview | subset | all), DIRS, FILTERS, TIE, FREE, TAG, TGROUPS, CONDS, MAXITER, SEED, SEEDFITS, RESUME.
# CHUNK is ignored (traces stay on the GPU). Pass --mem/--time on the sbatch line, sized from seff.
# 22105400 (pl5sj07 td4 subset, failed at the first kernel compile after build + freeing host traces): 24.0 GB MaxRSS, 1:52 -> 30G for a fit.
set -euo pipefail
source glusynapse_v2/env_v3.sh
V2=glusynapse_v2
if [ "$MODE" = preview ]; then
  D=$V2/extracted/ebner_preview,$V2/extracted/markram_delta-cooker
  PAIRS=180351-198084,181015-184976,182339-200396,182381-186941,182900-195902,183240-189246,183769-207759,186261-208490
else
  D=$V2/extracted/ebner_delta-prefire,$V2/extracted/markram_delta-prefire-tr
  PAIRS=$(cat $V2/subset24_pairs.txt)
fi
D=${DIRS:-$D}
EXTRA=()
if [ "$MODE" = all ]; then PAIRS=; fi
if [ -n "${SEEDFITS:-}" ]; then EXTRA+=(--seed-fits $SEEDFITS); fi
if [ -n "${RESUME:-}" ]; then EXTRA+=(--resume $RESUME); fi
FILTERS=${FILTERS:-'{"pre_drive": 1, "i_scale": 1e-5, "tau_T": 3.0}'}
TIE=${TIE:-'{"theta_NOi": "theta_Ti"}'}
FREE=${FREE:-theta_Ti,tau_NO,tau_Z,theta_Z}
TAG=${TAG:-}
echo "filters $FILTERS tie $TIE free $FREE"
nvidia-smi -L
echo "=== reduced fit GPU v3 ($MODE) $(date)"
python -u $V2/fit_v3.py --groups ${TGROUPS:-markram,nevian,ebner} --dirs $D ${PAIRS:+--pairs $PAIRS} "${EXTRA[@]}" --conditions ${CONDS:-control,mglu_block,post_nmdar,nmdar_block} --maxiter ${MAXITER:-80} --seed ${SEED:-1} \
    --filters "$FILTERS" --tie "$TIE" --free-filters $FREE \
    --x0 fit_results/delta-cooker.json --save $V2/results/reduced_gpu_$MODE$TAG
echo "=== done $(date)"
