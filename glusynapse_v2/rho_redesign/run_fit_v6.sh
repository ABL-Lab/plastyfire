#!/bin/bash
#SBATCH --job-name=v6_fit
#SBATCH --account=def-emuller
#SBATCH --gpus-per-node=nvidia_h100_80gb_hbm3_2g.20gb:1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=00:30:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/v6_fit_%j.out
# v6 = v5c with per-synapse thresholds in measured Ca quantities (EXP_PARAM.md; gpu_v6_rho / fit_v6). Copy of
# run_fit_v5.sh (not modified). Env as run_fit_v5.sh (TAG SEED FREE SET SEEDFITS SEEDSET JOINT FITGAMMA DROPT MAXITER
# EXTRA L5DIRS L5GROUPS L23DIRS WEIGHTS) plus, read by fit_v6.py from the env: CEXP_DIR (dir of <pathway>.csv, or
# STANDIN), SCALE_V / SCALE_E (json {col: coef}), THETA_PRE (json {"d": {...}, "p": {...}}), BOX (json), CEXP_MAP.
# STANDIN_SMOKE=1: after the main run, a second run (same data, CEXP_DIR=STANDIN, SCALE_V/E 48.3 x cpost_cai,
# THETA_PRE p = cpre_mg0, MAXITER 2) smoke-tests the join + scaled kernel; its numbers mean nothing.
# Sizing (measured v5, run_fit_v5.sh header): rescore 2:44 31.6 GB -> 40G 0:15 (x2 with STANDIN_SMOKE, still < 0:15);
#   fits 12-15 min 27-39 GB -> 50G 0:30.
# MEASURED S1_V6r 22155608 (repro bit-identical + STANDIN smoke, 2 runs): 5:30, 37.40 GB of 40G -> rescores 47G 0:15.
# MEASURED 2026-10-01 (2g, 39 targets): relabel 22156526 3:16 37.25 GB 73% CPU; fits E2 22156527/8 12:27/12:24 25.9 GB 99%; E1 22156531/2 13:57/13:03 32.8/26.7 GB 91-94%; E3 22156559/60 3:06/2:24 26.9/27.4 GB (stopped early, all-penalty DE) -> fits 41G 0:30, rescores 47G 0:15.
# MEASURED E3n (normalised ratio): rescore 22159388 2:46 28.8 GB; refit 22159389 4:04 37.2 GB (DE on penalty, E3 dropped).
set -euo pipefail
source glusynapse_v2/env_v3.sh
V2=glusynapse_v2; X=$V2/extracted
DIRS=${L5DIRS:-$X/ebner_delta-prefire-vca,$X/markram_delta-prefire-vca,$X/sj03_delta-prefire-vca,$X/sj03r50_delta-prefire-vca,$X/sj07_delta-prefire-vca}
PAIRS=$(cat $V2/subset24_pairs.txt)
FILTERS='{"pre_drive": 1, "i_scale": 1e-5, "t_drive": 4, "tau_E1": 100.0, "vamp_mode": 1, "v5_mode": 1, "rho_gamma": 1.0, "A_eCB": 1.0, "dpre_min": -0.29}'
nvidia-smi -L
echo "=== v6 fit $TAG seed $SEED $(date)"
E0="{}"; SET=${SET:-$E0}; SEEDSET=${SEEDSET:-$E0}
ARGS=(--free-filters ${FREE:-theta_V} --set "$SET")
[ -n "${SEEDFITS:-}" ] && ARGS+=(--seed-fits $SEEDFITS --seed-set "$SEEDSET")
[ "${JOINT:-0}" = 1 ] && ARGS+=(--joint)
[ -n "${DROPT:-}" ] && ARGS+=(--drop-targets "$DROPT")
[ "${FITGAMMA:-0}" = 1 ] && ARGS+=(--fit-gamma)
[ -n "${L23DIRS:-}" ] && ARGS+=(--l23-dirs $L23DIRS)
if [ -n "${EXTRA:-}" ]; then IFS=";" read -ra EXL <<< "$EXTRA"; for e in "${EXL[@]}"; do ARGS+=(--extra "$e"); done; fi
run() { python -u $V2/rho_redesign/fit_v6.py --dirs $DIRS --pairs $PAIRS --groups ${L5GROUPS:-paired_l5,sjostrom07} \
    --filters "$FILTERS" --x0 fit_results/delta-cooker.json --seed $SEED "${ARGS[@]}" "$@"; }
run --maxiter ${MAXITER:-300} --save $V2/rho_redesign/results/v6_$TAG
if [ "${STANDIN_SMOKE:-0}" = 1 ]; then
  echo "=== STAND-IN smoke (not a measurement) $(date)"
  CEXP_DIR=STANDIN SCALE_V='{"cpost_cai": 48.3}' SCALE_E='{"cpost_cai": 48.3}' THETA_PRE='{"p": {"cpre_mg0": 1.0}}' \
    run --maxiter 2 --save /scratch/dhuruva/split1/cexp_standin/v6_${TAG}_standin
fi
echo "=== done $(date)"
