#!/bin/bash
#SBATCH --job-name=w2_joint
# KERNEL (env, default gpu_v7_rho.py; gpu_v8_rho.py = gate_src options, A14).
# WAVE 2 joint fit (A9): copy of run_fit_v7.sh; env (L5DIRS L5GROUPS L23DIRS EXTRA JOINT FITGAMMA DROPT SET FREE) is passed by submit_joint_w2.sh. Outputs results/v7_W2_<tag>.json (TAG=W2_*).
# MEASURED (3g, 3 models, 66 targets): rescore 22159363 6:45 61.76 GB, 22159364 6:58 50.23 GB -> fits 78G 0:30 (MaxRSS 61.8+25%). Rescores exit 1 on REPRO DIFF (target set differs from seed) -> do not chain refits with afterok.
# MEASURED fits: 22161286 16:17 39.69 GB, 22161287 17:23 50.37 GB (78G ok; next 63G 0:30)
# Sizing: 3-model measured v4 22135207 (502+504+1527 rec) 11:41, MaxRSS 66.1 GB, GPU 7.1+9.1+13.5=29.7 GB (> 2g 20 GB -> 3g.40gb). Wave 2 ~2900 rec (x1.15): 76 GB +25% = 95G; fits 0:30, rescores 0:15.
#SBATCH --account=def-emuller
#SBATCH --gpus-per-node=nvidia_h100_80gb_hbm3_3g.40gb:1
#SBATCH --cpus-per-task=1
#SBATCH --mem=95G
#SBATCH --time=00:30:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/w2_joint_%j.out
# v7 = v6 (gpu_v6_rho / fit_v6, not modified) + LTD options in SET: veto_T (ms, bAP veto of the eCB step) and ecb_ref 1
# (theta_eCB in units of the own single-bAP VDCC pool jump); gpu_v7_rho.py runs fit_v6.main. Env as run_fit_v5.sh (TAG SEED FREE SET SEEDFITS SEEDSET JOINT FITGAMMA DROPT MAXITER
# EXTRA L5DIRS L5GROUPS L23DIRS WEIGHTS) plus, read by fit_v6.py from the env: CEXP_DIR (dir of <pathway>.csv, or
# STANDIN), SCALE_V / SCALE_E (json {col: coef}), THETA_PRE (json {"d": {...}, "p": {...}}), BOX (json), CEXP_MAP.
# STANDIN_SMOKE=1: after the main run, a second run (same data, CEXP_DIR=STANDIN, SCALE_V/E 48.3 x cpost_cai,
# THETA_PRE p = cpre_mg0, MAXITER 2) smoke-tests the join + scaled kernel; its numbers mean nothing.
# Sizing: rescore = A5 22155608 (v6 rescore + smoke) 5:30, MaxRSS 39.2 GB -> 50G 0:15; fits (v5 12-15 min, 27-39 GB) 50G 0:30.

# MEASURED S1_V6r 22155608 (repro bit-identical + STANDIN smoke, 2 runs): 5:30, 37.40 GB of 40G -> rescores 47G 0:15.
set -euo pipefail
source glusynapse_v2/env_v3.sh
V2=glusynapse_v2; X=$V2/extracted
DIRS=${L5DIRS:-$X/ebner_delta-prefire-vca,$X/markram_delta-prefire-vca,$X/sj03_delta-prefire-vca,$X/sj03r50_delta-prefire-vca,$X/sj07_delta-prefire-vca}
PAIRS=$(cat $V2/subset24_pairs.txt)
FILTERS='{"pre_drive": 1, "i_scale": 1e-5, "t_drive": 4, "tau_E1": 100.0, "vamp_mode": 1, "v5_mode": 1, "rho_gamma": 1.0, "A_eCB": 1.0, "dpre_min": -0.29}'
nvidia-smi -L
echo "=== v7 fit $TAG seed $SEED $(date)"
E0="{}"; SET=${SET:-$E0}; SEEDSET=${SEEDSET:-$E0}
ARGS=(--free-filters ${FREE:-theta_V} --set "$SET")
[ -n "${SEEDFITS:-}" ] && ARGS+=(--seed-fits $SEEDFITS --seed-set "$SEEDSET")
[ "${JOINT:-0}" = 1 ] && ARGS+=(--joint)
[ -n "${DROPT:-}" ] && ARGS+=(--drop-targets "$DROPT")
[ "${FITGAMMA:-0}" = 1 ] && ARGS+=(--fit-gamma)
[ -n "${L23DIRS:-}" ] && ARGS+=(--l23-dirs $L23DIRS)
if [ -n "${EXTRA:-}" ]; then IFS=";" read -ra EXL <<< "$EXTRA"; for e in "${EXL[@]}"; do ARGS+=(--extra "$e"); done; fi
run() { python -u $V2/rho_redesign/${KERNEL:-gpu_v7_rho.py} --dirs $DIRS --pairs $PAIRS --groups ${L5GROUPS:-paired_l5,sjostrom07} \
    --filters "$FILTERS" --x0 fit_results/delta-cooker.json --seed $SEED "${ARGS[@]}" "$@"; }
run --maxiter ${MAXITER:-300} --save $V2/rho_redesign/results/v7_$TAG
if [ "${STANDIN_SMOKE:-0}" = 1 ]; then
  echo "=== STAND-IN smoke (not a measurement) $(date)"
  CEXP_DIR=STANDIN SCALE_V='{"cpost_cai": 48.3}' SCALE_E='{"cpost_cai": 48.3}' THETA_PRE='{"p": {"cpre_mg0": 1.0}}' \
    run --maxiter 2 --save /scratch/dhuruva/split1/cexp_standin/v6_${TAG}_standin
fi
echo "=== done $(date)"
