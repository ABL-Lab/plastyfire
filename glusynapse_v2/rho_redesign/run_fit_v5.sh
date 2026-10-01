#!/bin/bash
#SBATCH --job-name=v5_fit
#SBATCH --account=def-emuller
#SBATCH --gpus-per-node=nvidia_h100_80gb_hbm3_2g.20gb:1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=00:45:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/v5_fit_%j.out
# v5 shared-pool rule (V5_DESIGN.md; gpu_v5_rho / fit_v5). Copy of run_fit_v4.sh (not modified). Env: TAG SEED
# FREE (free filters, default theta_V) SET (json merged into the v5 filters) SEEDFITS SEEDSET CHECK (1 = --check-v4)
# JOINT (1) FITGAMMA (1) DROPT MAXITER EXTRA ("NAME:GROUPS:DIRS:BASIS[:GEOM[:PAIRS]];...") L5DIRS L5GROUPS.
# Sizing basis (same data and loader as the v4 joint fits; the v5 kernel holds fewer registers, the DE vector is 7 not 18):
#   check (MAXITER=0): 22133806 (2g, JOINT 30 + 9 targets, MAXITER 0): 2:42, 39.1 GB MaxRSS -> 50G, 0:15.
#   fits: 22132568 / 22133630 (2g, JOINT, 39 targets, 18 free): 39.4-39.5 GB, 20-26 min -> 50G, 0:45.
#   v5 MAXITER 300 x 56 members ~ 2.3 x fewer evaluations per generation than v4 (144) at 150 gens x 2 -> same budget.
# GEOM_L23 (env, optional): alternative pair-geometry csv for the L2/3->L5 model (default ebner/pair_geometry_L23PC_L5TTPC.csv).
# MEASURED fits (2g, 39 targets, MAXITER 300): V5_s5 22135814 12:07 32.5 GB, V5_s6 22135815 7:14 27.2 GB, V5b_s5 22135953 13:03 27.2 GB, V5c_s5 22135954 13:57 39.1 GB (98-99% CPU) -> 50G, 0:30 (600 gens 0:45). 3-path check V5e2zchk 22136256 (3g): 5:17, 54.0 GB, 90% CPU.
# JOBS 2026-10-01: check 22135813 (V5chk), fits 22135814 (V5_s5 seeded) / 22135815 (V5_s6 unseeded) afterok. MEASURED V5chk 22135813: 2:44, 31.6 GB MaxRSS, 98% CPU (2-pathway check). 3-pathway (EXTRA l23l23): 3g, 79G (63.0 GB of 22135207), fits 0:45, check 0:15: 22136256-9.
set -euo pipefail
source glusynapse_v2/env_v3.sh
V2=glusynapse_v2; X=$V2/extracted
DIRS=${L5DIRS:-$X/ebner_delta-prefire-vca,$X/markram_delta-prefire-vca,$X/sj03_delta-prefire-vca,$X/sj03r50_delta-prefire-vca,$X/sj07_delta-prefire-vca}
PAIRS=$(cat $V2/subset24_pairs.txt)
# v5 constants: tau_E1 100 ms (E bracket), dpre_min -0.29 (E, Sjostrom 2003 CB1-agonist ceiling), A_eCB 1 (saturating
# limit), rho_gamma 1 (Chindemi), vamp_mode 1 (C1 gate form, also compiles the v4 kernel for --check-v4), t_drive 4
# (same loader as the v4 fits; the T chain is unused by the v5 kernel).
FILTERS='{"pre_drive": 1, "i_scale": 1e-5, "t_drive": 4, "tau_E1": 100.0, "vamp_mode": 1, "v5_mode": 1, "rho_gamma": 1.0, "A_eCB": 1.0, "dpre_min": -0.29}'
nvidia-smi -L
echo "=== v5 fit $TAG seed $SEED $(date)"
E0="{}"; SET=${SET:-$E0}; SEEDSET=${SEEDSET:-$E0}
ARGS=(--free-filters ${FREE:-theta_V} --set "$SET")
[ -n "${SEEDFITS:-}" ] && ARGS+=(--seed-fits $SEEDFITS --seed-set "$SEEDSET")
[ "${CHECK:-0}" = 1 ] && ARGS+=(--check-v4)
[ "${JOINT:-0}" = 1 ] && ARGS+=(--joint)
[ -n "${DROPT:-}" ] && ARGS+=(--drop-targets "$DROPT")
[ "${FITGAMMA:-0}" = 1 ] && ARGS+=(--fit-gamma)
[ -n "${L23DIRS:-}" ] && ARGS+=(--l23-dirs $L23DIRS)      # other emodel: L23DIRS + L23_BASIS_DIR + ANALYTICAL_BASIS_DIR (L5) env
if [ -n "${EXTRA:-}" ]; then IFS=";" read -ra EXL <<< "$EXTRA"; for e in "${EXL[@]}"; do ARGS+=(--extra "$e"); done; fi
python -u $V2/rho_redesign/fit_v5.py --dirs $DIRS --pairs $PAIRS --groups ${L5GROUPS:-paired_l5,sjostrom07} --filters "$FILTERS" \
    --x0 fit_results/delta-cooker.json --maxiter ${MAXITER:-300} --seed $SEED \
    --save $V2/rho_redesign/results/v5_$TAG "${ARGS[@]}"
echo "=== done $(date)"
