#!/bin/bash
#SBATCH --job-name=v4_fit
#SBATCH --account=def-emuller
#SBATCH --gpus-per-node=nvidia_h100_80gb_hbm3_1g.10gb:1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --time=00:30:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/v4_fit_%j.out
# v4 post-rule fits (RHO_REDESIGN.md) on the pl5sj07 -vca dirs (29 targets). Env: TAG SEED FREE (extra free filters appended
# to the td4 set) SET (json, fixed v4 options) SEEDFITS SEEDSET (json) FAST (1 = tau_fast kernel) CHECK (1 = --check-v3) JOINT (1) VGATE (1)
# FITGAMMA (1 = --fit-gamma: free gamma_d 20-250, gamma_p 100-600). FIX (json {name: value}: --fix-params, e.g. {"a01": 1.17}).
# Sizing: v3 fits 22107773-6 (same data): 21.6-26.3 GB MaxRSS, ~2 min load + 4.2-5.8 s/gen -> 33G, 0:30 for 150 gens.
# v4 fit 22108745 (L5 only): 24.0 GB MaxRSS, 11:44 on 1g -> 30G, 0:30.
# JOINT=1 sizing: pilot (MAXITER=2) on 2g.20gb, 1 CPU, 45G, 0:15 (basis: L5-only v4 fit 22108745 = 24.0 GB MaxRSS, 11:44 on 1g; L2/3 CPU eval 22108820 ~12 GB). Fill real-job sizing from the pilot seff.
# MEASURED 22128396/7 (2g, JOINT 30 L5 + 9 L2/3 targets, 502 + 504 records): 25:05 / 22:28, 32.97 GB MaxRSS = at the 33G limit -> 42G, 0:45 (+sj04/L23L23 records scale it further).
# MEASURED 22132567/8 (2g, JOINT 32 L5 incl. sj04 + 9 L2/3, 526 + 504 records): 22:30 / 20:58, 29.3 / 39.48 GB MaxRSS (unseeded peak) -> 50G, 0:45.
# MEASURED pilot 22113910 (2g, JOINT, 2 gens): 25.9 GB MaxRSS, 3:04 wall (~2 min load, 1.4-6.8 s/gen) -> 33G, 0:30 on 2g (L5+L2/3 traces ~11 GB GPU, too big for 1g).
# MAXITER=0 repro (fit_v4n, JOINT, seed C1jf): basis pilot 22113910 (2 gens, 25.9 GB, 3:04) and C1jf 22121225 (24.1 GB, 17:02, 97% CPU) -> 33G, 0:15 on 2g.
# MEASURED --fix-params repro MAXITER=0 22133806/8 (2g, JOINT 30 + 9 targets): 2:41 / 2:42, 39.05 / 39.11 GB MaxRSS, 93% CPU -> 50G, 0:15 (load dominates memory).
# LADDER (MINIMAL_LADDER.md) 22135018 check 0:15, 22135019-22 fits 0:45, all 2g 50G (basis 22132568 / 22133806); fill seff here.
# LADDER 3-pathway (FITPY=fit_v4n.py EXTRA l23l23): 22135203 repro 2g 50G 0:15; 22135204 check, 22135205-8 fits on 3g 79G (63.3 GB of 22134473 + 25%), 0:15 / 1:45 (26 min x 2533/1006 records + 50%); fill seff here.
# LADDER M2/M3 (39 targets): 22135509-13 on 2g 50G 0:45 (basis 22132568 39.5 GB 25 min); fill seff here.
# GEOM_L23 (env, optional): pair-geometry csv for the L2/3->L5 model (read by fit_*.py; default ebner/pair_geometry_L23PC_L5TTPC.csv; alternative ebner/pair_geometry_L23PC_L5TTPC_geo.csv, LETZKUS_LOCATION.md).
# MEASURED 22134473 (3g.40gb, 3 pathways 30 + 9 + 15 targets, 526+504+1527 records, MAXITER=2 pilot): 6:01, 63.28 GB MaxRSS -> full 79G, 1:45 (39-target fit 25 min x 2.5 records + 50%).
set -euo pipefail
source glusynapse_v2/env_v3.sh
V2=glusynapse_v2; X=$V2/extracted
DIRS=${L5DIRS:-$X/ebner_delta-prefire-vca,$X/markram_delta-prefire-vca,$X/sj03_delta-prefire-vca,$X/sj03r50_delta-prefire-vca,$X/sj07_delta-prefire-vca}   # L5DIRS overrides (e.g. + sj04_delta-prefire-vca)
PAIRS=$(cat $V2/subset24_pairs.txt)
FILTERS='{"pre_drive": 1, "i_scale": 1e-5, "t_drive": 4, "tau_E1": 100.0}'
nvidia-smi -L
echo "=== v4 fit $TAG seed $SEED $(date)"
E0="{}"; SET=${SET:-$E0}; SEEDSET=${SEEDSET:-$E0}
ARGS=(--free-filters theta_Te,tau_T,theta_Tg,dpre_min,theta_NOi,tau_NO,tau_Z,theta_Z${FREE:+,$FREE} --set "$SET")
[ -n "${SEEDFITS:-}" ] && ARGS+=(--seed-fits $SEEDFITS --seed-set "$SEEDSET")
[ "${FAST:-0}" = 1 ] && ARGS+=(--fast)
[ "${CHECK:-0}" = 1 ] && ARGS+=(--check-v3)
[ "${JOINT:-0}" = 1 ] && ARGS+=(--joint)
[ -n "${DROPT:-}" ] && ARGS+=(--drop-targets "$DROPT")      # validation-only targets, e.g. letzkus_3ap_200hz_dt-10ms@distal|control
[ "${VGATE:-0}" = 1 ] && ARGS+=(--vgate)
[ "${FITGAMMA:-0}" = 1 ] && ARGS+=(--fit-gamma)
[ -n "${FIX:-}" ] && ARGS+=(--fix-params "$FIX")      # json {name: value} held fixed, out of the DE vector (PARAM_REDUCTION.md)
# FITPY=fit_v4n.py: N-pathway fitter (PROTOCOLS.md); EXTRA="NAME:GROUPS:DIRS:BASIS[:GEOM[:PAIRS]];..." adds pathway models (fit_v4n only).
FITPY=${FITPY:-fit_v4.py}
[ -n "${L23DIRS:-}" ] && ARGS+=(--l23-dirs $L23DIRS)      # other emodel: L23DIRS + L23_BASIS_DIR + ANALYTICAL_BASIS_DIR (L5) env
if [ -n "${EXTRA:-}" ]; then IFS=";" read -ra EXL <<< "$EXTRA"; for e in "${EXL[@]}"; do ARGS+=(--extra "$e"); done; fi
python -u $V2/rho_redesign/$FITPY --dirs $DIRS --pairs $PAIRS --groups ${L5GROUPS:-paired_l5,sjostrom07} --filters "$FILTERS" \
    --tie '{"theta_NOi": "theta_Ti"}' --x0 fit_results/delta-cooker.json --maxiter ${MAXITER:-150} --seed $SEED \
    --save $V2/rho_redesign/results/v4_$TAG "${ARGS[@]}"
echo "=== done $(date)"
