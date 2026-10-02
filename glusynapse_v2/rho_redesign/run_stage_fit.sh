#!/bin/bash
# Stage 1 of the staged fit (CORE_TARGETS.md): the 7 L5 STDP-shape targets only, free a00 a01 a10 a11 gamma_d gamma_p (+ theta_eCB for C).
# Copy of run_r2_fit.sh (kernel gpu_v10_rho.py, split2 inputs, same L5 keep-list pairs). Usage (from plastyfire/, via sjob.sh -g):
# MEASURED stage-1 L5-only fits s1A 22266770/72: 1-1.6 min, 5.7-7.0 GB, so fits use 9G 0:15:00. Rescores with all three models stay 85G 0:20:00 (v10 tests 58-68 GiB).
# MEASURED stage-2 (11 targets) 22286005-17: fits 1:20-2:20, MaxRSS 7.1-8.5 GB, so fits now 11G 0:15:00; val rescores 7:29-8:23, MaxRSS 51-67 GB, so val 84G 0:15:00.
#   bash glusynapse_v2/rho_redesign/run_stage_fit.sh MODEL SEEDMODE [val]
#   MODEL  A pure Chindemi : SET {v5_mode 1, theta_V 0, A_eCB 0}, no gate_src / lic_src / ecb_src / no_mode / veto_T / ecb_ref -> lic -1, ecb 0, nom 0, theta_V <= 0 = no gate, A_eCB 0 = no eCB step = v4 M0 kernel path
#          B Chindemi + shaft licence (G8X gate_src 2, gate_win 100, gate_theta 0.15, v5_mode 2, t_exact 1), eCB off (A_eCB 0), NO off, no veto
#          C R1E rule : G8X SET unchanged (shaft licence, eCB on, ecb_ref 2, veto_T 25), free theta_eCB (default box), NO off
#   SEEDMODE s (DE seeded: A,B from s1_seed_cooker.json = delta-cooker a's + gammas; C from r1E_s.json; seed 5) | u (fully unseeded latin hypercube, no x0, seed 6)
#   val    rescore s1<M>_<S>.json with --maxiter 0 on ALL targets (L5 40 + l23 + l23l23, only the R2 DROPT dropped) -> s1<M>_<S>_val.*
#          (REPRO DIFF is expected there: the json was fitted on 7 targets; the script accepts exit 1 if the val json was written)
# Targets: the l5 target list = first two columns of r1D_s.csv minus the 7 core (KEEP below) is passed to --drop-targets; the log must say "over 7 targets".
# Output: glusynapse_v2/rho_redesign/s1<M>_<S>.{json,csv,_ckpt.npz} (+ _val.*). Grep: "pairs|chi2 L5|over 7 targets|chi2 total|REPRO|Traceback|v6 free|seed member".
# Sizing: UPPER BOUND, unmeasured for L5-only: fit 85G 0:45:00 (R1 all-model fits measured 72.9 GB, 27 min); val 85G 0:20:00 (v10 tests 58-68 GiB, 8-11 min).
set -euo pipefail
M=$1; S=$2; VAL=${3:-}
R=/project/rrg-emuller/dhuruva/plastyfire; V2=glusynapse_v2; RS=$V2/rho_redesign; W=/scratch/dhuruva/s2g0321; X=$W/extracted; T=delta-split2-ljp25g0321-prefire
source glusynapse_v2/env_v3.sh
export CEXP_DIR=$W/cexp ANALYTICAL_BASIS_DIR=$W/basis_l5l5 L23_BASIS_DIR=$W/basis_l23l5 PYTHONPATH=$R
DIRS=$X/ebner_$T-vca,$X/markram_$T-vca,$X/sj03_$T-vca,$X/sj03r50_$T-vca,$X/sj07_$T-vca,$X/l5extra_$T-vca
L23D=$X/ebner_l23l5_$T-vseg-rs
Z=$X/zilberter_l23l23_$T-vseg-rs,$X/l23l23extra_$T-vseg-rs
GEOM_L23=${GEOM_L23:-$R/ebner/pair_geometry_L23PC_L5TTPC.csv}
K5=$W/doublets/l5_keep_pairs.txt; K235=$W/doublets/l23l5_keep_pairs.txt; K2323=$W/doublets/l23l23_keep_pairs.txt
for f in $K5 $K235 $K2323; do [ -s $f ] || { echo "missing keep list $f"; exit 2; }; done
PAIRS=$(comm -12 <(tr ',' '\n' < $V2/subset24_pairs.txt | sed '/^$/d' | sort -u) <(tr ',' '\n' < $K5 | sed '/^$/d' | sort -u) | paste -sd,)
echo "pairs: l5 template $(tr ',' '\n' < $V2/subset24_pairs.txt | sed '/^$/d' | sort -u | wc -l) keep $(tr ',' '\n' < $K5 | sed '/^$/d' | wc -l) used $(echo $PAIRS | tr ',' '\n' | wc -l); l23l5 keep $(tr ',' '\n' < $K235 | sed '/^$/d' | wc -l); l23l23 keep $(tr ',' '\n' < $K2323 | sed '/^$/d' | wc -l)"
[ -n "$PAIRS" ] || { echo "empty L5 pair set"; exit 2; }
G8X='"v5_mode": 2, "veto_T": 25.0, "ecb_ref": 2, "gate_src": 2, "gate_win": 100.0, "gate_theta": 0.15, "theta_V": 3.5751446275588403, "t_exact": 1, "v10_counts": 1'
BOX=
case $M in
  A) SET='{"v5_mode": 1, "theta_V": 0.0, "A_eCB": 0.0}'; FREE=""; SEEDJ=$RS/s1_seed_cooker.json ;;
  B) SET='{"v5_mode": 2, "gate_src": 2, "gate_win": 100.0, "gate_theta": 0.15, "theta_V": 3.5751446275588403, "t_exact": 1, "veto_T": 25.0, "ecb_ref": 2, "v10_counts": 1, "A_eCB": 0.0}'; FREE=""; SEEDJ=$RS/s1_seed_cooker.json ;;
  C) SET="{$G8X}"; FREE=theta_eCB; SEEDJ=$RS/r1E_s.json ;;
  2C) SET="{$G8X, \"no_mode\": 3, \"tau_NO\": 6.7, \"A_NO\": 190.0, \"d_NO_max\": 0.5}"; FREE=theta_eCB,A_NO,d_NO_max; SEEDJ=$RS/s1C_s.json ;;
  2E) SET="{$G8X}"; FREE=theta_eCB; SEEDJ=$RS/s1C_s.json ;;
  2N) SET='{"v5_mode": 2, "veto_T": 25.0, "ecb_ref": 2, "theta_V": 0.0, "t_exact": 1, "v10_counts": 1, "no_mode": 3, "tau_NO": 6.7, "A_NO": 190.0, "d_NO_max": 0.5}'; FREE=theta_eCB,A_NO,d_NO_max; SEEDJ=$RS/s1A_s.json ;;
  *) echo "MODEL A|B|C (stage 1) or 2C|2E|2N (stage 2)" >&2; exit 1 ;;
esac
export BOX
FILTERS='{"pre_drive": 1, "i_scale": 1e-5, "t_drive": 4, "tau_E1": 100.0, "vamp_mode": 1, "v5_mode": 1, "rho_gamma": 1.0, "A_eCB": 1.0, "dpre_min": -0.29}'
DROPT="letzkus_3ap_200hz_dt-10ms@distal|control,letzkus_3ap_200hz_dt-10ms@distal|nmdar_block"
case $M in 2*) NAME=s$M; NT=11; KEEP2=",sjostrom_50hz_dt+10ms|nmdar_block,sjostrom_0.1hz_dt-10ms|mglu_block,sjostrom07_step200ms_pair|control,sjostrom07_step200ms_pair|no_block" ;; *) NAME=s1$M; NT=7; KEEP2="" ;; esac
OUT=$RS/${NAME}_$S
if [ "$VAL" = val ]; then
  # validation rescore: all targets, no stage-1 drops, extras on, maxiter 0 from the fitted json
  rm -f ${OUT}_val.json
  ARGS=(--free-filters "$FREE" --set "$SET" --fit-gamma --drop-targets "$DROPT" --l23-dirs $L23D
    --extra "l23:paired_l23l5:$L23D:$L23_BASIS_DIR:$GEOM_L23:$K235"
    --extra "l23l23:paired_l23l23,paired_l23l23_egger:$Z:$W/basis_l23l23::$K2323"
    --seed-fits $OUT.json --seed-set '{}' --seed 5 --maxiter 0 --save ${OUT}_val)
else
  KEEP='10Hz_10ms|control,10Hz_-10ms|control,sjostrom_0.1hz_dt+10ms|control,sjostrom_0.1hz_dt-10ms|control,sjostrom_20hz_dt+10ms|control,sjostrom_40hz_dt+10ms|control,sjostrom_50hz_dt-10ms|control'"$KEEP2"
  DROP=$(awk -F, -v keep="$KEEP" 'BEGIN{n=split(keep,k,","); for(i=1;i<=n;i++) K[k[i]]=1} NR>1{ if(NF!=6){print "bad csv row " NR > "/dev/stderr"; exit 3} t=$1"|"$2; if(!(t in K)) printf "%s%s", (c++?",":""), t}' $RS/r1D_s.csv)
  echo "kept $NT core targets; dropped $(echo $DROP | tr ',' '\n' | wc -l) other L5 targets"
  ARGS=(--free-filters "$FREE" --set "$SET" --fit-gamma --drop-targets "$DROPT,$DROP" --maxiter 300 --save $OUT)
  if [ "$S" = s ]; then ARGS+=(--seed 5 --seed-fits $SEEDJ --seed-set '{}'); else ARGS+=(--seed 6); fi
fi
echo "=== $NAME $S ${VAL:-fit} FREE '${FREE}'(+gamma_d,gamma_p) $(date)"; echo "SET $SET"; nvidia-smi -L
rc=0
python -u $RS/gpu_v10_rho.py --dirs $DIRS --pairs $PAIRS --groups paired_l5,sjostrom07,paired_l5_extra --filters "$FILTERS" \
  "${ARGS[@]}" || rc=$?
if [ $rc -ne 0 ]; then
  if [ "$VAL" = val ] && [ -s ${OUT}_val.json ]; then echo "val: exit $rc accepted (REPRO DIFF vs the 7-target json, outputs written)"; else exit $rc; fi
fi
echo "=== done $(date)"
