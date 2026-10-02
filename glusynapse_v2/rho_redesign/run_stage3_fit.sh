#!/bin/bash
# Stage 3 (STAGE3_DESIGN.md sections 2-3): 20 core targets on L5 (10), L2/3->L5 (4), L2/3->L2/3 (6); 3N adds sjostrom07_step200ms_pair|no_block (21).
# MEASURED 22291778-83 (3C s/u, 3N s, 20 cores): fits 8.6-11.7 min, 19.5-26.9 GB; vals 7.7-11.3 min, 49.4-70.4 GB (peak 70.4), so val now 88G 0:30:00 (84G was under +25%).
# MEASURED 22291785/86 (WL23 fits, 30 targets): 17.5-18.4 min, MaxRSS 66.0 GB (at the 63G limit), so fits now 83G 0:30:00.
# MEASURED all-pathway fit 22289873/74 (S1C rule, 30 targets): 16.5-17 min, 50.1-50.4 GB, so fits 63G 0:30:00.
# Based on run_stage_fit.sh (inputs, FILTERS, val branch) and run_s1c_l23fit.sh (--extra specs). Usage (from plastyfire/, via sjob.sh -g):
#   bash glusynapse_v2/rho_redesign/run_stage3_fit.sh MODEL SEEDMODE [val]
#   MODEL 3C : S1C rule (G8X shaft licence + eCB), free theta_eCB (+ a00 a01 a10 a11 gamma_d gamma_p); seeds s1C_s.json,s2E_s.json
#         3N : 3C + NO mode 3 (tau_NO 6.7), free A_NO, d_NO_max; seeds s2C_s.json,s1C_s.json (A_NO/d_NO_max not in SET)
#         4N : 3N rule on the 11 stage-2 L5 cores + all L2/3->L5 + L2/3->L2/3 (no Egger); seeds s2C_s.json, s1c_l23fit/s1CL_s.json
#         5N : same on all 40 L5 + all L2/3 targets, SEM-floor weights (w=(sem/0.05)^2)
#   SEEDMODE s (seed 5, --seed-fits) | u (seed 6, unseeded)
#   val : rescore s3<M>_<S>.json on ALL targets, unweighted, --maxiter 0 (FREE keys stripped from SET).
#   env KERNEL (default RS/gpu_v10_rho.py): fit kernel, so v11 models can be added as new MODEL cases.
# Weights (DE objective only; csv chi2 unweighted): Markram +10 x2; SEM floor 0.05 for 0.1 Hz +10 (0.64), Letzkus 1AP +10 and 3AP -10 prox (0.36);
#   20 Hz -10 mglu_block 0.5; zilberter_1ap_dt+10ms and dt-10ms 0.5 (known uniformity conflict).
# Log must say "over 20 targets" (21 for 3N), l23 4 targets, l23l23 6 targets. Grep: "pairs|chi2 L5|chi2 l23|over .* targets|chi2 total|REPRO|Traceback|weights|seed member".
# Sizing: MEASURED R1 3-pathway fits 22235493-502 max 72.9 GB 27:23 -> fit 91G 0:45:00; val 84G 0:15:00 (stage-2 rescores 51-67 GB, 8 min).
set -euo pipefail
M=$1; S=$2; VAL=${3:-}
R=/project/rrg-emuller/dhuruva/plastyfire; V2=glusynapse_v2; RS=$V2/rho_redesign; W=/scratch/dhuruva/s2g0321; X=$W/extracted; T=delta-split2-ljp25g0321-prefire
KERNEL=${KERNEL:-$RS/gpu_v10_rho.py}
source glusynapse_v2/env_v3.sh
export CEXP_DIR=$W/cexp ANALYTICAL_BASIS_DIR=$W/basis_l5l5 L23_BASIS_DIR=$W/basis_l23l5 PYTHONPATH=$R
DIRS=$X/ebner_$T-vca,$X/markram_$T-vca,$X/sj03_$T-vca,$X/sj03r50_$T-vca,$X/sj07_$T-vca,$X/l5extra_$T-vca
L23D=$X/ebner_l23l5_$T-vseg-rs
Z=$X/zilberter_l23l23_$T-vseg-rs,$X/l23l23extra_$T-vseg-rs
GEOM_L23=${GEOM_L23:-$R/ebner/pair_geometry_L23PC_L5TTPC.csv}
K5=$W/doublets/l5_keep_pairs.txt; K235=$W/doublets/l23l5_keep_pairs.txt; K2323=$W/doublets/l23l23_keep_pairs.txt
for f in $K5 $K235 $K2323; do [ -s $f ] || { echo "missing keep list $f"; exit 2; }; done
PAIRS=$(comm -12 <(tr ',' '\n' < $V2/subset24_pairs.txt | sed '/^$/d' | sort -u) <(tr ',' '\n' < $K5 | sed '/^$/d' | sort -u) | paste -sd,)
[ -n "$PAIRS" ] || { echo "empty L5 pair set"; exit 2; }
echo "pairs: l5 used $(echo $PAIRS | tr ',' '\n' | wc -l); l23l5 keep $(tr ',' '\n' < $K235 | sed '/^$/d' | wc -l); l23l23 keep $(tr ',' '\n' < $K2323 | sed '/^$/d' | wc -l); kernel $KERNEL"
G8X='"v5_mode": 2, "veto_T": 25.0, "ecb_ref": 2, "gate_src": 2, "gate_win": 100.0, "gate_theta": 0.15, "theta_V": 3.5751446275588403, "t_exact": 1, "v10_counts": 1'
case $M in
  3C) SET="{$G8X}"; FREE=theta_eCB; SEEDJ=$RS/s1C_s.json,$RS/s2E_s.json ;;
  3N) SET="{$G8X, \"no_mode\": 3, \"tau_NO\": 6.7}"; FREE=theta_eCB,A_NO,d_NO_max; SEEDJ=$RS/s2C_s.json,$RS/s1C_s.json ;;
  3V|3W) KERNEL=$RS/gpu_v11_rho.py; G8V=${G8X/\"v5_mode\": 2/\"v5_mode\": 1}; PK=1; [ $M = 3W ] && PK=0; SET="{$G8V, \"veto_t0\": 2, \"veto_Tv\": 17, \"veto_peak\": $PK, \"veto_peak_k\": 0.5}"; FREE=theta_eCB; SEEDJ=$RS/s1C_s.json,$RS/s2E_s.json ;;
  4V) KERNEL=$RS/gpu_v11_rho.py; G8V=${G8X/\"v5_mode\": 2/\"v5_mode\": 1}; SET="{$G8V, \"no_mode\": 3, \"tau_NO\": 6.7, \"veto_t0\": 2, \"veto_Tv\": 17, \"veto_peak\": 1, \"veto_peak_k\": 0.5}"; FREE=theta_eCB,A_NO,d_NO_max; SEEDJ=$RS/s2C_s.json,/scratch/dhuruva/s1c_l23fit/s1CL_s.json ;;
  4N|5N) SET="{$G8X, \"no_mode\": 3, \"tau_NO\": 6.7}"; FREE=theta_eCB,A_NO,d_NO_max; SEEDJ=$RS/s2C_s.json,/scratch/dhuruva/s1c_l23fit/s1CL_s.json ;;
  *) echo "MODEL 3C|3N|4N|5N|3V|3W|4V (stage 3)" >&2; exit 1 ;;
esac
export BOX=
FILTERS='{"pre_drive": 1, "i_scale": 1e-5, "t_drive": 4, "tau_E1": 100.0, "vamp_mode": 1, "v5_mode": 1, "rho_gamma": 1.0, "A_eCB": 1.0, "dpre_min": -0.29}'
DROPT="letzkus_3ap_200hz_dt-10ms@distal|control,letzkus_3ap_200hz_dt-10ms@distal|nmdar_block"
NAME=s$M; OUT=$RS/${NAME}_$S
EXTRA=(--l23-dirs $L23D
  --extra "l23:paired_l23l5:$L23D:$L23_BASIS_DIR:$GEOM_L23:$K235"
  --extra "l23l23:paired_l23l23,paired_l23l23_egger:$Z:$W/basis_l23l23::$K2323")
if [ "$VAL" = val ]; then
  rm -f ${OUT}_val.json
  # free params must come from the fitted json, not SET (SET values override the seed in _pack10)
  for k in ${FREE//,/ }; do SET=$(echo "$SET" | sed -E "s/, *\"$k\": *[-0-9.e+]+//; s/\"$k\": *[-0-9.e+]+, *//"); done
  ARGS=(--free-filters "$FREE" --set "$SET" --fit-gamma --drop-targets "$DROPT" "${EXTRA[@]}"
    --seed-fits $OUT.json --seed-set '{}' --seed 5 --maxiter 0 --save ${OUT}_val)
else
  KEEP5='10Hz_10ms|control,10Hz_-10ms|control,sjostrom_0.1hz_dt+10ms|control,sjostrom_0.1hz_dt-10ms|control,sjostrom_20hz_dt+10ms|control,sjostrom_40hz_dt+10ms|control,sjostrom_50hz_dt-10ms|control,sjostrom_20hz_dt-10ms|control,sjostrom_20hz_dt-10ms|mglu_block,sjostrom07_step200ms_pair|control'
  NT=20; [ "$M" = 3N ] && { KEEP5="$KEEP5,sjostrom07_step200ms_pair|no_block"; NT=21; }
  KEEP23='letzkus_1ap_dt+10ms|control,letzkus_3ap_200hz_dt+10ms@proximal|control,letzkus_3ap_200hz_dt+10ms@distal|control,letzkus_3ap_200hz_dt-10ms@proximal|control'
  KEEP2323='zilberter_1ap_dt+10ms|control,zilberter_1ap_dt-10ms|control,zilberter_train10_50hz_dt+4ms_last|control,zilberter_train10_50hz_dt-10ms_last|control,zilberter_train10_50hz_dt-10ms_last|mglu_block,zilberter_5ap_10hz_dt+10ms|control'
  dropof() { awk -F, -v keep="$1" -v pre="$3" 'BEGIN{n=split(keep,k,","); for(i=1;i<=n;i++) K[k[i]]=1} NR>1{ if(NF!=6){print "bad csv row " NR > "/dev/stderr"; exit 3} t=$1"|"$2; if(!(t in K)) printf "%s%s%s", (c++?",":""), pre, t}' $2; }
  DROP=$(dropof "$KEEP5" $RS/r1D_s.csv "l5/")
  DROP23=$(dropof "$KEEP23" $RS/s1C_s_val_l23.csv "l23/")
  DROP2323=$(dropof "$KEEP2323" $RS/s1C_s_val_l23l23.csv "l23l23/")
  WEIGHTS='{"l5/10Hz_10ms|control": 2.0, "l5/sjostrom_0.1hz_dt+10ms|control": 0.64, "l5/sjostrom_20hz_dt-10ms|mglu_block": 0.5, "l23/letzkus_1ap_dt+10ms|control": 0.36, "l23/letzkus_3ap_200hz_dt-10ms@proximal|control": 0.36, "l23l23/zilberter_1ap_dt+10ms|control": 0.5, "l23l23/zilberter_1ap_dt-10ms|control": 0.5}'
  if [ "$M" = 4N ] || [ "$M" = 4V ] || [ "$M" = 5N ]; then
    # 4N: 11 stage-2 L5 cores (CORE_TARGETS 1-11) + all L2/3->L5 (9 in csv) + all L2/3->L2/3 (15 in csv); 5N: all 40 L5 + same L2/3. Egger dropped (validation only).
    K11='10Hz_10ms|control,10Hz_-10ms|control,sjostrom_0.1hz_dt+10ms|control,sjostrom_0.1hz_dt-10ms|control,sjostrom_20hz_dt+10ms|control,sjostrom_40hz_dt+10ms|control,sjostrom_50hz_dt-10ms|control,sjostrom_50hz_dt+10ms|nmdar_block,sjostrom_0.1hz_dt-10ms|mglu_block,sjostrom07_step200ms_pair|control,sjostrom07_step200ms_pair|no_block'
    DROP=""; NT=40; [ "$M" != 5N ] && { DROP=$(dropof "$K11" $RS/r1D_s.csv "l5/"); NT=11; }
    DROP23=""; DROP2323="l23l23/egger1999_5ap_20hz_dt+10ms|control"
    WEIGHTS='{"l5/10Hz_10ms|control": 2.0, "l5/sjostrom_0.1hz_dt+10ms|control": 0.64, "l23/letzkus_1ap_dt+10ms|control": 0.36, "l23/letzkus_3ap_200hz_dt-10ms@proximal|control": 0.36, "l23l23/zilberter_1ap_dt+10ms|control": 0.5, "l23l23/zilberter_1ap_dt-10ms|control": 0.5}'
    if [ "$M" = 5N ]; then   # SEM floor 0.05 via w=(sem/0.05)^2 for tight-SEM targets (Markram exempt); 20 Hz -10 mglu_block 0.5 as in 3N
      WEIGHTS='{"l5/10Hz_10ms|control": 2.0, "l5/sjostrom_0.1hz_dt+10ms|control": 0.64, "l5/sjostrom_0.1hz_dt-10ms|mglu_block": 0.64, "l5/sjostrom_0.1hz_dt-10ms|nmdar_block": 0.64, "l5/2Hz_5ms|control": 0.64, "l5/sjostrom_40hz_dt0ms|control": 0.591, "l5/sjostrom_20hz_dt-10ms|mglu_block": 0.5, "l23/letzkus_1ap_dt+10ms|control": 0.36, "l23/letzkus_3ap_200hz_dt-10ms@proximal|control": 0.36, "l23l23/zilberter_1ap_dt+10ms|control": 0.5, "l23l23/zilberter_1ap_dt-10ms|control": 0.5, "l23l23/zilberter_train10_50hz_post_only|control": 0.64}'
    fi
  fi
  DROPS=$(printf '%s\n' "$DROP" "$DROP23" "$DROP2323" | sed '/^$/d' | paste -sd,)
  echo "kept $NT core targets; dropped l5 $(echo $DROP | tr ',' '\n' | wc -l), l23 $(echo $DROP23 | tr ',' '\n' | wc -l), l23l23 $(echo $DROP2323 | tr ',' '\n' | wc -l)"
  ARGS=(--free-filters "$FREE" --set "$SET" --fit-gamma --drop-targets "$DROPT,$DROPS" "${EXTRA[@]}" --weights "$WEIGHTS" --maxiter 300 --save $OUT)
  if [ "$S" = s ]; then ARGS+=(--seed 5 --seed-fits $SEEDJ --seed-set '{}'); else ARGS+=(--seed 6); fi
fi
echo "=== $NAME $S ${VAL:-fit} FREE '${FREE}'(+gamma_d,gamma_p) $(date)"; echo "SET $SET"; nvidia-smi -L
rc=0
python -u $KERNEL --dirs $DIRS --pairs $PAIRS --groups paired_l5,sjostrom07,paired_l5_extra --filters "$FILTERS" "${ARGS[@]}" || rc=$?
if [ $rc -ne 0 ]; then
  if [ "$VAL" = val ] && [ -s ${OUT}_val.json ]; then echo "val: exit $rc accepted (REPRO DIFF vs the fitted json, outputs written)"; else exit $rc; fi
fi
echo "=== done $(date)"
