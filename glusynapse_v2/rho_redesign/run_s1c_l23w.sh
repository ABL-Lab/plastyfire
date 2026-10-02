#!/bin/bash
# L5 vs L2/3 trade-off map (STAGE3_DESIGN step 3c): copy of run_s1c_l23fit.sh (same S1C rule, 7 L5 shape + 8 L2/3->L5 + 15
# L2/3->L2/3 targets) with all L2/3 target z^2 weights multiplied by WL23 (--weights "l23/<target>|<cond>": WL23 and
# "l23l23/...": WL23; the reported chi2 stays unweighted). Seeded from SEEDJ only (default s1C_s.json), seed 5.
# Usage (from plastyfire/, via sjob.sh -g): bash glusynapse_v2/rho_redesign/run_s1c_l23w.sh TAG WL23 [SEEDJ]
#   e.g. TAG=w025 WL23=0.25 ; TAG=w4 WL23=4 ; TAG=a12 WL23=1 SEEDJ=/scratch/dhuruva/s1c_l23w/s1C_s_a12.json (a00 1.2, a01 2.0)
# a bounds: --a-lo 0 --a-hi 5 (fit_v6 default), so a00 1.2 / a01 2.0 are inside.
# Output: /scratch/dhuruva/s1c_l23w/s1CW_<TAG>.{json,csv,_l23.csv,_l23l23.csv}. Grep: "weights|chi2 L5|chi2 l23|chi2 total|Traceback".
# Sizing: MEASURED r1E_s 22235501 (same rule, 65 targets) 22:29, 72.9 GB -> 91G 0:45.
set -euo pipefail
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
[ -n "$PAIRS" ] || { echo "empty L5 pair set"; exit 2; }
G8X='"v5_mode": 2, "veto_T": 25.0, "ecb_ref": 2, "gate_src": 2, "gate_win": 100.0, "gate_theta": 0.15, "theta_V": 3.5751446275588403, "t_exact": 1, "v10_counts": 1'
SET="{$G8X}"; FREE=theta_eCB; export BOX=
FILTERS='{"pre_drive": 1, "i_scale": 1e-5, "t_drive": 4, "tau_E1": 100.0, "vamp_mode": 1, "v5_mode": 1, "rho_gamma": 1.0, "A_eCB": 1.0, "dpre_min": -0.29}'
DROPT="letzkus_3ap_200hz_dt-10ms@distal|control,letzkus_3ap_200hz_dt-10ms@distal|nmdar_block,egger1999_5ap_20hz_dt+10ms|control"
KEEP='10Hz_10ms|control,10Hz_-10ms|control,sjostrom_0.1hz_dt+10ms|control,sjostrom_0.1hz_dt-10ms|control,sjostrom_20hz_dt+10ms|control,sjostrom_40hz_dt+10ms|control,sjostrom_50hz_dt-10ms|control'
DROP=$(awk -F, -v keep="$KEEP" 'BEGIN{n=split(keep,k,","); for(i=1;i<=n;i++) K[k[i]]=1} NR>1{ if(NF!=6){print "bad csv row " NR > "/dev/stderr"; exit 3} t=$1"|"$2; if(!(t in K)) printf "%s%s", (c++?",":""), t}' $RS/r1D_s.csv)
TAG=$1; WL23=$2; SEEDJ=${3:-$RS/s1C_s.json}
WJ=$(awk -F, -v w="$WL23" 'NR>1{printf "%s\"l23/%s|%s\": %s", (c++?", ":""), $1, $2, w}' $RS/s1C_s_val_l23.csv)
WJ2=$(awk -F, -v w="$WL23" 'NR>1{printf ", \"l23l23/%s|%s\": %s", $1, $2, w}' $RS/s1C_s_val_l23l23.csv)
WEIGHTS_JSON="{$WJ$WJ2}"
OD=/scratch/dhuruva/s1c_l23w; mkdir -p $OD; OUT=$OD/s1CW_$TAG
ARGS=(--free-filters "$FREE" --set "$SET" --fit-gamma --drop-targets "$DROPT,$DROP" --l23-dirs $L23D
  --extra "l23:paired_l23l5:$L23D:$L23_BASIS_DIR:$GEOM_L23:$K235"
  --extra "l23l23:paired_l23l23,paired_l23l23_egger:$Z:$W/basis_l23l23::$K2323"
  --weights "$WEIGHTS_JSON" --maxiter 300 --save $OUT --seed 5 --seed-fits $SEEDJ --seed-set '{}')
echo "=== s1CW $TAG WL23=$WL23 seed=$SEEDJ FREE '${FREE}'(+gamma_d,gamma_p) $(date)"; echo "SET $SET"; echo "WEIGHTS $WEIGHTS_JSON"; nvidia-smi -L
python -u $RS/gpu_v10_rho.py --dirs $DIRS --pairs $PAIRS --groups paired_l5,sjostrom07,paired_l5_extra --filters "$FILTERS" "${ARGS[@]}"
echo "=== done $(date)"
