#!/bin/bash
# Stage 4 driver: the stage-3 fits (run_stage3_fit.sh, same inputs, targets, weights and MODEL cases) with the fit-procedure
# fixes the user approved on 2026-10-02 (FIT_META_DIAG.md section 6; code spec for the fitter side in SPEC_FITPROC.md).
# Driver-only here (no kernel edit):
#   fix 2 multi-start : one job per START (s<N> seeded, u<N> unseeded LHS, N = DE seed), popsize POP 15, maxiter 600;
#                       "pick" mode selects the best start by the fitted objective and checks basin agreement (df <= 2).
#                       Admissible-only initial populations need the fitter flag --init-admissible (SPEC 2; env INITADM).
#   fix 3 anchors     : FIXANC=1 default (theta_eCB 0.535, A_NO 190, d_NO_max 0.57 in SET; run_stage3_fit.sh block).
#   fix 5 smoothing   : phase A with rho_sigma SIGMA (default 0.1, the existing ndtr readout of gpu_v4_rho._ratios, set
#                       through --set), then phase B: --resume from phase A's final population at rho_sigma 0 for MAXB
#                       generations, so the reported fit and all csvs are the binary readout. SIGMA=0: one phase.
#   fix 6 pair error  : "lopo" mode = leave-one-L5-pair-out / leave-one-group-out (L2/3 pathways, LOPO_G groups)
#                       rescores at --maxiter 0; "lopoagg" turns them into SE_pair per target (RS/<TAG>_sepair.csv);
#                       env SEPAIR=<that csv> multiplies each target's DE weight by SEM^2 / (SEM^2 + SE_pair^2), which is
#                       exactly chi2 with SEM_eff^2 = SEM^2 + SE_pair^2 in the objective (csv z stays on the data SEM).
#   fix 4 hinge       : needs the fitter flag --hinge (SPEC 4). Env HINGE (json) is passed through only when set.
#   user 2026-10-02   : DEMOTE=1 default (data-conflict rows to validation); depression-band seeds: every seed json with
#                       a00 < A00 (default 1.2) also enters as a copy with a00 = A00 and a10 shifted by the same amount,
#                       so theta_p - theta_d is unchanged at every synapse (--seed-set cannot do this: fit_v6 merges it
#                       UNDER the json's pre, and a00 lives in the json's "a" block, so the copy is made with jq).
# Usage (from plastyfire/, via glusynapse_v2/sjob.sh -g; pick / lopoagg / names are plain bash, no GPU):
#   bash glusynapse_v2/rho_redesign/run_stage4_fit.sh MODEL START [fit|val|lopo [IDX]|lopoagg|pick|name]
#   MODEL  3C 3N 4N 5N 3V 3W 4V, as run_stage3_fit.sh.   START  s5 s6 (seeded) u7 u8 u9 (unseeded); plain s / u = s5 / u6.
#   Recommended census: s5 s6 u7 u8 (2 seeded + 2 unseeded), then: run_stage4_fit.sh MODEL - pick
#   val     rescore <NAME>_<START>.json on all targets, unweighted, rho_sigma 0, --maxiter 0 (as stage 3).
#   lopo    IDX (or SLURM_ARRAY_TASK_ID): 0..NL5-1 drop L5 pair IDX (L5-only rescore, no extras); NL5..NL5+G-1 drop
#           L2/3->L5 group; then L2/3->L2/3 groups (G = LOPO_G, default 10; pairs assigned round-robin in sorted order).
#           Rescored json: env LOPOJ (default <NAME>_<START>.json); outputs /scratch/dhuruva/stage4/lopo/<TAG>/u<IDX>*.csv,
#           TAG = basename of LOPOJ. To rescore an old stage-3 json: FIXANC=0 DEMOTE=0 LOPOJ=$RS/s4N_u.json ... 4N u lopo.
#   lopoagg SE_pair = sqrt((n-1)/n sum (p_i - mean p)^2) per target -> RS/<TAG>_sepair.csv (key,pathway,target,condition,
#           n,sem,se_pair,ratio,factor). Pure awk.
# Env: KERNEL (default RS/gpu_v10_rho.py; 3V/3W/4V force gpu_v11), KERNELX (overrides the kernel after the MODEL case:
#      gpu_v12 -> RS/gpu_v12_rho.py, or a path; S6MECH_DESIGN.md arms run on 3V/4V with KERNELX=gpu_v12 + SETX), DEMOTE (1), FIXANC (1), MW (Markram +10 weight, 2.0),
#      POP (15), MAXITER (600, total generations: phase A = MAXITER - MAXB), MAXB (150), SIGMA (0.1), A00 (1.2; 0 = no
#      depression-band seeds), SEEDX (seed jsons), SEPAIR, HINGE, INITADM, STRATEGY, FITTER (the last four need FITTER=fit_v7 or fit_v8,
#      run via fit_launch.py; left unset they add no flag and the kernel uses fit_v6).
#      fit_v8 also: PHASEB (rescore default | de = v7 phase B), PHASEBK (30), POLISHK (3), POLISHN (200 evals per polished point); phase B = rescore + polish (fit_v8.py doc).
# Name: f<MODEL>[m<MW>][d][a][g<SIGMA>][e][h]_<START>; e.g. f4Nda_s5. Phase A outputs in /scratch/dhuruva/stage4.
# MEASURED 22307548 smoke 4N s5 (POP 15, 5+2 gen): 18:14, 57.7 GB, generation ~1.3 s (setup dominates); 600 gen adds ~13 min, so a fit is ~35 min -> 73G 1:00:00. LOPO: L5 unit 22307549 2:31 29.4 GB -> 37G 0:15; L2/3->L5 unit 22307550 5:00 38.2 GB -> 48G 0:15.
# Sizing (MEASURED, stage 3 at popsize 8 / maxiter 300): fits 8.6-27 min, 19.5-72.9 GB; val 70.4 GB 11:15 -> 88G 0:30.
#   Popsize 15 / 600 generations / two builds is about 4x the DE work: pilot the smoke test (POP 15, MAXITER 5, MAXB 2)
#   before sizing the census. LOPO L5-only rescores: MEASURED 22292757 3 rescores 7:31 37.6 GB -> 47G 0:15 each;
#   L2/3 group rescores: as val, 88G 0:30.
set -euo pipefail
M=$1; S=$2; MODE=${3:-fit}
R=/project/rrg-emuller/dhuruva/plastyfire; V2=glusynapse_v2; RS=$V2/rho_redesign; W=/scratch/dhuruva/s2g0321; X=$W/extracted; T=delta-split2-ljp25g0321-prefire
SCR=/scratch/dhuruva/stage4
KERNEL=${KERNEL:-$RS/gpu_v10_rho.py}
DEMOTE=${DEMOTE:-1}; FIXANC=${FIXANC:-1}; POP=${POP:-15}; MAXITER=${MAXITER:-600}; MAXB=${MAXB:-150}; SIGMA=${SIGMA:-0.1}; A00=${A00:-1.2}
case $S in s) S=s5;; u) S=u6;; esac
G8X='"v5_mode": 2, "veto_T": 25.0, "ecb_ref": 2, "gate_src": 2, "gate_win": 100.0, "gate_theta": 0.15, "theta_V": 3.5751446275588403, "t_exact": 1, "v10_counts": 1'
case $M in
  3C) SET="{$G8X}"; FREE=theta_eCB; SEEDJ=$RS/s1C_s.json,$RS/s2E_s.json ;;
  3N) SET="{$G8X, \"no_mode\": 3, \"tau_NO\": 6.7}"; FREE=theta_eCB,A_NO,d_NO_max; SEEDJ=$RS/s2C_s.json,$RS/s1C_s.json ;;
  3V|3W) KERNEL=$RS/gpu_v11_rho.py; G8V=${G8X/\"v5_mode\": 2/\"v5_mode\": 1}; PK=1; [ $M = 3W ] && PK=0; SET="{$G8V, \"veto_t0\": 2, \"veto_Tv\": 17, \"veto_peak\": $PK, \"veto_peak_k\": 0.5}"; FREE=theta_eCB; SEEDJ=$RS/s1C_s.json,$RS/s2E_s.json ;;
  4V) KERNEL=$RS/gpu_v11_rho.py; G8V=${G8X/\"v5_mode\": 2/\"v5_mode\": 1}; SET="{$G8V, \"no_mode\": 3, \"tau_NO\": 6.7, \"veto_t0\": 2, \"veto_Tv\": 17, \"veto_peak\": 1, \"veto_peak_k\": 0.5}"; FREE=theta_eCB,A_NO,d_NO_max; SEEDJ=$RS/s2C_s.json,/scratch/dhuruva/s1c_l23fit/s1CL_s.json ;;
  4N|5N) SET="{$G8X, \"no_mode\": 3, \"tau_NO\": 6.7}"; FREE=theta_eCB,A_NO,d_NO_max; SEEDJ=$RS/s2C_s.json,/scratch/dhuruva/s1c_l23fit/s1CL_s.json ;;
  *) echo "MODEL 3C|3N|4N|5N|3V|3W|4V" >&2; exit 1 ;;
esac
# KERNELX: kernel override after the MODEL case (s6mech 2026-10-03; 3V/3W/4V force gpu_v11 above). gpu_v12 or gpu_v12_rho.py ->
# $RS/gpu_v12_rho.py; a value with a / is a path (relative to the repo root, or absolute).
if [ -n "${KERNELX:-}" ]; then case $KERNELX in */*) KERNEL=$KERNELX ;; *) KERNEL=$RS/${KERNELX%_rho.py}_rho.py ;; esac; fi
if [ "$FIXANC" = 1 ]; then   # anchors and their sources: run_stage3_fit.sh FIXANC block (same values, same overrides)
  ANC_ECB=${ANC_ECB:-0.535}; ANC_ANO=${ANC_ANO:-190.0}; ANC_DNO=${ANC_DNO:-0.57}
  ADD="\"theta_eCB\": $ANC_ECB"
  case ",$FREE," in *",A_NO,"*) ADD="$ADD, \"A_NO\": $ANC_ANO";; esac
  case ",$FREE," in *",d_NO_max,"*) ADD="$ADD, \"d_NO_max\": $ANC_DNO";; esac
  SET="${SET%\}}, $ADD}"; FREE=$(echo "$FREE" | tr ',' '\n' | { grep -vxE 'theta_eCB|A_NO|d_NO_max' || true; } | paste -sd,)
fi
[ -n "${SEEDX:-}" ] && SEEDJ=$SEEDX
NAME=f$M; [ -n "${MW:-}" ] && NAME=${NAME}m${MW/./}
[ "$DEMOTE" = 1 ] && NAME=${NAME}d
[ "$FIXANC" = 1 ] && NAME=${NAME}a
[ "$SIGMA" != 0.1 ] && NAME=${NAME}g${SIGMA/./}
[ -n "${SEPAIR:-}" ] && NAME=${NAME}e
[ -n "${HINGE:-}" ] && NAME=${NAME}h
[ -n "${NAMEX:-}" ] && NAME=${NAME}${NAMEX}   # env NAMEX: free tag for smoke tests / sensitivity runs (kept out of the census pick)
OUT=$RS/${NAME}_$S

# ---- pure-bash modes (no inputs, no GPU) --------------------------------------------------------------------------
if [ "$MODE" = name ]; then echo "$OUT"; exit 0; fi
if [ "$MODE" = pick ]; then
  # best start by the fitted (weighted, binary-readout) objective de_fun; basin accepted if >= 2 starts within df <= 2
  for j in $RS/${NAME}_[su][0-9]*.json; do
    [ -s "$j" ] || continue; case $j in *_val.json) continue;; esac
    printf "%s %s %s %s\n" "$(jq -r '.de_fun' $j)" "$(jq -r '.chi2_total' $j)" "$(jq -r '.nfev' $j)" "$(basename $j .json)"
  done | sort -g | awk 'NR==1{b=$1; bn=$4} {n++; a+=($1-b<=2); printf "%-24s de_fun %9.3f  chi2_total %9.3f  gens %s  df %6.2f\n", $4, $1, $2, $3, $1-b}
      END{ if(!n){print "no fitted starts"; exit} printf "best %s; %d of %d starts within df 2 -> basin %s\n", bn, a, n, (a>=2?"ACCEPTED":"NOT CONFIRMED (run more starts)") }' \
    > $RS/${NAME}_pick.txt
  cat $RS/${NAME}_pick.txt; exit 0
fi
LOPOJ=${LOPOJ:-$OUT.json}; TAG=$(basename $LOPOJ .json); LD=$SCR/lopo/$TAG; LOPO_G=${LOPO_G:-10}
if [ "$MODE" = lopoagg ]; then
  # jackknife SE per target over the leave-one-out units of each pathway (u<idx>.csv L5; u<idx>_l23.csv; u<idx>_l23l23.csv)
  n5=$(cat $LD/n_l5 2>/dev/null || echo 0)
  [ $n5 -gt 0 ] || { echo "no $LD/n_l5: run the lopo array first"; exit 2; }
  for f in $LD/u*.csv; do b=$(basename $f .csv); i=${b#u}; i=${i%%_*}; pw=
    case $b in
      *_l23l23) if [ $i -ge $((n5 + LOPO_G)) ]; then pw=l23l23; fi ;;
      *_l23) if [ $i -ge $n5 ] && [ $i -lt $((n5 + LOPO_G)) ]; then pw=l23; fi ;;
      *) if [ $i -lt $n5 ]; then pw=l5; fi ;;
    esac
    if [ -n "$pw" ]; then awk -F, -v pw=$pw 'NR>1{print pw","$1","$2","$4","$5}' $f; fi
  done |
    awk -F, 'BEGIN{print "key,pathway,target,condition,n,sem,se_pair,ratio,factor"}
      { k=$1"/"$2"|"$3; n[k]++; s[k]+=$5; ss[k]+=$5*$5; sem[k]=$4; pw[k]=$1; t[k]=$2; c[k]=$3 }
      END{ for(k in n){ m=s[k]/n[k]; v=ss[k]/n[k]-m*m; if(v<0)v=0; se=sqrt((n[k]-1)*v)   # (n-1)/n sum (p-m)^2 = (n-1) var_pop
           printf "%s,%s,%s,%s,%d,%.6g,%.6g,%.4f,%.6f\n", k, pw[k], t[k], c[k], n[k], sem[k], se, se/sem[k], sem[k]^2/(sem[k]^2+se^2) } }' \
    > $RS/${TAG}_sepair.csv
  echo "$(($(wc -l < $RS/${TAG}_sepair.csv) - 1)) targets -> $RS/${TAG}_sepair.csv; largest SE_pair/SEM:"; tail -n +2 $RS/${TAG}_sepair.csv | sort -t, -k8 -gr | head -5
  exit 0
fi

# ---- GPU modes ---------------------------------------------------------------------------------------------------
[ "${DRY:-0}" = 1 ] || source glusynapse_v2/env_v3.sh   # DRY=1: print the kernel commands, run nothing (login-node safe)
export CEXP_DIR=$W/cexp ANALYTICAL_BASIS_DIR=$W/basis_l5l5 L23_BASIS_DIR=$W/basis_l23l5 PYTHONPATH=$R
[ -n "${FITTER:-}" ] && export FITTER   # SPEC_FITPROC.md part 2: RS/fit_launch.py aliases fit_v6 -> $FITTER (e.g. fit_v7), so the kernel file is not edited
DIRS=$X/ebner_$T-vca,$X/markram_$T-vca,$X/sj03_$T-vca,$X/sj03r50_$T-vca,$X/sj07_$T-vca,$X/l5extra_$T-vca
L23D=$X/ebner_l23l5_$T-vseg-rs
Z=$X/zilberter_l23l23_$T-vseg-rs,$X/l23l23extra_$T-vseg-rs
GEOM_L23=${GEOM_L23:-$R/ebner/pair_geometry_L23PC_L5TTPC.csv}
K5=$W/doublets/l5_keep_pairs.txt; K235=$W/doublets/l23l5_keep_pairs.txt; K2323=$W/doublets/l23l23_keep_pairs.txt
for f in $K5 $K235 $K2323; do [ -s $f ] || { echo "missing keep list $f"; exit 2; }; done
PAIRS=$(comm -12 <(tr ',' '\n' < $V2/subset24_pairs.txt | sed '/^$/d' | sort -u) <(tr ',' '\n' < $K5 | sed '/^$/d' | sort -u) | paste -sd,)
[ -n "$PAIRS" ] || { echo "empty L5 pair set"; exit 2; }
echo "pairs: l5 used $(echo $PAIRS | tr ',' '\n' | wc -l); l23l5 keep $(tr ',' '\n' < $K235 | sed '/^$/d' | wc -l); l23l23 keep $(tr ',' '\n' < $K2323 | sed '/^$/d' | wc -l); kernel $KERNEL"
export BOX=
FILTERS='{"pre_drive": 1, "i_scale": 1e-5, "t_drive": 4, "tau_E1": 100.0, "vamp_mode": 1, "v5_mode": 1, "rho_gamma": 1.0, "A_eCB": 1.0, "dpre_min": -0.29}'
DROPT="letzkus_3ap_200hz_dt-10ms@distal|control,letzkus_3ap_200hz_dt-10ms@distal|nmdar_block"
X23="l23:paired_l23l5:$L23D:$L23_BASIS_DIR:$GEOM_L23"; X2323="l23l23:paired_l23l23,paired_l23l23_egger:$Z:$W/basis_l23l23:"
EXTRA=(--l23-dirs $L23D --extra "$X23:$K235" --extra "$X2323:$K2323")
SET0=$SET   # binary readout (rho_sigma 0 = the default): val, lopo, phase B
run() { rc=0; if [ "${DRY:-0}" = 1 ]; then printf "DRY:"; printf " %q" "$@"; echo; return; fi; echo "=== $NAME $S $MODE $* $(date)"; nvidia-smi -L; python -u ${FITTER:+$RS/fit_launch.py} $KERNEL --dirs $DIRS --groups paired_l5,sjostrom07,paired_l5_extra --filters "$FILTERS" "$@" || rc=$?; }
strip_free() { for k in ${FREE//,/ }; do SET0=$(echo "$SET0" | sed -E "s/, *\"$k\": *[-0-9.e+]+//; s/\"$k\": *[-0-9.e+]+, *//"); done; }

if [ "$MODE" = val ]; then
  rm -f ${OUT}_val.json; strip_free
  run --pairs $PAIRS --free-filters "$FREE" --set "$SET0" --fit-gamma --drop-targets "$DROPT" "${EXTRA[@]}" \
    --seed-fits $OUT.json --seed-set '{}' --seed 5 --maxiter 0 --save ${OUT}_val
  if [ $rc -ne 0 ]; then [ -s ${OUT}_val.json ] && echo "val: exit $rc accepted (REPRO DIFF vs the fitted json, outputs written)" || exit $rc; fi
  echo "=== done $(date)"; exit 0
fi

if [ "$MODE" = lopo ]; then
  IDX=${4:-${SLURM_ARRAY_TASK_ID:?lopo needs IDX or an array task}}
  [ -s $LOPOJ ] || { echo "missing $LOPOJ"; exit 2; }
  mkdir -p $LD/keep; strip_free
  L5L=($(echo $PAIRS | tr ',' '\n' | sort)); NL5=${#L5L[@]}; echo $NL5 > $LD/n_l5
  grp() { tr ',' '\n' < $1 | sed '/^$/d' | sort | awk -v g=$2 -v G=$LOPO_G '(NR-1)%G!=g' | paste -sd,; }
  if [ $IDX -lt $NL5 ]; then
    P=$(printf '%s\n' "${L5L[@]}" | grep -vx "${L5L[$IDX]}" | paste -sd,); XA=(); echo "lopo $IDX: L5 without ${L5L[$IDX]} ($((NL5 - 1)) pairs), L5 only"
  elif [ $IDX -lt $((NL5 + LOPO_G)) ]; then
    g=$((IDX - NL5)); grp $K235 $g > $LD/keep/l23l5_g$g.txt; P=$PAIRS
    XA=(--l23-dirs $L23D --extra "$X23:$LD/keep/l23l5_g$g.txt"); echo "lopo $IDX: L2/3->L5 without group $g of $LOPO_G"
  elif [ $IDX -lt $((NL5 + 2 * LOPO_G)) ]; then
    g=$((IDX - NL5 - LOPO_G)); grp $K2323 $g > $LD/keep/l23l23_g$g.txt; P=$PAIRS
    XA=(--extra "$X2323:$LD/keep/l23l23_g$g.txt"); echo "lopo $IDX: L2/3->L2/3 without group $g of $LOPO_G"
  else echo "lopo IDX $IDX out of range (0..$((NL5 + 2 * LOPO_G - 1)))"; exit 2; fi
  rm -f $LD/u$IDX.json
  run --pairs $P --free-filters "$FREE" --set "$SET0" --fit-gamma --drop-targets "$DROPT" "${XA[@]}" \
    --seed-fits $LOPOJ --seed-set '{}' --seed 5 --maxiter 0 --save $LD/u$IDX
  if [ $rc -ne 0 ]; then [ -s $LD/u$IDX.json ] && echo "lopo: exit $rc accepted (REPRO DIFF expected: fewer pairs)" || exit $rc; fi
  echo "=== done $(date)"; exit 0
fi

[ "$MODE" = fit ] || { echo "MODE fit|val|lopo|lopoagg|pick|name" >&2; exit 1; }
# ---- fit: targets and weights exactly as run_stage3_fit.sh ----------------------------------------------------------
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
  K11='10Hz_10ms|control,10Hz_-10ms|control,sjostrom_0.1hz_dt+10ms|control,sjostrom_0.1hz_dt-10ms|control,sjostrom_20hz_dt+10ms|control,sjostrom_40hz_dt+10ms|control,sjostrom_50hz_dt-10ms|control,sjostrom_50hz_dt+10ms|nmdar_block,sjostrom_0.1hz_dt-10ms|mglu_block,sjostrom07_step200ms_pair|control,sjostrom07_step200ms_pair|no_block'
  DROP=""; NT=40; [ "$M" != 5N ] && { DROP=$(dropof "$K11" $RS/r1D_s.csv "l5/"); NT=11; }
  DROP23=""; DROP2323="l23l23/egger1999_5ap_20hz_dt+10ms|control"
  WEIGHTS='{"l5/10Hz_10ms|control": 2.0, "l5/sjostrom_0.1hz_dt+10ms|control": 0.64, "l23/letzkus_1ap_dt+10ms|control": 0.36, "l23/letzkus_3ap_200hz_dt-10ms@proximal|control": 0.36, "l23l23/zilberter_1ap_dt+10ms|control": 0.5, "l23l23/zilberter_1ap_dt-10ms|control": 0.5}'
  if [ "$M" = 5N ]; then
    WEIGHTS='{"l5/10Hz_10ms|control": 2.0, "l5/sjostrom_0.1hz_dt+10ms|control": 0.64, "l5/sjostrom_0.1hz_dt-10ms|mglu_block": 0.64, "l5/sjostrom_0.1hz_dt-10ms|nmdar_block": 0.64, "l5/2Hz_5ms|control": 0.64, "l5/sjostrom_40hz_dt0ms|control": 0.591, "l5/sjostrom_20hz_dt-10ms|mglu_block": 0.5, "l23/letzkus_1ap_dt+10ms|control": 0.36, "l23/letzkus_3ap_200hz_dt-10ms@proximal|control": 0.36, "l23l23/zilberter_1ap_dt+10ms|control": 0.5, "l23l23/zilberter_1ap_dt-10ms|control": 0.5, "l23l23/zilberter_train10_50hz_post_only|control": 0.64}'
  fi
fi
if [ "$DEMOTE" = 1 ]; then   # data conflicts under a uniform synapse-local rule (HARD_TARGETS.md), user 2026-10-02
  DEMOTED='l23l23/zilberter_1ap_dt+10ms|control,l23l23/zilberter_1ap_dt-10ms|control,l23l23/zilberter_train10_50hz_dt-10ms_last|mglu_block,l5/sjostrom_10hz_dt-10ms|control'
  for k in ${DEMOTED//,/ }; do
    case ",$DROP,$DROP23,$DROP2323," in *",$k,"*) ;; *) case $k in l5/*) DROP=${DROP:+$DROP,}$k ;; l23l23/*) DROP2323=${DROP2323:+$DROP2323,}$k ;; esac ;; esac
    e=$(printf "%s" "$k" | sed "s/[][\\.*^$+?(){}|#]/\\\\&/g"); WEIGHTS=$(echo "$WEIGHTS" | sed -E "s#, *\"$e\": *[0-9.]+##; s#\"$e\": *[0-9.]+, *##")
  done
fi
DROPS=$(printf '%s\n' "$DROP" "$DROP23" "$DROP2323" | sed '/^$/d' | paste -sd,)
cnt() { echo "$1" | tr ',' '\n' | sed '/^$/d' | wc -l; }; N5=$(( $(wc -l < $RS/r1D_s.csv) - 1 - $(cnt "$DROP") )); N23=$(( $(wc -l < $RS/s1C_s_val_l23.csv) - 1 - $(cnt "$DROP23") )); N2323=$(( $(wc -l < $RS/s1C_s_val_l23l23.csv) - 1 - $(cnt "$DROP2323") ))
echo "kept core targets: l5 $N5, l23 $N23, l23l23 $N2323 (total $((N5 + N23 + N2323))); dropped l5 $(cnt "$DROP"), l23 $(cnt "$DROP23"), l23l23 $(cnt "$DROP2323")"
[ -n "${MW:-}" ] && WEIGHTS=$(echo "$WEIGHTS" | sed -E "s#\"l5/10Hz_10ms\|control\": *[0-9.]+#\"l5/10Hz_10ms|control\": $MW#")
if [ -n "${SEPAIR:-}" ]; then   # fix 6: w *= SEM^2 / (SEM^2 + SE_pair^2) for every target with an SE (key = NAME/target|cond)
  [ -s "$SEPAIR" ] || { echo "missing SEPAIR $SEPAIR"; exit 2; }
  FAC=$(awk -F, 'NR>1{printf "%s\"%s\": %s", (c++?", ":"{"), $1, $9} END{print (c?"}":"{}")}' "$SEPAIR")
  WEIGHTS=$(jq -c -n --argjson w "$WEIGHTS" --argjson f "$FAC" 'reduce ($f|keys[]) as $k ($w; .[$k] = (($w[$k] // 1) * $f[$k]))')
  echo "SEPAIR $SEPAIR: $(echo "$FAC" | jq 'length') target weights scaled by SEM^2/(SEM^2+SE_pair^2)"
fi
echo "weights $WEIGHTS"
mkdir -p $SCR/seeds
XF=(); [ -n "${HINGE:-}" ] && XF+=(--hinge "$HINGE"); [ -n "${INITADM:-}" ] && XF+=(--init-admissible "$INITADM"); [ -n "${STRATEGY:-}" ] && XF+=(--strategy "$STRATEGY")
if [ "${FITTER:-}" = fit_v8 ]; then [ -n "${PHASEB:-}" ] && XF+=(--phaseb "$PHASEB"); [ -n "${PHASEBK:-}" ] && XF+=(--phaseb-k "$PHASEBK"); [ -n "${POLISHK:-}" ] && XF+=(--polish-k "$POLISHK"); [ -n "${POLISHN:-}" ] && XF+=(--polish-nfev "$POLISHN"); fi
SD=${S:1}; SEEDARGS=(--seed $SD)
if [ "${S:0:1}" = s ]; then
  SJ=$SEEDJ
  if awk -v v=$A00 'BEGIN{exit !(v > 0)}'; then   # depression-band copies (user 2026-10-02): a00 -> A00, a10 shifted by the same amount
    for j in ${SEEDJ//,/ }; do
      [ "$(jq --argjson A $A00 '.a.a00 < $A' $j)" = true ] || continue
      o=$SCR/seeds/$(basename $j .json)_a00h${A00/./}.json
      jq --argjson A $A00 '.a as $a | .a.a00 = $A | .a.a20 = $A | .a.a10 = ($a.a10 + $A - $a.a00) | .a.a30 = .a.a10' $j > $o
      SJ=$SJ,$o; echo "depression-band seed $o: a00 $(jq '.a.a00' $j) -> $A00, a10 $(jq '.a.a10' $j) -> $(jq '.a.a10' $o)"
    done
  fi
  SEEDARGS+=(--seed-fits $SJ --seed-set '{}')
fi
COMMON=(--pairs $PAIRS --free-filters "$FREE" --fit-gamma --drop-targets "$DROPT,$DROPS" "${EXTRA[@]}" --weights "$WEIGHTS" --popsize $POP "${XF[@]}")
echo "FREE '${FREE}' (+gamma_d, gamma_p); SET $SET; POP $POP MAXITER $MAXITER (phase B $MAXB at rho_sigma 0) SIGMA $SIGMA"
if awk -v v=$SIGMA 'BEGIN{exit !(v > 0)}'; then
  MAXA=$((MAXITER - MAXB)); [ $MAXA -gt 0 ] || { echo "MAXITER must exceed MAXB"; exit 2; }
  SA=$SCR/${NAME}_${S}_sA; rm -f ${SA}_ckpt.npz
  SETA="${SET%\}}, \"rho_sigma\": $SIGMA}"
  run "${COMMON[@]}" --set "$SETA" "${SEEDARGS[@]}" --maxiter $MAXA --save $SA; [ $rc -eq 0 ] || exit $rc
  [ "${DRY:-0}" = 1 ] || [ -s ${SA}_ckpt.npz ] || { echo "phase A wrote no ${SA}_ckpt.npz"; exit 2; }
  run "${COMMON[@]}" --set "$SET0" --seed $SD --resume ${SA}_ckpt.npz --maxiter $MAXB --save $OUT; [ $rc -eq 0 ] || exit $rc
  [ "${DRY:-0}" = 1 ] || echo "phase A (rho_sigma $SIGMA) de_fun $(jq .de_fun $SA.json) -> phase B (binary) de_fun $(jq .de_fun $OUT.json)"
else
  run "${COMMON[@]}" --set "$SET0" "${SEEDARGS[@]}" --maxiter $MAXITER --save $OUT; [ $rc -eq 0 ] || exit $rc
fi
echo "=== done $(date)"
