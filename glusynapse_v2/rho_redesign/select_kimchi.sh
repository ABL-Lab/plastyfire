#!/bin/bash
# Apply the pinned selection rule (select_stage4.sh) to kimchi results: each result json carries .results.val_rows
# ([pathway, target, condition, mean, sem, pred, z, role]) and .results.n_free. Pure jq/awk, no compute.
# Usage (from anywhere): bash select_kimchi.sh SWEEP [extra rho_redesign RUN prefixes ...]
#   Converts <ledger>/results/SWEEP/*.json into RUN_val{,_l23,_l23l23}.csv + RUN.json in /scratch/dhuruva/kimchi_select/SWEEP,
#   links the core-set csvs and SE_pair csv there, then runs select_stage4.sh on all of them (OUT = rho_redesign/select_SWEEP.csv).
#   Extra args are rorqual runs in rho_redesign (e.g. f4Vdah_u8) to rank alongside.
#   Also writes basin_SWEEP.csv: per model, best de_fun and how many starts lie within 2 of it.
set -euo pipefail
SW=$1; shift
RS=/project/rrg-emuller/dhuruva/plastyfire/glusynapse_v2/rho_redesign
LED=/lustre09/project/6070394/dhuruva/kimchi_ledger_glusyn/results/$SW
W=/scratch/dhuruva/kimchi_select/$SW; mkdir -p $W
for f in $RS/f4Ndah_s6.csv $RS/f4Ndah_s6_l23.csv $RS/f4Ndah_s6_l23l23.csv $RS/s4N_u_sepair.csv; do ln -sf $f $W/; done
runs=()
for j in $LED/*.json; do
  r=$(basename $j .json)
  [ "$(jq -r '.results.val_rows | length' $j 2>/dev/null || echo 0)" -gt 0 ] || { echo "skip $r (no val_rows)" >&2; continue; }
  for p in "l5:" "l23:_l23" "l23l23:_l23l23"; do pw=${p%%:*}; s=${p#*:}
    { echo "target,condition,target_mean,target_sem,pred,z"
      jq -r --arg pw $pw '.results.val_rows[] | select(.[0]==$pw) | [.[1],.[2],.[3],.[4],.[5],.[6]] | @csv' $j | tr -d '"'; } > $W/${r}_val$s.csv
  done
  jq '{n_free: .results.n_free, de_fun: .results.de_fun, model: .results.model, cluster: .cluster}' $j > $W/$r.json
  runs+=($r)
done
for r in "$@"; do for s in "" _l23 _l23l23; do ln -sf $RS/${r}_val$s.csv $W/; done; ln -sf $RS/$r.json $W/; runs+=($r); done
cd $W && OUT=$RS/select_$SW.csv bash $RS/select_stage4.sh "${runs[@]}"
# basin census per model (kimchi runs only): best de_fun and starts within 2
for r in $(ls $LED/*.json); do jq -r '[(.run | split("_")[1]), .run, .results.de_fun // ""] | @csv' $r; done | tr -d '"' |
  awk -F, '$3!=""{m=$1; if(!(m in b) || $3<b[m]) b[m]=$3; all[NR]=$0}
    END{print "model,best_de_fun,n_runs,n_within2"; for(i in all){split(all[i],a,","); n[a[1]]++; if(a[3]-b[a[1]]<=2) w[a[1]]++}
        for(m in b) printf "%s,%.2f,%d,%d\n", m, b[m], n[m], w[m]}' | sort | tee $RS/basin_$SW.csv | column -t -s,
