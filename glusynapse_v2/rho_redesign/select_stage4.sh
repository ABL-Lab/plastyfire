#!/bin/bash
# Pinned screen selection rule (DECISIONS 2026-10-02), pure awk on small csvs, no compute.
#   1. Markram HIT: 10Hz_10ms and 10Hz_-10ms (L5 control) both |z| <= 1 on the data SEM.
#   2. Among hits, rank by chi2_eff/n over the CORE rows, with SEM_eff^2 = SEM^2 + SE_pair^2 (SE_pair from LOPO, SEPAIR csv).
#   3. Within 10 % of the best chi2_eff/n counts as tied; ties go to fewer free parameters (n_free from the run json).
#   4. (user 2026-10-02) L5 slope gate: OLS slope of pred on data over all 40 L5 rows must be >= 0.6; hits that pass rank
#      first; if none passes, the hit with the highest slope leads. This is the one place validation rows enter.
#   Validation rows are never used: CORE = the fit targets of the CORE_FROM run (one common set for every model).
# Usage (from rho_redesign/): bash select_stage4.sh RUN [RUN ...]
#   RUN = prefix of RUN_val{,_l23,_l23l23}.csv and RUN.json.
#   env CORE_FROM (default f4Ndah_s6: 11 L5 + 9 L2/3->L5 + 12 L2/3->L2/3), SEPAIR (default s4N_u_sepair.csv), OUT (select_stage4.csv)
set -euo pipefail
CORE_FROM=${CORE_FROM:-f4Ndah_s6}; SEPAIR=${SEPAIR:-s4N_u_sepair.csv}; OUT=${OUT:-select_stage4.csv}
core=$(for p in "l5:" "l23:_l23" "l23l23:_l23l23"; do pw=${p%%:*}; s=${p#*:}
  [ -s $CORE_FROM$s.csv ] && awk -F, -v pw=$pw 'NR>1{print pw"/"$1"|"$2}' $CORE_FROM$s.csv; done | paste -sd' ')
TMP=$(mktemp -p /scratch/dhuruva); trap "rm -f $TMP" EXIT
for R in "$@"; do
  nf=$(jq -r '[.. | objects | .n_free? // empty] | first // empty' $R.json 2>/dev/null)
  [ -n "$nf" ] || nf=$(jq -r '(.x // .params // {}) | length' $R.json 2>/dev/null || echo NA)
  for p in "l5:" "l23:_l23" "l23l23:_l23l23"; do pw=${p%%:*}; s=${p#*:}
    awk -F, -v pw=$pw 'NR>1{print pw"/"$1"|"$2","$3","$4","$5","$6}' ${R}_val$s.csv; done |
  awk -F, -v core="$core" -v run=$R -v nf=$nf -v sp=$SEPAIR '
    BEGIN{n=split(core,c," "); for(i=1;i<=n;i++) isc[c[i]]=1
          while((getline l < sp)>0){split(l,a,","); if(a[1]!="key") se[a[1]]=a[7]}}
    {k=$1; m=$2; s=$3; pr=$4; z=$5
     if(k=="l5/10Hz_10ms|control"){p10=pr; z10=z} if(k=="l5/10Hz_-10ms|control"){m10=pr; zm10=z}
     if(k in isc){e=(k in se)?se[k]:0; chi+=(pr-m)^2/(s^2+e^2); nn++}
     if(k ~ /^l5\//){q++; sx+=m; sy+=pr; sxx+=m*m; sxy+=m*pr}}
    END{hit=(z10*z10<=1 && zm10*zm10<=1)?"HIT":"MISS"; sl=(q*sxy-sx*sy)/(q*sxx-sx*sx); g=(sl>=0.6)?"PASS":"low"
        key=(g=="PASS")?chi/nn:1000-sl
        printf "%s,%s,%s,%.2f,%.3f,%.3f,%.2f,%d,%.3f,%s,%.4f\n", run,hit,g,sl,p10,m10,chi,nn,chi/nn,nf,key}'
done > $TMP
# rank: hits first; among hits, L5 slope gate PASS (>= 0.6) first ranked by chi2eff_n, then gate-low ranked by highest slope;
# mark ties within 10 % of the best PASS hit (ties go to fewer free params)
sort -t, -k2,2 -k3,3 -k11,11g $TMP | cut -d, -f1-10 | awk -F, 'BEGIN{print "run,hit,slope_gate,slopeL5,mk_p10,mk_m10,chi2eff,n,chi2eff_n,n_free,tied_best"}
  {t=""; if($2=="HIT" && $3=="PASS"){ if(best=="") best=$9; t=($9<=1.1*best)?"TIED":"" } print $0","t}' | tee $OUT | column -t -s,
