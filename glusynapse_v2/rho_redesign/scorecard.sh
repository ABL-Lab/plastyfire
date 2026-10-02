#!/bin/bash
# Fair scorecard for staged fits (FIT_META_DIAG.md fix 1, user-approved 2026-10-02). Pure awk on small csvs, no compute.
# Usage (from anywhere): bash glusynapse_v2/rho_redesign/scorecard.sh [MODEL[:FITCSV] ...]
#   MODEL   prefix of the val rescore csvs RS/MODEL_val{,_l23,_l23l23}.csv (columns target,condition,target_mean,target_sem,pred,z)
#   FITCSV  prefix of the fit csvs (rows = fitted targets); default RS/MODEL; absolute path allowed. A missing pathway
#           fit csv means nothing was fitted on that pathway (all its val rows are held out).
#   env EXCL (default sjostrom_40hz_dt0ms|control): one L5 row also reported excluded (it has SEM 0.038 and dominates).
#   env OUT  (default RS/scorecard.csv)
# Metrics (val csv; z is unweighted):
#   hoL5/hoL23/hoL2323 = chi2 over HELD-OUT rows only (not in the fit csv), "chi2/n"; hoX = total held-out chi2 without EXCL.
#   slope = OLS slope of pred on data over all val rows (PASS >= 0.6, STAGE3_DESIGN.md section 4); slopeHO over held-out rows.
#   sign  = non-null rows (|mean - 1| > 0.1) with pred on the same side of 1 as the data, all rows / held-out rows.
#   Markram = 10Hz_10ms and 10Hz_-10ms (L5 control) pred and z; HIT when both |z| <= 1 (hard pass), else MISS.
# Default model list = the 2026-10-02 comparison set (s1C_w fit csvs live in /scratch/dhuruva/s1c_l23w/s1CW_w025).
set -euo pipefail
RS=/project/rrg-emuller/dhuruva/plastyfire/glusynapse_v2/rho_redesign
EXCL=${EXCL:-sjostrom_40hz_dt0ms|control}
OUT=${OUT:-$RS/scorecard.csv}
SPECS=("$@")
[ ${#SPECS[@]} -gt 0 ] || SPECS=(s1C_s s3C_s s3C_u s3N_s s3N_u s4N_s s4N_u s5N_s s3V_s s3W_s s4V_s s4V_u s1C_w:/scratch/dhuruva/s1c_l23w/s1CW_w025)
echo "model,hoL5_chi2,hoL5_n,hoL5x_chi2,hoL5x_n,hoL23_chi2,hoL23_n,hoL2323_chi2,hoL2323_n,ho_chi2,ho_n,hoX_chi2,hoX_n,slope,slopeHO,slope_pass,sign_hit,sign_n,signHO_hit,signHO_n,m10_pred,m10_z,mm10_pred,mm10_z,markram" > $OUT
printf "%-8s %-15s %-15s %-12s %-12s %-15s %-6s %-7s %-7s %-7s %s\n" model "hoL5 chi2/n" "hoL5-dt0" "hoL23L5" "hoL23L23" "hoTot-dt0" slope slopeHO sign signHO "Markram +10 / -10"
for spec in "${SPECS[@]}"; do
  m=${spec%%:*}; fc=$m; [ "$spec" != "$m" ] && fc=${spec#*:}
  case $fc in /*) ;; *) fc=$RS/$fc ;; esac
  [ -s $RS/${m}_val.csv ] || { echo "skip $m (no $RS/${m}_val.csv)"; continue; }
  for s in "" _l23 _l23l23; do
    pw=L5; [ "$s" = _l23 ] && pw=L23; [ "$s" = _l23l23 ] && pw=L2323
    [ -s $RS/${m}_val$s.csv ] || continue
    awk -F, -v pw=$pw 'FNR==NR{ if(FNR>1) F[$1"|"$2]=1; next } FNR>1{ print pw","(($1"|"$2) in F ? 1 : 0)","$1"|"$2","$3","$4","$5","$6 }' \
      <( [ -s $fc$s.csv ] && cat $fc$s.csv || echo header ) $RS/${m}_val$s.csv
  done | awk -F, -v m=$m -v excl="$EXCL" -v out=$OUT '
    function ols(sx,sy,sxx,sxy,n){ return (n*sxx-sx*sx)!=0 ? (n*sxy-sx*sy)/(n*sxx-sx*sx) : 0 }
    { pw=$1; f=$2; k=$3; y=$4; p=$6; z=$7; z2=z*z
      sx+=y; sy+=p; sxx+=y*y; sxy+=y*p; n++
      if(!f){ hc[pw]+=z2; hn[pw]++; hx+=y; hy+=p; hxx+=y*y; hxy+=y*p; hnn++
              if(!(pw=="L5" && k==excl)){ hxc[pw]+=z2; hxn[pw]++ } }
      if((y-1>0.1)||(1-y>0.1)){ sn++; ok=((y>1)==(p>1)); sh+=ok; if(!f){ shn++; shh+=ok } }
      if(pw=="L5" && k=="10Hz_10ms|control"){ mp=p; mz=z } if(pw=="L5" && k=="10Hz_-10ms|control"){ mmp=p; mmz=z } }
    END{ sl=ols(sx,sy,sxx,sxy,n); slh=ols(hx,hy,hxx,hxy,hnn)
      tot=hc["L5"]+hc["L23"]+hc["L2323"]; totn=hn["L5"]+hn["L23"]+hn["L2323"]
      tx=hxc["L5"]+hxc["L23"]+hxc["L2323"]; txn=hxn["L5"]+hxn["L23"]+hxn["L2323"]
      mk=((mz<=1&&mz>=-1)&&(mmz<=1&&mmz>=-1))?"HIT":"MISS"
      printf "%s,%.2f,%d,%.2f,%d,%.2f,%d,%.2f,%d,%.2f,%d,%.2f,%d,%.3f,%.3f,%s,%d,%d,%d,%d,%.3f,%.2f,%.3f,%.2f,%s\n", m, hc["L5"],hn["L5"],hxc["L5"],hxn["L5"],hc["L23"],hn["L23"],hc["L2323"],hn["L2323"],tot,totn,tx,txn,sl,slh,(sl>=0.6?"PASS":"FAIL"),sh,sn,shh,shn,mp,mz,mmp,mmz,mk >> out
      printf "%-8s %-15s %-15s %-12s %-12s %-15s %-6s %-7s %-7s %-7s %s %.3f(z%+.1f) / %.3f(z%+.1f)\n", m, sprintf("%.1f/%d",hc["L5"],hn["L5"]), sprintf("%.1f/%d",hxc["L5"],hxn["L5"]),
        sprintf("%.1f/%d",hc["L23"],hn["L23"]), sprintf("%.1f/%d",hc["L2323"],hn["L2323"]), sprintf("%.1f/%d",tx,txn),
        sprintf("%.2f%s",sl,(sl>=0.6?"":"x")), sprintf("%.2f",slh), sprintf("%d/%d",sh,sn), sprintf("%d/%d",shh,shn), mk, mp, mz, mmp, mmz }'
done
echo "wrote $OUT"
