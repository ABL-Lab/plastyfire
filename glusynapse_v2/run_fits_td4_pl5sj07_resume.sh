#!/bin/bash
# Resume the pl5sj07 fits (run_fits_td4_pl5sj07.sh, 22098444-7) that hit their 3:30 limit, from <save>_ckpt.npz, for the remaining
# 150 - nit generations. --resume starts DE from the saved population with a fresh maxiter.
# Sizing (measured on 22098444-7, 2026-09-30): sstat MaxRSS 32-35 GB -> 44G. Objective calls grew to 108-147 s on 1g by nit 85;
# on 2g.20gb the pilots ran 0.65x the 1g call time -> time = rem x 150 s x 0.7 x 1.5, rounded up to 15 min. 2g + CHUNK 10 (pilot 22098229).
# 2026-09-30 17:26: switched to the v3 kernel (run_reduced_gpu_v3.sh, gpu_v3.py; == v2 to 1e-15 on td4, 22107148). Pilot 22107523 1g: 4.4-4.9 s/gen,
# 26.3 GB MaxRSS, 1:52 wall for 3 gens -> ~40 gens = ~4 min + load -> 33G, 0:30.
# Usage: bash run_fits_td4_pl5sj07_resume.sh <jobid>...   (only jobs whose log has no final "chi2" line are resumed)
V2=glusynapse_v2; X=$V2/extracted; LOGS=/project/rrg-emuller/dhuruva/plastyfire/logs
DIRS=$X/ebner_delta-prefire-vca,$X/markram_delta-prefire-vca,$X/sj03_delta-prefire-vca,$X/sj03r50_delta-prefire-vca,$X/sj07_delta-prefire-vca
FREE=theta_Te,tau_T,theta_Tg,dpre_min,theta_NOi,tau_NO,tau_Z,theta_Z
cd /project/rrg-emuller/dhuruva/plastyfire
for j in "$@"; do
  f=$(ls $LOGS/fit_pl5sj07_*_$j.out) || continue
  grep -q '^chi2' $f && { echo "$j finished, no resume"; continue; }
  tag=$(basename $f .out); tag=${tag#fit_pl5sj07_}; tag=${tag%_$j}          # e.g. td4_s1
  td=${tag%%_*}; td=${td#td}; seed=${tag##*_s}
  nit=$(grep -o 'ckpt nit [0-9]*' $f | tail -1 | awk '{print $3}'); rem=$((150 - nit))
  [ $rem -le 0 ] && { echo "$j at nit $nit, nothing left"; continue; }
  min=$(( (rem*150*7*15/100/60 + 14) / 15 * 15 )); tl=$(printf '%d:%02d:00' $((min/60)) $((min%60)))
  ck=$V2/results/reduced_gpu_subset_pl5sj07_${tag}_ckpt.npz
  cp $ck ${ck%.npz}_nit$nit.npz                                               # keep the checkpoint the resume starts from
  id=$(env MODE=subset TGROUPS=paired_l5,sjostrom07 CONDS=control,mglu_block,post_nmdar,nmdar_block,no_block MAXITER=$rem SEED=$seed \
    TAG=_pl5sj07_$tag DIRS=$DIRS FREE=$FREE RESUME=${ck%.npz}_nit$nit.npz CHUNK=10 \
    FILTERS="{\"pre_drive\": 1, \"i_scale\": 1e-5, \"t_drive\": $td, \"tau_E1\": 100.0}" \
    sbatch --parsable --account=def-emuller --gpus-per-node=nvidia_h100_80gb_hbm3_1g.10gb:1 --cpus-per-task=1 --mem=33G --time=00:30:00 \
      --job-name=fit_${tag}_r --output=$LOGS/fit_pl5sj07_${tag}_r_%j.out $V2/run_reduced_gpu_v3.sh)
  echo "$j $tag nit $nit -> resume $id ($rem gens, $tl)"
done
