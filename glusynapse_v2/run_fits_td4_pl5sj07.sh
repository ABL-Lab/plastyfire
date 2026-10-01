#!/bin/bash
# Submits the t_drive 4 fits (3 seeds) and the t_drive 3 reference on the 5 -vca dirs, paired_l5 + sjostrom07 (29 targets).
# Sizing: td3 fits on 4 dirs measured 24-27 GB, 33-72 min, 97% CPU (DECISIONS.md) -> 33G; sj07 adds 0.93 of 6.6 GB compressed
# (+14%, ~+3.7 GB in memory) -> (27 + 3.7) x 1.25 = 39G; time 1:45 as td3 (72 min + 50%). seff unavailable when submitted (slurmdb down).
# 2026-09-30 moved to 1g.10gb (3g queue 681 pending / 16 running; our GPU fairshare 0.10). Pilots 22098228 (1g, CHUNK 4) and 22098229 (2g, CHUNK 10):
# both ok, objective call 13-15 s on 1g, 9-10 s on 2g, vs 7 s on a full H100 (td3_s5 ref 64 min) -> 1g ~2.1x = ~135 min, +50% -> 3:30. Checkpointed (RESUME).
V2=glusynapse_v2; X=$V2/extracted
DIRS=$X/ebner_delta-prefire-vca,$X/markram_delta-prefire-vca,$X/sj03_delta-prefire-vca,$X/sj03r50_delta-prefire-vca,$X/sj07_delta-prefire-vca
FREE=theta_Te,tau_T,theta_Tg,dpre_min,theta_NOi,tau_NO,tau_Z,theta_Z
S5=$V2/results/reduced_gpu_subset_pl5r50_td3_s5.json
run() {  # td seed seedfits
  env MODE=subset TGROUPS=paired_l5,sjostrom07 CONDS=control,mglu_block,post_nmdar,nmdar_block,no_block MAXITER=150 SEED=$2 \
    TAG=_pl5sj07_td$1_s$2 DIRS=$DIRS FREE=$FREE SEEDFITS=${3:-} CHUNK=4 \
    FILTERS="{\"pre_drive\": 1, \"i_scale\": 1e-5, \"t_drive\": $1, \"tau_E1\": 100.0}" \
    sbatch --parsable --account=def-emuller --gpus-per-node=nvidia_h100_80gb_hbm3_1g.10gb:1 --cpus-per-task=1 --mem=39G --time=03:30:00 \
      --job-name=fit_td$1_s$2 --output=/project/rrg-emuller/dhuruva/plastyfire/logs/fit_pl5sj07_td$1_s$2_%j.out $V2/run_reduced_gpu.sh
}
cd /project/rrg-emuller/dhuruva/plastyfire
echo "td4 s1 $(run 4 1 $S5)"; echo "td4 s2 $(run 4 2 $S5)"; echo "td4 s3 $(run 4 3)"; echo "td3 s1 $(run 3 1 $S5)"
