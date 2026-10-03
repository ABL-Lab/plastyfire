#!/bin/bash
# Stage-4 census (2026-10-02): models 3W, 4V, 4N x starts s5 s6 (seeded) u7 u8 (unseeded), with fit_v7 options:
# admissible init, Markram +/-10 hinge (zero inside 1 SEM, lambda 100), DEMOTE=1, FIXANC=1, a00 1.2 seeds, POP 15, 600 gen.
# Each fit is followed by its validation rescore (afterok). Pick the basin afterwards: run_stage4_fit.sh MODEL - pick.
# Usage (from plastyfire/): bash glusynapse_v2/rho_redesign/submit_stage4_census.sh
# Sizing:
#   fits: memory = stage-3 MaxRSS of the same model x 1.25 (3W 23.0 GB -> 29G; 4V 50.2 GB -> 63G; 4N 65.5 GB -> 82G;
#         smoke 22308544 4N 58.7 GB fits inside). Time: smoke 12:01 for 7 generations + 600 x 1.3 s, so ~25 min -> 1:00:00
#         (extrapolated, not measured; tighten after the first finishes).
#   vals: MEASURED peak 72.9 GB 9:52 (22300527) -> 92G 0:30:00.
# MEASURED census 22314719-42: fits 3W 12:34-15:42 peak 19.0 GB; 4V 27:41-36:08 peak 56.6 GB; 4N 25:31-38:15 peak 58.7 GB;
#   vals 7:25-8:59 peak 62.0 GB. Next: fits 3W 24G 0:30, 4V 71G 1:00, 4N 74G 1:00; vals 78G 0:15; kimchi fit+val 78G 1:15.
set -euo pipefail
RS=glusynapse_v2/rho_redesign
# options live in stage4_env.sh so that the batch jobs see them regardless of sbatch --export
for spec in "3W:29G:$RS/s3W_s.json,$RS/s3C_u.json" "4V:63G:$RS/s4V_u.json,$RS/s4V_s.json" "4N:82G:$RS/s4N_u.json,$RS/s4Nm8_s.json"; do
  M=${spec%%:*}; rest=${spec#*:}; MEM=${rest%%:*}; SEEDS=${rest#*:}
  for S in s5 s6 u7 u8; do
    f=$(glusynapse_v2/sjob.sh -g -n f4_${M}_$S -m $MEM -t 01:00:00 -b "stage-3 $M MaxRSS x1.25; smoke 22308544 12:01 + 600 gen x 1.3 s" -- \
        "source $RS/stage4_env.sh; SEEDX=$SEEDS bash $RS/run_stage4_fit.sh $M $S fit")
    v=$(glusynapse_v2/sjob.sh -g -n v4_${M}_$S -m 92G -t 00:30:00 -d afterok:$f -b 'MEASURED 22300527 val 72.9 GB 9:52' -- \
        "source $RS/stage4_env.sh; SEEDX=$SEEDS bash $RS/run_stage4_fit.sh $M $S val")
    echo "$M $S fit $f val $v"
  done
done
