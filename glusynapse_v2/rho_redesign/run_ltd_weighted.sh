#!/bin/bash
# A7: LTD-weighted refits of v5c on delta-split1 (39 targets). WEIGHTS env (fit_v5.py): z^2 x w for every target whose
# data mean < 1.0 in the DE objective only; json/csv report unweighted chi2, chi2_ltd/chi2_ltp, AIC/BIC unweighted.
# Rescore S1_V5cW_r (MAXITER 0, WEIGHTS ltd 2, seed s6_39) must print REPRO OK and unweighted chi2 87.84; fits afterok on it.
# Sizing (measured, run_split1_refits.sh / run_fit_v5.sh): v5 check 2:44 31.6 GB -> 40G 0:15; v5 fits 12-15 min 27-39 GB -> 50G 0:30.
# MEASURED: rescore 22154961 4:43, 37.25 GB (of 40G, 52% CPU); fits 22154962-5 11:43-13:15, 26.2-37.2 GB (of 50G, 95-99% CPU).
set -euo pipefail
R=/lustre09/project/6070394/dhuruva/plastyfire; V2=$R/glusynapse_v2; X=$V2/extracted; RS=$V2/rho_redesign; RES=$RS/results
COMMON=(L5DIRS=$X/ebner_delta-split1-prefire-vca,$X/markram_delta-split1-prefire-vca,$X/sj03_delta-split1-prefire-vca,$X/sj03r50_delta-split1-prefire-vca,$X/sj07_delta-split1-prefire-vca
  L23DIRS=$X/ebner_l23l5_delta-split1-prefire-vseg-rs ANALYTICAL_BASIS_DIR=/scratch/dhuruva/split1/basis_l5l5
  L23_BASIS_DIR=/scratch/dhuruva/split1/basis_l23l5 JOINT=1 FITGAMMA=1
  "DROPT=letzkus_3ap_200hz_dt-10ms@distal|control,letzkus_3ap_200hz_dt-10ms@distal|nmdar_block")
V5=("${COMMON[@]}" 'SET={"v5_mode": 2}' FREE=theta_V,theta_eCB)
G2=(--account=def-emuller --gpus-per-node=nvidia_h100_80gb_hbm3_2g.20gb:1)
SEEDS=$RES/v5_S1_V5c_s6_39.json,$RES/v5_S1_V5c_s5_39.json
j0=$(env "${V5[@]}" 'WEIGHTS={"ltd": 2.0}' TAG=S1_V5cW_r SEED=1 MAXITER=0 SEEDFITS=$RES/v5_S1_V5c_s6_39.json \
     sbatch --parsable "${G2[@]}" -J S1_V5cW_r --mem=40G --time=00:15:00 $RS/run_fit_v5.sh)
out="S1_V5cW_r $j0"
for w in 2 4; do
  W="WEIGHTS={\"ltd\": $w.0}"
  a=$(env "${V5[@]}" "$W" TAG=S1_V5cW${w}_s5 SEED=5 SEEDFITS=$SEEDS sbatch --parsable --dependency=afterok:$j0 "${G2[@]}" -J S1_V5cW${w}_s5 --mem=50G --time=00:30:00 $RS/run_fit_v5.sh)
  b=$(env "${V5[@]}" "$W" TAG=S1_V5cW${w}_s6 SEED=6 sbatch --parsable --dependency=afterok:$j0 "${G2[@]}" -J S1_V5cW${w}_s6 --mem=50G --time=00:30:00 $RS/run_fit_v5.sh)
  out="$out S1_V5cW${w}_s5 $a S1_V5cW${w}_s6 $b"
done
echo "$out"
