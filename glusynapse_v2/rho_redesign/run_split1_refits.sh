#!/bin/bash
# delta-split1 refits (login node: sbatch only). v5c (best on BIC) + v4 C1 (best chi2), seeded from V5c_s7 / C1Ajn_s5,
# + unseeded basin check, + MAXITER 0 rescore of the old fits on split1 data, + 54-target v5c after ext_l23l23.
# Data: extracted/*_delta-split1-prefire-* (sj03 dir includes r50 bursts), bases /scratch/dhuruva/split1/basis_*.
# Sizing (measured, run_fit_v5.sh / run_fit_v4.sh headers): v5 check 2:44 31.6 GB -> 40G 0:15; v5 fits 12-15 min
#   27-39 GB -> 50G 0:30; v4 C1 joint fits 20-26 min 39.5 GB -> 50G 0:45; 3-path v5c 15:53 63.0 GB -> 3g 79G 0:30.
set -euo pipefail
R=/lustre09/project/6070394/dhuruva/plastyfire; V2=$R/glusynapse_v2; X=$V2/extracted; RS=$V2/rho_redesign; RES=$RS/results
DEP39=afterok:${EXT_L5:-22136800}:${EXT_L23L5:-22136807}; DEP54=$DEP39:${EXT_L23L23:-22136812}   # 22136800/807 TIMEOUT (MaxRSS 26.1/34.6 GB); reruns EXT_L5/EXT_L23L5 env; DEP54=$DEP39:${EXT_L23L23:-22136812}
COMMON=(L5DIRS=$X/ebner_delta-split1-prefire-vca,$X/markram_delta-split1-prefire-vca,$X/sj03_delta-split1-prefire-vca,$X/sj03r50_delta-split1-prefire-vca,$X/sj07_delta-split1-prefire-vca
  L23DIRS=$X/ebner_l23l5_delta-split1-prefire-vseg-rs ANALYTICAL_BASIS_DIR=/scratch/dhuruva/split1/basis_l5l5
  L23_BASIS_DIR=/scratch/dhuruva/split1/basis_l23l5 JOINT=1 FITGAMMA=1
  "DROPT=letzkus_3ap_200hz_dt-10ms@distal|control,letzkus_3ap_200hz_dt-10ms@distal|nmdar_block")
V5=("${COMMON[@]}" 'SET={"v5_mode": 2}' FREE=theta_V,theta_eCB)
V4=("${COMMON[@]}" 'SET={"vamp_mode": 1}' FREE=theta_Te,tau_T,theta_Tg,dpre_min,theta_NOi,tau_NO,tau_Z,theta_Z,theta_V,rho_gamma)
G2=(--account=def-emuller --gpus-per-node=nvidia_h100_80gb_hbm3_2g.20gb:1)
sub() { local dep=$1 opts=$2 script=$3; shift 3; env "$@" sbatch --parsable --dependency=$dep "${G2[@]}" $opts $RS/$script; }
s=${SEEDV5:-$RES/v5_V5c_s7.json}; c=${SEEDC1:-$RES/v4_C1Ajn_s5.json}; T=${TAGSFX:-}   # seeds / tag suffix overridable
j1=$(sub $DEP39 "-J S1_V5c_r --mem=40G --time=00:15:00" run_fit_v5.sh "${V5[@]}" TAG=S1_V5c_r$T SEED=1 MAXITER=0 SEEDFITS=$s)
j2=$(sub $DEP39 "-J S1_V5c_s5 --mem=50G --time=00:30:00" run_fit_v5.sh "${V5[@]}" TAG=S1_V5c_s5$T SEED=5 SEEDFITS=$s,$c)
j3=$(sub $DEP39 "-J S1_V5c_s6 --mem=50G --time=00:30:00" run_fit_v5.sh "${V5[@]}" TAG=S1_V5c_s6$T SEED=6)
j4=$(sub $DEP39 "-J S1_C1_r --mem=40G --time=00:15:00" run_fit_v4.sh "${V4[@]}" TAG=S1_C1_r$T SEED=1 MAXITER=0 SEEDFITS=$c)
j5=$(sub $DEP39 "-J S1_C1_s5 --mem=50G --time=00:45:00" run_fit_v4.sh "${V4[@]}" TAG=S1_C1_s5$T SEED=5 SEEDFITS=$c)
Z="l23l23:paired_l23l23:$X/zilberter_l23l23_delta-split1-prefire-vseg-rs:/scratch/dhuruva/split1/basis_l23l23"
j6=$(env "${V5[@]}" TAG=S1_V5cz_s5$T SEED=5 SEEDFITS=$RES/v5_V5cz_s5.json,$s EXTRA="$Z" sbatch --parsable --dependency=$DEP54 \
     --account=def-emuller --gpus-per-node=nvidia_h100_80gb_hbm3_3g.40gb:1 -J S1_V5cz_s5 --mem=79G --time=00:30:00 $RS/run_fit_v5.sh)
echo "S1_V5c_r $j1 S1_V5c_s5 $j2 S1_V5c_s6 $j3 S1_C1_r $j4 S1_C1_s5 $j5 S1_V5cz_s5 $j6"
