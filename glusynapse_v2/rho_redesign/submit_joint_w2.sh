#!/bin/bash
# Wave-2 joint fit submit (login node: sbatch only). Jobs J0 (v5c rescore, v7 off) J1 (V7v_s5 rescore) -> J2 seeded, J3/J4 unseeded afterok J1.
set -euo pipefail
R=/lustre09/project/6070394/dhuruva/plastyfire; V2=$R/glusynapse_v2; X=$V2/extracted; RS=$V2/rho_redesign; RES=$RS/results
COMMON=(L5DIRS=$X/ebner_delta-split1-prefire-vca,$X/markram_delta-split1-prefire-vca,$X/sj03_delta-split1-prefire-vca,$X/sj03r50_delta-split1-prefire-vca,$X/sj07_delta-split1-prefire-vca,$X/l5extra_delta-split1-prefire-vca
  L5GROUPS=paired_l5,sjostrom07,paired_l5_extra
  L23DIRS=$X/ebner_l23l5_delta-split1-prefire-vseg-rs ANALYTICAL_BASIS_DIR=/scratch/dhuruva/split1/basis_l5l5
  L23_BASIS_DIR=/scratch/dhuruva/split1/basis_l23l5 JOINT=1 FITGAMMA=1 FREE=theta_V,theta_eCB
  "DROPT=letzkus_3ap_200hz_dt-10ms@distal|control,letzkus_3ap_200hz_dt-10ms@distal|nmdar_block"
  "EXTRA=l23l23:paired_l23l23,paired_l23l23_egger:$X/zilberter_l23l23_delta-split1-prefire-vseg-rs,$X/l23l23extra_delta-split1-prefire-vseg-rs:/scratch/dhuruva/split1/basis_l23l23")
OFF='SET={"v5_mode": 2}'; ON='SET={"v5_mode": 2, "veto_T": 25.0}'
V5=$RES/v5_S1_V5c_s6_39.json; V7=$RES/v7_S1_V7v_s5.json
G=(--account=def-emuller --gpus-per-node=nvidia_h100_80gb_hbm3_3g.40gb:1)
sub() { local dep=$1 name=$2 time=$3; shift 3; local d=(); [ "$dep" != none ] && d=(--dependency=$dep)
  env "${COMMON[@]}" "$@" sbatch --parsable "${d[@]}" "${G[@]}" -J $name --mem=${MEM:-78G} --time=$time $RS/run_joint_w2.sh; }
j0=$(sub none W2_J0 00:15:00 "$OFF" TAG=W2_J0 SEED=1 MAXITER=0 SEEDFITS=$V5)
j1=$(sub none W2_J1 00:15:00 "$ON" TAG=W2_J1 SEED=1 MAXITER=0 SEEDFITS=$V7)
D=afterok:$j1
j2=$(sub $D W2_J2 00:30:00 "$ON" TAG=W2_J2 SEED=5 SEEDFITS=$V7)
j3=$(sub $D W2_J3 00:30:00 "$ON" TAG=W2_J3 SEED=6)
j4=$(sub $D W2_J4 00:30:00 "$ON" TAG=W2_J4 SEED=7)
echo "J0 $j0 J1 $j1 J2 $j2 J3 $j3 J4 $j4"
