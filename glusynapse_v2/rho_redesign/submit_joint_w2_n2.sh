#!/bin/bash
# Wave-2 joint V7vn2 (ecb_ref 2; CEXP_DIR+SCALE_E set by gpu_v7_rho.py; CEXP_MAP default l23l23 -> L23L23). No dependencies.
set -euo pipefail
R=/lustre09/project/6070394/dhuruva/plastyfire; V2=$R/glusynapse_v2; X=$V2/extracted; RS=$V2/rho_redesign; RES=$RS/results
COMMON=(L5DIRS=$X/ebner_delta-split1-prefire-vca,$X/markram_delta-split1-prefire-vca,$X/sj03_delta-split1-prefire-vca,$X/sj03r50_delta-split1-prefire-vca,$X/sj07_delta-split1-prefire-vca,$X/l5extra_delta-split1-prefire-vca
  L5GROUPS=paired_l5,sjostrom07,paired_l5_extra
  L23DIRS=$X/ebner_l23l5_delta-split1-prefire-vseg-rs ANALYTICAL_BASIS_DIR=/scratch/dhuruva/split1/basis_l5l5
  L23_BASIS_DIR=/scratch/dhuruva/split1/basis_l23l5 JOINT=1 FITGAMMA=1 FREE=theta_V,theta_eCB
  "DROPT=letzkus_3ap_200hz_dt-10ms@distal|control,letzkus_3ap_200hz_dt-10ms@distal|nmdar_block"
  "EXTRA=l23l23:paired_l23l23,paired_l23l23_egger:$X/zilberter_l23l23_delta-split1-prefire-vseg-rs,$X/l23l23extra_delta-split1-prefire-vseg-rs:/scratch/dhuruva/split1/basis_l23l23")
ON='SET={"v5_mode": 2, "veto_T": 25.0, "ecb_ref": 2}'
VN=$RES/v7_S1_V7vn2_s5.json
G=(--account=def-emuller --gpus-per-node=nvidia_h100_80gb_hbm3_3g.40gb:1)
sub() { env "${COMMON[@]}" "${@:3}" sbatch --parsable "${G[@]}" -J $1 --mem=${MEM:-78G} --time=$2 $RS/run_joint_w2.sh; }
n2=$(sub W2_N2 00:30:00 "$ON" TAG=W2_N2 SEED=5 SEEDFITS=$VN)
n3=$(sub W2_N3 00:30:00 "$ON" TAG=W2_N3 SEED=6)
echo "N2 $n2 N3 $n3"
