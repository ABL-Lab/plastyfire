#!/bin/bash
# A14 wave 2: G1 shaft gate (gpu_v8_rho.py) on the v7x arrival timing (SET t_exact 1, gpu_v7x_rho). Env = submit_gate_v8.sh.
# Sizing: measured G8 src-2 runs: rescore 22177166 68.0 GB, fits 16-23 min at 40-68 GB -> 87G; fits 0:35, rescore 0:20.
# G1a_r0 = t_exact 0 MAXITER 0 on v7_G8_G1a_s5 (log must print REPRO OK; exit 0 either way, so nothing depends on it).
set -euo pipefail
R=/lustre09/project/6070394/dhuruva/plastyfire; V2=$R/glusynapse_v2; X=$V2/extracted; RS=$V2/rho_redesign; RES=$RS/results
COMMON=(KERNEL=gpu_v8_rho.py
  L5DIRS=$X/ebner_delta-split1-prefire-vca,$X/markram_delta-split1-prefire-vca,$X/sj03_delta-split1-prefire-vca,$X/sj03r50_delta-split1-prefire-vca,$X/sj07_delta-split1-prefire-vca,$X/l5extra_delta-split1-prefire-vca
  L5GROUPS=paired_l5,sjostrom07,paired_l5_extra
  L23DIRS=$X/ebner_l23l5_delta-split1-prefire-vseg-rs ANALYTICAL_BASIS_DIR=/scratch/dhuruva/split1/basis_l5l5
  L23_BASIS_DIR=/scratch/dhuruva/split1/basis_l23l5 JOINT=1 FITGAMMA=1
  "DROPT=letzkus_3ap_200hz_dt-10ms@distal|control,letzkus_3ap_200hz_dt-10ms@distal|nmdar_block"
  "EXTRA=l23l23:paired_l23l23,paired_l23l23_egger:$X/zilberter_l23l23_delta-split1-prefire-vseg-rs,$X/l23l23extra_delta-split1-prefire-vseg-rs:/scratch/dhuruva/split1/basis_l23l23")
B='"v5_mode": 2, "veto_T": 25.0, "ecb_ref": 2, "gate_src": 2, "gate_win": 100.0'
A='"gate_theta": 0.15, "theta_V": 3.5751446275588403'
G=(--account=def-emuller --gpus-per-node=nvidia_h100_80gb_hbm3_3g.40gb:1)
sub() { local j=$1 t=$2; shift 2; env "${COMMON[@]}" "$@" sbatch --parsable "${G[@]}" -J $j --mem=87G --time=$t $RS/run_joint_w2.sh; }
r=$(sub G8_G1a_r0 00:20:00 "SET={$B, $A, \"t_exact\": 0}" TAG=G8_G1a_r0 SEED=5 SEEDFITS=$RES/v7_G8_G1a_s5.json MAXITER=0 FREE=theta_eCB)
a=$(sub G8X_G1a_s5 00:35:00 "SET={$B, $A, \"t_exact\": 1}" TAG=G8X_G1a_s5 SEED=5 SEEDFITS=$RES/v7_G8_G1a_s5.json FREE=theta_eCB)
u=$(sub G8X_G1a_u6 00:35:00 "SET={$B, $A, \"t_exact\": 1}" TAG=G8X_G1a_u6 SEED=6 FREE=theta_eCB)
f=$(sub G8X_G1_s5 00:35:00 "SET={$B, \"t_exact\": 1}" TAG=G8X_G1_s5 SEED=5 SEEDFITS=$RES/v7_G8_G1_s5.json FREE=theta_V,theta_eCB)
echo "G1a_r0 $r G8X_G1a_s5 $a G8X_G1a_u6 $u G8X_G1_s5 $f"
