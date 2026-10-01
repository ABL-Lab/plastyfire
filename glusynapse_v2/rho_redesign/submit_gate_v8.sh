#!/bin/bash
# A14 gate_src options (gpu_v8_rho.py, GATE_IMPL.md) on the W2_N3 65-target joint env (copy of submit_joint_w2_n2.sh COMMON).
# Sizing: repro = measured W2 rescores 6:45-6:58, 50-62 GB -> 78G 0:15. gate_src 2 adds the shaft trace on the host during
# each build (L2/3->L2/3 1356M steps x 4 B = 5.4 GB + uint16 copy 2.7 GB) -> (61.6 + 8) x 1.25 = 87G; GPU +6.9 GB uint16
# (31 -> 38 GB of the 3g.40gb slice). Fits: measured 17-19 min + ~3 min host calibration (_open_q) + slower src-2 kernel -> 0:40
# (src 2), 0:35 (src 1), 0:30 (src 3). Refits wait afterok on the G1 rescore (gate_src != 0 exits 0; only R0 checks REPRO OK).
set -euo pipefail
R=/lustre09/project/6070394/dhuruva/plastyfire; V2=$R/glusynapse_v2; X=$V2/extracted; RS=$V2/rho_redesign; RES=$RS/results
COMMON=(KERNEL=gpu_v8_rho.py
  L5DIRS=$X/ebner_delta-split1-prefire-vca,$X/markram_delta-split1-prefire-vca,$X/sj03_delta-split1-prefire-vca,$X/sj03r50_delta-split1-prefire-vca,$X/sj07_delta-split1-prefire-vca,$X/l5extra_delta-split1-prefire-vca
  L5GROUPS=paired_l5,sjostrom07,paired_l5_extra
  L23DIRS=$X/ebner_l23l5_delta-split1-prefire-vseg-rs ANALYTICAL_BASIS_DIR=/scratch/dhuruva/split1/basis_l5l5
  L23_BASIS_DIR=/scratch/dhuruva/split1/basis_l23l5 JOINT=1 FITGAMMA=1
  "DROPT=letzkus_3ap_200hz_dt-10ms@distal|control,letzkus_3ap_200hz_dt-10ms@distal|nmdar_block"
  "EXTRA=l23l23:paired_l23l23,paired_l23l23_egger:$X/zilberter_l23l23_delta-split1-prefire-vseg-rs,$X/l23l23extra_delta-split1-prefire-vseg-rs:/scratch/dhuruva/split1/basis_l23l23")
B='"v5_mode": 2, "veto_T": 25.0, "ecb_ref": 2'
N3=$RES/v7_W2_N3.json; TV='"theta_V": 3.5751446275588403'
G=(--account=def-emuller --gpus-per-node=nvidia_h100_80gb_hbm3_3g.40gb:1)
sub() { local j=$1 t=$2 m=$3; shift 3; env "${COMMON[@]}" "$@" sbatch --parsable "${G[@]}" ${DEP:-} -J $j --mem=$m --time=$t $RS/run_joint_w2.sh; }
F2=FREE=theta_V,theta_eCB
r0=$(sub G8_R0 00:15:00 78G "SET={$B, \"gate_src\": 0}" TAG=G8_R0 SEED=6 SEEDFITS=$N3 MAXITER=0 $F2)
g1r=$(sub G8_G1r015 00:20:00 87G "SET={$B, \"gate_src\": 2, \"gate_win\": 100.0, \"gate_theta\": 0.15}" TAG=G8_G1r015 SEED=6 SEEDFITS=$N3 MAXITER=0 $F2)
g5r=$(sub G8_G5r050 00:20:00 87G "SET={$B, \"gate_src\": 2, \"gate_win\": 100.0, \"gate_theta\": 0.5}" TAG=G8_G5r050 SEED=6 SEEDFITS=$N3 MAXITER=0 $F2)
DEP="--dependency=afterok:$g1r"
ids=""
for spec in "G1:2:100.0:00:40:00:87G" "G2:2:0.0:00:40:00:87G" "S1:1:100.0:00:35:00:78G" "C3:3:0.0:00:30:00:78G"; do
  IFS=: read -r nm src win hh mm ss mem <<< "$spec"
  S="SET={$B, \"gate_src\": $src, \"gate_win\": $win}"
  a=$(sub G8_${nm}_s5 $hh:$mm:$ss $mem "$S" TAG=G8_${nm}_s5 SEED=5 SEEDFITS=$N3 $F2)
  b=$(sub G8_${nm}_u6 $hh:$mm:$ss $mem "$S" TAG=G8_${nm}_u6 SEED=6 $F2)
  ids="$ids $nm:$a,$b"
done
ga=$(sub G8_G1a_s5 00:40:00 87G "SET={$B, \"gate_src\": 2, \"gate_win\": 100.0, \"gate_theta\": 0.15, $TV}" TAG=G8_G1a_s5 SEED=5 SEEDFITS=$N3 FREE=theta_eCB)
echo "R0 $r0 G1r015 $g1r G5r050 $g5r |$ids G1a:$ga"
