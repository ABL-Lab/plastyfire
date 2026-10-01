#!/bin/bash
# Full prefire BCL validation of a 3-pathway v8 joint fit (G1 shaft gate, live GluSynapseV8; login node: sbatch + shell
# only). Copy of submit_valid_v7.sh with run_valid_v8.sh / run_valid_v8_cmp.sh. Env: FIT ECB_REF TAG.
#   1 gate smoke (STAGE=smoke, L2/3->L2/3 pair 10212-10960: Zilberter train10 50 Hz +4 + Egger) -> all shards afterok
#   -> compare afterok (per-pathway chi2 live / offline / fit and shaft-licence openings / open time live vs offline).
#   Shards (STAGE=shard, 32 workers): L5 3 x 8 pairs, L2/3->L5 2 x ~56 GEOM all_protocols pairs, L2/3->L2/3 8 x 15 pairs.
# Sizing basis (MEASURED v7 W2_N3, submit_valid_v7.sh header; +~10% for the shaft trace and the gate events):
#   shards L5 24-29 min 15.4 GB -> 20G 0:45; L23L5 20-22 min 18.5 GB -> 23G 0:35; L23L23 13-16 min 12.3 GB -> 16G 0:30
#   (32 CPU, 84-95% eff). Requested 3 x 32 x 0.75 + 2 x 32 x 0.583 + 8 x 32 x 0.5 = 237 CPU.h (used ~150).
#   Compare: v7 W2 full compare MEASURED 12:26, 15.75 GB -> 20G 0:30 (v8 keeps one more trace per record).
#   MEASURED G8X_G1a_u6: gate 22195075 3:37 1.33 GB; shards L5 22:19 14.3 GB 93%, L23L5 22:07 20.2 GB 92%, L23L23 15:30 12.0 GB 81%; compare 22195089 OOM at 20G; 22197386 at 40G: 27:58, 21.8 GB -> next 28G 0:45.
#   Gate smoke: v8 smoke 22193267 MEASURED 2:15, 1.44 GB (2 CPU) -> 1800M 0:15.
set -euo pipefail
: "${FIT:?}" "${ECB_REF:?}" "${TAG:?}"
R=/project/rrg-emuller/dhuruva/plastyfire; BV=$R/glusynapse_v2/bcl_validation; X=$R/glusynapse_v2/extracted
D=/scratch/dhuruva/bcl_valid_v8_${TAG}; mkdir -p $D
F=$FIT; [ "${F:0:1}" = / ] || F=$R/$F
split() { local n=$1; shift; local a=("$@") k i; for ((k = 0; k < n; k++)); do local s=(); for ((i = k; i < ${#a[@]}; i += n)); do s+=("${a[i]}"); done; (IFS=,; echo "${s[*]}"); done; }
mapfile -t P5 < <(jq -r '.args.pairs' $F | tr , '\n')
mapfile -t P23 < <(awk -F, 'NR > 1 && $3 == "True" {print $1 "-" $2}' $R/ebner/pair_geometry_L23PC_L5TTPC.csv | sort -u)
mapfile -t P2323 < <(ls $X/zilberter_l23l23_delta-split1-prefire-vseg-rs $X/l23l23extra_delta-split1-prefire-vseg-rs | grep npz | sed 's/__.*//' | sort -u)
echo "pairs: L5 ${#P5[@]}, L23 ${#P23[@]}, L23L23 ${#P2323[@]}"
E="ALL,FIT=$FIT,ECB_REF=$ECB_REF,TAG=$TAG"
ST='[{"path": "L23L23", "pair": "10212-10960", "proto": "zilberter_train10_50hz_dt+4ms_last", "cond": "control", "phase": "prefire", "express": "cooker"}, {"path": "L23L23", "pair": "10212-10960", "proto": "egger1999_5ap_20hz_dt+10ms", "cond": "control", "phase": "prefire", "express": "cooker"}]'
g=$(SMOKE_TASKS="$ST" sbatch --parsable -J valid_v8_gate --mem=1800M --time=00:15:00 --export=$E,STAGE=smoke,SMOKE_TASKS $BV/run_valid_v8.sh)
ids=(); k=0
sub() { ids+=($(PL5=$1 PL23=$2 PL2323=$3 sbatch --parsable -J valid_v8_sh$k --dependency=afterok:$g --cpus-per-task=32 --mem=$5 \
    --time=$4 --export=$E,STAGE=shard,PL5,PL23,PL2323,BCL_OUT=$D/shard_$k.jsonl $BV/run_valid_v8.sh)); k=$((k + 1)); }
for s in $(split 3 "${P5[@]}"); do sub "$s" none none 00:45:00 20G; done
for s in $(split 2 "${P23[@]}"); do sub none "$s" none 00:35:00 23G; done
for s in $(split 8 "${P2323[@]}"); do sub none none "$s" 00:30:00 16G; done
dep=$(IFS=:; echo "${ids[*]}")
c=$(sbatch --parsable --dependency=afterok:$dep --mem=28G --time=00:45:00 --export=$E,STAGE=full $BV/run_valid_v8_cmp.sh)
echo "gate $g shards ${ids[*]} compare $c"
