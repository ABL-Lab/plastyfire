#!/bin/bash
# Full prefire BCL validation of a 3-pathway v7 joint fit (login node: sbatch + shell only). Env: FIT ECB_REF TAG.
#   1 gate smoke (run_valid_v7.sh STAGE=smoke, L2/3->L2/3 pair 10212-10960: Zilberter train10 50 Hz +4 + Egger, the
#     new pathway and the ecb_ref 2 cexp path; exits 1 if a task fails) -> all shards afterok -> compare afterok.
#   Shards (run_valid_v7.sh STAGE=shard, 32 workers): L5 3 x 8 pairs, L2/3->L5 2 x ~56 GEOM all_protocols pairs,
#   L2/3->L2/3 8 x 15 pairs.
# Sizing basis (per-task wall from S1_V5c merged.jsonl: L5 mean 171 s (max 655), L2/3->L5 141 s; L2/3->L2/3 not
#   measured: 171 s x 202/152 s prefire tstop = 227 s, an upper estimate for the smaller L2/3 cell; tasks L5 729,
#   L2/3->L5 486, L2/3->L2/3 ~2000 incl. mglu_block and Banerjee): work 34.6 + 19.0 + 126 = 180 task-h, at the V5c
#   shards' 80% CPU efficiency ~225 CPU.h. Memory: V5c shard peak 15.7 GB at 32 workers x 1.25 -> 20G.
#   Time = per-shard work / 32 / 0.8 x 1.5: L5 27 min -> 0:45, L2/3->L5 22 min -> 0:45, L2/3->L2/3 37 min -> 1:00.
#   A timed-out shard resumes from its jsonl when resubmitted. Gate smoke: 22166461 (2:12, 1.24 GB) -> 2 CPU 1536M 0:15.
set -euo pipefail
: "${FIT:?}" "${ECB_REF:?}" "${TAG:?}"
R=/project/rrg-emuller/dhuruva/plastyfire; BV=$R/glusynapse_v2/bcl_validation; X=$R/glusynapse_v2/extracted
D=/scratch/dhuruva/bcl_valid_v7_${TAG}; mkdir -p $D
F=$FIT; [ "${F:0:1}" = / ] || F=$R/$F
split() { local n=$1; shift; local a=("$@") k i; for ((k = 0; k < n; k++)); do local s=(); for ((i = k; i < ${#a[@]}; i += n)); do s+=("${a[i]}"); done; (IFS=,; echo "${s[*]}"); done; }
mapfile -t P5 < <(jq -r '.args.pairs' $F | tr , '\n')
mapfile -t P23 < <(awk -F, 'NR > 1 && $3 == "True" {print $1 "-" $2}' $R/ebner/pair_geometry_L23PC_L5TTPC.csv | sort -u)
mapfile -t P2323 < <(ls $X/zilberter_l23l23_delta-split1-prefire-vseg-rs $X/l23l23extra_delta-split1-prefire-vseg-rs | grep npz | sed 's/__.*//' | sort -u)
echo "pairs: L5 ${#P5[@]}, L23 ${#P23[@]}, L23L23 ${#P2323[@]}"
E="ALL,FIT=$FIT,ECB_REF=$ECB_REF,TAG=$TAG"
ST='[{"path": "L23L23", "pair": "10212-10960", "proto": "zilberter_train10_50hz_dt+4ms_last", "cond": "control", "phase": "prefire", "express": "cooker"}, {"path": "L23L23", "pair": "10212-10960", "proto": "egger1999_5ap_20hz_dt+10ms", "cond": "control", "phase": "prefire", "express": "cooker"}]'
g=$(SMOKE_TASKS="$ST" sbatch --parsable -J valid_v7_gate --export=$E,STAGE=smoke,SMOKE_TASKS $BV/run_valid_v7.sh)
ids=(); k=0
sub() { ids+=($(PL5=$1 PL23=$2 PL2323=$3 sbatch --parsable -J valid_v7_sh$k --dependency=afterok:$g --cpus-per-task=32 --mem=20G \
    --time=$4 --export=$E,STAGE=shard,PL5,PL23,PL2323,BCL_OUT=$D/shard_$k.jsonl $BV/run_valid_v7.sh)); k=$((k + 1)); }
for s in $(split 3 "${P5[@]}"); do sub "$s" none none 00:45:00; done
for s in $(split 2 "${P23[@]}"); do sub none "$s" none 00:45:00; done
for s in $(split 8 "${P2323[@]}"); do sub none none "$s" 01:00:00; done
dep=$(IFS=:; echo "${ids[*]}")
c=$(sbatch --parsable --dependency=afterok:$dep --mem=38G --time=00:15:00 --export=$E,STAGE=full $BV/run_valid_v7_cmp.sh)
# MEASURED W2_N3: shards L5 24-29 min 14.5-15.4 GB, L23L5 20-22 min 17.8-18.5 GB, L23L23 13-16 min 11.5-12.3 GB (32 CPU, 84-95% eff). Compare 22174649 TIMED OUT at 0:15 (15.7 GB at kill; 36% CPU) -> resubmitted 25G 1:00.
echo "gate $g shards ${ids[*]} compare $c"
