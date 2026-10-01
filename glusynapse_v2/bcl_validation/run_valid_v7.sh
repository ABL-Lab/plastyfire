#!/bin/bash
# Prefire BCL validation of a v7 fit (gpu_v7_rho via fit_v6: veto_T, ecb_ref 0 / 2) with live GluSynapseV7
# (live_v7.py, mod_build_v7 from compile_v7.sh), split1 data / sims / circuits as run_valid_S1_V5c.sh, all pathways of
# the fit: L5 (+ l5extra), L2/3->L5 (--joint), L2/3->L2/3 (--extra l23l23: Zilberter + Egger; Banerjee validation only).
# Parameters (env): FIT = results json (relative to plastyfire/ or absolute), ECB_REF = 0 | 2 (must equal the fit's),
#   TAG = output name (/scratch/dhuruva/bcl_valid_v7_${TAG}/), STAGE = smoke | pilot | full | shard.
#   shard: PL5 / PL23 / PL2323 pair lists ("none" = no pair of that pathway), BCL_OUT = $D/shard_<k>.jsonl
#   (submit_valid_v7.sh builds them). Done tasks in BCL_OUT are skipped, so a timed-out shard is resubmitted as is.
# STAGE=smoke: SMOKE_TASKS (json list; default L5 180351-198084 x {sjostrom_20hz_dt+10ms, sjostrom_40hz_dt+10ms}
#   control, the +10 trains where the veto acts), then compare_prefire_v7.py --show-syn inline; exit 1 if a task failed.
# Other stages: compare with run_valid_v7_cmp.sh.
# Sizing basis (smoke, the #SBATCH defaults): smoke 22166461 MEASURED 2:12, 1.24 GB, 76% CPU (2 workers, 2 L5 tasks +
#   compare) -> 1536M (1.24 x 1.25 = 1.55), 2 CPU, 0:15.
# Shards (sbatch overrides, submit_valid_v7.sh): V5c shards 22147501-07 MEASURED 13.1-15.7 GB at 32 workers (0.49 GB
#   per worker on the peak shard), 4:20-9:33, 80% CPU -> 32 CPU, 20G (15.7 x 1.25); time from the per-task wall of
#   S1_V5c merged.jsonl (L5 mean 171 s, L2/3->L5 141 s), see submit_valid_v7.sh.
#SBATCH --job-name=valid_v7
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=2
#SBATCH --mem=1536M
#SBATCH --time=00:15:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/valid_v7_%j.out
set -euo pipefail
: "${FIT:?FIT=<results json>}" "${ECB_REF:?ECB_REF=0|2}" "${TAG:?TAG=<name>}"
export ECB_REF
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1
S1=/scratch/dhuruva/split1
SR=$S1/refitting_results/fitting/n120/seed20262009
export BCL_EMODEL=delta-split1
export BCL_L5_SIMS=$SR/Sabrina_L5TTPC_L5TTPC_STDP/simulations
export BCL_L5_CIRCUIT=$ROOT/data/dhuruva_delta-split1_circuit_config.json
export BCL_L23_SIMS=$SR/Ebner2019_L23PC_L5TTPC/simulations
export BCL_L23_CIRCUIT=$ROOT/data/dhuruva_split1_l23l5_circuit_config.json
export BCL_L23L23_SIMS=$SR/Zilberter2009_L23PC_L23PC/simulations:$SR/L23L23extra/simulations
export BCL_L23L23_CIRCUIT=$ROOT/data/dhuruva_split1_l23l23_circuit_config.json
export ANALYTICAL_BASIS_DIR=$S1/basis_l5l5 L23_BASIS_DIR=$S1/basis_l23l5 L23L23_BASIS_DIR=$S1/basis_l23l23 CEXP_DIR=$S1/cexp
D=/scratch/dhuruva/bcl_valid_v7_${TAG}; mkdir -p $D
STAGE=${STAGE:-smoke}
OUT=${BCL_OUT:-$D/live_${STAGE}.jsonl}
LIVE="python -u glusynapse_v2/bcl_validation/live_v7.py --fit $FIT --out $OUT --workers ${SLURM_CPUS_PER_TASK:-1} --timeout 3600"
case $STAGE in
  smoke)
    T=$D/tasks_smoke.json
    DEF='[{"path": "L5", "pair": "180351-198084", "proto": "sjostrom_20hz_dt+10ms", "cond": "control", "phase": "prefire", "express": "cooker"}, {"path": "L5", "pair": "180351-198084", "proto": "sjostrom_40hz_dt+10ms", "cond": "control", "phase": "prefire", "express": "cooker"}]'
    if [ -n "${SMOKE_TASKS:-}" ]; then echo "$SMOKE_TASKS" > $T; else echo "$DEF" > $T; fi
    $LIVE --tasks $T
    python -u glusynapse_v2/bcl_validation/compare_prefire_v7.py --fit $FIT --results $OUT --save $D/cmp_smoke --show-syn
    jq -e -s 'all(.ok)' $OUT > /dev/null || { echo "smoke: failed live tasks"; exit 1; } ;;
  pilot) export BCL_PAIRS_L23L23=${PL2323:-10212-10960}; $LIVE --from-fit --pairs-l5 180351-198084 --pairs-l23 10149-186264 --no-loc-filter ;;
  shard) export BCL_SHARD=1 BCL_PAIRS_L23L23=$PL2323; $LIVE --from-fit --pairs-l5 $PL5 --pairs-l23 $PL23 ;;
  full)  $LIVE --from-fit ;;
  *) echo "unknown STAGE $STAGE"; exit 1 ;;
esac
