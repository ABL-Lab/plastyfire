#!/bin/bash
# Prefire BCL validation of a v8 fit (gpu_v8_rho via fit_v6: G1 shaft gate gate_src 2, or gate_src 0 = v7) with live
# GluSynapseV8 (live_v8.py, mod_build_v8 from compile_v8.sh); split1 data / sims / circuits and stages as run_valid_v7.sh.
# Parameters (env): FIT = results json, ECB_REF = 0 | 2 (must equal the fit's), TAG = output name
#   (/scratch/dhuruva/bcl_valid_v8_${TAG}/), STAGE = smoke | pilot | full | shard (PL5 / PL23 / PL2323, BCL_OUT as v7).
# STAGE=smoke (G1 pilot): SMOKE_TASKS (json list; default L5 192879-186028, 1 basal + 5 apical synapses at 447-724 um
#   where the shaft gate decides, x {sjostrom_50hz_dt+10ms, sjostrom_20hz_dt+10ms} control, Sjostrom 2001 LTP
#   protocols), then compare_prefire_v8.py --show-syn inline (live vs the CPU port of the gpu_v8 G1 kernel, per-synapse
#   rest / licence openings / open time / rho); exit 1 if a task failed.
# STAGE=pilot: --from-fit on L5 192879-186028, L2/3->L5 10149-186264 (Letzkus distal, sh_distal), L2/3->L2/3 10212-10960.
# Other stages: compare with run_valid_v8_cmp.sh.
# Sizing basis (smoke, the #SBATCH defaults): V7 smoke 22166461 MEASURED 2:12, 1.24 GB, 76% CPU (2 workers, 2 L5 Sj01
#   tasks + compare) -> same work (2 L5 tasks; compare loads 2 records with one more trace) -> 1536M, 2 CPU, 0:15.
#SBATCH --job-name=valid_v8
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=2
#SBATCH --mem=1536M
#SBATCH --time=00:15:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/valid_v8_%j.out
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
D=/scratch/dhuruva/bcl_valid_v8_${TAG}; mkdir -p $D
STAGE=${STAGE:-smoke}
OUT=${BCL_OUT:-$D/live_${STAGE}.jsonl}
LIVE="python -u glusynapse_v2/bcl_validation/live_v8.py --fit $FIT --out $OUT --workers ${SLURM_CPUS_PER_TASK:-1} --timeout 3600"
case $STAGE in
  smoke)
    T=$D/tasks_smoke.json
    DEF='[{"path": "L5", "pair": "192879-186028", "proto": "sjostrom_50hz_dt+10ms", "cond": "control", "phase": "prefire", "express": "cooker"}, {"path": "L5", "pair": "192879-186028", "proto": "sjostrom_20hz_dt+10ms", "cond": "control", "phase": "prefire", "express": "cooker"}]'
    if [ -n "${SMOKE_TASKS:-}" ]; then echo "$SMOKE_TASKS" > $T; else echo "$DEF" > $T; fi
    $LIVE --tasks $T
    python -u glusynapse_v2/bcl_validation/compare_prefire_v8.py --fit $FIT --results $OUT --save $D/cmp_smoke --show-syn
    jq -e -s 'all(.ok)' $OUT > /dev/null || { echo "smoke: failed live tasks"; exit 1; } ;;
  pilot) export BCL_PAIRS_L23L23=${PL2323:-10212-10960}; $LIVE --from-fit --pairs-l5 192879-186028 --pairs-l23 10149-186264 --no-loc-filter ;;
  shard) export BCL_SHARD=1 BCL_PAIRS_L23L23=$PL2323; $LIVE --from-fit --pairs-l5 $PL5 --pairs-l23 $PL23 ;;
  full)  $LIVE --from-fit ;;
  *) echo "unknown STAGE $STAGE"; exit 1 ;;
esac
