#!/bin/bash
# Prefire BCL validation of fit v5_S1_V5c_s6_39 (v5c, emodel delta-split1; 39 targets) with live GluSynapseV5 (continuous
# MEASURED pilot 22144694: 8:54, 5.24 GB (12 workers), 80% CPU; compare 22144695 0:20, 0.93 GB. Live = offline exactly (0/258 rho disagreements). FULL: og full 22139631 1:45:37, 21.99 GB at 32 workers; split1 pilot uses 1.26x og memory -> 32 CPU, 28G, 2:45.
# rule), as run_valid_V5c_s7.sh. Data / sims / circuits / bases are the split1 ones, set by the BCL_* env of live_v5.py
# (BCL_EMODEL tag -> extracted/*_delta-split1-prefire-*); mechanism = mod_build_v5 + the shared DEES lib (split1 channels).
# STAGE=pilot: L5 180351-198084 + L2/3 10149-186264 (--no-loc-filter), the V5c_s7 pilot pairs (both exist in split1).
# STAGE=full : all records; the pilot's results in the same jsonl are skipped. Compare: run_valid_S1_V5c_cmp.sh.
# Sizing basis: V5c_s7 pilot 22138636 MEASURED 9:05, MaxRSS 4.16 GB, CPU eff 86% (12 workers, 38 tasks)
#   -> 12 CPU, 6G (4.16 x 1.25 = 5.2), 0:15 (9:05 x 1.5 = 13:38). Full: resize from this pilot's seff (V5c_s7 full: 32 CPU, 22G).
#SBATCH --job-name=valid_S1_V5c
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=12
#SBATCH --mem=6G
#SBATCH --time=00:15:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/valid_S1_V5c_%j.out
set -euo pipefail
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
FIT=glusynapse_v2/rho_redesign/results/v5_S1_V5c_s6_39.json
OUT=${BCL_OUT:-/scratch/dhuruva/bcl_valid_S1_V5c/live.jsonl}   # STAGE=shard: BCL_OUT per shard, PL5 / PL23 pair lists
mkdir -p /scratch/dhuruva/bcl_valid_S1_V5c
SEL=""
if [ "${STAGE:-pilot}" = "shard" ]; then SEL="--pairs-l5 $PL5 --pairs-l23 $PL23"; export BCL_SHARD=1; fi
if [ "${STAGE:-pilot}" = "pilot" ]; then
    SEL="--pairs-l5 180351-198084 --pairs-l23 10149-186264 --no-loc-filter"
fi
python -u glusynapse_v2/bcl_validation/live_v5.py --fit $FIT --from-fit $SEL --out $OUT \
    --workers ${SLURM_CPUS_PER_TASK:-1} --timeout 3600
