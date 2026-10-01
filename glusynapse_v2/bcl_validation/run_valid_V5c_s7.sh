#!/bin/bash
# Prefire BCL validation of fit v5_V5c_s7 (v5c, v5_mode 2): live GluSynapseV5 (continuous rule) on every record behind
# the fit's targets (live_v5.py --from-fit: prefire runs, cooker expression = the fit's own traces), then
# run_valid_V5c_s7_cmp.sh (afterok) compares per record and per target. BCL_V5c_s7.md.
# STAGE=pilot: L5 180351-198084 (30 records) + L2/3 10149-186264 (distal, 8 records, --no-loc-filter) = 38 tasks,
#   the design of run_valid_C1Ajn_s5.sh STAGE=pilot.
# STAGE=full : all 1295 records (716 L5 + 579 L2/3); the pilot's results in the same jsonl are skipped.
# Sizing basis: v4 pilot 22135321 (same 38 tasks, 12 workers) MEASURED 8:34, MaxRSS 4.97 GB, CPU eff 84% ->
#   12 CPU, 7G (4.97 x 1.25 = 6.2), 0:15 (8:34 x 1.5 = 12:51). GluSynapseV5 adds one state (W) to V4.
# STAGE=full: resize from this pilot's seff (v4 full plan: 32 workers, 2:30, 17G, ~80 CPU-h requested).
# MEASURED pilot 22138636: 9:05, MaxRSS 4.16 GB, CPU eff 86% (12 workers). Pilot verdict: live vs offline same records chi2 410.85 vs 406.84, rho disagreements 0/258; one L5 eCB edge event (10Hz_10ms, dpre_maxdiff 0.29).
# FULL plan: v4 full 22135837 at 1:18 had 736/1295 records and MaxRSS 17.0 GB at 32 workers (its 17G limit) -> 32 CPU, 22G (17.0 x 1.25), 3:30 (~2:20 x 1.5).
# MEASURED full 22139631: 1:45:37, MaxRSS 21.99 GB of 22G, CPU eff 93% (32 workers, 1257 new records).
#SBATCH --job-name=valid_V5c_s7
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=12
#SBATCH --mem=7G
#SBATCH --time=00:15:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/valid_V5c_s7_%j.out
set -euo pipefail
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1
FIT=glusynapse_v2/rho_redesign/results/v5_V5c_s7.json
OUT=/scratch/dhuruva/bcl_valid_V5c_s7/live.jsonl
mkdir -p /scratch/dhuruva/bcl_valid_V5c_s7
SEL=""
if [ "${STAGE:-pilot}" = "pilot" ]; then
    SEL="--pairs-l5 180351-198084 --pairs-l23 10149-186264 --no-loc-filter"
fi
python -u glusynapse_v2/bcl_validation/live_v5.py --fit $FIT --from-fit $SEL --out $OUT \
    --workers ${SLURM_CPUS_PER_TASK:-1} --timeout 3600
