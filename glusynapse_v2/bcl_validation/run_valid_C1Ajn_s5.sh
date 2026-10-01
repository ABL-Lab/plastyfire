#!/bin/bash
# Prefire BCL validation of fit v4_C1Ajn_s5: live GluSynapseV4 (continuous rule, V4_BIN=0) on every record behind the
# MEASURED full 22135837 (32 workers, 1295 records): 1:58:16, MaxRSS 16.99 GB of 17G (at the limit; next time 22G), CPU eff 89%.
# fit's targets (live_v4.py --from-fit: prefire runs, cooker expression = the fit's own traces), then
# run_valid_C1Ajn_s5_cmp.sh (afterok) compares per record and per target. BCL_C1Ajn_s5.md.
# STAGE=pilot: L5 180351-198084 (30 records) + L2/3 10149-186264 (distal, 8 records, --no-loc-filter) = 38 tasks.
# STAGE=full : all 1295 records (716 L5 + 579 L2/3); the pilot's results in the same jsonl are skipped.
# Sizing basis (pilot), equiv 22131696 MEASURED (non-bin, 11 workers): MaxRSS 4.66 GB (0.42 GB/worker), task walls
#   prefire 10Hz 48 s, sj50 98 s, sj07 222-232 s, L2/3 sj50 88 s, letzkus 3AP 167 s; 0.1 Hz (502 s sim) unmeasured,
#   ~380 s by s/sim-s. Pilot ~7000 s of tasks on 12 workers -> makespan ~10-12 min, +50% -> 0:30.
#   Mem 12 x 0.5 GB (bin-run 7.53/15) x 1.25 -> 8G. ~2 CPU-h used, 6 requested.
# MEASURED pilot 22135321 (12 workers, 38 tasks): 8:34, MaxRSS 4.97 GB (0.41 GB/worker), CPU eff 84%, 1:42:48 core;
#   walls: 10Hz 53 s, 0.1 Hz 196 s, r50 burst 297 s (max 301), 152 s-sim 97 s, sj07 pair 265 / pre 156 / post 134 s,
#   L2/3 sj50 77, 1ap 96, nopost 63, 3AP 135-148 s. Full: 1295 records ~46 CPU-h of task wall, ~44.7 left after the pilot.
#   STAGE=full: 32 workers -> 44.7/32/0.84 = 1.66 h, +50% -> 2:30 (> max task 301 s); mem 32 x 0.41 x 1.25 + 0.5 -> 17G.
#   sbatch --cpus-per-task=32 --mem=17G --time=02:30:00 --export=ALL,STAGE=full run_valid_C1Ajn_s5.sh (~80 CPU-h requested).
#SBATCH --job-name=valid_C1Ajn
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=12
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/valid_C1Ajn_%j.out
set -euo pipefail
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 V4_BIN=0
FIT=glusynapse_v2/rho_redesign/results/v4_C1Ajn_s5.json
OUT=/scratch/dhuruva/bcl_valid_C1Ajn_s5/live.jsonl
mkdir -p /scratch/dhuruva/bcl_valid_C1Ajn_s5
SEL=""
if [ "${STAGE:-pilot}" = "pilot" ]; then
    SEL="--pairs-l5 180351-198084 --pairs-l23 10149-186264 --no-loc-filter"
fi
python -u glusynapse_v2/bcl_validation/live_v4.py --fit $FIT --from-fit $SEL --out $OUT \
    --workers ${SLURM_CPUS_PER_TASK:-1} --timeout 3600
