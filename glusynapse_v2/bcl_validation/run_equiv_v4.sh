#!/bin/bash
# Equivalence test of the live GluSynapseV4 rule vs the offline v4 kernel, fit v4_C1Ajd_s5 (tasks_equiv.json: L5 pair
# 180351-198084 x sjostrom_50hz_dt-10ms, 10Hz_-10ms, sjostrom07_step200ms_pair (+ mglu_block, no_block, nmdar_block);
# L2/3 proximal pair 7471-187420 x letzkus_3ap_200hz_dt+10ms, sjostrom_50hz_dt+10ms). 22 runs: 7 prefire (cooker
# expression), 15 full (cooker or v4 expression), then compare_v4.py --mode equiv.
# Sizing: MEASURED 22131696 (11 workers, 22 tasks): 10:31, MaxRSS 4.66 GB, CPU eff 84%, max task wall 412 s (full).
#   22132496 FAILED in 0:10 (2.7 GB): float32 dt_ms 0.2499962 failed the bin assert -> compare with round(dt,4) as batch_v2 does.
#   Rerun (2026-10-01) with GluSynapseV4 bin mode: a self-event every 0.25 ms per synapse during the induction (CVODE
#   restarts), task walls not yet measured, assume <= 3x. 15 workers = all 15 full runs in one wave, the 7 prefire runs
#   in the freed slots: 4.66 GB x 15/11 x 1.25 = 8G, + ~0.1 GB/worker for the extra CVODE steps recorded -> 10G; ~3 x 7 min + compare 3 min ~ 25 min, +50% -> 0:45; 15 CPUs,
#   ~11 CPU-h.
#SBATCH --job-name=equiv_v4
#SBATCH --account=rrg-emuller
# MEASURED 22132674 (bin mode, 15 workers): TIMEOUT at 0:45, 7.53 GB MaxRSS, CPU 98%; 5/22 done, task walls 863-2340 s (~5.7x the
#   continuous path). Resume (skips ok tasks): 17 left -> 17 workers one wave, ~2400-3000 s + compare -> 1:30; 7.53 x 17/15 x 1.25 -> 11G; ~25 CPU-h.
# MEASURED 22134590 (resume, 17 tasks, 17 workers): 1:23:50, 7.76 GB MaxRSS, CPU 73%. Result 20/22 pass.
#SBATCH --cpus-per-task=17
#SBATCH --mem=11G
#SBATCH --time=01:30:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/equiv_v4_%j.out
set -euo pipefail
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1
FIT=glusynapse_v2/rho_redesign/results/v4_C1Ajd_s5.json
OUT=/scratch/dhuruva/bcl_validation/equiv_C1Ajd_s5.jsonl
python -u glusynapse_v2/bcl_validation/live_v4.py --fit $FIT --tasks glusynapse_v2/bcl_validation/tasks_equiv.json \
    --out $OUT --workers 17
python -u glusynapse_v2/bcl_validation/compare_v4.py --fit $FIT --results $OUT --mode equiv \
    --save glusynapse_v2/rho_redesign/results/bcl_equiv_C1Ajd_s5
