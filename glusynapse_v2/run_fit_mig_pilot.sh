#!/bin/bash
# Pilot: can the pl5sj07 t_drive 4 fit run on a smaller MIG slice (1g.10gb / 2g.20gb)? The 3g.40gb queue had 681 pending / 16 running
# (2026-09-30). MAXITER=2 only, to measure GPU memory peak and objective-call time. Traces are streamed (CHUNK GB) so they fit the slice.
# Usage: SLICE=1g.10gb CHUNK=4 sbatch --gpus-per-node=nvidia_h100_80gb_hbm3_$SLICE:1 run_fit_mig_pilot.sh
# Host memory 39G is as run_fits_td4_pl5sj07.sh (td3 fits measured 24-27 GB + sj07).
#SBATCH --job-name=fit_mig_pilot
#SBATCH --account=def-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=39G
#SBATCH --time=00:30:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/fit_mig_pilot_%j.out
set -euo pipefail
V2=glusynapse_v2; X=$V2/extracted
export XLA_PYTHON_CLIENT_PREALLOCATE=false   # so the sampled memory is the real peak
( while true; do nvidia-smi --query-compute-apps=used_memory --format=csv,noheader 2>/dev/null | head -1; sleep 5; done ) > logs/fit_mig_pilot_${SLURM_JOB_ID}_gpumem.txt &
SMI=$!
env MODE=subset TGROUPS=paired_l5,sjostrom07 CONDS=control,mglu_block,post_nmdar,nmdar_block,no_block MAXITER=2 SEED=1 \
  TAG=_migpilot_${SLURM_JOB_ID} CHUNK=${CHUNK:-4} \
  DIRS=$X/ebner_delta-prefire-vca,$X/markram_delta-prefire-vca,$X/sj03_delta-prefire-vca,$X/sj03r50_delta-prefire-vca,$X/sj07_delta-prefire-vca \
  FREE=theta_Te,tau_T,theta_Tg,dpre_min,theta_NOi,tau_NO,tau_Z,theta_Z \
  SEEDFITS=$V2/results/reduced_gpu_subset_pl5r50_td3_s5.json \
  FILTERS='{"pre_drive": 1, "i_scale": 1e-5, "t_drive": 4, "tau_E1": 100.0}' \
  bash $V2/run_reduced_gpu.sh
kill $SMI || true
echo "GPU mem peak (MiB): $(sort -n logs/fit_mig_pilot_${SLURM_JOB_ID}_gpumem.txt | tr -dc '0-9\n' | sort -n | tail -1)"
