#!/bin/bash
#SBATCH --job-name=v2_postspk
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000M   # measured MaxRSS 804 MB
#SBATCH --time=03:00:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/v2_postspk_%A_%a.out
# Somatic spikes vs pulses in every induction protocol, all pairs (count_post_spikes.py); array over pair chunks.
set -euo pipefail
source .venv/bin/activate
SIMS=refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
PAIRS=$(ls $SIMS | awk -v i=$SLURM_ARRAY_TASK_ID -v n=${NCHUNK:-12} 'NR % n == i' | paste -sd,)
python -u glusynapse_v2/count_post_spikes.py --pairs $PAIRS ${PROTOS:+--protos $PROTOS} --out glusynapse_v2/results/post_spikes${OUTTAG:-}_chunk$SLURM_ARRAY_TASK_ID.csv
