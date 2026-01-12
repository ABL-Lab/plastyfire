#!/bin/bash
#SBATCH --job-name=plasty_opt
#SBATCH --account=ctb-emuller
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/plasty_opt_%j.log
#SBATCH --mail-type=ALL

# Production script to run plastyfire optimization using SLURM for distributed execution
# This uses the fixed SLURM pipeline with proper job waiting

# Ensure we are in the script directory
cd "$(dirname "$0")"

# Create logs directory if it doesn't exist
mkdir -p logs

# Ensure environment is loaded
source /home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/setupenv.sh

python modelfitter.py \
    --use-slurm \
    --gen=100 \
    --pop_size=50 \
    --sample_size=20 \
    --seed=19091998 \
    --max-jobs=800 \
    --fitness-schedule=fitness_schedule_mrk97_08.yaml \
    -v
