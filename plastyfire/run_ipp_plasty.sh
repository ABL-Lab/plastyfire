#!/bin/bash
#SBATCH --job-name=plasty_ipp
#SBATCH --account=ctb-emuller
#SBATCH --time=24:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=63
#SBATCH --mem=150G
#SBATCH --output=logs/plasty_ipp_%j.log
#SBATCH --mail-type=ALL

# IPyParallel execution script for plastyfire
# This creates a persistent cluster within the Slurm allocation

# Ensure we are in the script directory
# Ensure we are in the correct directory
cd /home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/plastyfire

# Create logs directory if it doesn't exist
mkdir -p logs

# Ensure environment is loaded
source /home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/setupenv.sh

# Profile configuration
PROFILE=job_${SLURM_JOB_ID}
export IPP_PROFILE_DIR=${HOME}/.ipython/profile_${PROFILE}

echo "Creating IPP profile ${PROFILE}..."
ipython profile create --parallel --profile=${PROFILE}

echo "Starting IPController..."
ipcontroller --ip="*" --profile=${PROFILE} &
sleep 10

echo "Starting IPEngines..."
srun --account=ctb-emuller ipengine --profile=${PROFILE} &
sleep 30

echo "Starting Optimization..."
python modelfitter.py \
    --use-ipp \
    --gen=100 \
    --pop_size=50 \
    --sample_size=20 \
    --seed=19091998 \
    --recipe-path=/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/biodata/recipe_mod.csv \
    --profile=${PROFILE} \
    -v

echo "Optimization finished. Cleaning up..."
# Cleanup is handled by Slurm killing the job, but we can be polite
pkill -P $$
rm -rf ${IPP_PROFILE_DIR}
