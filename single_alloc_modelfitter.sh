#!/bin/sh
#SBATCH --job-name=plast_opt
#SBATCH --account=proj96
#SBATCH --partition=prod
#SBATCH --qos=longjob
#SBATCH --time=72:00:00
#SBATCH --nodes=64
#SBATCH --constraint=cpu
#SBATCH --cpus-per-task=2
#SBATCH --no-requeue
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --chdir=/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/refitting_resultsfitting/n100/seed19091997/
#SBATCH --output=opt-%j.log

# Set environment
source /home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/setupenv.sh
krenew -b -K 10
set -x
set -e
echo "TMPDIR:" $TMPDIR
export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=ipyparallel.${SLURM_JOBID}

echo "Launching controller"
ipcontroller --init --ip='*' --sqlitedb --ping=30000 --profile=${IPYTHON_PROFILE} &
sleep 1m

echo "Launching engines"
srun ipengine --timeout=500 --profile=${IPYTHON_PROFILE} &
sleep 5m

# Set next job
cp /home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/single_alloc_modelfitter.sh .
sbatch --dependency=afterany:${SLURM_JOBID} single_alloc_modelfitter.sh

# Run
python /home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/plastyfire/modelfitter.py --gen=100 --sample_size=100 --seed=19091997 --ipp_id=${SLURM_JOBID} -v

