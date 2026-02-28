#!/bin/bash

# Parameters
#SBATCH --array=0-99%20
#SBATCH --cpus-per-task=20
#SBATCH --error=/lustre06/project/6077694/dhuruva/plastyfire/basis_results/logs/%A_%a_0_log.err
#SBATCH --job-name=submitit
#SBATCH --mem=64GB
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --output=/lustre06/project/6077694/dhuruva/plastyfire/basis_results/logs/%A_%a_0_log.out
#SBATCH --partition=prod
#SBATCH --signal=USR2@90
#SBATCH --time=120
#SBATCH --wckey=submitit

# command
export SUBMITIT_EXECUTOR=slurm
srun --unbuffered --output /lustre06/project/6077694/dhuruva/plastyfire/basis_results/logs/%A_%a_%t_log.out --error /lustre06/project/6077694/dhuruva/plastyfire/basis_results/logs/%A_%a_%t_log.err /cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/python/3.11.5/bin/python -u -m submitit.core._submit /lustre06/project/6077694/dhuruva/plastyfire/basis_results/logs
