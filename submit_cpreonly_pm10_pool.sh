#!/bin/bash
# One Slurm job running all 700 c_pre-only BCL STDP sims in a 60-way process pool.
#
# c_pre-only refit scoring ONLY dt=-10 and +10 (fit.py --fit-dt=-10,10, no tail
# anchors). Offline weighted err 1.736 over those 2 points. +5 and the tails
# are UNSCORED predictions in this fit, not fitted values.
#
# theta_d = a00*c_pre, theta_p = a10*c_pre, with the c_post
# coefficients pinned to zero. Four free params (a00, a10, gamma_d, gamma_p)
# instead of six. Params + provenance live in run_de_fit2_pool.py::CPRE_ONLY_PM10_ARGS.
#
#   a00 = 1.014763  a10 = 1.809655  gamma_d = 161.8697  gamma_p = 281.6517
#   param hash 21f3bec952eb  ->  simulation_edges_21f3bec952eb.pkl per workdir
#
# Does NOT touch DE fit #2 (b8c7ff3ecf0a) or the first c_pre-only run (bdbf06f915d0).
#
# Sizing is the measured DE-fit-2 envelope, unchanged:
#   elapsed  : ~1m10s typical, 16m55s max (cache-miss pairs do a live threshold search)
#   MaxRSS   : 4.0 GB peak per simulation  => 60 * 4 GB = 240 GB (one Narval node)
#   walltime : 700 / 60 * ~1.2 min ~= 15 min; 2h leaves room for the slow tail.
#
# bluecellulab computes theta at runtime from the c_pre/c_post cache + the
# a-params, so edges.h5 theta_d/theta_p are NOT read on this path and no
# threshold injection is needed here.
#
# Usage:  sbatch submit_cpreonly_pm10_pool.sh
#         sbatch submit_cpreonly_pm10_pool.sh --skip-existing

#SBATCH --job-name=cpreonly_pm10_pool
#SBATCH --account=ctb-emuller
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=60
#SBATCH --mem=240G
#SBATCH --time=02:00:00
#SBATCH --chdir=/lustre06/project/6077694/dhuruva/plastyfire
#SBATCH --output=/lustre06/project/6077694/dhuruva/plastyfire/cpreonly_pm10_pool_%j.log

unset PATH PYTHONPATH LD_LIBRARY_PATH CMAKE_PREFIX_PATH
unset LIBRARY_PATH JPKGDIR NEURODAMUS_DIR HOC_LIBRARY_PATH
hash -r

module --force purge
module load StdEnv/2023 scipy-stack/2024a gcc/12.3 openmpi/4.1.5 \
    hdf5-mpi/1.14.2 cmake/3.31.0 mpi4py/4.0.3 python/3.11.5

export JPKGDIR=/project/def-emuller/opt/jupyterhub-pkgs/build-07.04.2026-py311
export PYTHONPATH=/cvmfs/soft.computecanada.ca/easybuild/python/site-packages
export PYTHONPATH=$PYTHONPATH:/cvmfs/soft.computecanada.ca/custom/python/site-packages
export PYTHONPATH=$PYTHONPATH:$JPKGDIR/lib/python3.11/site-packages/
export PYTHONPATH=$PYTHONPATH:$JPKGDIR/lib/python
export PYTHONPATH=$PYTHONPATH:$JPKGDIR/lib/python3.11/site-packages/neurodamus/core/hoc
export PYTHONPATH=$PYTHONPATH:$HOME/.local/lib/python3.11/site-packages
export HOC_LIBRARY_PATH=$JPKGDIR/lib/python3.11/site-packages/neurodamus/data/hoc:$JPKGDIR/share/neurodamus_neocortex/hoc
export LD_LIBRARY_PATH=$JPKGDIR/lib
export PATH=$PATH:/opt/software/slurm/bin/:$JPKGDIR/bin
export PYTHONPATH=$PYTHONPATH:/lustre06/project/6077694/dhuruva/plastyfire

# NEURON/numpy must not each grab all 60 cores — the pool provides the parallelism.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

python run_de_fit2_pool.py --cpre-only-pm10 --workers 60 "$@"
