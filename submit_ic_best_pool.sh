#!/bin/bash
# One Slurm job running all 700 IC-best BCL STDP sims in a 60-way process pool.
# Replaces submit_ic_best_700.sh (700 separate sbatch jobs, each waiting ~2h in
# the queue to do ~1m10s of work).
#
# Sizing from sacct on the 470 completed one-per-job runs:
#   elapsed  : ~1m10s typical, 16m55s max (cache-miss pairs do a live threshold search)
#   MaxRSS   : 4.0 GB peak per simulation
# => 60 workers * 4 GB = 240 GB. Narval standard nodes have 249 GB usable, so
#    this fills one node. Drop --workers/--mem together if you want a smaller
#    allocation (e.g. --workers 32 with --mem 130G).
#
# Walltime: 700 sims / 60 workers * ~1.2 min ~= 15 min, plus a margin for the
# slow cache-miss tail. 2h is generous; the job exits as soon as it finishes.
#
# Usage:  sbatch submit_ic_best_pool.sh
#         sbatch submit_ic_best_pool.sh --skip-existing

#SBATCH --job-name=icb_pool
#SBATCH --account=ctb-emuller
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=60
#SBATCH --mem=240G
#SBATCH --time=02:00:00
#SBATCH --chdir=/lustre06/project/6077694/dhuruva/plastyfire
#SBATCH --output=/lustre06/project/6077694/dhuruva/plastyfire/icb_pool_%j.log

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

python run_ic_best_pool.py --workers 60 "$@"
