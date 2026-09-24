#!/bin/bash
#SBATCH --job-name=fit_cpreonly_s1
#SBATCH --account=ctb-emuller
#SBATCH --cpus-per-task=60
#SBATCH --ntasks=1
#SBATCH --mem=26G
#SBATCH --time=01:00:00
#SBATCH --chdir=/lustre06/project/6077694/dhuruva/plastyfire/analytical_method
#SBATCH --output=/lustre06/project/6077694/dhuruva/plastyfire/analytical_method/fit_cpreonly_s1_%j.log

unset PATH PYTHONPATH LD_LIBRARY_PATH CMAKE_PREFIX_PATH
unset LIBRARY_PATH JPKGDIR NEURODAMUS_DIR HOC_LIBRARY_PATH
hash -r
module --force purge
module load StdEnv/2023 scipy-stack/2024a gcc/12.3 openmpi/4.1.5 hdf5-mpi/1.14.2 cmake/3.31.0 mpi4py/4.0.3 python/3.11.5
export JPKGDIR=/project/def-emuller/opt/jupyterhub-pkgs/build-07.04.2026-py311
export PYTHONPATH=/cvmfs/soft.computecanada.ca/easybuild/python/site-packages
export PYTHONPATH=$PYTHONPATH:/cvmfs/soft.computecanada.ca/custom/python/site-packages
export PYTHONPATH=$PYTHONPATH:$JPKGDIR/lib/python3.11/site-packages/
export PYTHONPATH=$PYTHONPATH:$JPKGDIR/lib/python
export PYTHONPATH=$PYTHONPATH:$HOME/.local/lib/python3.11/site-packages
export LD_LIBRARY_PATH=$JPKGDIR/lib
export PATH=$PATH:/opt/software/slurm/bin/:$JPKGDIR/bin
export PYTHONPATH=$PYTHONPATH:/lustre06/project/6077694/dhuruva/plastyfire

# stride 1 = the full 0.25 ms grid, no decimation. Coarser strides
# over-credit potentiation (a brief excursion above theta_p is credited for a
# whole sample interval), which biases the fit toward the in-vitro targets and
# makes the error look better than it is: 10 ms -> 4.47, 1 ms -> 5.68,
# 0.5 ms -> 8.22, with +10ms falling 1.126 -> 1.112 -> 1.103 as the grid
# refines. BCL measures 1.041. stride 1 removes that bias.
#
# Memory: effcai (4729 x 168001 float32) plus the transposed copy = 6.4 GB,
# shared across workers by fork copy-on-write. 96G leaves headroom for the
# transpose peak and the per-worker (N,) temporaries.
# Time: ~10.3 s/eval at 0.25 ms; popsize 24 = 96 individuals = 2 waves of 60
# -> ~21 s/gen. 6 h covers ~1000 generations.

echo "host: $(hostname)  cpus: $SLURM_CPUS_PER_TASK  mem: ${SLURM_MEM_PER_NODE}M"
echo "seed=${SEED:=3}  popsize=${POPSIZE:=24}  maxiter=${MAXITER:=400}  tol=${TOL:=0.001}"

time python fit.py --cpre-only --fit-gamma --de \
    --seed "$SEED" --popsize "$POPSIZE" --tol "$TOL" \
    --workers 60 --maxiter "$MAXITER" \
    --stride 1 --fit-dt=-10,5,10 --anchor-tails "$@"
