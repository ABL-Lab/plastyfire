#!/bin/bash
# Plot the c_pre-only STDP curve and, for comparison, DE fit #2 — both from the
# same 700 pair/protocol workdirs, same EPSP basis, same script.
#
# Whole-node allocation on purpose: plot_stdp_ic_best.py uses a bare
# ProcessPoolExecutor(), which sizes itself from os.cpu_count() (the physical
# node), not from --cpus-per-task. A smaller --mem would risk OOM.

#SBATCH --job-name=cpreonly_plot
#SBATCH --account=ctb-emuller
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=60
#SBATCH --mem=240G
#SBATCH --time=00:40:00
#SBATCH --chdir=/lustre06/project/6077694/dhuruva/plastyfire
#SBATCH --output=/lustre06/project/6077694/dhuruva/plastyfire/cpreonly_plot_%j.log

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

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "################ c_pre-only (bdbf06f915d0) ################"
python plot_stdp_ic_best.py --param-hash bdbf06f915d0 \
    --label "c_pre-only (4 param)" --out stdp_cpreonly_bcl.png

echo
echo "################ DE fit #2 (b8c7ff3ecf0a) ################"
python plot_stdp_ic_best.py --param-hash b8c7ff3ecf0a \
    --label "DE fit #2 (6 param)" --out stdp_defit2_bcl.png
