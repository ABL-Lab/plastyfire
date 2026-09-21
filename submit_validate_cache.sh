#!/bin/bash
#SBATCH --job-name=det_test
#SBATCH --account=ctb-emuller
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=01:30:00
#SBATCH --chdir=/lustre06/project/6077694/dhuruva/plastyfire
#SBATCH --output=/lustre06/project/6077694/dhuruva/plastyfire/logs/det_test_%j.log
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
export PYTHONPATH=$PYTHONPATH:$JPKGDIR/lib/python3.11/site-packages/neurodamus/core/hoc
export PYTHONPATH=$PYTHONPATH:$HOME/.local/lib/python3.11/site-packages
export HOC_LIBRARY_PATH=$JPKGDIR/lib/python3.11/site-packages/neurodamus/data/hoc:$JPKGDIR/share/neurodamus_neocortex/hoc
export LD_LIBRARY_PATH=$JPKGDIR/lib
export PATH=$PATH:/opt/software/slurm/bin/:$JPKGDIR/bin
export PYTHONPATH=$PYTHONPATH:/lustre06/project/6077694/dhuruva/plastyfire
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

python validate_cpost_determinism.py \
    --cache cpre_cpost_cache/ion_channels_tau278.pkl \
    --circuit-config data/dhuruva_modified_ion_channels_circuit_config.json \
    --pairs 5 --seed 0 --workers 10
