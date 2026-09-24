#!/bin/bash
#SBATCH --job-name=analytical_extract
#SBATCH --account=ctb-emuller
#SBATCH --cpus-per-task=60
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --chdir=/lustre06/project/6077694/dhuruva/plastyfire/analytical_method
#SBATCH --output=/lustre06/project/6077694/dhuruva/plastyfire/analytical_method/extract_%j.log

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

echo "host: $(hostname)  cpus: $SLURM_CPUS_PER_TASK  mem: ${SLURM_MEM_PER_NODE}M"
time python extract.py --workers 60
echo "--- output ---"
ls extracted | wc -l
du -sh extracted
