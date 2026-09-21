#!/bin/bash
# Step 1 for the b78bbff arm: c_pre/c_post on the pre-ec0d47c emodel.
#
# b78bbff (2025-10-17) is the emodel the Jun 30 STDP figure ran on — apical
# gNaTgbar_NaTg = 0.04 (not the halved 0.019856), geom_nseg_fixed(20) (not 40),
# apical Ca 0.0025 (not 0.003333). Everything downstream of the emodel is held at
# today's values, so this isolates the ec0d47c conductance refit.
#
# Output: cpre_cpost_cache/b78bbff_tau278.pkl   (new name; nothing overwritten)

#SBATCH --job-name=b78_cpre
#SBATCH --account=ctb-emuller
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=192G
#SBATCH --time=06:00:00
#SBATCH --chdir=/lustre06/project/6077694/dhuruva/plastyfire
#SBATCH --output=/lustre06/project/6077694/dhuruva/plastyfire/logs/b78_cpre_%j.log

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

python precompute_cpre_cpost.py \
    --params chindemi_aparams \
    --results-dir /project/ctb-emuller/dhuruva/plastyfire/refitting_results \
    --circuit-config data/dhuruva_b78bbff_circuit_config.json \
    --output cpre_cpost_cache/b78bbff_tau278.pkl \
    --workers 48
