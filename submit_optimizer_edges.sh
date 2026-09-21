#!/bin/bash
#SBATCH --job-name=edges_optimizer
#SBATCH --account=ctb-emuller
#SBATCH --cpus-per-task=64
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

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

cd /project/ctb-emuller/dhuruva/plastyfire

echo "=== GA optimizer: sample_size=30, pop_size=16, gen=100, seeded with chindemi ==="
python plastyfire/modelfitter_edges.py \
    --sample_size 30  \
    --pop_size    16  \
    --gen         100 \
    --max-workers 62  \
    --cpre-cpost-cache /project/ctb-emuller/dhuruva/plastyfire/cpre_cpost_cache/chindemi_tau278.pkl \
    --basis-dir   basis_results_edges_mini/ \
    --seed        1234 \
    --seed-individual "101.5,216.2,1.002,1.954,1.159,2.483,1.127,2.456,5.236,1.782" \
    --log-file    optimizer_edges_n30_seeded \
    -v
echo "Exit code: $?"
