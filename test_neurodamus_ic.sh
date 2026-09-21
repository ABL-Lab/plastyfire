#!/bin/bash
# Interactive test for neurodamus ion-channels run on ONE pair/protocol.
# Run this directly on a compute node to diagnose failures before submitting 700 jobs.
#
# Usage (on compute node):
#   bash test_neurodamus_ic.sh
#   bash test_neurodamus_ic.sh 180164-197248 10Hz_10ms

PAIR="${1:-180164-197248}"
PROTO="${2:-10Hz_10ms}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="$SCRIPT_DIR/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations/$PAIR/$PROTO"
SIM_CFG="$WORKDIR/prefire_simulation_config_ic.json"

if [[ ! -f "$SIM_CFG" ]]; then
    echo "ERROR: $SIM_CFG not found. Run submit_neurodamus_ic.py --dry-run first."
    exit 1
fi

DEES_SPECIAL="/project/ctb-emuller/dhuruva/DEES_cell_packages/x86_64/special"
DEES_LIB="/project/ctb-emuller/dhuruva/DEES_cell_packages/x86_64"

# Try JPKGDIR candidates in order — DEES special was compiled against 20.03.2025
JPKGDIR_CANDIDATES=(
    "/project/def-emuller/opt/jupyterhub-pkgs/build-20.03.2025-py311"
    "/project/def-emuller/opt/jupyterhub-pkgs/build-07.04.2026-py311"
    "/project/def-emuller/opt/jupyterhub-pkgs/build-12.02.2026-py311"
)

JPKGDIR=""
for candidate in "${JPKGDIR_CANDIDATES[@]}"; do
    if [[ -f "$candidate/bin/neurodamus_init.py" ]]; then
        JPKGDIR="$candidate"
        echo "Found JPKGDIR: $JPKGDIR"
        break
    fi
done

if [[ -z "$JPKGDIR" ]]; then
    echo "ERROR: No valid JPKGDIR found."
    exit 1
fi

echo "=== Test setup ==="
echo "Workdir  : $WORKDIR"
echo "Config   : $SIM_CFG"
echo "Special  : $DEES_SPECIAL"
echo "JPKGDIR  : $JPKGDIR"
echo ""

# ---------- environment ----------
unset PATH PYTHONPATH LD_LIBRARY_PATH CMAKE_PREFIX_PATH
unset LIBRARY_PATH JPKGDIR_ENV NEURODAMUS_DIR HOC_LIBRARY_PATH
hash -r

module --force purge
module load StdEnv/2023 scipy-stack/2025a gcc/12.3 openmpi/4.1.5 hdf5-mpi/1.14.4 \
    cmake/3.31.0 cuda/12.9 mpi4py/4.0.3 pytest/8.2.2 rust/1.91.0 boost/1.85.0
module load python/3.11.5

export JPKGDIR="$JPKGDIR"
export PYTHONPATH=/cvmfs/soft.computecanada.ca/easybuild/python/site-packages:/cvmfs/soft.computecanada.ca/custom/python/site-packages
export PYTHONPATH=$PYTHONPATH:$JPKGDIR/lib/python3.11/site-packages/
export PATH=$PATH:/opt/software/slurm/bin/:$JPKGDIR/bin
export PYTHONPATH=$PYTHONPATH:$JPKGDIR/lib/python
export PYTHONPATH=$PYTHONPATH:$JPKGDIR/lib/python3.11/site-packages/neurodamus/core/hoc
export HOC_LIBRARY_PATH=$JPKGDIR/lib/python3.11/site-packages/neurodamus/data/hoc:$JPKGDIR/share/neurodamus_neocortex/hoc
# DEES special must come first so its libnrnmech.so (with cad/Ka_kampa/kBK) is found
# JPKGDIR/lib must provide the libnrniv.so that matches the build of DEES special
export LD_LIBRARY_PATH=$DEES_LIB:$JPKGDIR/lib

echo "=== libnrniv.so being used ==="
ldd "$DEES_SPECIAL" 2>/dev/null | grep nrn
echo ""
echo "=== NEURON version ==="
python -c "import neuron; print('neuron', neuron.__version__)" 2>/dev/null || echo "(neuron import failed)"
echo ""

echo "=== Running neurodamus (single rank, no srun) ==="
cd "$WORKDIR"
"$DEES_SPECIAL" -mpi \
    -python "$JPKGDIR/bin/neurodamus_init.py" \
    --configFile="$SIM_CFG" --lb-mode=WholeCell --verbose 2>&1 | head -80

echo ""
echo "Exit code: $?"
