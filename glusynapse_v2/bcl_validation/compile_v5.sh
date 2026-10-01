#!/bin/bash
# Compile glusynapse_v2/mod/GluSynapseV5.mod alone into glusynapse_v2/mod_build_v5/x86_64 (dedicated directory, loaded
# next to the DEES library with h.nrn_load_dll; the shared DEES_cell_packages/x86_64 is not touched), then smoke_v5.py.
# Sizing: as compile_v4.sh, MEASURED 22131695: 0:07, 224 MB -> 512M, 0:15 floor.
#SBATCH --job-name=compile_v5
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=512M
#SBATCH --time=00:15:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/compile_v5_%j.out
set -euo pipefail
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate
B=$ROOT/glusynapse_v2/mod_build_v5
mkdir -p $B/mod && cp $ROOT/glusynapse_v2/mod/GluSynapseV5.mod $B/mod/
cd $B && rm -rf x86_64 && nrnivmodl mod > nrnivmodl.log 2>&1 || { tail -40 nrnivmodl.log; exit 1; }
tail -3 nrnivmodl.log; ls -la x86_64/libnrnmech.so
cd $ROOT && python -u glusynapse_v2/bcl_validation/smoke_v5.py
