#!/bin/bash
# Compile glusynapse_v2/mod/GluSynapseV7.mod alone into glusynapse_v2/mod_build_v7/x86_64 (dedicated directory, loaded
# next to the DEES library with h.nrn_load_dll; the shared DEES_cell_packages/x86_64 and mod_build_v5 are not touched),
# then smoke_v7.py (single-compartment veto / uE checks). Copy of compile_v5.sh.
# Sizing basis: compile_v5 22138635 MEASURED 0:11, MaxRSS 215 MB of 512M (42%) -> 215 x 1.25 = 270 MB -> 320M;
#   0:11 x 1.5 -> 0:15 floor; 1 CPU (nrnivmodl + smoke are serial).
#SBATCH --job-name=compile_v7
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=320M
#SBATCH --time=00:15:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/compile_v7_%j.out
set -euo pipefail
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate
B=$ROOT/glusynapse_v2/mod_build_v7
mkdir -p $B/mod && cp $ROOT/glusynapse_v2/mod/GluSynapseV7.mod $B/mod/
cd $B && rm -rf x86_64 && nrnivmodl mod > nrnivmodl.log 2>&1 || { tail -40 nrnivmodl.log; exit 1; }
tail -3 nrnivmodl.log; ls -la x86_64/libnrnmech.so
cd $ROOT && python -u glusynapse_v2/bcl_validation/smoke_v7.py
