#!/bin/bash
# Compile glusynapse_v2/mod/GluSynapseV4.mod alone into glusynapse_v2/mod_build_v4/x86_64 (loaded next to the DEES
# library with h.nrn_load_dll; no duplicate mechanism names), then a 200 ms single-compartment smoke test.
# Sizing: one-mod nrnivmodl + a one-compartment NEURON run; minimum request 1 CPU, 2G, 0:15.
# MEASURED 22131695: 0:07, 224 MB, smoke OK -> 512M (224 MB + 25%, rounded up), 0:15 floor.
#SBATCH --job-name=compile_v4
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=512M
#SBATCH --time=00:15:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/compile_v4_%j.out
set -euo pipefail
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate
B=$ROOT/glusynapse_v2/mod_build_v4
mkdir -p $B/mod && cp $ROOT/glusynapse_v2/mod/GluSynapseV4.mod $B/mod/
cd $B && rm -rf x86_64 && nrnivmodl mod > nrnivmodl.log 2>&1 || { tail -40 nrnivmodl.log; exit 1; }
tail -3 nrnivmodl.log; ls -la x86_64/libnrnmech.so
cd $ROOT && python -u glusynapse_v2/bcl_validation/smoke_v4.py
