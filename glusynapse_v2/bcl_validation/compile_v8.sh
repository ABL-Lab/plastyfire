#!/bin/bash
# Compile glusynapse_v2/mod/GluSynapseV8.mod alone into glusynapse_v2/mod_build_v8/x86_64 (dedicated directory, loaded
# next to the DEES library and mod_build_v7 with h.nrn_load_dll; DEES_cell_packages/x86_64, mod_build_v5 / v7 are not
# touched), then smoke_v8.py (single compartment: gate_src 0 == V7 bit for bit, shaft-licence timeline, licence
# replaces Vg). Copy of compile_v7.sh.
# Sizing basis: compile_v7 22166460 MEASURED 0:14, MaxRSS 213 MB -> 213 x 1.25 = 266 MB -> 320M; smoke_v8 runs 5 short
#   single-compartment runs (V7's 4) -> 0:14 x 1.5 -> 0:15 floor; 1 CPU (nrnivmodl + smoke are serial).
#SBATCH --job-name=compile_v8
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=320M
#SBATCH --time=00:15:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/compile_v8_%j.out
set -euo pipefail
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate
B=$ROOT/glusynapse_v2/mod_build_v8
mkdir -p $B/mod && cp $ROOT/glusynapse_v2/mod/GluSynapseV8.mod $B/mod/
cd $B && rm -rf x86_64 && nrnivmodl mod > nrnivmodl.log 2>&1 || { tail -40 nrnivmodl.log; exit 1; }
tail -3 nrnivmodl.log; ls -la x86_64/libnrnmech.so
cd $ROOT && python -u glusynapse_v2/bcl_validation/smoke_v8.py
