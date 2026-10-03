#!/bin/bash
# Per-synapse Cpre / Cpre_APV / Cpre_Mg0 / Cpost (effcai + free cai) and VDCC/NMDA charge, delta-split1. See measure_cexp.py.
# Pilot: PATHNAME=L5 START=0 STOP=1 (one pair; run for L5, L23L5, L23L23).  Array: START=$((SLURM_ARRAY_TASK_ID*CHUNK)).
# Reference cpre/cpost = /scratch/dhuruva/split1/cache/*.pkl (defit2 globals, same code path as precompute_cpre_cpost.py).
# Sizing basis: split1 cache jobs 22136795 (24 L5 pairs, 8 CPU, 2:04, 5.6 GB), 22136803 (120 L23L5, 16 CPU, 4:44, 13.9 GB).
# MEASURED pilots (1 pair, 1 CPU, 8G req): L5 22155299 1:07, 588 MB (7 syn, 56 s); L23L5 22155300 1:03, 1.06 GB (3 syn, 52 s); L23L23 22155301 0:56, 1.31 GB (6 syn, 42 s).
# MEASURED arrays (CHUNK=12, 2G req): elapsed 1.5-8:44 per task (~40 s per pair), CPU eff ~97%; MaxRSS L5 0.45-1.94 GB, L23L5 up to 2.09 GB (AT the 2G limit, several tasks; the 6 failures were the no-single-AP post cells, not OOM), L23L23 0.69-1.72 GB -> next time 2.7G (2.09 x 1.25) and keep 0:30. Total ~1.5 CPU-h.
# MEASURED E1 ljp25 arrays (2700M req): 22177266 L5 8:40 1.58 GB; 22177268 L23L5 peak 8:34 2.76 GB (above 2.7G req: next 3.5G); 22177270 L23L23 4:48 1.43 GB.
# MEASURED A17 CEXP_DIAG=1 3-pair jobs (START=0 STOP=3, 3G req): 22184271 L5 E1 5:00, 647 MB, CPU eff 63%; 22184272 L5 E1+JS 2:39, 489 MB, 93%; 22184273 L23L23 E1 1:50, 1.53 GB, 77% -> 3-pair diag next time: 1 CPU, 2G, 0:15.
# MEASURED CEXP_SABATINI pilots (1 pair, 7 syn, 10 trials): 22199261 ljp0 3:48 934 MB; 22199262 ljp25 3:51 936 MB (~3.9x old per-pair) -> arrays CHUNK 12: old max 8:44 x3.9 = 34 min +50% -> 1:00, 2G (1.58 GB old peak +25%). Arrays 22199566 / 22199567, merges 22199568 / 22199569.
# MEASURED CEXP_SABATINI arrays (CHUNK 12, 2G req): 22199566 ljp0 37:09/36:35, 1.36/2.09 GB; 22199567 ljp25 38:13/42:13, 1.54/2.09 GB (AT the 2G limit) -> next time 2.7G 1:00.
# MEASURED T32 arrays (CHUNK 12, 2700M req, L5): 22218430 s25 37:34, 0.95 GB; 22218434 s12_g110 38:15, 1.99 GB (peak) -> L5-only next time 2.5G 1:00.
#SBATCH --job-name=cexp
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/cexp_%x_%A_%a.out
set -euo pipefail
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
export OMP_NUM_THREADS=1
mkdir -p ${CEXP_OUT:-/scratch/dhuruva/split1/cexp}
CHUNK=${CHUNK:-1}
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then S=$((SLURM_ARRAY_TASK_ID*CHUNK)); E=$((S+CHUNK)); else S=${START:-0}; E=${STOP:-1}; fi
if [ "${MERGE:-0}" = 1 ]; then python -u glusynapse_v2/cexp/measure_cexp.py --path $PATHNAME --merge
else python -u glusynapse_v2/cexp/measure_cexp.py --path $PATHNAME --start $S --stop $E; fi
