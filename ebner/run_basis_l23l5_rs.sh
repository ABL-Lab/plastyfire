#!/bin/bash
# T29: EPSP basis on the re-split L23PC->L5TTPC edges, og-delta, one sbatch per pair (12 CPUs = 12 concurrent trial procs).
# Sizing from 22066328-22066447 (same 120 pairs, old edges): max MaxRSS 12.25 GB, max elapsed 14:49, mean ~6 min,
# CPU eff 70-87% -> 15g (max+25%), 00:30:00 (max+50% rounded up). ~170 CPU-h allocated (logged in things_to_note.md).
# MEASURED: 120 jobs 22086449-22086568 all completed, MaxRSS max 11.6 GiB (of 15g), CPU eff ~72%.
set -euo pipefail
ROOT=/home/dhuruva/ws/plastyfire
cd $ROOT && source .venv/bin/activate
python -u submit_basis_edges.py --stdp-name Ebner2019_L23PC_L5TTPC --output-dir $ROOT/basis_results_edges_ebner_l23l5_delta_rs \
    --circuit-config $ROOT/data/dhuruva_delta_l23l5_circuit_config.json --mem 15g --walltime 00:30:00 --skip-existing
