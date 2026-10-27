#!/bin/bash
# EPSP basis (all-0, singletons, all-1 rho configs; C01 test pulses) for the 120 L23PC->L5TTPC pairs, og-delta:
# one sbatch per pair. Sizing: Sabrina delta basis jobs max MaxRSS 10.7 GB (12 workers) -> 14g; L23->L5 pairs
# have up to 24 synapses (L5-L5 runs took <= 22 min) -> 2 h.
set -euo pipefail
ROOT=/home/dhuruva/ws/plastyfire
cd $ROOT && source .venv/bin/activate
python -u submit_basis_edges.py --stdp-name Ebner2019_L23PC_L5TTPC --output-dir $ROOT/basis_results_edges_ebner_l23l5_delta \
    --circuit-config $ROOT/data/dhuruva_delta_circuit_config.json --mem 14g --walltime 02:00:00 --skip-existing
