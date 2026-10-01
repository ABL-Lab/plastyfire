#!/bin/bash
# MEASURED 22128322: 120/120 ok, max task 1:39, MaxRSS max 5.33 GB, CPU eff 74% -> next 7G, 0:15 (L2/3 posts are far cheaper than L5).
# EPSP basis, L2/3 PC -> L2/3 PC (Zilberter) on the new l23l23 edges, og-delta: one array task per pair (task i = i-th pair
# of the index csv), 12 CPUs = 12 concurrent trial procs, 5 trials (as submit_basis_edges.py, but logs in logs/ and no
# per-pair batch files). Sizing from basis_l23l5_rs 22086449-22086568 (120 pairs): MaxRSS max 11.6 GiB, max elapsed 14:49,
# CPU eff ~72%, ~90 CPU-h used -> 15G (+25%), 0:30 (+50% rounded up). Tasks beyond the pair count exit 0.
#SBATCH --job-name=basis_l23l23
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=12
#SBATCH --mem=15G
#SBATCH --time=00:30:00
#SBATCH --array=0-119
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/basis_l23l23_%A_%a.out
set -euo pipefail
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT && source .venv/bin/activate && export PYTHONPATH=$ROOT
S=$ROOT/refitting_results/fitting/n120/seed20262009
SIMS=$S/Zilberter2009_L23PC_L23PC/simulations
OUT=$ROOT/basis_results_edges_zilberter_l23l23_delta_rs
PAIR=$(tail -n +2 $S/index_Zilberter2009_L23PC_L23PC.csv | cut -d, -f1 | awk '!s[$0]++' | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")
[ -z "$PAIR" ] && { echo "no pair for task $SLURM_ARRAY_TASK_ID"; exit 0; }
PRE=${PAIR%-*}; POST=${PAIR#*-}; CSV=$OUT/basis_${PRE}_${POST}.csv
[ -f "$CSV" ] && { echo "exists $CSV"; exit 0; }
CFG=""
for d in $(ls $SIMS/$PAIR | sort); do
    [ -f $SIMS/$PAIR/$d/simulation_config.json ] && [ -f $SIMS/$PAIR/$d/prespikes.h5 ] && { CFG=$SIMS/$PAIR/$d/simulation_config.json; break; }
done
[ -z "$CFG" ] && { echo "no simulation_config.json + prespikes.h5 for $PAIR"; exit 1; }
mkdir -p $OUT
python -u run_basis_pair_edges.py --pre-gid $PRE --post-gid $POST --sim-config $CFG --output-csv $CSV \
    --num-trials 5 --workers 12 --circuit-config $ROOT/data/dhuruva_delta_l23l23_circuit_config.json
