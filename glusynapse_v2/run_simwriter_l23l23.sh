#!/bin/bash
# MEASURED 22128319: 18:58, MaxRSS 19.8 GB, CPU eff 59% on 16; 120 pairs, 1555 workdirs -> next 25G, 0:30.
# L2/3 PC -> L2/3 PC (Zilberter 2009): find_pairs (120 pairs, seed 20262009) + workdirs for configs/Zilberter2009_L23PC_L23PC.yaml.
# Sizing from simwriter_l23l5 22043887 (find + calibrate + write, 120 pairs x 5 protocols, 30 pool workers): 17:56,
# MaxRSS 34.7 GiB, TotalCPU 1:47 h on 32 CPUs (19% eff). Memory = 30 find_pairs workers (hard-coded) -> 34.7 + 25% = 44G.
# Work here: 13 protocols, ~46 calibration sims per post cell vs ~28 (x1.6), no cached cells (x1.5) -> ~43 min + 50% -> 1:15.
# CPUs 16 (measured ~6 busy on average at 32; 2.4x the CPU work = ~4.3 CPU-h still fits).
#SBATCH --job-name=simwriter_l23l23
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=16
#SBATCH --mem=44G
#SBATCH --time=01:15:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/simwriter_l23l23_%j.out
set -euo pipefail
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT/plastyfire && source $ROOT/.venv/bin/activate && export PYTHONPATH=$ROOT
python -u simwriter.py --config $ROOT/configs/Zilberter2009_L23PC_L23PC.yaml --workers 16
IDX=$ROOT/refitting_results/fitting/n120/seed20262009/index_Zilberter2009_L23PC_L23PC.csv
echo "index rows $(($(wc -l < $IDX) - 1)), pairs $(tail -n +2 $IDX | cut -d, -f1 | sort -u | wc -l)"
