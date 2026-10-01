#!/bin/bash
# T25 r50: workdirs for configs/Sjostrom2003b_L5TTPC_L5TTPC.yaml (2 burst protocols, nreps 50, 120 pairs; copy of run_simwriter_sjostrom2003.sh).
# pairs_from the existing index). Post amplitudes reuse the cached single_cells pkls (same stimulus keys as
# the Ebner sjostrom_* protocols). Memory: the L2/3->L5 run peaked at 25.2 GB with 30 workers, incl. find_pairs.
#SBATCH --job-name=simwriter_sj03r50
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:15:00
# seff 22068339: 39 s, 456 MB, 0.7% CPU eff on 32 cores -> 1 CPU / 1G / 0:15.
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/simwriter_sj03r50_%j.out
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT/plastyfire && source $ROOT/.venv/bin/activate && export PYTHONPATH=$ROOT
python -u simwriter.py --config $ROOT/configs/Sjostrom2003b_L5TTPC_L5TTPC.yaml
