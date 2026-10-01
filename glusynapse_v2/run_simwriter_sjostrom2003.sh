#!/bin/bash
# T25 step 6: workdirs for configs/Sjostrom2003_L5TTPC_L5TTPC.yaml (11 protocols on the 120 Sabrina pairs,
# pairs_from the existing index). Post amplitudes reuse the cached single_cells pkls (same stimulus keys as
# the Ebner sjostrom_* protocols). Memory: the L2/3->L5 run peaked at 25.2 GB with 30 workers, incl. find_pairs.
#SBATCH --job-name=simwriter_sj03
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/simwriter_sj03_%j.out
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT/plastyfire && source $ROOT/.venv/bin/activate && export PYTHONPATH=$ROOT
python -u simwriter.py --config $ROOT/configs/Sjostrom2003_L5TTPC_L5TTPC.yaml
