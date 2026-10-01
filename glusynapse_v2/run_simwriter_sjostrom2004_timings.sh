#!/bin/bash
# Re-run of run_simwriter_sjostrom2004.sh for the 4 added Fig 4 timing ids (existing id is skipped/unchanged by design: same step key).
# MEASURED 22133660 (4 new ids): 0:13, 331 MB MaxRSS, 4% CPU (all step keys cached) -> oversized; the pre+100c rerun needs a new key (run_simwriter_sj04_pre100c.sh).
# MEASURED 22127910 (1 protocol): 1:39, MaxRSS 1.43 GB, 81% CPU on 8. The step-threshold search is shared by the 4 ids (same step key),
# so workload ~ same; keep 2G, 0:15, 8 CPU = workers.
#SBATCH --job-name=simwriter_sj04_tim
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=8
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/simwriter_sj04_tim_%j.out
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT/plastyfire && source $ROOT/.venv/bin/activate && export PYTHONPATH=$ROOT
python -u simwriter.py --config $ROOT/configs/Sjostrom2004_L5TTPC_L5TTPC.yaml --workers 8
