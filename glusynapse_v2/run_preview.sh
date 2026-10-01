#!/bin/bash
#SBATCH --job-name=v2_preview
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=16
#SBATCH --mem=180G
#SBATCH --time=03:00:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/v2_preview_%j.out
# Early look while ebner T5 runs: 8 pairs whose 26 non-Letzkus protocols are finished.
set -euo pipefail
source .venv/bin/activate
V2=glusynapse_v2; OUT=$V2/extracted/ebner_preview
SIMS=refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
PAIRS=180351-198084,181015-184976,182339-200396,182381-186941,182900-195902,183240-189246,183769-207759,186261-208490
PROTOS=$(python -c "import csv;print(','.join(dict.fromkeys(r['protocol_id'] for r in csv.DictReader(open('$SIMS/../../index_Ebner2019_L5TTPC_L5TTPC.csv')) if not r['protocol_id'].startswith('letzkus'))))")
python -u $V2/extract_v2.py --param-hash delta-prefire --protocols $PROTOS --pairs $PAIRS --window --workers 16 --skip-existing --out $OUT
python -u $V2/validate_v2.py --dirs $OUT --a fit_results/delta-cooker.json --pairs 8 --same-run
python -u $V2/preview_eval.py --dirs $OUT
