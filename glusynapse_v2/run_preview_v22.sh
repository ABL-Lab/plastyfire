#!/bin/bash
#SBATCH --job-name=v22_prevfit
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=32
#SBATCH --mem=180G
#SBATCH --time=08:00:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/v22_prevfit_%j.out
# v2.2 preview: as v2.1 (pre_drive 1) plus the pre-rate threshold theta_Z on the NO pathway.
set -euo pipefail
source .venv/bin/activate
V2=glusynapse_v2
SIMS=refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations
PAIRS=180351-198084,181015-184976,182339-200396,182381-186941,182900-195902,183240-189246,183769-207759,186261-208490
PROTOS=$(python -c "import csv;print(','.join(dict.fromkeys(r['protocol_id'] for r in csv.DictReader(open('$SIMS/../../index_Ebner2019_L5TTPC_L5TTPC.csv')) if not r['protocol_id'].startswith('letzkus'))))")
python -u $V2/tests/test_units.py
python -u $V2/extract_v2.py --param-hash delta-prefire --protocols $PROTOS --pairs $PAIRS --window --workers 16 --skip-existing --out $V2/extracted/ebner_preview
python -u $V2/extract_v2.py --param-hash delta-cooker --pairs $PAIRS --workers 16 --skip-existing --out $V2/extracted/markram_delta-cooker
D=$V2/extracted/ebner_preview,$V2/extracted/markram_delta-cooker
echo "=== v2.2 fit (pre_drive 1) $(date)"
python -u $V2/fit_v2.py --dirs $D --pairs $PAIRS --model v2 --workers 32 --maxiter 80 \
    --filters '{"pre_drive": 1, "i_scale": 1e-5}' \
    --free-filters theta_Ti,tau_T,theta_NOi,tau_NO,tau_Z,theta_Z --filter-steps 6 \
    --x0 fit_results/delta-cooker.json --save $V2/results/preview_fit_v22
echo "=== done $(date)"
