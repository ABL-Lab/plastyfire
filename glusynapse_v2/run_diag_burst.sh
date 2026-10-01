#!/bin/bash
#SBATCH --job-name=v2_diagburst
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=780M   # measured MaxRSS 618 MB (syn1)
#SBATCH --time=02:00:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/lustre09/project/6070394/dhuruva/plastyfire/logs/v2_diagburst_%j.out
# Single-cell burst dissection (diag_burst.py). env: PAIR, VARIANTS, PROTOS.
set -euo pipefail
source .venv/bin/activate
python -u glusynapse_v2/diag_burst.py --pair $PAIR ${VARIANTS:+--variants $VARIANTS} ${PROTOS:+--protos $PROTOS} ${BDC:+--bd-cands $BDC} ${SYNC:+--syn-cands $SYNC} ${WIDTH:+--width $WIDTH} ${AMP:+--amp $AMP} \
    --out glusynapse_v2/results/diag_burst_${PAIR}${TAG:-}.json
