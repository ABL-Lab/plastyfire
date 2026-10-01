#!/bin/bash
# MEASURED 22128318: 0:03, MaxRSS 43 MB; 13 ids, 15 targets, missing none.
# Parse check: targets.py paired_l23l23 + Zilberter yaml protocols (ids in PAIRED_L23L23 must exist in the yaml).
# Sizing from simwriter_l23l5_nulls 22123227 (python start-up + yaml): 20 s, MaxRSS 315 MB -> 1 CPU, 1G, 0:15.
#SBATCH --job-name=check_l23l23
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:15:00
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/check_l23l23_%j.out
set -euo pipefail
ROOT=/project/rrg-emuller/dhuruva/plastyfire
cd $ROOT/glusynapse_v2 && source $ROOT/.venv/bin/activate && export PYTHONPATH=$ROOT:$ROOT/plastyfire
python -u targets.py paired_l23l23
python -u - <<'EOF'
import yaml, targets
ids = {p["id"] for p in yaml.safe_load(open("/project/rrg-emuller/dhuruva/plastyfire/configs/Zilberter2009_L23PC_L23PC.yaml"))["protocols"]}
miss = {k for k, _ in targets.PAIRED_L23L23} - ids
print(len(ids), "yaml ids;", len(targets.PAIRED_L23L23), "targets; missing ids:", miss or "none")
assert not miss
EOF
