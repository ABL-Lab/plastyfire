#!/bin/bash
#SBATCH --job-name=validate_cpre_cpost
#SBATCH --account=ctb-emuller
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

cd /project/ctb-emuller/dhuruva/plastyfire
source setupenv.sh

echo "=== validate_cpre_cpost.py ==="
python validate_cpre_cpost.py
echo "Exit code: $?"
