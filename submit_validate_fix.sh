#!/bin/bash
#SBATCH --job-name=validate_fix
#SBATCH --account=ctb-emuller
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

cd /home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire
source setupenv.sh

echo "=== validate_fix.py ==="
python validate_fix.py
echo "Exit code: $?"
