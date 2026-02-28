#!/bin/bash
#SBATCH --job-name=getenv
#SBATCH --account=ctb-emuller
#SBATCH --time=00:05:00
#SBATCH --mem=1G
#SBATCH --output=env.txt
#SBATCH --error=env.err

env
which python
python -c "import sys; print(sys.executable); import numpy; print(numpy.__file__)"
