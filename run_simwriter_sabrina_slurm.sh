#!/bin/bash
#SBATCH --job-name=simwriter_sabrina
#SBATCH --account=rrg-emuller
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=60
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/simwriter_sabrina_%j.out
#SBATCH --error=/project/rrg-emuller/dhuruva/plastyfire/logs/simwriter_sabrina_%j.err

cd /project/rrg-emuller/dhuruva/plastyfire
source /project/rrg-emuller/dhuruva/plastyfire/.venv/bin/activate

python -c "
from plastyfire.simwriter import OptSimWriter

print('Generating Sabrina_L5TTPC_L5TTPC_STDP protocols (npairs=120)...')
writer = OptSimWriter('configs/Sabrina_L5TTPC_L5TTPC_STDP.yaml')
pairs = writer.find_pairs()
writer.write_sim_files(pairs)
print('Done! Simulation folders generated successfully.')
"
