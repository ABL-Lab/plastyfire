#!/bin/bash
source /home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/setupenv.sh

python -c "
from plastyfire.simwriter import OptSimWriter

# 1. L5TTPC to L5TTPC
print('\n[1/4] Generating L5TTPC_L5TTPC protocols...')
writer = OptSimWriter('configs/L5TTPC_L5TTPC_STDP.yaml')
pairs = writer.find_pairs()
writer.write_sim_files(pairs)

# 2. L23PC to L5TTPC
print('\n[2/4] Generating L23PC_L5TTPC protocols...')
writer = OptSimWriter('configs/L23PC_L5TTPC_STDP.yaml')
pairs = writer.find_pairs()
writer.write_sim_files(pairs)

# 3. L23PC to L23PC
print('\n[3/4] Generating L23PC_L23PC protocols...')
writer = OptSimWriter('configs/L23PC_L23PC_STDP.yaml')
pairs = writer.find_pairs()
writer.write_sim_files(pairs)

# 4. L4PC to L23PC
print('\n[4/4] Generating L4PC_L23PC protocols...')
writer = OptSimWriter('configs/L4PC_L23PC_STDP.yaml')
pairs = writer.find_pairs()
writer.write_sim_files(pairs)

print('\nDone! Simulation folders generated successfully.')
"
