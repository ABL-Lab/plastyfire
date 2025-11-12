#!/bin/bash
# Test script to run plastyfire optimization using SLURM for distributed execution
# This uses the fixed SLURM pipeline with proper job waiting

python modelfitter.py \
    --use-slurm \
    --gen=4 \
    --pop_size=4 \
    --sample_size=4 \
    --seed=165 \
    --max-jobs=50 \
    -v
