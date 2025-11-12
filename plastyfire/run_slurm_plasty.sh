#!/bin/bash
# Production script to run plastyfire optimization using SLURM for distributed execution
# This uses the fixed SLURM pipeline with proper job waiting

python modelfitter.py \
    --use-slurm \
    --gen=100 \
    --pop_size=50 \
    --sample_size=10 \
    --seed=19091998 \
    --max-jobs=500 \
    -v
