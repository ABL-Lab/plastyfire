#!/bin/bash

# Define paths
BASE_RESULTS_DIR="/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/trace_results/DHURUVA_PARAMS_V5"
OUT_DIR="$BASE_RESULTS_DIR/figures"
PYTHON_SCRIPT="/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/plot_detailed_traces.py"

# Create output directory for figures
mkdir -p "$OUT_DIR"

echo "Using output directory: $OUT_DIR"

# Loop through all items in the base directory
for pair_dir in "$BASE_RESULTS_DIR"/*; do
    # Check if it is a directory
    if [ -d "$pair_dir" ]; then
        pair_name=$(basename "$pair_dir")
        
        # Skip the figures directory or logs
        if [ "$pair_name" == "figures" ] || [ "$pair_name" == "logs" ]; then
            continue
        fi
        
        echo "Processing pair: $pair_name"
        
        # Run the python script for this specific folder
        python3 "$PYTHON_SCRIPT" --pair-dir "$pair_dir" --out-base-dir "$OUT_DIR"
    fi
done

echo "Done! All figures have been generated in $OUT_DIR"
