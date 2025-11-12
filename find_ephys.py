#!/usr/bin/env python3
"""Find ephys file with specific EPSP value"""
import os
import pickle
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'plastyfire'))
from plastyfire.evaluator import compute_epsp_prefire

# Target EPSP values
target_before = 1.5654242228629773
target_after = 17.960802743181993

ephys_dir = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/refitting_results/fitting/n100/seed19091997/ephys_data"

# Get all ephys files
ephys_files = [f for f in os.listdir(ephys_dir) if f.endswith('.pkl')]

print(f"Searching through {len(ephys_files)} ephys files...")

matches = []
for filename in ephys_files:
    filepath = os.path.join(ephys_dir, filename)
    try:
        epsp = compute_epsp_prefire(filepath)

        # Check if this matches either the before or after value
        if abs(epsp - target_before) < 0.0001:
            print(f"FOUND 'before' match: {filename}")
            print(f"  EPSP value: {epsp}")
            matches.append(('before', filename, filepath))
        elif abs(epsp - target_after) < 0.0001:
            print(f"FOUND 'after' match: {filename}")
            print(f"  EPSP value: {epsp}")
            matches.append(('after', filename, filepath))
    except Exception as e:
        # Skip files that cause errors
        pass

print(f"\nTotal matches found: {len(matches)}")

# Print details
for match_type, filename, filepath in matches:
    print(f"\n{match_type.upper()} file: {filename}")
    print(f"Full path: {filepath}")

    # Extract gids and rho from filename
    parts = filename.replace('ephys_data_', '').replace('.pkl', '').split('_')
    print(f"Pre GID: {parts[0]}")
    print(f"Post GID: {parts[1]}")
    print(f"Rho config: {parts[2:]}")
