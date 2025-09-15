#!/usr/bin/env python3
"""
Recover all cached results that are missing from the CSV file
"""

import os
import pickle
import pandas as pd
import hashlib
from pathlib import Path

def recover_cached_results():
    """Recover all cached results and add them to the CSV"""
    
    # Paths
    cache_dir = Path("plastyfire/.cache")
    csv_path = Path("plastyfire/simulation_epsp_df.csv")
    
    # Load existing CSV
    if csv_path.exists():
        existing_df = pd.read_csv(csv_path)
        existing_hashes = set(existing_df['param_hash'].unique())
        print(f"Existing CSV has {len(existing_df)} entries with {len(existing_hashes)} unique param hashes")
    else:
        existing_df = pd.DataFrame()
        existing_hashes = set()
        print("No existing CSV found")
    
    # Process all cache files
    cache_files = list(cache_dir.glob("*.pkl"))
    print(f"Found {len(cache_files)} cache files")
    
    all_cached_data = []
    processed_hashes = set()
    
    for cache_file in cache_files:
        try:
            with open(cache_file, 'rb') as f:
                cache_data = pickle.load(f)
            
            # Extract param hash from filename
            param_hash = cache_file.stem[:12]  # First 12 chars
            
            if param_hash in processed_hashes:
                continue  # Skip duplicates
            processed_hashes.add(param_hash)
            
            if "resdb" in cache_data and isinstance(cache_data["resdb"], pd.DataFrame):
                res_db = cache_data["resdb"]
                
                # Create DataFrame with cached data
                for _, row in res_db.iterrows():
                    all_cached_data.append({
                        'pkl_file': f"cached_{param_hash}_{row['protocol_id']}.pkl",
                        'protocol_id': row['protocol_id'],
                        'param_hash': param_hash,
                        'epsp_ratio': row['epsp_ratio']
                    })
                
                print(f"Processed cache file {param_hash}: {len(res_db)} protocols")
            else:
                print(f"Cache file {param_hash} has no resdb or invalid format")
                
        except Exception as e:
            print(f"Error processing {cache_file}: {e}")
    
    # Create DataFrame from all cached data
    if all_cached_data:
        cached_df = pd.DataFrame(all_cached_data)
        print(f"\nTotal cached entries: {len(cached_df)}")
        print(f"Unique param hashes in cache: {cached_df['param_hash'].nunique()}")
        print(f"Unique protocols: {cached_df['protocol_id'].nunique()}")
        
        # Find missing entries
        cached_hashes = set(cached_df['param_hash'].unique())
        missing_hashes = cached_hashes - existing_hashes
        print(f"Missing param hashes: {len(missing_hashes)}")
        
        if missing_hashes:
            print("Missing hashes:", sorted(missing_hashes))
            
            # Add missing entries to CSV
            missing_df = cached_df[cached_df['param_hash'].isin(missing_hashes)]
            print(f"Adding {len(missing_df)} missing entries to CSV")
            
            # Append to existing CSV
            if existing_df.empty:
                missing_df.to_csv(csv_path, index=False)
            else:
                missing_df.to_csv(csv_path, mode='a', header=False, index=False)
            
            print(f"Updated CSV saved to {csv_path}")
        else:
            print("No missing entries found")
    else:
        print("No cached data found")

if __name__ == "__main__":
    recover_cached_results()
