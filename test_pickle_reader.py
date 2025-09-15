#!/usr/bin/env python3
"""
Simple test script to verify pickle files can be read correctly
"""
import pickle
import os
import glob

def test_pickle_files(data_dir="ephys_data"):
    """Test reading pickle files from the ephys_data directory"""
    
    # Check if directory exists
    if not os.path.exists(data_dir):
        print(f"Directory {data_dir} does not exist yet. Run the main script first.")
        return
    
    # Find all pickle files
    pickle_files = glob.glob(os.path.join(data_dir, "*.pkl"))
    
    if not pickle_files:
        print(f"No pickle files found in {data_dir}")
        return
    
    print(f"Found {len(pickle_files)} pickle files:")
    
    for pkl_file in pickle_files:
        print(f"\nTesting file: {pkl_file}")
        try:
            with open(pkl_file, 'rb') as f:
                data = pickle.load(f)
            
            print(f"  Keys: {list(data.keys())}")
            print(f"  Pre GID: {data.get('pre_gid', 'N/A')}")
            print(f"  Post GID: {data.get('post_gid', 'N/A')}")
            print(f"  Rho config: {data.get('rho_config', 'N/A')}")
            print(f"  Time array shape: {data['t'].shape if 't' in data else 'N/A'}")
            print(f"  Voltage array shape: {data['v'].shape if 'v' in data else 'N/A'}")
            print(f"  Pre spikes type: {type(data['pre_spikes']) if 'pre_spikes' in data else 'N/A'}")
            
        except Exception as e:
            print(f"  ERROR reading file: {e}")

if __name__ == "__main__":
    test_pickle_files()
