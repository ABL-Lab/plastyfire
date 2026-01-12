
import os
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt
import glob
import logging

# Add plastyfire to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "plastyfire"))
from plastyfire.evaluator import compute_epsp_prefire

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

RESULTS_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/evaluation_results/Chindemi_params/sim_results"
PROTOCOL = "mrk97_08"

def get_epsp_ratio(sim_file):
    try:
        with open(sim_file, "rb") as f:
            raw_results = pickle.load(f)
        
        # Extract info from filename
        # Format: simulation_{param_hash}_{protocol}_{pre}_{post}.pkl
        basename = os.path.basename(sim_file)
        parts = basename.split("_")
        # simulation, hash, protocol (can be multiple parts), pre, post.pkl
        # Assuming protocol is mrk97_08 (2 parts)
        # simulation, hash, mrk97, 08, pre, post.pkl
        pre_gid = parts[-2]
        post_gid = parts[-1].replace(".pkl", "")
        
        rho_data = raw_results["rho_GB"]
        if len(rho_data) > 100:
            rho_data = np.transpose(rho_data)

        initial_rho = [0 if k[0] < 0.5 else 1 for k in rho_data]
        final_rho = [0 if k[-1] < 0.5 else 1 for k in rho_data]
        
        initial_rho_str = "_".join(map(str, initial_rho))
        final_rho_str = "_".join(map(str, final_rho))
        
        ephys_before_path = os.path.join(RESULTS_DIR, f"ephys_data_{pre_gid}_{post_gid}_{initial_rho_str}.pkl")
        ephys_after_path = os.path.join(RESULTS_DIR, f"ephys_data_{pre_gid}_{post_gid}_{final_rho_str}.pkl")
        
        if not os.path.exists(ephys_before_path):
            logger.warning(f"Missing ephys file: {ephys_before_path}")
            return None
        if not os.path.exists(ephys_after_path):
            logger.warning(f"Missing ephys file: {ephys_after_path}")
            return None
            
        epsp_before = compute_epsp_prefire(ephys_before_path)
        epsp_after = compute_epsp_prefire(ephys_after_path)
        
        if epsp_before == 0:
            return None
            
        return epsp_after / epsp_before
        
    except Exception as e:
        logger.error(f"Error processing {sim_file}: {e}")
        return None

def main():
    sim_files = glob.glob(os.path.join(RESULTS_DIR, f"simulation_*_{PROTOCOL}_*.pkl"))
    logger.info(f"Found {len(sim_files)} simulation files for {PROTOCOL}")
    
    ratios = []
    for sim_file in sim_files:
        ratio = get_epsp_ratio(sim_file)
        if ratio is not None:
            ratios.append(ratio)
            
    logger.info(f"Computed {len(ratios)} ratios")
    
    if not ratios:
        logger.error("No ratios computed. Exiting.")
        return

    # Plot Histogram
    plt.figure(figsize=(10, 6))
    plt.hist(ratios, bins=20, edgecolor='black', alpha=0.7)
    plt.axvline(np.mean(ratios), color='red', linestyle='dashed', linewidth=1, label=f'Mean: {np.mean(ratios):.4f}')
    plt.xlabel('EPSP Ratio')
    plt.ylabel('Count')
    plt.title(f'EPSP Ratio Histogram for {PROTOCOL}')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    output_file = f"{PROTOCOL}_histogram.png"
    plt.savefig(output_file, dpi=300)
    logger.info(f"Histogram saved to {output_file}")
    print(f"Mean EPSP Ratio: {np.mean(ratios):.4f}")
    print(f"SEM EPSP Ratio: {np.std(ratios, ddof=1) / np.sqrt(len(ratios)):.4f}")

if __name__ == "__main__":
    main()
