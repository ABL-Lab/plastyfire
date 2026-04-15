import sys
import pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def main(pkl_path, out_path="trace_plot.png"):
    print(f"Loading {pkl_path}...")
    try:
        with open(pkl_path, "rb") as f:
            data = pickle.load(f)
    except Exception as e:
        print(f"Error loading {pkl_path}: {e}")
        return
        
    syn_ids = list(data["c_pre"].keys())
    if len(syn_ids) == 0:
        print("No synapses found in data.")
        return
        
    syn_id = syn_ids[0]
    syn_index = 0 # first synapse in syn_props lists
    
    t = data["t"]
    cai_CR = data["cai_CR"][syn_id]
    effcai_GB = data["effcai_GB"][syn_id]
    
    cpre = data["c_pre"][syn_id]
    cpost = data["c_post"][syn_id]
    
    theta_d = data["syn_props"]["theta_d_GB"][syn_index]
    theta_p = data["syn_props"]["theta_p_GB"][syn_index]
    
    # Handle rho if available
    rho = None
    if "rho_GB" in data:
        rho = data["rho_GB"][syn_id]
    
    # 3 subplots if rho is included, else 2
    n_plots = 3 if rho is not None else 2
    fig, axs = plt.subplots(n_plots, 1, figsize=(10, 4 * n_plots), sharex=True)
    if n_plots == 2:
        axs = [axs[0], axs[1]]  # Ensure list indexing works same way
    
    # Plot cai_CR
    axs[0].plot(t, cai_CR, label="cai_CR", color="blue")
    axs[0].set_ylabel("Ca2+ concentration (mM)")
    axs[0].set_title(f"Synapse {syn_id} - Calcium Traces")
    axs[0].legend()
    
    # Plot effcai_GB
    axs[1].plot(t, effcai_GB, label="effcai_GB", color="orange")
    axs[1].axhline(theta_d, color="red", linestyle="--", label=f"theta_d ({theta_d:.4g})")
    axs[1].axhline(theta_p, color="green", linestyle="--", label=f"theta_p ({theta_p:.4g})")
    axs[1].axhline(cpre, color="purple", linestyle=":", label=f"cpre ({cpre:.4g})")
    axs[1].axhline(cpost, color="brown", linestyle=":", label=f"cpost ({cpost:.4g})")
    axs[1].set_ylabel("Effective Ca2+")
    axs[1].legend(loc='upper right')
    
    if rho is not None:
        axs[2].plot(t, rho, label="rho_GB", color="purple")
        axs[2].set_ylabel("Rho")
        axs[2].set_ylim([-0.1, 1.1])
        axs[2].axhline(0.5, color="k", linestyle="--", alpha=0.5)
        axs[2].legend()

    axs[-1].set_xlabel("Time (ms)")
    
    # Also plot the prespikes and postspikes as vertical lines
    for ax in axs:
        for p_spike in data.get("prespikes", []):
            ax.axvline(p_spike, color="gray", linestyle="-", alpha=0.3, label="Pre spike" if p_spike == data["prespikes"][0] else "")
        for p_spike in data.get("postspikes", []):
            ax.axvline(p_spike, color="black", linestyle="-", alpha=0.3, label="Post spike" if p_spike == data["postspikes"][0] else "")
            
    # Add axvspan for induction period (commonly fastforward point is 1000, but we can just use the spikes)
    # Deduplicate legend
    handles, labels = axs[0].get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    axs[0].legend(by_label.values(), by_label.keys(), loc='upper right')
        
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    print(f"Saved plot to {out_path}")
    print(f"  Max effcai_GB: {np.max(effcai_GB):.4g}")
    print(f"  cpre: {cpre:.4g}, cpost: {cpost:.4g}")
    print(f"  theta_d: {theta_d:.4g}, theta_p: {theta_p:.4g}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("pkl_path", help="Path to simulation_traces.pkl")
    parser.add_argument("--out", default="trace_plot.png", help="Output plot filename")
    args = parser.parse_args()
    main(args.pkl_path, args.out)
