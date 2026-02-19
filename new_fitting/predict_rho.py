"""
Predict rho traces from effcai traces using the Graupner-Brunel plasticity model.

The GluSynapse.mod file defines the relationship between effcai and rho as:
    rho_GB' = ( - rho_GB*(1 - rho_GB)*(rho_star_GB - rho_GB)
             + pot_GB*gamma_p_GB*(1 - rho_GB)
             - dep_GB*gamma_d_GB*rho_GB ) / ((1e3)*tau_ind_GB)
             
Where:
    pot_GB = 1 if effcai_GB > theta_p_GB else 0
    dep_GB = 1 if effcai_GB > theta_d_GB else 0
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt
import os

# Parameters from GluSynapse.mod (GLOBAL)
tau_ind_GB = 70.0  # s (time constant for induction)
rho_star_GB = 0.5
tau_effca_GB = 278.3177658387  # ms

# Parameters provided by the user (Chindemi_params)
CHINDEMI_PARAMS = {
    "gamma_d_GB_GluSynapse": 101.5,
    "gamma_p_GB_GluSynapse": 216.2,
    "a00": 1.002,
    "a01": 1.954,
    "a10": 1.159,
    "a11": 2.483,
    "a20": 1.127,
    "a21": 2.456,
    "a30": 5.236,
    "a31": 1.782,
}

def calculate_thresholds(c_pre, c_post, loc, params):
    """
    Calculate theta_d and theta_p based on synapse properties and fit parameters.
    Based on logic in simulator.py _set_local_params
    """
    # Initialize with default values (0.0 seems appropriate if keys missing, 
    # though in simulator.py they are set conditionally)
    theta_d = 0.0
    theta_p = 0.0
    
    if loc == "basal":
        if "a00" in params and "a01" in params:
            theta_d = params["a00"] * c_pre + params["a01"] * c_post
        if "a10" in params and "a11" in params:
            theta_p = params["a10"] * c_pre + params["a11"] * c_post
    elif loc == "apical":
        if "a20" in params and "a21" in params:
            theta_d = params["a20"] * c_pre + params["a21"] * c_post
        if "a30" in params and "a31" in params:
            theta_p = params["a30"] * c_pre + params["a31"] * c_post
            
    return theta_d, theta_p

def compute_rho_euler(effcai_trace, t, theta_d, theta_p, params, rho0=0.0):
    """
    Compute rho using Euler method.
    
    rho_GB' = ( - rho_GB*(1 - rho_GB)*(rho_star_GB - rho_GB)
             + pot_GB*gamma_p_GB*(1 - rho_GB)
             - dep_GB*gamma_d_GB*rho_GB ) / ((1e3)*tau_ind_GB)
    """
    n_points = len(t)
    rho = np.zeros(n_points)
    rho[0] = rho0
    
    gamma_d = params["gamma_d_GB_GluSynapse"]
    gamma_p = params["gamma_p_GB_GluSynapse"]
    tau_ind_ms = tau_ind_GB * 1000.0 # Convert to ms
    
    # Pre-calculate constants
    inv_tau = 1.0 / tau_ind_ms
    
    for i in range(n_points - 1):
        dt = t[i + 1] - t[i]
        
        # Determine states
        effcai = effcai_trace[i]
        # Logic from GluSynapse:
        # WATCH (effcai_GB > theta_d_GB) 2 -> dep_GB = 1
        # WATCH (effcai_GB < theta_d_GB) 3 -> dep_GB = 0
        dep = 1.0 if effcai > theta_d else 0.0
        pot = 1.0 if effcai > theta_p else 0.0
        
        curr_rho = rho[i]
        
        # Cubic term (bistability)
        cubic_term = -curr_rho * (1.0 - curr_rho) * (rho_star_GB - curr_rho)
        
        # Potentiation term
        pot_term = pot * gamma_p * (1.0 - curr_rho)
        
        # Depression term
        dep_term = -dep * gamma_d * curr_rho
        
        # Derivative
        drho = (cubic_term + pot_term + dep_term) * inv_tau
        
        # Update
        rho[i + 1] = curr_rho + dt * drho
        
        # Clip to [0, 1] as rho is a probability/fraction
        # Although theoretically it should stay within bounds if ODE is stable and dt small enough,
        # numerical errors might push it out.
        if rho[i+1] > 1.0: rho[i+1] = 1.0
        if rho[i+1] < 0.0: rho[i+1] = 0.0
        
    return rho

if __name__ == "__main__":
    # Load data
    data_path = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/trace_results/Chindemi_params/180164-197248/10Hz_5ms/simulation_traces.pkl"
    if not os.path.exists(data_path):
        print(f"Error: Data file not found at {data_path}")
        # Try to find recent pickle if specific one fails, or just list dir
        base_dir = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/trace_results/Chindemi_params"
        print(f"Checking {base_dir}...")
        try:
             # Find a valid file if hardcoded one doesn't exist
             for root, dirs, files in os.walk(base_dir):
                 for file in files:
                     if file.endswith("simulation_traces.pkl") or file.endswith(".pkl"):
                         data_path = os.path.join(root, file)
                         print(f"Found alternative: {data_path}")
                         break
                 if os.path.exists(data_path): break
        except Exception as e:
            print(f"Error searching for files: {e}")
            
    if not os.path.exists(data_path):
        print("Could not find any simulation traces file.")
        exit(1)
        
    print(f"Loading data from {data_path}")
    with open(data_path, "rb") as f:
        data = pickle.load(f)
    
    # Get time vector
    # Some files might have 't' key, others 'time'. simulation_traces.pkl usually has 't' or it's inferred.
    if "t" in data:
        t = data["t"]
    elif "time" in data:
        t = data["time"]
    else:
        # Infer from shape assuming dt=0.1 or similar if not found? 
        # But 'analytical_effcai.py' used: t = data.get("t", np.arange(data["effcai_GB"].shape[-1]))
         # Let's hope it's there or consistent with analytical_effcai.py
         # If shape is [synapses, time], t length should be time
         if "effcai_GB" in data:
             shape = np.shape(data["effcai_GB"])
             n_steps = shape[-1] if len(shape) > 1 else shape[0]
             # Check if t is consistent
             t = np.arange(n_steps) * 0.025 # Default dt simulation usually
             print("Warning: 't' not found in data, created dummy time vector with dt=0.025ms")

    
    # Get traces from simulation
    # Data shape is [synapses, timepoints]
    effcai_sim = np.asarray(data["effcai_GB"])
    rho_sim = np.asarray(data["rho_GB"])
    
    # Transpose if needed (expecting [synapses, timepoints])
    # T vector length check
    if len(t) != effcai_sim.shape[-1]:
        if len(t) == effcai_sim.shape[0]:
            effcai_sim = effcai_sim.T
            rho_sim = rho_sim.T
        else:
            print(f"Error: Time vector length {len(t)} mismatch with data shape {effcai_sim.shape}")
            # Try to fix t if it's just a variable in the pickle that was saved differently
            pass
            
    # Get synapse properties
    # In some pickles, synprop might be a dict inside data, or keys in data
    synprop = data.get("synprop", {})
    if not synprop and "synapseID" in data:
        # Maybe flat structure
        synprop = data
        
    syn_ids = synprop.get("synapseID", [])
    if not isinstance(syn_ids, list) and not isinstance(syn_ids, np.ndarray):
         # Try to find list of synapses another way or create dummy
         n_syns = effcai_sim.shape[0]
         syn_ids = np.arange(n_syns)
         print(f"Warning: synapseID not found, using indices 0..{n_syns-1}")
         
    c_pres = synprop.get("Cpre", np.zeros(len(syn_ids)))
    c_posts = synprop.get("Cpost", np.zeros(len(syn_ids)))
    locs = synprop.get("loc", ["basal"] * len(syn_ids)) # Default to basal if not specified? 
    
    print(f"Found {len(syn_ids)} synapses.")
    
    # Filter to valid indices
    valid_indices = range(min(len(syn_ids), effcai_sim.shape[0]))
    
    # Setup plotting
    # Select a few synapses to plot, preferably with different behaviors if possible
    # Or just first 3
    plot_indices = list(valid_indices)[:3]
    
    n_plots = len(plot_indices)
    if n_plots == 0:
        print("No valid synapses to plot.")
        exit(0)
        
    fig, axes = plt.subplots(n_plots, 1, figsize=(12, 4*n_plots), sharex=True)
    if n_plots == 1:
        axes = [axes]
        
    mse_list = []
    
    print("\nProcessing synapses...")
        
    for idx_i, idx_data in enumerate(plot_indices):
        syn_id = syn_ids[idx_data]
        
        # Handle if properties are missing/scalar/array
        c_pre = c_pres[idx_data] if len(c_pres) > idx_data else 0.0
        c_post = c_posts[idx_data] if len(c_posts) > idx_data else 0.0
        loc = locs[idx_data] if len(locs) > idx_data else "basal"
        
        # Calculate thresholds
        theta_d, theta_p = calculate_thresholds(c_pre, c_post, loc, CHINDEMI_PARAMS)
        
        print(f"\nSynapse {idx_data} (ID: {syn_id}, Loc: {loc}):")
        print(f"  Cpre: {c_pre:.4f}, Cpost: {c_post:.4f}")
        print(f"  theta_d: {theta_d:.6f}, theta_p: {theta_p:.6f}")
        
        # Get traces
        effcai_trace = effcai_sim[idx_data]
        rho_trace_sim = rho_sim[idx_data]
        rho0 = rho_trace_sim[0] # Initial condition from sim
        
        print(f"  Rho range sim: [{np.min(rho_trace_sim):.4f}, {np.max(rho_trace_sim):.4f}]")
        
        # Predict rho
        rho_predicted = compute_rho_euler(effcai_trace, t, theta_d, theta_p, CHINDEMI_PARAMS, rho0=rho0)
        
        # Calculate error
        mse = np.mean((rho_trace_sim - rho_predicted)**2)
        rmse = np.sqrt(mse)
        mse_list.append(rmse)
        print(f"  RMSE: {rmse:.6f}")
        
        # Plot
        ax = axes[idx_i]
        
        # Plot rho
        ax.plot(t/1000, rho_trace_sim, 'k-', label='Simulated rho', linewidth=2, alpha=0.6)
        ax.plot(t/1000, rho_predicted, 'r--', label='Predicted rho', linewidth=1.5, alpha=0.9)
        
        # Add thresholds info to title
        ax.set_title(f"Synapse {syn_id} ({loc}) - RMSE: {rmse:.6f}\n" 
                     f"theta_d={theta_d:.4f}, theta_p={theta_p:.4f}")
        ax.set_ylabel("rho")
        if idx_i == 0:
            ax.legend(loc="upper right")
        ax.grid(True, alpha=0.3)
        
        # Optional: Plot effcai (scaled) or logical states?
        # Let's plot 'pot' and 'dep' regions
        # Scale effcai to fit? No, use twinx
        ax2 = ax.twinx()
        ax2.plot(t/1000, effcai_trace, 'g-', label='effcai', alpha=0.2, linewidth=0.5)
        # Add threshold lines
        ax2.axhline(theta_d, color='blue', linestyle=':', alpha=0.4, label='theta_d')
        ax2.axhline(theta_p, color='orange', linestyle=':', alpha=0.4, label='theta_p')
        ax2.set_ylabel("effcai (mM)", color='green')
        ax2.tick_params(axis='y', labelcolor='green')

    axes[-1].set_xlabel("Time (s)")
    plt.tight_layout()
    output_png = "rho_comparison.png"
    plt.savefig(output_png, dpi=150)
    print(f"\nSaved {output_png}")
    print(f"Average RMSE: {np.mean(mse_list):.6f}")
