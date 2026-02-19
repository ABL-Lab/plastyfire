"""
Analytical solver for effective calcium (effcai_GB) from calcium traces (cai_CR).

The GluSynapse.mod file defines the relationship as:
    effcai_GB' = -effcai_GB/tau_effca_GB + (cai_CR - min_ca_CR)

This is a first-order linear ODE that can be solved analytically using an exponential integrator.

Parameters from GluSynapse.mod:
    tau_effca_GB = 200 ms  (time constant for effective calcium)
    min_ca_CR = 70e-6 mM = 70 nM  (baseline calcium concentration)
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt


def compute_effcai_analytical(cai_trace, t, tau_effca=200.0, min_ca=70e-6, effcai0=0.0):
    """
    Compute effective calcium using an analytical exponential integrator.
    
    The ODE is: d(effcai)/dt = -effcai/tau + (cai - min_ca)
    
    Using the exponential integrator (exact solution between timesteps):
        effcai[n+1] = effcai[n] * exp(-dt/tau) + (cai[n] - min_ca) * tau * (1 - exp(-dt/tau))
    
    Parameters
    ----------
    cai_trace : np.ndarray
        Calcium concentration trace (mM)
    t : np.ndarray
        Time vector (ms)
    tau_effca : float
        Time constant for effective calcium (ms), default 200 ms
    min_ca : float
        Baseline calcium concentration (mM), default 70e-6 mM = 70 nM
    effcai0 : float
        Initial effective calcium (mM), default 0
    
    Returns
    -------
    effcai : np.ndarray
        Effective calcium trace (same units as input after integration)
    """
    n_points = len(t)
    effcai = np.zeros(n_points)
    effcai[0] = effcai0
    
    for i in range(n_points - 1):
        dt = t[i + 1] - t[i]
        decay = np.exp(-dt / tau_effca)
        driving = cai_trace[i] - min_ca
        # Exponential integrator: exact solution for constant driving between steps
        effcai[i + 1] = effcai[i] * decay + driving * tau_effca * (1 - decay)
    
    return effcai


def compute_effcai_euler(cai_trace, t, tau_effca=200.0, min_ca=70e-6, effcai0=0.0):
    """
    Compute effective calcium using Euler method (as in NEURON with METHOD euler).
    
    The ODE is: d(effcai)/dt = -effcai/tau + (cai - min_ca)
    
    Euler update:
        effcai[n+1] = effcai[n] + dt * (-effcai[n]/tau + (cai[n] - min_ca))
    
    Parameters
    ----------
    cai_trace : np.ndarray
        Calcium concentration trace (mM)
    t : np.ndarray
        Time vector (ms)
    tau_effca : float
        Time constant for effective calcium (ms), default 200 ms
    min_ca : float
        Baseline calcium concentration (mM), default 70e-6 mM = 70 nM
    effcai0 : float
        Initial effective calcium (mM), default 0
    
    Returns
    -------
    effcai : np.ndarray
        Effective calcium trace
    """
    n_points = len(t)
    effcai = np.zeros(n_points)
    effcai[0] = effcai0
    
    for i in range(n_points - 1):
        dt = t[i + 1] - t[i]
        driving = cai_trace[i] - min_ca
        deriv = -effcai[i] / tau_effca + driving
        effcai[i + 1] = effcai[i] + dt * deriv
    
    return effcai


if __name__ == "__main__":
    # Load data
    data_path = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/trace_results/Chindemi_params/180164-197248/10Hz_5ms/simulation_traces.pkl"
    with open(data_path, "rb") as f:
        data = pickle.load(f)
    
    # Get time vector
    t = data.get("t", np.arange(data["effcai_GB"].shape[-1]))
    
    # Get cai and effcai from simulation
    cai = np.asarray(data["cai_CR"])
    effcai_sim = np.asarray(data["effcai_GB"])
    
    # Ensure shape is (n_synapses, n_timepoints)
    if cai.shape[0] == len(t) and cai.shape[1] != len(t):
        cai = cai.T
    if effcai_sim.shape[0] == len(t) and effcai_sim.shape[1] != len(t):
        effcai_sim = effcai_sim.T
    
    # Analyze time steps
    dts = np.diff(t)
    print("\n" + "="*60)
    print("Time Step Analysis")
    print("="*60)
    print(f"Min dt: {np.min(dts):.4f} ms")
    print(f"Max dt: {np.max(dts):.4f} ms")
    print(f"Mean dt: {np.mean(dts):.4f} ms")
    print(f"Number of large gaps (> 100 ms): {np.sum(dts > 100)}")
    
    # Parameters from GluSynapse.mod
    tau_effca_GB = 278.3177658387  # ms
    min_ca_CR = 70e-6     # mM (70 nM)
    
    # Compute analytical solution for first synapse
    synapse_idx = 0
    cai_trace = cai[synapse_idx]
    effcai_trace_sim = effcai_sim[synapse_idx]
    
    # Use Euler method (as in NEURON with METHOD euler)
    # Make a more robust version that handles gaps by resetting or skipping
    def compute_effcai_euler_robust(cai_trace, t, tau_effca=200.0, min_ca=70e-6, effcai0=0.0):
        n_points = len(t)
        effcai = np.zeros(n_points)
        effcai[0] = effcai0
        
        for i in range(n_points - 1):
            dt = t[i + 1] - t[i]
            
            # If gap is too large, reset or assume steady state?
            # Standard Euler becomes unstable if dt > 2*tau
            # Here tau=200, so dt > 400 is unstable.
            if dt > 100.0: 
                # Reset to steady state or initial condition for the purpose of comparison
                # Or just skip update (hold previous value)
                # But physically, if gap is huge, it should decay to 0 (if driving is 0)
                # Let's approximate decay over large dt
                decay = np.exp(-dt / tau_effca)
                driving = cai_trace[i] - min_ca
                # Use analytical step for large gaps to avoid explosion
                effcai[i + 1] = effcai[i] * decay + driving * tau_effca * (1 - decay)
            else:
                driving = cai_trace[i] - min_ca
                deriv = -effcai[i] / tau_effca + driving
                effcai[i + 1] = effcai[i] + dt * deriv
        return effcai

    effcai_euler = compute_effcai_euler_robust(cai_trace, t, tau_effca_GB, min_ca_CR)
    
    # Also compute exponential integrator for comparison

    effcai_analytical = compute_effcai_analytical(cai_trace, t, tau_effca_GB, min_ca_CR)
    
    # Focus on first 50 seconds (50000 ms) where activity occurs
    mask = t < 50000
    t_plot = t[mask] / 1000  # convert to seconds
    
    # --- Comparison 1: NEURON vs Analytical ---
    fig1, ax1 = plt.subplots(figsize=(10, 5))
    ax1.plot(t_plot, effcai_trace_sim[mask], 'k-', label='NEURON simulation', linewidth=1.5, alpha=0.6)
    ax1.plot(t_plot, effcai_analytical[mask], 'r--', label='Analytical (Exp. Integrator)', linewidth=1.5, alpha=0.8)
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("effcai_GB (mM)")
    ax1.set_title("NEURON vs Analytical Solution")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("effcai_comparison_analytical.png", dpi=150)
    print("Saved effcai_comparison_analytical.png")
    
    # --- Comparison 2: NEURON vs Euler ---
    fig2, ax2 = plt.subplots(figsize=(10, 5))
    ax2.plot(t_plot, effcai_trace_sim[mask], 'k-', label='NEURON simulation', linewidth=1.5, alpha=0.6)
    ax2.plot(t_plot, effcai_euler[mask], 'b--', label='Python Euler', linewidth=1.5, alpha=0.8)
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel("effcai_GB (mM)")
    ax2.set_title("NEURON vs Python Euler Implementation")
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("effcai_comparison_euler.png", dpi=150)
    print("Saved effcai_comparison_euler.png")
    
    # statistics
    diff_analytical = effcai_trace_sim - effcai_analytical
    diff_euler = effcai_trace_sim - effcai_euler

    print("\n" + "="*60)
    print("Comparison Statistics")
    print("="*60)
    print(f"Analytical - Max absolute error: {np.max(np.abs(diff_analytical)):.6e}")
    print(f"Analytical - RMS error:          {np.sqrt(np.mean(diff_analytical**2)):.6e}")
    print("-" * 60)
    print(f"Euler      - Max absolute error: {np.max(np.abs(diff_euler)):.6e}")
    print(f"Euler      - RMS error:          {np.sqrt(np.mean(diff_euler**2)):.6e}")
    print("="*60)

