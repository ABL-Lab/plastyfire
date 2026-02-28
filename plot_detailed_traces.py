import os
import sys
import pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Fallback threshold params (used only when theta_d/theta_p not saved in pkl)
# Update this to the active param set when switching versions
DHURUVA_PARAMS_FALLBACK = {
    "a00": 1.0003490121729528,
    "a01": 1.0877438002708641,
    "a10": 1.3917207001380245,
    "a11": 4.129505808304705,
    "a20": 2.3079274077684175,
    "a21": 1.9619384872398133,
    "a30": 4.456560780131431,
    "a31": 4.80485148600644,
}

def load_data(pkl_path):
    with open(pkl_path, "rb") as f:
        return pickle.load(f)

def calculate_mean_thresholds(synprops):
    """Compute mean theta_d and theta_p from synprops, avoiding hardcoded params where possible."""
    if not synprops:
        return None, None
        
    # If the thresholds were saved directly by the simulation, use them!
    if "theta_d_GB" in synprops and "theta_p_GB" in synprops:
        return np.mean(synprops["theta_d_GB"]), np.mean(synprops["theta_p_GB"])
        
    # Fallback for older data without these keys
    if "Cpre" not in synprops:
        return None, None
    Cpre = np.asarray(synprops["Cpre"])
    Cpost = np.asarray(synprops["Cpost"])
    locs = synprops.get("loc", ["basal"] * len(Cpre))
    
    tds = []
    tps = []
    for cp, cpo, loc in zip(Cpre, Cpost, locs):
        if loc == "basal":
            td = DHURUVA_PARAMS_FALLBACK["a00"] * cp + DHURUVA_PARAMS_FALLBACK["a01"] * cpo
            tp = DHURUVA_PARAMS_FALLBACK["a10"] * cp + DHURUVA_PARAMS_FALLBACK["a11"] * cpo
        else:
            td = DHURUVA_PARAMS_FALLBACK["a20"] * cp + DHURUVA_PARAMS_FALLBACK["a21"] * cpo
            tp = DHURUVA_PARAMS_FALLBACK["a30"] * cp + DHURUVA_PARAMS_FALLBACK["a31"] * cpo
        tds.append(td)
        tps.append(tp)
    return np.mean(tds), np.mean(tps)

def plot_single_protocol(ax_cai, ax_eff, ax_rho, ax_prime, ax_cacicr, data, title, is_left=True, x_lim=None):
    t = np.array(data['t']) / 1000.0
    
    pre_spikes = np.array(data.get('prespikes', [])) / 1000.0
    post_spikes = np.array(data.get('postspikes', [])) / 1000.0
    
    # --- 1. Calcium Panel ---
    if 'cai_CR' in data:
        cai = data['cai_CR']
        # Convert mM to uM
        cai_uM = cai * 1000.0
        cai_mean = np.mean(cai_uM, axis=1) if cai_uM.ndim > 1 else cai_uM
        ax_cai.plot(t, cai_mean, 'k-', lw=1.5)
        
    y_max = ax_cai.get_ylim()[1]
    if y_max < 0.1: y_max = 2.0
        
    y_spike_pre = y_max * 1.05
    y_spike_post = y_max * 0.95
    for s in pre_spikes:
        ax_cai.plot([s, s], [y_spike_pre, y_spike_pre + y_max*0.1], color='#66b3e6', lw=3)
    for s in post_spikes:
        ax_cai.plot([s, s], [y_spike_post, y_spike_post + y_max*0.1], color='#ffaa55', lw=3)
        
    ax_cai.set_title(title, fontsize=14)
    if is_left:
        ax_cai.set_ylabel("[Ca2+]i ($\mu$M)", fontsize=12)

    # --- 2. Effcai (c*) Panel ---
    if 'effcai_GB' in data:
        eff = data['effcai_GB']
        eff_mean = np.mean(eff, axis=1) if eff.ndim > 1 else eff
        ax_eff.plot(t, eff_mean, 'k-', lw=1.5)
        
        # Calculate thresholds dynamically
        td, tp = calculate_mean_thresholds(data.get("synprop", {}))
        if td is not None and tp is not None:
            ax_eff.axhline(tp, color='#ff6666', linestyle='--', alpha=0.8, lw=1)
            ax_eff.axhline(td, color='#66cc66', linestyle='--', alpha=0.8, lw=1)
            if is_left:
                ax_eff.text(t[0], tp, r" $\theta_p$", color='#ff6666', fontsize=12, va='bottom')
                ax_eff.text(t[0], td, r" $\theta_d$", color='#66cc66', fontsize=12, va='bottom')
    
    if is_left:
        ax_eff.set_ylabel("$c^*$", fontsize=12)

    # --- 3. Rho Panel ---
    if 'rho_GB' in data:
        rho = data['rho_GB']
        rho_mean = np.mean(rho, axis=1) if rho.ndim > 1 else rho
        ax_rho.plot(t, rho_mean, 'k-', lw=1.5)
        ax_rho.set_ylim(-0.05, 1.05)
        ax_rho.axhline(0.5, color='gray', linestyle='--', alpha=0.5)
    
    if is_left:
        ax_rho.set_ylabel(r"$\rho$", fontsize=12)

    # --- 4. Prime CICR Panel ---
    if 'prime_CICR' in data:
        prime = data['prime_CICR']
        prime_mean = np.mean(prime, axis=1) if prime.ndim > 1 else prime
        ax_prime.plot(t, prime_mean, 'k-', lw=1.5)
    
    if is_left:
        ax_prime.set_ylabel("Prime CICR", fontsize=12)

    # --- 5. Ca CICR Panel ---
    if 'ca_cicr' in data:
        ca_c = data['ca_cicr']
        # Convert mM to uM
        ca_c_uM = ca_c * 1000.0
        ca_c_mean = np.mean(ca_c_uM, axis=1) if ca_c_uM.ndim > 1 else ca_c_uM
        ax_cacicr.plot(t, ca_c_mean, 'k-', lw=1.5)
        
    if is_left:
        ax_cacicr.set_ylabel("Ca CICR ($\mu$M)", fontsize=12)
    ax_cacicr.set_xlabel("Time (s)", fontsize=12)
    
    for ax in [ax_cai, ax_eff, ax_rho, ax_prime, ax_cacicr]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
    # X Limits
    if x_lim is not None:
        ax_cai.set_xlim(x_lim)
    else:
        # Auto-Zoom
        all_spikes = list(pre_spikes) + list(post_spikes)
        if len(all_spikes) > 0:
            x_min = max(0, min(all_spikes) - 0.05)
            x_max = max(all_spikes) + 0.6
            ax_cai.set_xlim(x_min, x_max)
        else:
            ax_cai.set_xlim(0, max(t) if len(t) > 0 else 0.8)
        
    ax_eff.set_xlim(ax_cai.get_xlim())
    ax_rho.set_xlim(ax_cai.get_xlim())
    ax_prime.set_xlim(ax_cai.get_xlim())
    ax_cacicr.set_xlim(ax_cai.get_xlim())

def main():
    import argparse
    import glob
    parser = argparse.ArgumentParser(description="Plot detailed traces")
    parser.add_argument("--pair-dir", required=True, help="Pair directory containing protocol traces")
    parser.add_argument("--out-base-dir", required=True, help="Base directory to save figures")
    args = parser.parse_args()
    
    pair_name = os.path.basename(os.path.normpath(args.pair_dir))
    
    pkl_files = glob.glob(os.path.join(args.pair_dir, "*", "simulation_traces.pkl"))
    if not pkl_files:
        print(f"Cannot find traces in {args.pair_dir}")
        return
        
    for pkl_path in pkl_files:
        protocol_dir = os.path.basename(os.path.dirname(pkl_path))
        d = load_data(pkl_path)
        
        # protocol_dir e.g. '10Hz_10ms'
        if "_" in protocol_dir:
            freq, delay = protocol_dir.split("_", 1)
        else:
            freq, delay = protocol_dir, "unknown"
            
        out_dir = os.path.join(args.out_base_dir, freq, delay)
        os.makedirs(out_dir, exist_ok=True)
        
        out_file = os.path.join(out_dir, f"{pair_name}_traces.png")
        zoomed_file = os.path.join(out_dir, f"{pair_name}_traces_zoomed.png")
        
        # Plot 1: Standard view
        fig, axes = plt.subplots(5, 1, figsize=(6, 10), sharex=False)
        plot_single_protocol(axes[0], axes[1], axes[2], axes[3], axes[4], d, f"{pair_name} - {protocol_dir}", is_left=True)
        plt.tight_layout()
        plt.savefig(out_file, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved {out_file}")
        
        # Plot 2: Fixed zoom 0-3 seconds
        fig2, axes2 = plt.subplots(5, 1, figsize=(6, 10), sharex=False)
        plot_single_protocol(axes2[0], axes2[1], axes2[2], axes2[3], axes2[4], d, f"{pair_name} - {protocol_dir}", is_left=True, x_lim=(0, 3.0))
        plt.tight_layout()
        plt.savefig(zoomed_file, dpi=300, bbox_inches='tight')
        plt.close(fig2)
        print(f"Saved {zoomed_file}")

if __name__ == '__main__':
    main()
