import os
import sys
import pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import glob
from collections import defaultdict

# Use the same fallback thresholds
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

def calculate_first_synapse_thresholds(synprops):
    if not synprops:
        return None, None
    if "theta_d_GB" in synprops and "theta_p_GB" in synprops:
        return synprops["theta_d_GB"][0], synprops["theta_p_GB"][0]
        
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
    return tds[0], tps[0]

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Plot stacked detailed traces by frequency")
    parser.add_argument("--pair-dir", required=True, help="Pair directory containing protocol traces")
    parser.add_argument("--out-dir", required=True, help="Output directory to save figures")
    parser.add_argument("--x-lim", type=float, nargs=2, default=[-2.0, 50.0], help="X-axis limits [min, max]")
    parser.add_argument("--zoomed-x-lim", type=float, nargs=2, default=[0.0, 3.0], help="Zoomed X-axis limits")
    args = parser.parse_args()
    
    pair_name = os.path.basename(os.path.normpath(args.pair_dir))
    pkl_files = glob.glob(os.path.join(args.pair_dir, "*", "simulation_traces.pkl"))
    if not pkl_files:
        print(f"Cannot find traces in {args.pair_dir}")
        return
        
    # Group by frequency
    grouped_data = defaultdict(list)
    for pkl_path in pkl_files:
        protocol_dir = os.path.basename(os.path.dirname(pkl_path))
        if "_" in protocol_dir:
            freq, delay = protocol_dir.split("_", 1)
        else:
            freq, delay = protocol_dir, "unknown"
        grouped_data[freq].append((delay, pkl_path))
        
    os.makedirs(args.out_dir, exist_ok=True)
    
    # Colormap
    cmap = plt.get_cmap("coolwarm")
    
    for freq, protocols in grouped_data.items():
        print(f"Processing freq {freq} with {len(protocols)} delays...")
        # Sort delays numerically if possible
        def extract_delay_num(d):
            import re
            m = re.match(r"(-?\d+)", d)
            return float(m.group(1)) if m else float('inf')
        protocols.sort(key=lambda x: extract_delay_num(x[0]))
        
        for view_name, x_lim in [("full", args.x_lim), ("zoomed", args.zoomed_x_lim)]:
            fig, axes = plt.subplots(5, 1, figsize=(8, 12), sharex=True)
            ax_cai, ax_eff, ax_rho, ax_prime, ax_cacicr = axes
            
            for ax in axes:
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                
            has_thresholds = False
            num_protocols = len(protocols)
            for idx, (delay, pkl_path) in enumerate(protocols):
                color = cmap(float(idx) / max(1, num_protocols - 1))
                data = load_data(pkl_path)
                t = np.array(data['t']) / 1000.0
                
                # Mask to restrict drawn data to the view limits so Matplotlib correctly scales Y
                mask = (t >= x_lim[0] - 0.5) & (t <= x_lim[1] + 0.5)
                t_plot = t[mask]
                
                # Plot traces for the FIRST synapse
                if 'cai_CR' in data:
                    cai = data['cai_CR'] * 1000.0
                    cai_plot = cai[:, 0] if cai.ndim > 1 else cai
                    ax_cai.plot(t_plot, cai_plot[mask], color=color, lw=1.5, label=delay)
                
                if 'effcai_GB' in data:
                    eff = data['effcai_GB']
                    eff_plot = eff[:, 0] if eff.ndim > 1 else eff
                    ax_eff.plot(t_plot, eff_plot[mask], color=color, lw=1.5)
                    
                    if not has_thresholds:
                        td, tp = calculate_first_synapse_thresholds(data.get("synprop", {}))
                        if td is not None and tp is not None:
                            ax_eff.axhline(tp, color='#ff6666', linestyle='--', alpha=0.8, lw=1)
                            ax_eff.axhline(td, color='#66cc66', linestyle='--', alpha=0.8, lw=1)
                            ax_eff.text(max(0, x_lim[0]), tp, r" $\theta_p$", color='#ff6666', fontsize=12, va='bottom')
                            ax_eff.text(max(0, x_lim[0]), td, r" $\theta_d$", color='#66cc66', fontsize=12, va='bottom')
                            has_thresholds = True

                if 'rho_GB' in data:
                    rho = data['rho_GB']
                    rho_plot = rho[:, 0] if rho.ndim > 1 else rho
                    ax_rho.plot(t_plot, rho_plot[mask], color=color, lw=1.5)
                
                if 'P_CICR' in data:
                    prime = data['P_CICR']
                    prime_plot = prime[:, 0] if prime.ndim > 1 else prime
                    ax_prime.plot(t_plot, prime_plot[mask], color=color, lw=1.5)
                    
                if 'ca_ryr_CICR' in data:
                    ca_c = data['ca_ryr_CICR'] * 1000.0
                    ca_c_plot = ca_c[:, 0] if ca_c.ndim > 1 else ca_c
                    ax_cacicr.plot(t_plot, ca_c_plot[mask], color=color, lw=1.5)
            
            # Formatting
            ax_cai.set_ylabel(r"[Ca2+]i ($\mu$M)", fontsize=12)
            # Make sure legend is drawn outside the plot
            ax_cai.legend(title=r"$\Delta t$", loc='center left', bbox_to_anchor=(1.05, 0.5))
            ax_cai.set_title(f"{pair_name} - {freq}", fontsize=14)
            ax_cai.set_ylim(bottom=-0.05)
            
            ax_eff.set_ylabel("$c^*$", fontsize=12)
            ax_eff.set_ylim(bottom=-0.01)
            
            ax_rho.set_ylabel(r"$\rho$", fontsize=12)
            ax_rho.set_ylim(-0.05, 1.05)
            ax_rho.axhline(0.5, color='gray', linestyle='--', alpha=0.5)
            
            ax_prime.set_ylabel("P CICR", fontsize=12)
            
            ax_cacicr.set_ylabel(r"Ca RyR CICR ($\mu$M)", fontsize=12)
            ax_cacicr.set_ylim(bottom=-0.05)
            ax_cacicr.set_xlabel("Time (s)", fontsize=12)
            
            ax_cai.set_xlim(x_lim)
            
            # Align labels
            fig.align_ylabels(axes)
            plt.tight_layout()
            
            out_file = os.path.join(args.out_dir, f"{pair_name}_{freq}_stacked_{view_name}.png")
            plt.savefig(out_file, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved {out_file}")

if __name__ == '__main__':
    main()
