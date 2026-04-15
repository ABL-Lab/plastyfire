import os
import sys
import pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Add plastyfire to path so we can import its modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

from bluepysnap import Simulation as BpySimulation
from plastyfire.simulator import c_pre_finder, c_post_finder, ParamsGenerator, EXTRA_RECIPE_PATH


def load_data(pkl_path):
    with open(pkl_path, "rb") as f:
        return pickle.load(f)


def compute_thresholds_from_calibration(sim_config_path, fit_params,
                                        node_pop="S1nonbarrel_neurons",
                                        edge_pop="S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"):
    """
    Computes per-synapse theta_d and theta_p thresholds using dedicated calibration
    simulations, exactly as runconnectedpair() does it.

    Returns
    -------
    c_pre  : dict  {global_syn_id -> float}  peak effcai from pre-only spike
    c_post : dict  {global_syn_id -> float}  peak effcai from post-only spike
    """
    sim = BpySimulation(sim_config_path)
    pre_gid  = sim.node_sets.content["precell"]["node_id"][0]
    post_gid = sim.node_sets.content["postcell"]["node_id"][0]

    # Build the same stimulus that runconnectedpair uses for c_post_finder
    width = sim.config["inputs"]["pulse0"]["width"]
    amp   = sim.config["inputs"]["pulse0"]["amp_start"]
    stimulus = {"nspikes": 1, "freq": 0.1, "width": width, "offset": 1000, "amp": amp}

    # Generate per-synapse extra params (loc, Use, Dep, … from the recipe)
    pgen = ParamsGenerator(sim.circuit, node_pop, edge_pop, EXTRA_RECIPE_PATH)
    syn_extra_params = pgen.generate_params(pre_gid, post_gid)

    # Run calibration sims with the same baseline params (no a-values needed here)
    baseline_params = fit_params.copy() if fit_params is not None else {}

    print("  [calibration] Running c_pre_finder …")
    c_pre = c_pre_finder(sim_config_path, baseline_params, syn_extra_params,
                         pre_gid, post_gid,
                         node_pop=node_pop, edge_pop=edge_pop, fixhp=True)

    print("  [calibration] Running c_post_finder …")
    c_post = c_post_finder(sim_config_path, baseline_params, syn_extra_params,
                           pre_gid, post_gid, stimulus,
                           node_pop=node_pop, edge_pop=edge_pop, fixhp=True)

    return c_pre, c_post


def derive_thresholds(c_pre, c_post, fit_params):
    """
    Apply the linear threshold formula to each synapse.

    fit_params must contain a00/a01/a10/a11 (basal) and optionally a20/a21/a30/a31 (apical).
    syn_locs is a dict {global_syn_id -> 'basal'|'apical'}.

    Returns
    -------
    theta_d : dict  {global_syn_id -> float}
    theta_p : dict  {global_syn_id -> float}
    """
    theta_d, theta_p = {}, {}
    for syn_id in c_pre:
        cp  = c_pre[syn_id]
        cpo = c_post[syn_id]
        # Default to basal coefficients
        td = fit_params["a00"] * cp + fit_params["a01"] * cpo
        tp = fit_params["a10"] * cp + fit_params["a11"] * cpo
        theta_d[syn_id] = td
        theta_p[syn_id] = tp
    return theta_d, theta_p


def plot_single_protocol(ax_cai, ax_eff, ax_rho, data, title,
                         theta_d_first=None, theta_p_first=None,
                         is_left=True, x_lim=None):
    t = np.array(data['t']) / 1000.0

    pre_spikes  = np.array(data.get('prespikes',  [])) / 1000.0
    post_spikes = np.array(data.get('postspikes', [])) / 1000.0

    # Helper: dict-of-arrays → 1-D array for the first synapse
    def extract_array(val):
        if isinstance(val, dict):
            arr = np.array(list(val.values())[0])
        else:
            arr = np.array(val)
        return arr[:, 0] if arr.ndim > 1 else arr

    # --- 1. Calcium Panel ---
    if 'cai_CR' in data:
        cai = extract_array(data['cai_CR'])
        cai_uM = cai * 1000.0          # mM → µM
        ax_cai.plot(t, cai_uM, 'k-', lw=1.5)

    y_max = ax_cai.get_ylim()[1]
    if y_max < 0.1:
        y_max = 2.0

    y_spike_pre  = y_max * 1.05
    y_spike_post = y_max * 0.95
    for s in pre_spikes:
        ax_cai.plot([s, s], [y_spike_pre,  y_spike_pre  + y_max * 0.1], color='#66b3e6', lw=3)
    for s in post_spikes:
        ax_cai.plot([s, s], [y_spike_post, y_spike_post + y_max * 0.1], color='#ffaa55', lw=3)

    ax_cai.set_title(title, fontsize=14)
    if is_left:
        ax_cai.set_ylabel("[Ca2+]i ($\\mu$M)", fontsize=12)

    # --- 2. Effcai (c*) Panel ---
    if 'effcai_GB' in data:
        eff_plot = extract_array(data['effcai_GB'])
        ax_eff.plot(t, eff_plot, 'k-', lw=1.5)

        # Use thresholds passed in from calibration runs
        td, tp = theta_d_first, theta_p_first

        # Final fallback: read theta_d/theta_p saved directly in synprops by the simulator
        if td is None or tp is None:
            for key in ("syn_props", "synprop"):
                synprops = data.get(key, {})
                if "theta_d_GB" in synprops and "theta_p_GB" in synprops:
                    raw_td = synprops["theta_d_GB"]
                    raw_tp = synprops["theta_p_GB"]
                    td = raw_td[0] if isinstance(raw_td, (list, np.ndarray)) else raw_td
                    tp = raw_tp[0] if isinstance(raw_tp, (list, np.ndarray)) else raw_tp
                    break

        if td is not None and tp is not None:
            ax_eff.axhline(tp, color='#ff6666', linestyle='--', alpha=0.8, lw=1)
            ax_eff.axhline(td, color='#66cc66', linestyle='--', alpha=0.8, lw=1)
            if is_left:
                ax_eff.text(t[0], tp, r" $\theta_p$" + f" ({tp:.4f})", color='#ff6666', fontsize=12, va='bottom')
                ax_eff.text(t[0], td, r" $\theta_d$" + f" ({td:.4f})", color='#66cc66', fontsize=12, va='bottom')

    if is_left:
        ax_eff.set_ylabel("$c^*$", fontsize=12)

    # --- 3. Rho Panel ---
    if 'rho_GB' in data:
        rho_plot = extract_array(data['rho_GB'])
        ax_rho.plot(t, rho_plot, 'k-', lw=1.5)
        ax_rho.set_ylim(-0.05, 1.05)
        ax_rho.axhline(0.5, color='gray', linestyle='--', alpha=0.5)

    if is_left:
        ax_rho.set_ylabel(r"$\rho$", fontsize=12)

    for ax in [ax_cai, ax_eff, ax_rho]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    if x_lim is not None:
        ax_cai.set_xlim(x_lim)
    else:
        ax_cai.set_xlim(-2.0, 50.0)

    ax_eff.set_xlim(ax_cai.get_xlim())
    ax_rho.set_xlim(ax_cai.get_xlim())


def plot_cicr_protocol(ax_caer, ax_ip3, ax_S, ax_rgate, ax_Jserca, ax_Jleak,
                       ax_Jcicr, ax_h, ax_minf, ax_ninf, ax_nves,
                       data, title, is_left=True, x_lim=None):
    t = np.array(data['t']) / 1000.0

    def extract_array(val):
        if isinstance(val, dict):
            arr = np.array(list(val.values())[0])
        else:
            arr = np.array(val)
        return arr[:, 0] if arr.ndim > 1 else arr

    panels = [
        ('Ca_ER_CICR',    ax_caer,   'b-',  "Ca_ER"),
        ('IP3_CICR',      ax_ip3,    'r-',  "IP3"),
        ('S_CICR',        ax_S,      'g-',  "S (SERCA)"),
        ('ryr_gate_CICR', ax_rgate,  'm-',  "ReLU gate"),
        ('J_serca_CICR',  ax_Jserca, 'g--', "J_serca"),
        ('J_leak_CICR',   ax_Jleak,  'k--', "J_leak"),
        ('J_cicr_CICR',   ax_Jcicr,  'm--', "J_cicr"),
        ('h_CICR',        ax_h,      'b--', "h (inact)"),
        ('m_inf_CICR',    ax_minf,   '-',   "m_inf"),
        ('n_inf_CICR',    ax_ninf,   '-',   "n_inf"),
        ('Nves_CICR',     ax_nves,   '-',   "Nves (frac)"),
    ]
    colors_override = {'m_inf_CICR': 'orange', 'n_inf_CICR': 'cyan', 'Nves_CICR': 'brown'}

    for key, ax, style, ylabel in panels:
        if key in data:
            y = extract_array(data[key])
            color = colors_override.get(key)
            if color:
                ax.plot(t, y, style, color=color, lw=1.5)
            else:
                ax.plot(t, y, style, lw=1.5)
        if is_left:
            ax.set_ylabel(ylabel, fontsize=10)

    ax_caer.set_title(f"CICR Diagnostics: {title}", fontsize=12)
    ax_nves.set_xlabel("Time (s)", fontsize=12)

    all_axes = [ax_caer, ax_ip3, ax_S, ax_rgate, ax_Jserca, ax_Jleak,
                ax_Jcicr, ax_h, ax_minf, ax_ninf, ax_nves]
    for ax in all_axes:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlim(x_lim if x_lim is not None else (-2.0, 50.0))


def main():
    import argparse
    import glob

    parser = argparse.ArgumentParser(description="Plot detailed traces")
    parser.add_argument("--pair-dir",      required=True,
                        help="Pair directory containing protocol sub-directories")
    parser.add_argument("--out-base-dir",  required=True,
                        help="Base directory to save figures")
    parser.add_argument("--sim-config",    default=None,
                        help="Path to prefire_simulation_config.json for calibration. "
                             "If omitted, auto-detected from --pair-dir.")
    # Threshold params — only needed when the pkl lacks theta_d/theta_p AND
    # no --sim-config is available for a calibration run.
    parser.add_argument("--a00", type=float, default=None)
    parser.add_argument("--a01", type=float, default=None)
    parser.add_argument("--a10", type=float, default=None)
    parser.add_argument("--a11", type=float, default=None)
    parser.add_argument("--tau-effca", type=float, default=None,
                        dest="tau_effca",
                        help="tau_effca_GB_GluSynapse value used in the calibration run")
    args = parser.parse_args()

    pair_name = os.path.basename(os.path.normpath(args.pair_dir))

    # --- Locate calibration sim config ---
    sim_config = args.sim_config
    if sim_config is None:
        # Heuristic: look one level up for prefire_simulation_config.json
        candidate = os.path.join(args.pair_dir, "prefire_simulation_config.json")
        if os.path.isfile(candidate):
            sim_config = candidate
        else:
            # Try the first protocol sub-directory
            sub = next(iter(glob.glob(os.path.join(args.pair_dir, "*",
                                                    "prefire_simulation_config.json"))), None)
            if sub:
                sim_config = sub

    # --- Build fit_params for calibration ---
    fit_params = {}
    if args.tau_effca is not None:
        fit_params["tau_effca_GB_GluSynapse"] = args.tau_effca
    # Pass a-values through so c_post_finder's _set_local_params can set thresholds
    for name, val in [("a00", args.a00), ("a01", args.a01),
                      ("a10", args.a10), ("a11", args.a11)]:
        if val is not None:
            fit_params[name] = val

    # --- Run calibration once per pair (expensive) ---
    c_pre_dict  = None
    c_post_dict = None
    theta_d_dict = None
    theta_p_dict = None

    if sim_config is not None and os.path.isfile(sim_config):
        print(f"Running cpre/cpost calibration from: {sim_config}")
        try:
            c_pre_dict, c_post_dict = compute_thresholds_from_calibration(
                sim_config, fit_params if fit_params else None)

            # Derive thresholds if a-values provided
            if all(k in fit_params for k in ("a00", "a01", "a10", "a11")):
                theta_d_dict, theta_p_dict = derive_thresholds(
                    c_pre_dict, c_post_dict, fit_params)
                print(f"  Calibrated thresholds (first syn): "
                      f"theta_d={list(theta_d_dict.values())[0]:.6f}  "
                      f"theta_p={list(theta_p_dict.values())[0]:.6f}")
        except Exception as exc:
            print(f"  WARNING: calibration failed ({exc}). "
                  f"Falling back to values in pkl.")
    else:
        print("No sim-config found — thresholds will be read from pkl (theta_d_GB/theta_p_GB).")

    # First-synapse scalar thresholds for plotting
    theta_d_first = list(theta_d_dict.values())[0] if theta_d_dict else None
    theta_p_first = list(theta_p_dict.values())[0] if theta_p_dict else None

    # --- Find and plot every protocol pkl ---
    pkl_files = glob.glob(os.path.join(args.pair_dir, "*", "simulation_traces.pkl"))
    if not pkl_files:
        print(f"Cannot find traces in {args.pair_dir}")
        return

    for pkl_path in pkl_files:
        protocol_dir = os.path.basename(os.path.dirname(pkl_path))
        d = load_data(pkl_path)

        freq, delay = (protocol_dir.split("_", 1) + ["unknown"])[:2] \
            if "_" in protocol_dir else (protocol_dir, "unknown")

        out_dir = os.path.join(args.out_base_dir, freq, delay)
        os.makedirs(out_dir, exist_ok=True)

        base      = os.path.join(out_dir, f"{pair_name}_traces")
        plot_title = f"{pair_name} - {protocol_dir}"

        # Plot 1: full view
        fig, axes = plt.subplots(3, 1, figsize=(6, 6), sharex=False)
        plot_single_protocol(axes[0], axes[1], axes[2], d, plot_title,
                             theta_d_first=theta_d_first,
                             theta_p_first=theta_p_first,
                             is_left=True)
        plt.tight_layout()
        plt.savefig(f"{base}.png", dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved {base}.png")

        # Plot 2: zoomed 0–3 s
        fig2, axes2 = plt.subplots(3, 1, figsize=(6, 6), sharex=False)
        plot_single_protocol(axes2[0], axes2[1], axes2[2], d, plot_title,
                             theta_d_first=theta_d_first,
                             theta_p_first=theta_p_first,
                             is_left=True, x_lim=(0, 3.0))
        plt.tight_layout()
        plt.savefig(f"{base}_zoomed.png", dpi=300, bbox_inches='tight')
        plt.close(fig2)
        print(f"Saved {base}_zoomed.png")


if __name__ == '__main__':
    main()
