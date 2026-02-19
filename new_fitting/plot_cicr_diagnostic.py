#!/usr/bin/env python
"""
Publication-ready CICR diagnostic: 10Hz_10ms vs 10Hz_-10ms side-by-side.

Plots 5 rows × 2 columns:
  Rows: cai, effcai, rho, mGluR bucket, ER state
  Cols: 10Hz_10ms | 10Hz_-10ms

Usage:
    # Default params:
    python plot_cicr_diagnostic.py

    # Optimized params from JSON:
    python plot_cicr_diagnostic.py --json result_nn_2proto.json

    # Manually override specific params:
    python plot_cicr_diagnostic.py --json result_nn_2proto.json --pair-idx 5 --syn-idx 2
"""
import argparse
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numba import njit
from fit_params_CICR import (
    preload_all_data, _compute_effcai_from_cai,
    DEFAULT_PARAMS, PARAM_NAMES, TAU_EFFCA_GB, MIN_CA_CR,
)

# ── Publication style ──
plt.rcParams.update({
    'font.size': 14,
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 10,
    'lines.linewidth': 1.2,
    'axes.linewidth': 1.2,
    'xtick.major.width': 1.0,
    'ytick.major.width': 1.0,
    'figure.dpi': 150,
    'savefig.dpi': 200,
    'font.family': 'sans-serif',
})

COLORS = {
    'raw': '#4C72B0',       # steel blue
    'cicr': '#C44E52',      # muted red
    'theta_d': '#DD8452',   # orange
    'theta_p': '#55A868',   # green
    'mGluR': '#8172B3',     # purple
    'thresh': '#937860',    # tan
    'er': '#64B5CD',        # light blue
    'cicr_ca': '#C44E52',   # red
}


@njit(cache=True)
def _compute_rho_trace(effcai, t, theta_d, theta_p, gamma_d, gamma_p, rho0):
    """Euler-integrate rho ODE for one synapse, return FULL trace."""
    n = len(t)
    rho = np.zeros(n, dtype=np.float64)
    rho[0] = rho0
    inv_tau = 1.0 / (70.0 * 1000.0)
    rho_star = 0.5
    r = rho0
    for i in range(n - 1):
        dt = t[i + 1] - t[i]
        dep = 1.0 if effcai[i] > theta_d else 0.0
        pot = 1.0 if effcai[i] > theta_p else 0.0
        drho = (-r * (1.0 - r) * (rho_star - r)
                + pot * gamma_p * (1.0 - r)
                - dep * gamma_d * r) * inv_tau
        r_new = r + dt * drho
        if r_new < 0.0:
            r_new = 0.0
        elif r_new > 1.0:
            r_new = 1.0
        r = r_new
        rho[i + 1] = r
    return rho


@njit(cache=True)
def _apply_cicr_debug(cai_1syn, t, c_pre_val, c_post_val,
                      tau_IP3, IP3_threshold, a_evt, b_evt,
                      a_er, b_er, g_serca, g_release,
                      sat_cap, sigmoid_slope, K_act):
    """Single-synapse Primed CICR with full debug traces.

    Mirrors _apply_cicr_batch from fit_params_CICR.py but returns intermediate
    traces (Priming integrator, ER Ca, CICR Ca contribution) for visualization.
    """
    K_serca = 0.1
    fc = 0.83
    fe = 0.17
    caAvg = 0.0017
    MAX_DT_SUB = 0.1

    n_time = len(cai_1syn)
    cai_total = np.copy(cai_1syn)
    priming_trace = np.zeros(n_time, dtype=np.float64)
    ca_er_trace = np.zeros(n_time, dtype=np.float64)
    ca_cicr_trace = np.zeros(n_time, dtype=np.float64)

    ca_event_thresh = a_evt * c_pre_val + b_evt * c_post_val
    gamma_er = a_er * c_pre_val + b_er * c_post_val
    ca_er = (caAvg - fc * cai_1syn[0]) / fe
    if ca_er < 0.0:
        ca_er = 0.0
    Priming_val = 0.0
    ca_cicr = 0.0
    was_above = 0
    event_peak = 0.0  # peak cai while above ca_event_thresh

    priming_trace[0] = Priming_val
    ca_er_trace[0] = ca_er
    ca_cicr_trace[0] = ca_cicr

    for i in range(n_time - 1):
        dt_full = t[i + 1] - t[i]
        if dt_full <= 0.0:
            cai_total[i + 1] = cai_1syn[i + 1] + ca_cicr
            priming_trace[i + 1] = Priming_val
            ca_er_trace[i + 1] = ca_er
            ca_cicr_trace[i + 1] = ca_cicr
            continue

        # --- Option A: peak-tracking event detection ---
        is_above = 1 if cai_1syn[i] > ca_event_thresh else 0
        event_amp = 0.0
        if is_above == 1 and was_above == 0:
            event_peak = cai_1syn[i]
        elif is_above == 1:
            if cai_1syn[i] > event_peak:
                event_peak = cai_1syn[i]
        elif was_above == 1:
            raw_amp = event_peak - ca_event_thresh
            event_amp = (raw_amp / sat_cap) if raw_amp < sat_cap else 1.0
            event_peak = 0.0
        was_above = is_above

        # --- IP3 priming integrator ---
        Priming_val = Priming_val * np.exp(-dt_full / tau_IP3) + event_amp

        # --- Sigmoid coincidence gate ---
        arg = -sigmoid_slope * (Priming_val - IP3_threshold)
        if arg > 500.0:
            P_gate = 0.0
        elif arg < -500.0:
            P_gate = 1.0
        else:
            P_gate = 1.0 / (1.0 + np.exp(arg))

        n_sub = max(1, int(np.ceil(dt_full / MAX_DT_SUB)))
        dt_sub = dt_full / n_sub

        for _ in range(n_sub):
            ca_cyt = cai_1syn[i] + ca_cicr

            ca_cicr_um = 1000.0 * ca_cicr
            if ca_cicr > gamma_er:
                J_serca_raw = g_serca * (ca_cicr_um**2) / (K_serca**2 + ca_cicr_um**2)
                J_serca = min(J_serca_raw, ca_cicr / dt_sub)
            else:
                J_serca = 0.0

            if ca_er > ca_cyt:
                J_release = g_release * P_gate * (ca_cyt / (ca_cyt + K_act)) * (ca_er - ca_cyt)
            else:
                J_release = 0.0

            J_net = J_release - J_serca
            ca_er = ca_er - (fc / fe) * J_net * dt_sub
            if ca_er < 0.0:
                ca_er = 0.0
            ca_cicr = ca_cicr + J_net * dt_sub
            if ca_cicr < 0.0:
                ca_cicr = 0.0

        cai_total[i + 1] = cai_1syn[i + 1] + ca_cicr
        priming_trace[i + 1] = Priming_val
        ca_er_trace[i + 1] = ca_er
        ca_cicr_trace[i + 1] = ca_cicr

    return cai_total, priming_trace, ca_er_trace, ca_cicr_trace


def params_from_json(json_path):
    """Load parameter dict from result JSON file."""
    with open(json_path) as f:
        data = json.load(f)

    # Find the best result (any method)
    for method_name, method_data in data.get("results", {}).items():
        if "best_params" in method_data:
            param_vec = method_data["best_params"]
            dp = {}
            for i, name in enumerate(PARAM_NAMES):
                dp[name] = param_vec[i]
            print(f"Loaded optimized params from '{method_name}' in {json_path}")
            return dp

    raise ValueError(f"No best_params found in {json_path}")


def compute_single(pd_item, syn_idx, dp):
    """Run full diagnostic for one pair/synapse, return dict of traces."""
    t = pd_item["t"]
    cai = pd_item["cai"]
    s = syn_idx

    cp = pd_item["c_pre"][s]
    cq = pd_item["c_post"][s]

    gamma_d = dp["gamma_d_GB_GluSynapse"]
    gamma_p = dp["gamma_p_GB_GluSynapse"]

    if pd_item["is_apical"][s]:
        theta_d = dp["a20"] * cp + dp["a21"] * cq
        theta_p = dp["a30"] * cp + dp["a31"] * cq
    else:
        theta_d = dp["a00"] * cp + dp["a01"] * cq
        theta_p = dp["a10"] * cp + dp["a11"] * cq

    cai_1syn = cai[s]
    cai_total, priming_trace, ca_er_trace, ca_cicr_trace = _apply_cicr_debug(
        cai_1syn, t, cp, cq,
        dp["tau_IP3"], dp["IP3_threshold"], dp["a_evt"], dp["b_evt"],
        dp["a_er"], dp["b_er"], dp["g_serca_cicr"], dp["g_release_cicr"],
        dp["sat_cap"], dp["sigmoid_slope"], dp["K_act"],
    )

    effcai_no = _compute_effcai_from_cai(cai_1syn, t, TAU_EFFCA_GB, MIN_CA_CR)
    effcai_ci = _compute_effcai_from_cai(cai_total, t, TAU_EFFCA_GB, MIN_CA_CR)

    rho_no = _compute_rho_trace(effcai_no, t, theta_d, theta_p,
                                 gamma_d, gamma_p, pd_item["rho0"][s])
    rho_ci = _compute_rho_trace(effcai_ci, t, theta_d, theta_p,
                                 gamma_d, gamma_p, pd_item["rho0"][s])

    return {
        "t": t, "cai_raw": cai_1syn, "cai_total": cai_total,
        "effcai_no": effcai_no, "effcai_ci": effcai_ci,
        "rho_no": rho_no, "rho_ci": rho_ci,
        "priming": priming_trace, "ca_er": ca_er_trace, "ca_cicr": ca_cicr_trace,
        "theta_d": theta_d, "theta_p": theta_p,
        "cp": cp, "cq": cq,
        "ca_event_thresh": dp["a_evt"] * cp + dp["b_evt"] * cq,
        "IP3_threshold": dp["IP3_threshold"],
        "gamma_er": dp["a_er"] * cp + dp["b_er"] * cq,
    }


def print_stats(label, d):
    """Print diagnostic stats for one protocol."""
    cai_diff = d["cai_total"] - d["cai_raw"]
    print(f"\n  --- {label} ---")
    print(f"  n_timesteps:         {len(d['t'])}")
    print(f"  cai_CR range:        [{d['cai_raw'].min():.6f}, {d['cai_raw'].max():.6f}] mM")
    print(f"  CICR contribution:   [{cai_diff.min():.6f}, {cai_diff.max():.6f}] mM")
    n_above = np.sum(d["cai_raw"] > d["ca_event_thresh"])
    print(f"  # cai > thresh:      {n_above}/{len(d['t'])}")
    print(f"  Priming max:         {d['priming'].max():.6f}  (thresh={d['IP3_threshold']})")
    cicr_active = d["priming"].max() > d["IP3_threshold"]
    print(f"  CICR gate activated: {'YES' if cicr_active else 'NO'}")
    print(f"  ca_er range:         [{d['ca_er'].min():.6f}, {d['ca_er'].max():.6f}] mM")
    print(f"  ca_cicr max:         {d['ca_cicr'].max():.6f} mM")
    print(f"  effcai max (no/ci):  {d['effcai_no'].max():.6f} / {d['effcai_ci'].max():.6f}")
    print(f"  theta_d={d['theta_d']:.4f}, theta_p={d['theta_p']:.4f}")
    print(f"  rho final (no/ci):   {d['rho_no'][-1]:.6f} / {d['rho_ci'][-1]:.6f}")


def _pick_random_pair_syn(protocol_data, protocols, seed=None):
    """Pick a random pair/synapse shared across all requested protocols.

    Selects a pair index that exists in ALL protocols (same index, not
    necessarily same pair), then a random synapse from that pair.
    Prefers synapses whose cai_CR peak exceeds 1 uM (0.001 mM) so the
    diagnostic shows non-trivial dynamics.
    """
    rng = np.random.default_rng(seed)

    # Use the first protocol to enumerate valid pair indices
    ref_proto = protocols[0]
    ref_pairs = protocol_data[ref_proto]
    n_pairs = len(ref_pairs)

    # Filter to pairs that exist in all protocols at the same index
    valid_pair_idxs = [
        i for i in range(n_pairs)
        if all(i < len(protocol_data[p]) for p in protocols)
    ]
    if not valid_pair_idxs:
        raise RuntimeError("No pair index is valid across all protocols")

    pair_idx = int(rng.choice(valid_pair_idxs))
    pd_item = ref_pairs[pair_idx]

    # Among synapses in this pair, prefer ones with a healthy cai peak
    cai = pd_item["cai"]          # (n_syn, n_time)
    cai_peaks = cai.max(axis=1)   # (n_syn,)
    good_syns = np.where(cai_peaks > 0.001)[0]   # > 1 uM
    if len(good_syns) == 0:
        good_syns = np.arange(cai.shape[0])        # fall back: take any
    syn_idx = int(rng.choice(good_syns))

    return pair_idx, syn_idx


def main():
    parser = argparse.ArgumentParser(description="CICR diagnostic: 10Hz_10ms vs 10Hz_-10ms")
    parser.add_argument("--pair-idx", type=int, default=None)
    parser.add_argument("--syn-idx", type=int, default=None)
    parser.add_argument("--random-pair", action="store_true",
                        help="Randomly pick a pair/synapse (prefers cai peak > 1 uM)")
    parser.add_argument("--seed", type=int, default=None,
                        help="RNG seed for --random-pair (default: random)")
    parser.add_argument("--json", type=str, default=None,
                        help="Load optimized params from result JSON file")
    parser.add_argument("--output", type=str, default="cicr_diagnostic.png")
    args = parser.parse_args()

    protocols = ["10Hz_10ms", "10Hz_-10ms"]

    # Load params
    if args.json:
        dp = params_from_json(args.json)
        param_label = f"Optimized ({args.json})"
    else:
        dp = dict(DEFAULT_PARAMS)
        param_label = "Default params"

    print(f"Using: {param_label}")
    print(f"  Key CICR params: tau_IP3={dp['tau_IP3']:.1f}, "
          f"IP3_threshold={dp['IP3_threshold']:.2f}, "
          f"sat_cap={dp['sat_cap']:.4f}, sigmoid_slope={dp['sigmoid_slope']:.2f}, "
          f"K_act={dp['K_act']:.5f}")
    print(f"  ca_event_thresh = a_evt*Cpre + b_evt*Cpost: a_evt={dp['a_evt']:.5f}, b_evt={dp['b_evt']:.5f}")
    print(f"  gamma_er         = a_er*Cpre  + b_er*Cpost:  a_er={dp['a_er']:.5f},  b_er={dp['b_er']:.5f}")
    print(f"  g_serca={dp['g_serca_cicr']:.4f}, g_release={dp['g_release_cicr']:.4f}")

    print("Loading data...")
    protocol_data = preload_all_data()

    # ── Resolve pair/synapse index ──
    if args.random_pair:
        pair_idx, syn_idx = _pick_random_pair_syn(protocol_data, protocols, seed=args.seed)
        print(f"Randomly selected: pair_idx={pair_idx}, syn_idx={syn_idx}  "
              f"(seed={args.seed}, rerun with --pair-idx {pair_idx} --syn-idx {syn_idx})")
    else:
        pair_idx = args.pair_idx if args.pair_idx is not None else 0
        syn_idx  = args.syn_idx  if args.syn_idx  is not None else 0

    results = {}
    for proto in protocols:
        pairs = protocol_data.get(proto, [])
        if not pairs:
            print(f"No pairs for {proto}")
            return
        if pair_idx >= len(pairs):
            print(f"pair-idx {pair_idx} out of range for {proto} (max {len(pairs)-1})")
            return
        pd_item = pairs[pair_idx]
        n_syn = pd_item["cai"].shape[0]
        if syn_idx >= n_syn:
            print(f"syn-idx {syn_idx} out of range (max {n_syn-1})")
            return

        print(f"Running {proto}...")
        results[proto] = compute_single(pd_item, syn_idx, dp)
        print_stats(proto, results[proto])

    # ── Plot: 5 rows × 1 column ──
    fig, axes = plt.subplots(5, 1, figsize=(10, 16), sharex=True)
    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])

    row_labels = [
        r'Ca$^{2+}$ ($\mu$M)',
        'effcai',
        r'$\rho$',
        'IP3 Priming Level',
        r'ER Ca$^{2+}$ ($\mu$M)',
    ]

    # Auto-detect stimulation window from first protocol's cai
    first_d = results[protocols[0]]
    t_s = first_d["t"] / 1000.0
    # Find last time cai is significantly above baseline
    cai_above = first_d["cai_raw"] > first_d["cai_raw"][0] * 2
    if np.any(cai_above):
        last_event_s = t_s[np.where(cai_above)[0][-1]]
        xlim_max = min(last_event_s * 1.5, t_s[-1])  # 50% padding
    else:
        xlim_max = t_s[-1]

    linestyles = ['-', '--']
    
    # Specific colors for each protocol
    # Proto 1 (Blue/Red theme)
    # Proto 2 (Green/Purple theme)
    PROTO_COLORS = [
        {
            'raw': '#1f77b4',       # Dark Blue
            'cicr': '#d62728',      # Dark Red
            'mGluR': '#9467bd',     # Dark Purple
            'er': '#17becf',        # Dark Cyan
            'cicr_ca': '#d62728',   # Dark Red
        },
        {
            'raw': '#2ca02c',       # Green
            'cicr': '#e377c2',      # Pink
            'mGluR': '#8c564b',     # Brown
            'er': '#bcbd22',        # Olive
            'cicr_ca': '#e377c2',   # Pink
        }
    ]
    
    for i, proto in enumerate(protocols):
        d = results[proto]
        ls = linestyles[i % len(linestyles)]
        # Cycle through color sets if more than 2 protocols
        p_colors = PROTO_COLORS[i % len(PROTO_COLORS)]
        
        t_s = d["t"] / 1000.0
        suffix = f" ({proto})"

        # Row 0: cai
        ax = axes[0]
        ax.plot(t_s, d["cai_total"] * 1000, color=p_colors['cicr'], ls=ls, alpha=0.8,
                linewidth=1.2, label=proto)
        if i == 0:
            ax.axhline(d["ca_event_thresh"] * 1000, color=COLORS['thresh'], ls=':',
                       lw=1.5, alpha=0.7, label='Event thresh')
            ax.set_title("CICR Protocols Comparison", fontsize=18, fontweight='bold')
            ax.set_ylabel(row_labels[0])

        # Row 1: effcai
        ax = axes[1]
        ax.plot(t_s, d["effcai_ci"], color=p_colors['cicr'], ls=ls, alpha=0.8,
                linewidth=1.2, label=proto)
        if i == 0:
            ax.axhline(d["theta_d"], color=COLORS['theta_d'], ls=':', lw=1.5,
                       alpha=0.8, label=fr'$\theta_d$ = {d["theta_d"]:.3f}')
            ax.axhline(d["theta_p"], color=COLORS['theta_p'], ls=':', lw=1.5,
                       alpha=0.8, label=fr'$\theta_p$ = {d["theta_p"]:.3f}')
            ax.set_ylabel(row_labels[1])

        # Row 2: rho
        ax = axes[2]
        ax.plot(t_s, d["rho_ci"], color=p_colors['cicr'], ls=ls, alpha=0.9,
                linewidth=1.5, label=proto)
        if i == 0:
            ax.axhline(0.5, color='gray', ls=':', lw=1.0, alpha=0.5, label=r'$\rho^*$=0.5')
            ax.set_ylim(-0.05, 1.05)
            ax.set_ylabel(row_labels[2])

        # Row 3: IP3 priming bucket
        ax = axes[3]
        ax.plot(t_s, d["priming"], color=p_colors['mGluR'], ls=ls, alpha=0.9,
                linewidth=1.2, label='IP3 priming' + suffix)
        if i == 0:
            ax.axhline(d["IP3_threshold"], color=COLORS['cicr'], ls=':', lw=1.5,
                       alpha=0.8, label=f'IP3 threshold = {d["IP3_threshold"]:.2f}')
            ax.set_ylabel(row_labels[3])

        # Row 4: ER state
        ax = axes[4]
        ax.plot(t_s, d["ca_er"] * 1000, color=p_colors['er'], ls=ls, alpha=0.9,
                linewidth=1.2, label=r'ca$_{\mathrm{ER}}$' + suffix)
        ax.plot(t_s, d["ca_cicr"] * 1000, color=p_colors['cicr_ca'], ls=ls, alpha=0.9,
                linewidth=1.2, label=r'ca$_{\mathrm{CICR}}$' + suffix)
        if i == 0:
            ax.axhline(d["gamma_er"] * 1000, color=COLORS['theta_d'], ls=':', lw=1.5,
                       alpha=0.7, label=fr'$\gamma_{{ER}}$ = {d["gamma_er"]*1000:.2f} µM')
            ax.set_ylabel(row_labels[4])
            ax.set_xlabel('Time (s)', fontsize=14)

    # Clean up
    for ax in axes:
        ax.set_xlim(0, xlim_max)
        ax.tick_params(direction='in', top=True, right=True)
        ax.spines['top'].set_visible(True)
        ax.spines['right'].set_visible(True)
        # Place legend outside or neatly inside
        ax.legend(loc='upper right', fontsize=8, ncol=2)

    plt.tight_layout(h_pad=0.8, w_pad=1.5)
    fig.suptitle(f'CICR Diagnostic — {param_label}\nPair #{pair_idx}, Synapse #{syn_idx}',
                 fontsize=18, fontweight='bold', y=1.02)
    plt.savefig(args.output, dpi=200, bbox_inches='tight', facecolor='white')
    print(f"\nSaved to {args.output}")


if __name__ == "__main__":
    main()
