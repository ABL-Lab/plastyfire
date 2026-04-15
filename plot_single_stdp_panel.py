#!/usr/bin/env python3
"""Plot a single STDP experiment summary panel from one simulation pickle.

The output is styled after the reference panel supplied by the user:
  - top-left: induction protocol schematic (pre/post spike timing)
  - top-right: baseline vs long-term average EPSP traces with a scale bar
  - middle: EPSP amplitudes across time, plus baseline / induction / long-term markers
  - bottom: rho_GB traces across the protocol for all synapses and their mean

The script uses the existing Experiment helper to extract EPSPs from the raw
voltage trace, so the measured amplitudes are consistent with the rest of the repo.
"""

import argparse
import os
import pickle
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np

REPO_DIR = os.path.dirname(os.path.abspath(__file__))
if REPO_DIR not in sys.path:
    sys.path.insert(0, REPO_DIR)

from plastyfire.ephysutils import Experiment


DEFAULT_SIM_PATH = os.path.join(
    REPO_DIR,
    "full_trace_results",
    "CHINDEMI_PARAMS",
    "184033-189853",
    "10Hz_10ms",
    "simulation_traces.pkl",
)
DEFAULT_OUTPUT = os.path.join(REPO_DIR, "figures", "stdp_panel_184033_189853_20Hz_10ms.png")
DEFAULT_REFERENCE_OUTPUT = os.path.join(
    REPO_DIR,
    "figures",
    "stdp_panel_184033_189853_20Hz_10ms_reference_layout.png",
)


def infer_protocol_from_path(sim_path):
    freq_dt = os.path.basename(os.path.dirname(sim_path))
    freq_hz_str, dt_ms_str = freq_dt.split("_")
    freq_hz = float(freq_hz_str.replace("Hz", ""))
    dt_ms = float(dt_ms_str.replace("ms", ""))
    pair_name = os.path.basename(os.path.dirname(os.path.dirname(sim_path)))
    return pair_name, freq_hz, dt_ms


def make_spike_train_trace(freq_hz, dt_ms, n_pairs=5):
    period_ms = 1000.0 / freq_hz
    t = np.linspace(0.0, n_pairs * period_ms + 40.0, 2000)
    pre_spikes = np.arange(n_pairs) * period_ms + 10.0
    post_spikes = pre_spikes + dt_ms

    def spike_shape(time_axis, spike_times):
        trace = np.zeros_like(time_axis)
        for spike_time in spike_times:
            rise = np.exp(-((time_axis - spike_time) / 0.9) ** 2)
            reset = np.exp(-((time_axis - (spike_time + 2.5)) / 2.5) ** 2)
            trace += 1.0 * rise - 0.18 * reset
        return trace

    pre_trace = spike_shape(t, pre_spikes)
    post_trace = spike_shape(t, post_spikes)
    return t, pre_trace, post_trace


def add_scale_bar(ax, x0, y0, dx, dy, x_label, y_label, color="k", lw=2.0):
    ax.plot([x0, x0 + dx], [y0, y0], color=color, lw=lw, solid_capstyle="butt")
    ax.plot([x0, x0], [y0, y0 + dy], color=color, lw=lw, solid_capstyle="butt")
    ax.text(x0 + dx / 2.0, y0 - 0.06 * dy, x_label, ha="center", va="top", fontsize=8)
    ax.text(x0 - 0.06 * dx, y0 + dy / 2.0, y_label, ha="right", va="center", fontsize=8)


def extract_synapse_matrix(trace_data, n_timepoints):
    if isinstance(trace_data, dict):
        matrix = np.array(list(trace_data.values()), dtype=float)
    else:
        matrix = np.asarray(trace_data, dtype=float)
        if matrix.ndim == 1:
            matrix = matrix.reshape(1, -1)
        elif matrix.ndim == 2:
            if matrix.shape[0] == n_timepoints and matrix.shape[1] != n_timepoints:
                matrix = matrix.T
            elif matrix.shape[1] != n_timepoints and matrix.shape[0] > matrix.shape[1]:
                matrix = matrix.T
    return matrix


def map_time_for_display(t, induction_end):
    mapped = np.asarray(t, dtype=float).copy()
    baseline_mask = mapped <= 0.0
    induction_mask = (mapped > 0.0) & (mapped <= induction_end)
    post_mask = mapped > induction_end

    if np.any(baseline_mask):
        baseline_span = max(abs(float(np.min(mapped[baseline_mask]))), 1e-9)
        mapped[baseline_mask] = 10.0 * mapped[baseline_mask] / baseline_span

    induction_span = max(float(induction_end), 1e-9)
    if np.any(induction_mask):
        mapped[induction_mask] = 2.0 * mapped[induction_mask] / induction_span

    if np.any(post_mask):
        post_end = float(np.max(mapped[post_mask]))
        post_span = max(post_end - 2.0, 1e-9)
        mapped[post_mask] = 2.0 + 38.0 * (mapped[post_mask] - 2.0) / post_span

    return mapped


def build_panel(sim_path, output_path, trace_average_n=10, summary_n=10, reference_layout=False):
    with open(sim_path, "rb") as handle:
        sim_data = pickle.load(handle)

    _, freq_hz, dt_ms = infer_protocol_from_path(sim_path)
    experiment = Experiment(data=sim_data, c01duration=4.0, c02duration=4.0, period=4.0)

    c01_spikes = np.asarray(experiment.cxspikes["C01"], dtype=float)
    c02_spikes = np.asarray(experiment.cxspikes["C02"], dtype=float)
    c01_times = experiment.normalize_time(c01_spikes)
    c02_times = experiment.normalize_time(c02_spikes)
    c01_epsp = np.asarray(experiment.epsp["C01"], dtype=float)
    c02_epsp = np.asarray(experiment.epsp["C02"], dtype=float)

    induction_spikes = np.asarray(sim_data["prespikes"][len(c01_spikes):-len(c02_spikes)], dtype=float)
    induction_start_real = float(experiment.normalize_time(induction_spikes[0])) if induction_spikes.size else 0.0
    induction_end_real = float(experiment.normalize_time(induction_spikes[-1])) if induction_spikes.size else 0.0

    rho_t = np.asarray(sim_data["t"], dtype=float)
    rho_traces = extract_synapse_matrix(sim_data["rho_GB"], len(rho_t))
    rho_time = experiment.normalize_time(rho_t)

    if reference_layout:
        c01_times = np.linspace(-10.0, 0.0, len(c01_epsp), endpoint=False) + 10.0 / (2 * len(c01_epsp))
        c02_times = np.linspace(2.0, 40.0, len(c02_epsp), endpoint=False) + 38.0 / (2 * len(c02_epsp))
        induction_left = 0.0
        induction_right = 2.0
        rho_time_plot = map_time_for_display(rho_time, induction_end_real)
    else:
        display_induction_width = max(induction_end_real - induction_start_real, 0.28)
        induction_center = 0.5 * (induction_start_real + induction_end_real)
        induction_left = induction_center - display_induction_width / 2.0
        induction_right = induction_center + display_induction_width / 2.0
        rho_time_plot = rho_time

    baseline_mean = float(np.mean(c01_epsp[-summary_n:]))
    long_term_mean = float(np.mean(c02_epsp[-summary_n:]))
    baseline_x = (float(c01_times[-summary_n]), float(c01_times[-1]))
    long_term_x = (float(c02_times[-summary_n]), float(c02_times[-1]))

    baseline_trace = np.mean(experiment.cxtrace["C01"][-trace_average_n:], axis=0)
    long_term_trace = np.mean(experiment.cxtrace["C02"][-trace_average_n:], axis=0)
    tdense = np.asarray(experiment.cxtrace["t"], dtype=float)
    trace_t = tdense - tdense[0]
    baseline_trace = baseline_trace - baseline_trace[0]
    long_term_trace = long_term_trace - long_term_trace[0]

    max_rho_points = 3000
    rho_stride = max(1, int(np.ceil(rho_time_plot.size / max_rho_points)))
    rho_time_plot = rho_time_plot[::rho_stride]
    rho_traces_plot = rho_traces[:, ::rho_stride]
    rho_mean = np.mean(rho_traces_plot, axis=0)

    if reference_layout:
        fig = plt.figure(figsize=(3.45, 4.05), dpi=220)
        grid = fig.add_gridspec(3, 1, height_ratios=[0.88, 2.08, 0.92], hspace=0.05)
    else:
        fig = plt.figure(figsize=(3.7, 4.3), dpi=200)
        grid = fig.add_gridspec(3, 1, height_ratios=[0.95, 2.25, 0.95], hspace=0.06)
    ax_top = fig.add_subplot(grid[0])
    ax_main = fig.add_subplot(grid[1])
    ax_rho = fig.add_subplot(grid[2], sharex=ax_main)

    fig.text(0.03 if reference_layout else 0.04, 0.972, "b", fontsize=15 if reference_layout else 14)

    ax_top.set_axis_off()

    ax_proto = ax_top.inset_axes([0.00, 0.10, 0.66, 0.88] if reference_layout else [0.00, 0.16, 0.63, 0.80])
    proto_t, proto_pre, proto_post = make_spike_train_trace(freq_hz, dt_ms)
    proto_pre_y = 1.33 if reference_layout else 1.32
    proto_post_y = 0.43 if reference_layout else 0.42
    proto_lw = 0.85 if reference_layout else 0.8
    proto_fs = 8.5 if reference_layout else 8
    ax_proto.plot(proto_t, proto_pre_y + 0.42 * proto_pre, color="k", lw=proto_lw)
    ax_proto.plot(proto_t, proto_post_y + 0.42 * proto_post, color="k", lw=proto_lw, ls="--")
    ax_proto.text(-20 if reference_layout else -18, proto_pre_y, "pre", ha="right", va="center", fontsize=proto_fs)
    ax_proto.text(-20 if reference_layout else -18, proto_post_y, "post", ha="right", va="center", fontsize=proto_fs)
    ax_proto.text(18, -0.06 if reference_layout else -0.18,
                  f"Frequency = {int(freq_hz)} Hz\n$\Delta t$ = {int(dt_ms)} ms",
                  ha="left", va="top", fontsize=proto_fs)
    ax_proto.set_xlim(-22 if reference_layout else -20, proto_t.max() + 8)
    ax_proto.set_ylim(-0.30 if reference_layout else -0.35, 1.95)
    ax_proto.axis("off")

    ax_trace = ax_top.inset_axes([0.69, 0.06, 0.29, 0.90] if reference_layout else [0.68, 0.08, 0.30, 0.86])
    if reference_layout:
        ax_trace.plot(trace_t, long_term_trace, color="#f39c35", lw=0.95)
        ax_trace.plot(trace_t, baseline_trace, color="#69a8e5", lw=0.95)
    else:
        ax_trace.plot(trace_t, baseline_trace, color="#69a8e5", lw=1.0)
        ax_trace.plot(trace_t, long_term_trace, color="#f39c35", lw=1.0)
    trace_peak = max(float(np.max(baseline_trace)), float(np.max(long_term_trace)))
    ax_trace.set_xlim(0, 60)
    ax_trace.set_ylim(-0.03 if reference_layout else -0.02, max(trace_peak * 1.08, 0.60 if reference_layout else 0.55))
    if reference_layout:
        y0 = ax_trace.get_ylim()[1] * 0.72
        add_scale_bar(ax_trace, 39, y0, 20, 0.5, "20 ms", "0.5 mV", lw=1.7)
    else:
        add_scale_bar(ax_trace, 39, ax_trace.get_ylim()[1] * 0.73, 20, 0.5, "20 ms", "0.5 mV")
    ax_trace.axis("off")

    scatter_kwargs = dict(s=5.2 if reference_layout else 4.5, color="#d3d3d3", edgecolors="none", alpha=0.95)
    ax_main.scatter(c01_times, c01_epsp, **scatter_kwargs)
    ax_main.scatter(c02_times, c02_epsp, **scatter_kwargs)
    ax_main.axvspan(induction_left, induction_right, color="#d9d9d9", zorder=0)
    ax_main.hlines(baseline_mean, baseline_x[0], baseline_x[1], color="#69a8e5", lw=2.0, zorder=3)
    ax_main.hlines(long_term_mean, long_term_x[0], long_term_x[1], color="#f39c35", lw=2.0, zorder=3)

    handles = [
        Line2D([0], [0], color="#69a8e5", lw=2.0, label="baseline"),
        Line2D([0], [0], color="#f39c35", lw=2.0, label="long term"),
        Rectangle((0, 0), 1, 1, facecolor="#d9d9d9", edgecolor="#d9d9d9", label="induction"),
    ]
    ax_main.legend(handles=handles, frameon=False, loc="lower right", fontsize=8.4 if reference_layout else 8)

    ax_main.set_ylabel("EPSP amplitude (mV)", fontsize=9.4 if reference_layout else 9)
    ax_main.tick_params(axis="y", labelsize=8.6 if reference_layout else 8, length=2.5, pad=1.2 if reference_layout else 1.5)
    ax_main.tick_params(axis="x", labelbottom=False, length=2.5)
    ax_main.spines["top"].set_visible(False)
    ax_main.spines["right"].set_visible(False)
    ymax = max(float(np.max(c01_epsp)), float(np.max(c02_epsp)), baseline_mean, long_term_mean)
    if reference_layout:
        xlim = (-12, 42)
        xticks = [-10, 0, 10, 20, 30, 40]
        ax_main.set_xlim(*xlim)
        ax_main.set_xticks(xticks)
        ax_main.set_ylim(0.0, max(3.05, ymax * 1.03))
    else:
        xlim = (min(float(c01_times[0]), -0.6), float(c02_times[-1]) + 0.25)
        xticks = None
        ax_main.set_xlim(*xlim)
        ax_main.set_ylim(0.0, ymax * 1.06)

    for trace in rho_traces_plot:
        ax_rho.plot(rho_time_plot, trace, color="#bdbdbd", lw=0.7, alpha=0.45, zorder=1)
    ax_rho.plot(rho_time_plot, rho_mean, color="#2b2b2b", lw=1.4, zorder=3)
    ax_rho.axhline(0.5, color="#666666", lw=0.9, ls="--", zorder=2)
    ax_rho.axvspan(induction_left, induction_right, color="#d9d9d9", zorder=0)
    ax_rho.set_ylabel("rho_GB", fontsize=9.4 if reference_layout else 9)
    ax_rho.set_xlabel("Time (min)", fontsize=9.4 if reference_layout else 9)
    ax_rho.set_ylim(-0.03, 1.03)
    ax_rho.set_yticks([0.0, 0.5, 1.0])
    ax_rho.tick_params(axis="both", labelsize=8.6 if reference_layout else 8, length=2.5, pad=1.2 if reference_layout else 1.5)
    ax_rho.spines["top"].set_visible(False)
    ax_rho.spines["right"].set_visible(False)
    ax_rho.set_xlim(*xlim)
    if xticks is not None:
        ax_rho.set_xticks(xticks)
    rho_handles = [
        Line2D([0], [0], color="#2b2b2b", lw=1.4, label="mean rho"),
        Line2D([0], [0], color="#666666", lw=0.9, ls="--", label="rho = 0.5"),
    ]
    ax_rho.legend(handles=rho_handles, frameon=False, loc="upper right", fontsize=7.8 if reference_layout else 7.5)
    ax_rho.text(0.02, 0.86, f"n = {rho_traces.shape[0]} syn", transform=ax_rho.transAxes,
                ha="left", va="top", fontsize=7.8 if reference_layout else 7.5)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Plot a single STDP summary panel from one simulation pickle.")
    parser.add_argument("--sim-path", default=DEFAULT_SIM_PATH, help="Path to simulation_traces.pkl")
    parser.add_argument("--output", default=None, help="Output PNG path")
    parser.add_argument("--trace-average-n", type=int, default=10,
                        help="Number of late EPSP traces to average for the top-right inset")
    parser.add_argument("--summary-n", type=int, default=10,
                        help="Number of late EPSPs used for the baseline / long-term bars")
    parser.add_argument("--reference-layout", action="store_true",
                        help="Remap the display timeline to match the supplied reference panel layout")
    args = parser.parse_args()
    output_path = args.output
    if output_path is None:
        output_path = DEFAULT_REFERENCE_OUTPUT if args.reference_layout else DEFAULT_OUTPUT
    build_panel(os.path.abspath(args.sim_path), os.path.abspath(output_path),
                trace_average_n=args.trace_average_n,
                summary_n=args.summary_n,
                reference_layout=args.reference_layout)
    print(f"Saved {os.path.abspath(output_path)}")


if __name__ == "__main__":
    main()
