"""
Plot effcai_GB traces of a bluecellulab (BCL) edges run together with the
per-synapse theta_d / theta_p thresholds that the same run used.

Context
-------
The current edges simulator (plastyfire/simulator_edges.py) writes
simulation_traces.pkl with cai_CR / shaft_cai / ica_NMDA / ica_VDCC / rho_GB,
but *not* effcai_GB.  effcai_GB is a pure leaky integral of cai_CR
(GluSynapse.mod line 361):

    effcai_GB' = -effcai_GB/tau_effca_GB + (cai_CR - min_ca_CR)

so it is reconstructed here exactly (exponential integrator on the uniform
grid), with the tau_effca_GB that the run actually used.

theta_d / theta_p are read back from the pool log the run wrote, which prints
one line per synapse:

    syn gid=299807469  loc=basal  cp=...  cq=...  theta_d=...  theta_p=...

Usage
-----
python plot_effcai_thetas.py \
    --pair-dir  .../simulations/180217-186276 \
    --out-dir   figures/effcai_thetas_180217-186276 \
    --tau-effca 278.3177658387
"""

import os
import re
import glob
import pickle
import argparse
from collections import OrderedDict

import numpy as np
from scipy.signal import lfilter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MIN_CA_CR_DEFAULT = 70e-6      # mM, GluSynapse.mod
THETA_D_COLOR = "#d62728"      # red
THETA_P_COLOR = "#2ca02c"      # green
DEFAULT_STYLE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "onerule.mplstyle")
TAU_EFFCA_DEFAULT = 200.0      # ms, GluSynapse.mod default (override per run!)

THETA_RE = re.compile(
    r"syn gid=(?P<gid>\d+)\s+loc=(?P<loc>\w+)\s+cp=(?P<cp>[-\d.eE+]+)\s+"
    r"cq=(?P<cq>[-\d.eE+]+)\s+theta_d=(?P<td>[-\d.eE+]+)\s+theta_p=(?P<tp>[-\d.eE+]+)")


def parse_thetas(proto_dir, log_glob, ref_mtime=None):
    """Return {gid: dict(loc, cp, cq, theta_d, theta_p)} from the run's pool log."""
    cands = []
    for path in glob.glob(os.path.join(proto_dir, log_glob)):
        with open(path, "r", errors="ignore") as fh:
            hits = THETA_RE.findall(fh.read())
        if hits:
            cands.append(path)
    if not cands:
        return {}, None
    # Prefer the log written closest in time to the traces file
    if ref_mtime is not None:
        cands.sort(key=lambda p: abs(os.path.getmtime(p) - ref_mtime))
    else:
        cands.sort(key=os.path.getmtime, reverse=True)
    log_path = cands[0]
    out = {}
    with open(log_path, "r", errors="ignore") as fh:
        for m in THETA_RE.finditer(fh.read()):
            out[int(m.group("gid"))] = dict(loc=m.group("loc"),
                                            cp=float(m.group("cp")),
                                            cq=float(m.group("cq")),
                                            theta_d=float(m.group("td")),
                                            theta_p=float(m.group("tp")))
    return out, log_path


def effcai_from_cai(cai, dt, tau_effca, min_ca):
    """Exact exponential-integrator solution of the GluSynapse effcai ODE."""
    decay = np.exp(-dt / tau_effca)
    b = [0.0, (1.0 - decay) * tau_effca]
    a = [1.0, -decay]
    return lfilter(b, a, np.asarray(cai, dtype=np.float64) - min_ca)


def delay_of(proto):
    m = re.search(r"_(-?\d+)ms", proto)
    return float(m.group(1)) if m else np.inf


def load_protocol(pkl_path, tau_effca, min_ca, step):
    with open(pkl_path, "rb") as fh:
        data = pickle.load(fh)
    t = np.asarray(data["t"], dtype=np.float64)
    dt = float(t[1] - t[0])
    gids = [int(g) for g in data["global_ids"]]
    eff, rho = OrderedDict(), OrderedDict()
    for gid in gids:
        eff[gid] = effcai_from_cai(data["cai_CR"][gid], dt, tau_effca, min_ca)[::step]
        rho[gid] = np.asarray(data["rho_GB"][gid], dtype=np.float64)[::step]
    return dict(t=t[::step] / 1000.0, gids=gids, effcai=eff, rho=rho,
                prespikes=np.asarray(data["prespikes"], dtype=np.float64) / 1000.0)


def _decorate(ax, th, xlim, fs=6, show_values=True):
    if th is None:
        return
    ax.axhline(th["theta_d"], color=THETA_D_COLOR, ls="--", lw=0.8, alpha=0.9)
    ax.axhline(th["theta_p"], color=THETA_P_COLOR, ls="--", lw=0.8, alpha=0.9)
    if show_values:
        ax.text(xlim[0], th["theta_d"], r"  $\theta_d$=%.3f" % th["theta_d"],
                color=THETA_D_COLOR, fontsize=fs, va="bottom", ha="left")
        ax.text(xlim[0], th["theta_p"], r"  $\theta_p$=%.3f" % th["theta_p"],
                color=THETA_P_COLOR, fontsize=fs, va="bottom", ha="left")


def figure_per_protocol(pair, proto, run, thetas, out_dir, xlim, with_rho):
    gids = run["gids"]
    ncol = 2
    nrow = int(np.ceil(len(gids) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(11, 2.0 * nrow), sharex=True)
    axes = np.atleast_1d(axes).ravel()
    t = run["t"]
    for ax, gid in zip(axes, gids):
        th = thetas.get(gid)
        ax.plot(t, run["effcai"][gid], color="#1f77b4", lw=0.8)
        _decorate(ax, th, xlim)
        if with_rho:
            axr = ax.twinx()
            axr.plot(t, run["rho"][gid], color="0.6", lw=0.8, alpha=0.7)
            axr.set_ylim(-0.05, 1.05)
            axr.set_ylabel(r"$\rho$", fontsize=8, color="0.5")
            axr.tick_params(labelsize=7, colors="0.5")
        loc = th["loc"] if th else "?"
        ax.set_title("gid %d  (%s)" % (gid, loc), fontsize=9)
        ax.set_ylabel("effcai_GB", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.set_xlim(xlim)
        ax.spines["top"].set_visible(False)
    for ax in axes[len(gids):]:
        ax.set_visible(False)
    for ax in axes[max(0, len(gids) - ncol):len(gids)]:
        ax.set_xlabel("time (s)", fontsize=8)
    fig.suptitle("%s  %s  -  reconstructed effcai_GB vs thresholds" % (pair, proto),
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = os.path.join(out_dir, "%s_%s_effcai_thetas.png" % (pair, proto))
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("saved", out)


def figure_delay_grid(pair, runs, thetas, out_dir, xlim, protos, rows_per_fig,
                      with_rho, tag):
    """rows = synapses, left column = protos[0], right column = protos[1]."""
    missing = [p for p in protos if p not in runs]
    if missing:
        raise SystemExit("protocol(s) not found in pair dir: %s" % ", ".join(missing))
    gids = runs[protos[0]]["gids"]
    chunks = [gids[i:i + rows_per_fig] for i in range(0, len(gids), rows_per_fig)]
    for ci, chunk in enumerate(chunks, start=1):
        nrow = len(chunk)
        fig, axes = plt.subplots(nrow, 2, figsize=(4.2, 1.25 * nrow),
                                 sharex=True, squeeze=False)
        for r, gid in enumerate(chunk):
            th = thetas.get(gid)
            row_max = max(float(np.max(runs[p]["effcai"][gid])) for p in protos)
            top = max(row_max, th["theta_p"] if th else 0.0) * 1.15
            for c, proto in enumerate(protos):
                ax = axes[r][c]
                run = runs[proto]
                ax.plot(run["t"], run["effcai"][gid], color="#68a8e0", lw=0.6)
                _decorate(ax, th, xlim, fs=5.5, show_values=(c == 0))
                if with_rho:
                    axr = ax.twinx()
                    axr.plot(run["t"], run["rho"][gid], color="0.6", lw=0.6, alpha=0.8)
                    axr.set_ylim(-0.05, 1.05)
                    axr.spines["right"].set_visible(True)
                    axr.tick_params(labelsize=5, colors="0.5", length=1.5)
                    if c == 1:
                        axr.set_ylabel(r"$\rho$", color="0.5")
                    else:
                        axr.set_yticklabels([])
                ax.set_xlim(xlim)
                ax.set_ylim(0, top)
                ax.tick_params(labelsize=6)
                if r == 0:
                    ax.set_title(r"$\Delta t$ = %s" % proto.split("_", 1)[1])
                if c == 0:
                    ax.set_ylabel("effcai_GB")
                    ax.text(-0.32, 0.5, "gid %d\n%s" % (gid, th["loc"] if th else "?"),
                            transform=ax.transAxes, rotation=90, va="center",
                            ha="center", fontsize=6)
                else:
                    ax.set_yticklabels([])
                if r == nrow - 1:
                    ax.set_xlabel("time (s)")
        fig.suptitle("%s  effcai_GB vs thresholds" % pair, fontsize=7)
        suffix = "" if len(chunks) == 1 else "_p%d" % ci
        out = os.path.join(out_dir, "%s_effcai_thetas_%s_vs_%s%s%s.png" % (
            pair, protos[0].split("_", 1)[1], protos[1].split("_", 1)[1], suffix, tag))
        fig.savefig(out, dpi=400, bbox_inches="tight")
        plt.close(fig)
        print("saved", out)


def figure_by_synapse(pair, runs, thetas, out_dir, xlim, tag):
    protos = list(runs)
    gids = runs[protos[0]]["gids"]
    # turbo (not coolwarm): the mid-range delay must not land on white
    cmap = plt.get_cmap("turbo")
    ncol = 2
    nrow = int(np.ceil(len(gids) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(11, 2.0 * nrow), sharex=True)
    axes = np.atleast_1d(axes).ravel()
    for ax, gid in zip(axes, gids):
        for i, proto in enumerate(protos):
            c = cmap(0.05 + 0.85 * i / max(1, len(protos) - 1))
            run = runs[proto]
            ax.plot(run["t"], run["effcai"][gid], color=c, lw=0.8,
                    label=proto.split("_", 1)[1])
        _decorate(ax, thetas.get(gid), xlim)
        loc = thetas[gid]["loc"] if gid in thetas else "?"
        ax.set_title("gid %d  (%s)" % (gid, loc), fontsize=9)
        ax.set_ylabel("effcai_GB", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.set_xlim(xlim)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    for ax in axes[len(gids):]:
        ax.set_visible(False)
    for ax in axes[max(0, len(gids) - ncol):len(gids)]:
        ax.set_xlabel("time (s)", fontsize=8)
    axes[0].legend(title=r"$\Delta t$", fontsize=7, title_fontsize=7,
                   loc="upper right", ncol=2, frameon=False)
    fig.suptitle("%s  -  effcai_GB across STDP delays vs per-synapse thresholds" % pair,
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = os.path.join(out_dir, "%s_effcai_thetas_by_synapse%s.png" % (pair, tag))
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("saved", out)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pair-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--traces-subdir", default="bluecellulab_results_optimizer",
                    help="subdir of each protocol dir holding simulation_traces.pkl "
                         "(use '.' if the pkl sits directly in the protocol dir)")
    ap.add_argument("--log-glob", default="*.log",
                    help="glob (inside the protocol dir) for the pool log carrying "
                         "the 'syn gid=... theta_d=... theta_p=...' lines")
    ap.add_argument("--tau-effca", type=float, default=TAU_EFFCA_DEFAULT)
    ap.add_argument("--min-ca", type=float, default=MIN_CA_CR_DEFAULT)
    ap.add_argument("--decimate", type=int, default=20,
                    help="keep every Nth sample when plotting (dt=0.025 ms -> 0.5 ms)")
    ap.add_argument("--xlim", type=float, nargs=2, default=None, help="seconds")
    ap.add_argument("--zoom-xlim", type=float, nargs=2, default=None,
                    help="if given, also emit a zoomed set of figures")
    ap.add_argument("--protocols", nargs="*", default=None)
    ap.add_argument("--no-rho", action="store_true")
    ap.add_argument("--style", default=DEFAULT_STYLE,
                    help="matplotlib style file (default: onerule.mplstyle next to this script)")
    ap.add_argument("--compare", nargs=2, default=["10Hz_-10ms", "10Hz_10ms"],
                    metavar=("LEFT", "RIGHT"),
                    help="the two protocols to put side by side (left, right column)")
    ap.add_argument("--rows-per-fig", type=int, default=4,
                    help="synapses per figure in the side-by-side layout")
    ap.add_argument("--extra-figures", action="store_true",
                    help="also emit the per-protocol 8-synapse grids and the "
                         "all-delays overlay")
    args = ap.parse_args()

    if args.style and os.path.exists(args.style):
        plt.style.use(args.style)
        print("style: %s" % args.style)
    elif args.style:
        print("WARNING: style file not found: %s" % args.style)

    pair = os.path.basename(os.path.normpath(args.pair_dir))
    os.makedirs(args.out_dir, exist_ok=True)

    pattern = os.path.join(args.pair_dir, "*", args.traces_subdir, "simulation_traces.pkl")
    pkls = sorted(glob.glob(pattern))
    if not pkls:
        raise SystemExit("no simulation_traces.pkl under %s" % pattern)

    runs, thetas = OrderedDict(), {}
    entries = []
    for pkl in pkls:
        proto_dir = os.path.dirname(os.path.dirname(pkl)) if args.traces_subdir != "." \
            else os.path.dirname(pkl)
        proto = os.path.basename(proto_dir)
        wanted = args.protocols if args.protocols else (
            None if args.extra_figures else args.compare)
        if wanted and proto not in wanted:
            continue
        entries.append((delay_of(proto), proto, proto_dir, pkl))
    entries.sort()

    for _, proto, proto_dir, pkl in entries:
        th, log_path = parse_thetas(proto_dir, args.log_glob, os.path.getmtime(pkl))
        if th and not thetas:
            thetas = th
            print("thresholds from %s" % log_path)
        print("loading %s ..." % proto)
        runs[proto] = load_protocol(pkl, args.tau_effca, args.min_ca, args.decimate)

    if not thetas:
        print("WARNING: no theta_d/theta_p found in logs (%s)" % args.log_glob)

    t0 = runs[next(iter(runs))]["t"]
    views = [("", args.xlim or (float(t0[0]), float(t0[-1])))]
    if args.zoom_xlim:
        views.append(("_zoom", tuple(args.zoom_xlim)))

    for tag, xlim in views:
        figure_delay_grid(pair, runs, thetas, args.out_dir, xlim, args.compare,
                          args.rows_per_fig, not args.no_rho, tag)
        if args.extra_figures:
            for proto, run in runs.items():
                figure_per_protocol(pair, proto + tag, run, thetas,
                                    args.out_dir, xlim, not args.no_rho)
            figure_by_synapse(pair, runs, thetas, args.out_dir, xlim, tag)

    for gid, th in thetas.items():
        peaks = {p: float(np.max(r["effcai"][gid])) for p, r in runs.items()}
        best = max(peaks.values())
        print("gid %d %-6s theta_d=%.4f theta_p=%.4f  max effcai over protocols=%.4f  "
              "(%s)" % (gid, th["loc"], th["theta_d"], th["theta_p"], best,
                        "crosses theta_p" if best > th["theta_p"] else
                        ("crosses theta_d only" if best > th["theta_d"] else "no crossing")))


if __name__ == "__main__":
    main()
