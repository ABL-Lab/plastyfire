"""
Compare rho transitions between Neurodamus (CoreNEURON) and bluecellulab (CVODE).

ND:  out/rho.h5                          — full time series, 420000 steps @ dt=0.1 ms
BCL: bluecellulab_results/rho.h5         — initial/final only (old format)
     bluecellulab_results/rho_timeseries.npy — full time series (new format, if present)

Global edge ID → ND local element_id mapping (pre=180164, post=197248):
  348477224 → 2173  rho0=1
  348477225 → 2174  rho0=0
  348477226 → 2175  rho0=1
  348477227 → 2176  rho0=0
"""

import os
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WORKDIR = os.path.dirname(os.path.abspath(__file__))

# Global ID → ND local element_id
GLOBAL_TO_LOCAL = {348477224: 2173, 348477225: 2174, 348477226: 2175, 348477227: 2176}
RHO0            = {348477224: 1.0,  348477225: 0.0,  348477226: 1.0,  348477227: 0.0}

# ── Load ND rho time series ───────────────────────────────────────────────────
nd_h5 = os.path.join(WORKDIR, "out", "rho.h5")
with h5py.File(nd_h5) as f:
    nd_data = f["report/S1nonbarrel_neurons/data"][()]     # (N_t, 4)
    nd_time_meta = f["report/S1nonbarrel_neurons/mapping/time"][()]  # [t0, t1, dt]
    nd_eids = f["report/S1nonbarrel_neurons/mapping/element_ids"][()]
nd_t = np.arange(nd_time_meta[0], nd_time_meta[1] + nd_time_meta[2], nd_time_meta[2])[:len(nd_data)]
# Build dict: local_eid → trajectory
nd_rho = {int(eid): nd_data[:, i] for i, eid in enumerate(nd_eids)}
print(f"ND:  {len(nd_t)} steps, t={nd_t[0]:.1f}..{nd_t[-1]:.1f} ms,  eids={nd_eids}")

# ── Load BCL rho ──────────────────────────────────────────────────────────────
bcl_ts_path = os.path.join(WORKDIR, "bluecellulab_results", "rho_timeseries.npy")
bcl_ids_path = os.path.join(WORKDIR, "bluecellulab_results", "rho_global_ids.npy")
bcl_h5_path  = os.path.join(WORKDIR, "bluecellulab_results", "rho.h5")

bcl_has_timeseries = os.path.exists(bcl_ts_path) and os.path.exists(bcl_ids_path)

if bcl_has_timeseries:
    bcl_ts  = np.load(bcl_ts_path)     # (N_t, 1+N_syns)
    bcl_gids = np.load(bcl_ids_path).tolist()
    bcl_t   = bcl_ts[:, 0]
    bcl_rho_by_gid = {int(gid): bcl_ts[:, i+1] for i, gid in enumerate(bcl_gids)}
    print(f"BCL (timeseries): {len(bcl_t)} steps, t={bcl_t[0]:.1f}..{bcl_t[-1]:.1f} ms")
else:
    print(f"BCL timeseries not found; using initial/final from rho.h5")
    with h5py.File(bcl_h5_path) as f:
        bcl_data = f["report/S1nonbarrel_neurons/data"][()]  # (2, 4)
        bcl_gids_raw = f["report/S1nonbarrel_neurons/mapping/element_ids"][()]
    bcl_init  = {int(g): float(bcl_data[0, i]) for i, g in enumerate(bcl_gids_raw)}
    bcl_final = {int(g): float(bcl_data[1, i]) for i, g in enumerate(bcl_gids_raw)}
    bcl_rho_by_gid = None
    print(f"BCL (initial/final): gids={list(bcl_init.keys())}")
    for g in sorted(bcl_init):
        print(f"  gid {g}  rho0={bcl_init[g]:.3f}  rho_final={bcl_final[g]:.3f}")

# ── Build matched synapse list ────────────────────────────────────────────────
synapses = [(gid, GLOBAL_TO_LOCAL[gid]) for gid in sorted(GLOBAL_TO_LOCAL)]

# ── Print comparison table ────────────────────────────────────────────────────
print("\n{:>12s} {:>8s} {:>8s} {:>10s} {:>10s} {:>10s}".format(
    "global_id", "local_id", "rho0", "BCL_final", "ND_final", "delta(BCL-ND)"))
print("-" * 70)
for gid, lid in synapses:
    rho0 = RHO0[gid]
    nd_final = float(nd_rho[lid][-1])
    if bcl_rho_by_gid is not None:
        bcl_fin = float(bcl_rho_by_gid[gid][-1])
    else:
        bcl_fin = bcl_final[gid]
    print(f"  {gid:>10d} {lid:>8d} {rho0:>8.1f} {bcl_fin:>10.4f} {nd_final:>10.4f} {bcl_fin-nd_final:>+11.4f}")

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(4, 1, figsize=(14, 10), sharex=True)
colors = {"ND": "steelblue", "BCL": "tomato"}

t_sec = nd_t / 1000.0  # ms → s

for ax, (gid, lid) in zip(axes, synapses):
    rho0 = RHO0[gid]
    nd_trace = nd_rho[lid]

    ax.plot(t_sec, nd_trace, color=colors["ND"], lw=1.0, label=f"ND (CoreNEURON)")

    if bcl_rho_by_gid is not None:
        bcl_t_sec = bcl_t / 1000.0
        bcl_trace = bcl_rho_by_gid[gid]
        # Downsample to ND dt for readability
        bcl_interp = np.interp(t_sec, bcl_t_sec, bcl_trace)
        ax.plot(t_sec, bcl_interp, color=colors["BCL"], lw=1.0, alpha=0.8,
                label="BCL (CVODE, fixed RNG)")
        delta = float(bcl_trace[-1]) - float(nd_trace[-1])
        ax.set_title(f"syn {gid} (ND:{lid})  rho₀={rho0:.0f}  "
                     f"BCL_final={bcl_trace[-1]:.3f}  ND_final={nd_trace[-1]:.3f}  "
                     f"Δ={delta:+.3f}", fontsize=9)
    else:
        # Old format: draw start/end as horizontal dashed lines
        nd_end_t = nd_t[-1] / 1000.0
        bcl_fin = bcl_final[gid]
        ax.hlines(rho0,    0, nd_end_t, colors=colors["BCL"], linestyles="--", lw=1.5,
                  label=f"BCL initial={rho0:.2f}")
        ax.hlines(bcl_fin, 0, nd_end_t, colors=colors["BCL"], linestyles="-",  lw=2.0,
                  label=f"BCL final={bcl_fin:.3f}  ⚠ STALE (pre SK-fix)")
        delta = bcl_fin - float(nd_trace[-1])
        ax.set_title(f"syn {gid} (ND:{lid})  rho₀={rho0:.0f}  "
                     f"BCL_final={bcl_fin:.3f} ⚠  ND_final={nd_trace[-1]:.3f}  "
                     f"Δ={delta:+.3f}", fontsize=9)

    ax.set_ylim(-0.05, 1.05)
    ax.set_ylabel("rho_GB")
    ax.axhline(0.5, color="gray", lw=0.5, ls=":")
    ax.legend(fontsize=8, loc="upper right")

axes[-1].set_xlabel("Time (s)")

status = "BCL timeseries (re-run with SK fix + fixed RNG)" if bcl_has_timeseries \
         else "⚠  BCL result is STALE (Apr 27, before SK fix in ND)\n   Re-run BCL with fixed simulator_edges.py to get matching timeseries"
fig.suptitle(f"ρ Transitions: BCL vs ND — 10Hz +5ms (pre→post)\n{status}", fontsize=10)
plt.tight_layout()

out_png = os.path.join(WORKDIR, "bcl_vs_nd_rho.png")
plt.savefig(out_png, dpi=150)
print(f"\nSaved: {out_png}")
