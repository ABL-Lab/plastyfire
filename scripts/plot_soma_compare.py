import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

WORKDIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations/180164-197248/10Hz_5ms"

# ── Load BCL ──────────────────────────────────────────────────────────────────
bcl = np.load(f"{WORKDIR}/bcl_test_soma.npy")
t_bcl, v_bcl = bcl[:, 0], bcl[:, 1]

# ── Load ND ───────────────────────────────────────────────────────────────────
with h5py.File(f"{WORKDIR}/out/soma.h5", "r") as f:
    tinfo = f["report/S1nonbarrel_neurons/mapping/time"][:]
    v_nd  = f["report/S1nonbarrel_neurons/data"][:, 0].astype(float)
t_start, t_end, dt = tinfo
t_nd = np.arange(len(v_nd)) * dt + t_start + dt  # ND reports at end of interval

def calc_vrest(t, v, window=(800, 980)):
    """Baseline = median of the quiet window just before first pre-spike (~1000 ms)."""
    mask = (t >= window[0]) & (t <= window[1])
    return np.median(v[mask]) if mask.any() else np.median(v)

v_rest_bcl = calc_vrest(t_bcl, v_bcl)
v_rest_nd  = calc_vrest(t_nd,  v_nd)
v_rest = v_rest_bcl  # used for title

# ── EPSP peak detection ───────────────────────────────────────────────────────
def find_epsp_peaks(t, v, rest, min_height_above_rest=0.05, min_distance_ms=50,
                    t_min=900.0):
    """Find EPSP peaks after t_min ms (ignores initial voltage drift artefact).

    Uses average dt for the distance calculation so that variable-step CVODE
    traces (where median dt is dominated by tiny steps during EPSPs) are
    handled correctly.
    """
    # Average dt: stable for both fixed-step and CVODE variable-step outputs
    avg_dt = (t[-1] - t[0]) / max(len(t) - 1, 1)
    min_dist_samples = max(1, int(min_distance_ms / avg_dt))
    peaks, props = find_peaks(v, height=rest + min_height_above_rest,
                              distance=min_dist_samples)
    # filter out pre-stimulus artefact
    peaks = peaks[t[peaks] >= t_min]
    return t[peaks], v[peaks]

pt_bcl, pv_bcl = find_epsp_peaks(t_bcl, v_bcl, v_rest_bcl)
pt_nd,  pv_nd  = find_epsp_peaks(t_nd,  v_nd,  v_rest_nd)

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(3, 1, figsize=(14, 11))

ax = axes[0]
ax.plot(t_nd,  v_nd,  lw=0.8, color="steelblue", label="Neurodamus (CoreNEURON)")
ax.plot(t_bcl, v_bcl, lw=0.8, color="tomato",    label="bluecellulab (CVODE)", alpha=0.85)
ax.scatter(pt_nd,  pv_nd,  s=20, color="steelblue", zorder=5)
ax.scatter(pt_bcl, pv_bcl, s=20, color="tomato",    zorder=5, marker="x")
ax.set_ylabel("Vm (mV)")
ax.set_title(f"Soma voltage: BCL vs Neurodamus — EPSP-only (no post-spikes), 1500 ms  |  Vrest≈{v_rest:.1f} mV")
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Zoom on first EPSP
ax2 = axes[1]
ax2.plot(t_nd,  v_nd,  lw=1.2, color="steelblue", label="ND")
ax2.plot(t_bcl, v_bcl, lw=1.2, color="tomato",    label="BCL", alpha=0.85)
ax2.set_xlim(950, 1150)
ax2.set_ylabel("Vm (mV)")
ax2.set_title("Zoom: first EPSP (950–1150 ms)")
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)

# EPSP amplitude & timing delta bar chart
ax3 = axes[2]
n = min(len(pt_bcl), len(pt_nd))
if n > 0:
    idx = np.arange(n)
    w = 0.35
    amp_nd  = pv_nd[:n]  - v_rest_nd
    amp_bcl = pv_bcl[:n] - v_rest_bcl
    ax3.bar(idx - w/2, amp_nd,  w, color="steelblue", label="ND  peak ΔV")
    ax3.bar(idx + w/2, amp_bcl, w, color="tomato",    label="BCL peak ΔV", alpha=0.85)
    ax3.set_xticks(idx)
    ax3.set_xticklabels([f"EPSP {i+1}" for i in idx])
    ax3.set_ylabel("Peak ΔV above rest (mV)")
    ax3.set_title("EPSP amplitude comparison")
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3, axis="y")

ax2.set_xlabel("Time (ms)")
plt.tight_layout()
out = f"{WORKDIR}/bcl_vs_nd_soma.png"
plt.savefig(out, dpi=150)
print(f"Saved: {out}")

# ── Numeric summary ───────────────────────────────────────────────────────────
print(f"\nVrest  BCL={v_rest_bcl:.3f} mV   ND={v_rest_nd:.3f} mV")
print(f"\nEPSP peaks BCL ({len(pt_bcl)}):  t={np.round(pt_bcl,3)}  amp={np.round(pv_bcl - v_rest_bcl, 4)} mV")
print(f"EPSP peaks ND  ({len(pt_nd)}):   t={np.round(pt_nd, 3)}  amp={np.round(pv_nd  - v_rest_nd,  4)} mV")

n = min(len(pt_bcl), len(pt_nd))
if n > 0:
    dt_peak = pt_bcl[:n] - pt_nd[:n]
    da_peak = (pv_bcl[:n] - v_rest_bcl) - (pv_nd[:n] - v_rest_nd)
    print(f"\nPeak timing delta BCL-ND (ms): {np.round(dt_peak, 4)}")
    print(f"  mean={dt_peak.mean():.4f}  std={dt_peak.std():.4f}")
    print(f"\nPeak amplitude delta BCL-ND (mV): {np.round(da_peak, 4)}")
    print(f"  mean={da_peak.mean():.4f}  std={da_peak.std():.4f}")

# ── Spike check ───────────────────────────────────────────────────────────────
def find_spikes(t, v, thr=-30):
    above = v > thr
    crossings = np.where(~above[:-1] & above[1:])[0]
    return t[crossings + 1]

spk_bcl = find_spikes(t_bcl, v_bcl)
spk_nd  = find_spikes(t_nd,  v_nd)
if len(spk_bcl) or len(spk_nd):
    print(f"\nBCL spikes ({len(spk_bcl)}): {np.round(spk_bcl, 3)}")
    print(f"ND  spikes ({len(spk_nd)}):  {np.round(spk_nd, 3)}")
    if len(spk_bcl) == len(spk_nd) and len(spk_bcl) > 0:
        dt_spk = spk_bcl - spk_nd
        print(f"Spike timing delta: mean={dt_spk.mean():.4f}  std={dt_spk.std():.4f}")
else:
    print("\nNo spikes detected (EPSP-only run).")
