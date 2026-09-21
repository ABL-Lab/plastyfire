"""
STDP curves across extracellular calcium levels (DE fit #2 params).

Reads the flat ca_scan/ output written by run_ca_scan_pool.py and overlays one
mean +/- SEM curve per Ca level. 2.0 mM is the level the model was FITTED at, so
it has no ca_scan/ directory: its curve comes from the fit-#2 BCL pool pkls in
the per-pair workdirs (simulation_edges_<hash>.pkl), which are 2.0 mM by
construction. Same rho keys, same ratio math -- only the file layout differs. Model curves only -- the in-vitro (Markram 1997)
points were dropped: they were measured at one calcium, so they are a reference
for the 2.0 mM fit, not for the swept levels.

Rendered with biodata/onerule.mplstyle (paper panel: 1.76 x 1.59 in, 7 pt).

Same EPSP basis and ratio formula as plot_stdp_ic_best.py. The basis was built
at 2.0 mM and is reused at every level on purpose: the reported quantity is a
ratio (after / before) measured at ONE calcium, so the Ca-dependent release
scaling is a common factor that divides out. (The *thresholds* are NOT reused --
those are rebuilt per level, see run_ca_scan_pool.py.)

Usage (compute node):
    python plot_ca_scan.py --workers 32
    python plot_ca_scan.py --ca 3.0 2.0 1.8 1.3 1.2 --out ca_scan.png
"""

import argparse
import multiprocessing
import os
import pickle

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from plastyfire.evaluator_edges import compute_epsp_from_basis, _load_basis

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASIS_DIR  = os.path.join(SCRIPT_DIR, "basis_results_edges_ion_channels")
OUT_ROOT   = os.path.join(SCRIPT_DIR, "ca_scan")
# 2.0 mM == the fit calcium. No sweep was run at it; the fit-#2 pool pkls already
# ARE the 2.0 mM measurement, so they are read in place from the workdirs.
FIT_CA          = 2.0
FIT_RESULTS_DIR = os.path.join(
    SCRIPT_DIR,
    "refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations")
FIT_PARAM_HASH  = "b8c7ff3ecf0a"
STYLE      = os.path.join(SCRIPT_DIR, "biodata", "onerule.mplstyle")

# Panel style. constrained_layout is on in there, so no tight_layout() below, and
# the 4-colour prop_cycle is what gives each Ca level its colour (in the plotted
# order: 3.0, 1.8, 1.3, 1.2 mM).
plt.style.use(STYLE)

PROTO_DT = {"10Hz_10ms": 10, "10Hz_5ms": 5, "10Hz_30ms": 30, "10Hz_50ms": 50,
            "10Hz_-10ms": -10, "10Hz_-30ms": -30, "10Hz_-50ms": -50}

CA_LEVELS = [3.0, FIT_CA, 1.8, 1.3, 1.2]


def ca_tag(ca):
    return f"ca_{ca:.2f}mM"


def _worker(args):
    """One pkl -> one EPSP ratio. Identical basis/ratio math to plot_stdp_ic_best."""
    ca, pkl, pair, dt = args
    try:
        with open(pkl, "rb") as f:
            raw = pickle.load(f)
        pre_gid, post_gid = pair.split("-")
        basis_df = _load_basis(BASIS_DIR, pre_gid, post_gid)
        ep_b_m, ep_b_s = compute_epsp_from_basis(basis_df, list(raw["initial_rho"]))
        ep_a_m, _      = compute_epsp_from_basis(basis_df, list(raw["final_rho"]))
        if ep_b_m == 0:
            return None
        cv2 = min((ep_b_s / ep_b_m) ** 2, 0.25)
        return {"ca": ca, "pair": pair, "dt": dt,
                "epsp_ratio": (ep_a_m / ep_b_m) * (1.0 + cv2)}
    except Exception as e:
        print(f"  Error {ca}mM {pair} dt={dt}: {e}")
        return None


def _fit_ca_tasks():
    """The 2.0 mM curve: fit-#2 pool pkls, one per pair/protocol workdir."""
    if not os.path.isdir(FIT_RESULTS_DIR):
        print(f"  ! {FIT_RESULTS_DIR} missing — skipping {FIT_CA} mM")
        return []
    out = []
    for pair in sorted(os.listdir(FIT_RESULTS_DIR)):
        if not os.path.isdir(os.path.join(FIT_RESULTS_DIR, pair)):
            continue
        for proto, dt in PROTO_DT.items():
            pkl = os.path.join(FIT_RESULTS_DIR, pair, proto,
                               f"simulation_edges_{FIT_PARAM_HASH}.pkl")
            if os.path.isfile(pkl):
                out.append((FIT_CA, pkl, pair, dt))
    if not out:
        print(f"  ! no simulation_edges_{FIT_PARAM_HASH}.pkl found — "
              f"skipping {FIT_CA} mM")
    return out


def collect(ca_levels, workers):
    tasks = []
    for ca in ca_levels:
        if ca == FIT_CA:
            tasks += _fit_ca_tasks()
            continue
        d = os.path.join(OUT_ROOT, ca_tag(ca))
        if not os.path.isdir(d):
            print(f"  ! {d} missing — skipping {ca} mM")
            continue
        for fn in sorted(os.listdir(d)):
            if not fn.endswith(".pkl"):
                continue
            pair, proto = fn[:-4].split("__")
            if proto in PROTO_DT:
                tasks.append((ca, os.path.join(d, fn), pair, PROTO_DT[proto]))
    if not tasks:
        raise SystemExit("No ca_scan pkls found — run run_ca_scan_pool.py first.")
    print(f"Reading {len(tasks)} pkls with {workers} workers …")
    with multiprocessing.Pool(workers) as pool:
        rows = pool.map(_worker, tasks)
    df = pd.DataFrame([r for r in rows if r])
    print(f"Collected {len(df)} EPSP ratios ({len(tasks) - len(df)} unusable)")
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ca", type=float, nargs="+", default=CA_LEVELS)
    ap.add_argument("--workers", type=int, default=32)
    ap.add_argument("--out", default="ca_scan_stdp.png")
    ap.add_argument("--csv", default="ca_scan_stdp.csv")
    args = ap.parse_args()

    df = collect(args.ca, args.workers)

    stats = (df.groupby(["ca", "dt"])["epsp_ratio"]
               .agg(["mean", "sem", "count"]).reset_index().sort_values(["ca", "dt"]))
    print("\n" + stats.to_string(index=False))
    stats.to_csv(args.csv, index=False)
    print(f"\nWrote {args.csv}")

    fig, ax = plt.subplots()
    # The style's prop_cycle is 4 colours and there are 5 Ca levels, so it would
    # wrap and give 1.2 mM the same blue as 3.0. Take the cycle and add one
    # neutral tail colour; the fit calcium (2.0) is drawn black so it reads as
    # the reference rather than as another sweep point.
    cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    levels = sorted(args.ca, reverse=True)
    colors, spare = {}, [c for c in cycle] + ["7b6ba8"]
    for ca in levels:
        colors[ca] = "black" if ca == FIT_CA else "#" + spare.pop(0).lstrip("#")
    for ca in levels:
        s = stats[stats["ca"] == ca]
        if s.empty:
            continue
        ax.errorbar(s["dt"], s["mean"], yerr=s["sem"], marker="o",
                    capsize=1.5, elinewidth=0.5, color=colors[ca],
                    zorder=3 if ca == FIT_CA else 2,
                    label=f"{ca} mM")
    ax.axhline(1.0, color="0.6", ls="--", lw=0.5, zorder=0)
    ax.axvline(0.0, color="0.6", ls="--", lw=0.5, zorder=0)
    ax.set_xlabel(r"$\Delta t$ (ms)")
    ax.set_ylabel("EPSP ratio")
    ax.legend()
    fig.savefig(args.out, dpi=300)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
