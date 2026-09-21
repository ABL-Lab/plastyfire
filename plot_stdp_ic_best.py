"""
Plot STDP curve for ion-channel optimizer best params vs in vitro (Markram 1997, 10 Hz).

Reads simulation_edges_0ce64fa83b85.pkl from every pair/protocol workdir,
computes EPSP ratio via ion-channel basis, plots mean ± SEM vs in vitro.

Panel styling follows biodata/onerule.mplstyle (same as plot_ca_scan.py), so the
figure comes out at the paper panel size rather than a screen-sized figure.

Usage:
    python plot_stdp_ic_best.py [--out stdp_ic_best.png]
"""

import argparse
import os
import pickle
from concurrent.futures import ProcessPoolExecutor, as_completed

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plastyfire.evaluator_edges import compute_epsp_from_basis, _load_basis

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
STYLE       = os.path.join(SCRIPT_DIR, "biodata", "onerule.mplstyle")
# Same paper-panel style as plot_ca_scan.py: 1.76x1.59 in, 7 pt, thin lines,
# constrained_layout on (so no tight_layout() below) and a 4-colour prop_cycle.
plt.style.use(STYLE)
RESULTS_DIR = os.path.join(
    SCRIPT_DIR,
    "refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations",
)
BASIS_DIR  = os.path.join(SCRIPT_DIR, "basis_results_edges_ion_channels")
PARAM_HASH = "0ce64fa83b85"  # ion-channel optimizer gen-4 best (untied apical/basal)
# The tied (--same-apical-basal) variant of the same gen-4 params hashes
# differently and is written by run_ic_best_pool_tied.py. Select it with
# --param-hash cdf3a1e1db98.
PARAM_HASH_TIED = "cdf3a1e1db98"

PROTO_TO_DT = {
    "10Hz_5ms":   5,
    "10Hz_10ms":  10,
    "10Hz_30ms":  30,
    "10Hz_50ms":  50,
    "10Hz_-10ms": -10,
    "10Hz_-30ms": -30,
    "10Hz_-50ms": -50,
}

INVITRO_DT   = [-10,    5,     10]
INVITRO_MEAN = [0.7922, 1.2038, 1.2013]
INVITRO_SEM  = [0.0259, 0.0644, 0.0626]


def _worker(args):
    pair, proto_dir, dt = args
    pre_gid, post_gid = pair.split("-")
    workdir = os.path.join(RESULTS_DIR, pair, proto_dir)
    pkl     = os.path.join(workdir, f"simulation_edges_{PARAM_HASH}.pkl")

    if not os.path.exists(pkl):
        return None

    try:
        with open(pkl, "rb") as f:
            raw = pickle.load(f)
        initial_rho = list(raw["initial_rho"])
        final_rho   = list(raw["final_rho"])
        basis_df    = _load_basis(BASIS_DIR, pre_gid, post_gid)
        ep_b_m, ep_b_s = compute_epsp_from_basis(basis_df, initial_rho)
        ep_a_m, _      = compute_epsp_from_basis(basis_df, final_rho)
        if ep_b_m == 0:
            return None
        cv2   = min((ep_b_s / ep_b_m) ** 2, 0.25)
        ratio = (ep_a_m / ep_b_m) * (1.0 + cv2)
        # Runs written with pairrunner --force fired the wrong number of post
        # spikes during induction. Keep them separable so a forced batch cannot
        # quietly become "the" STDP curve.
        forced = bool(raw.get("guardrail_forced", False))
        return {"pair": pair, "dt": dt, "epsp_ratio": ratio, "forced": forced,
                "n_post": raw.get("n_post_spikes"),
                "n_post_exp": raw.get("n_post_spikes_exp")}
    except Exception as e:
        print(f"  Error {pair} {proto_dir}: {e}")
        return None


def main():
    # _worker reads these as module globals; on fork they are inherited by the
    # pool, so assigning here before submitting is what makes the flags take effect.
    global PARAM_HASH, BASIS_DIR
    _default_hash, _default_basis = PARAM_HASH, BASIS_DIR

    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="stdp_ic_best.png")
    parser.add_argument("--param-hash", default=_default_hash,
                        help=f"param hash naming the pkl to read "
                             f"(default {_default_hash}; tied run = {PARAM_HASH_TIED})")
    parser.add_argument("--basis-dir", default=_default_basis)
    parser.add_argument("--drop-forced", action="store_true",
                        help="Exclude runs whose pkl has guardrail_forced=True "
                             "(written by pairrunner/pool --force, i.e. the post "
                             "cell missed the induction spike count)")
    parser.add_argument("--label", default=None,
                        help="legend label for the simulated curve (default 'In silico')")
    args = parser.parse_args()

    PARAM_HASH = args.param_hash
    BASIS_DIR  = args.basis_dir
    print(f"param_hash : {PARAM_HASH}")
    print(f"basis_dir  : {BASIS_DIR}")

    jobs = []
    for pair_dir in sorted(os.scandir(RESULTS_DIR), key=lambda e: e.name):
        if not pair_dir.is_dir():
            continue
        pair = pair_dir.name
        for proto_dir, dt in PROTO_TO_DT.items():
            workdir = os.path.join(RESULTS_DIR, pair, proto_dir)
            if os.path.isfile(os.path.join(workdir, "prefire_simulation_config.json")):
                jobs.append((pair, proto_dir, dt))

    print(f"Processing {len(jobs)} pair/protocol combinations …")
    rows, missing = [], 0
    with ProcessPoolExecutor() as ex:
        for r in as_completed(ex.submit(_worker, j) for j in jobs):
            result = r.result()
            if result is None:
                missing += 1
            else:
                rows.append(result)

    df = pd.DataFrame(rows)
    print(f"Collected {len(df)} EPSP ratios  ({missing} missing pkls)")

    if not df.empty and "forced" in df:
        n_forced = int(df["forced"].sum())
        if n_forced and args.drop_forced:
            print(f"Dropping {n_forced} forced run(s) (guardrail_forced=True)")
            df = df[~df["forced"]].reset_index(drop=True)
        elif n_forced:
            print(f"WARNING: {n_forced} of {len(df)} runs are FORCED "
                  f"(guardrail_forced=True — post spike count off target); "
                  f"pass --drop-forced to exclude them")
            for _, r in df[df["forced"]].sort_values("dt").iterrows():
                print(f"    forced: {r['pair']} dt={int(r['dt']):+4d} "
                      f"({r['n_post']}/{r['n_post_exp']} post spikes)")

    if df.empty:
        print(f"No data — no simulation_edges_{PARAM_HASH}.pkl found under {RESULTS_DIR}")
        return

    summary = (
        df.groupby("dt")["epsp_ratio"]
        .agg(["mean", "sem", "count"])
        .reset_index()
        .sort_values("dt")
    )
    print("\nSummary by dt:")
    for _, row in summary.iterrows():
        print(f"  dt={int(row['dt']):+4d} ms: mean={row['mean']:.4f} ± {row['sem']:.4f}  (n={int(row['count'])})")

    fig, ax = plt.subplots()

    # Colours, line widths and marker sizes all come from the style sheet; the
    # only thing set here is which cycle entry each curve gets (0 = blue in
    # silico, 1 = orange in vitro) and the thin capsize that matches it.
    cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    ax.errorbar(
        summary["dt"], summary["mean"], yerr=summary["sem"],
        marker="o", ls="-", color="#" + cycle[0].lstrip("#"),
        capsize=1.5, elinewidth=0.5, zorder=3,
        label=args.label or "In silico",
    )
    ax.errorbar(
        INVITRO_DT, INVITRO_MEAN, yerr=INVITRO_SEM,
        marker="^", ls="--", color="#" + cycle[1].lstrip("#"),
        capsize=1.5, elinewidth=0.5, zorder=2,
        label="In vitro (Markram 1997)",
    )

    ax.axhline(1.0, color="0.6", ls="--", lw=0.5, zorder=0)
    ax.axvline(0.0, color="0.6", ls="--", lw=0.5, zorder=0)
    ax.set_xlabel(r"$\Delta t$ (ms)")
    ax.set_ylabel("EPSP ratio")
    # Headroom above the +5/+10 ms peak so the top-right legend does not sit on
    # the data; without it the in-vitro label crosses the peak error bars.
    _top = max(summary["mean"] + summary["sem"]) 
    ax.set_ylim(top=_top + 0.30 * (_top - min(summary["mean"] - summary["sem"])))
    ax.legend(loc="upper right", handlelength=1.2, borderaxespad=0.2,
              labelspacing=0.3, handletextpad=0.4)

    out_path = os.path.join(SCRIPT_DIR, args.out)
    fig.savefig(out_path, dpi=300)
    print(f"\nSaved → {out_path}")


if __name__ == "__main__":
    main()
