#!/usr/bin/env python3
"""Pairwise (unitary) L5TTPC->L5TTPC EPSP: og-delta vs antic-delta emodel.

EPSP per pair = mean somatic EPSP amplitude (mV) at the initial rho0_GB from
dhuruva_modified_edges.h5, by linear superposition of the per-pair basis CSVs
(baseline all-depressed + singleton deltas); logic reused from
plot_basis_epsp_histogram.py (get_rho0_for_pair, load_basis, extrapolate_epsp).

Usage: epsp_hist_og_vs_antic.py [--og DIR] [--antic DIR]
"""
import argparse, os, sys
import h5py, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

ROOT = "/project/rrg-emuller/dhuruva/plastyfire"
sys.path.insert(0, ROOT)
from plot_basis_epsp_histogram import get_rho0_for_pair, load_basis, extrapolate_epsp  # noqa

EDGES_H5 = f"{ROOT}/data/dhuruva_modified_edges.h5"   # read-only
HERE = os.path.dirname(os.path.abspath(__file__))
EXP_MEAN, EXP_SD = 1.3, 1.1   # Markram et al. 1997 J Physiol 500:409 (literature value, not verified here)


def pairs_in(d):
    out = set()
    for f in os.listdir(d):
        if f.startswith("basis_") and f.endswith(".csv"):
            p = f[6:-4].split("_")
            if len(p) == 2 and all(x.isdigit() for x in p):
                out.add((int(p[0]), int(p[1])))
    return out


def epsp_for(basis_dir, pre, post, rho0):
    bm, bs, sm, ss, n = load_basis(pre, post, basis_dir=basis_dir)
    if len(rho0) != n:
        raise ValueError(f"rho0 len {len(rho0)} != n_syn {n}")
    return extrapolate_epsp(rho0, bm, bs, sm, ss, n_trials=0)[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--og", default=f"{ROOT}/basis_results_edges_sabrina_n120_delta")
    ap.add_argument("--antic", default=f"{ROOT}/basis_results_edges_sabrina_n120_antic-delta")
    a = ap.parse_args()

    pairs = sorted(pairs_in(a.og) & pairs_in(a.antic))
    rows = []
    with h5py.File(EDGES_H5, "r") as hf:
        for pre, post in pairs:
            rho0 = get_rho0_for_pair(hf, pre, post)
            if rho0 is None:
                print(f"skip {pre}_{post}: no edges"); continue
            try:
                rows.append((f"{pre}_{post}", epsp_for(a.og, pre, post, rho0),
                             epsp_for(a.antic, pre, post, rho0)))
            except (FileNotFoundError, ValueError) as e:
                print(f"skip {pre}_{post}: {e}")
    df = pd.DataFrame(rows, columns=["pair", "epsp_og", "epsp_antic"])
    df.to_csv(os.path.join(HERE, "epsp_og_vs_antic.csv"), index=False)

    print(f"{'':6s} {'n':>4s} {'mean':>6s} {'median':>6s} {'SD':>6s} {'min':>6s} {'max':>6s}")
    for c, lab in [("epsp_og", "og"), ("epsp_antic", "antic")]:
        v = df[c].values
        print(f"{lab:6s} {len(v):4d} {v.mean():6.2f} {np.median(v):6.2f} {v.std(ddof=1):6.2f} {v.min():6.2f} {v.max():6.2f}")
    ratio = df.epsp_antic / df.epsp_og
    p = wilcoxon(df.epsp_antic, df.epsp_og).pvalue
    print(f"median antic/og = {ratio.median():.3f} (IQR {ratio.quantile(.25):.3f}-{ratio.quantile(.75):.3f}); Wilcoxon p = {p:.2e}")
    print(f"frac within exp mean+-SD: og {((df.epsp_og-EXP_MEAN).abs()<=EXP_SD).mean():.2f}, "
          f"antic {((df.epsp_antic-EXP_MEAN).abs()<=EXP_SD).mean():.2f}")

    plt.style.use(f"{ROOT}/onerule.mplstyle")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.5, 3.3))
    hi = max(4.0, np.ceil(df[["epsp_og", "epsp_antic"]].values.max() / 0.2) * 0.2)
    bins = np.arange(0, hi + 1e-9, 0.2)
    ax.axvspan(EXP_MEAN - EXP_SD, EXP_MEAN + EXP_SD, color="0.85", zorder=0,
               label=f"exp. {EXP_MEAN}±{EXP_SD} mV\n(Markram 1997; lit. value,\nnot verified here)")
    ax.axvline(EXP_MEAN, color="0.4", lw=1, zorder=0)
    for c, lab, col in [("epsp_og", "og-delta", "C0"), ("epsp_antic", "antic-delta", "C1")]:
        v = df[c].values
        ax.hist(v, bins=bins, histtype="stepfilled", alpha=0.35, color=col)
        ax.hist(v, bins=bins, histtype="step", color=col, lw=1.2,
                label=f"{lab} (mean {v.mean():.2f}, med {np.median(v):.2f})")
        ax.axvline(v.mean(), color=col, ls="--", lw=1)
        ax.axvline(np.median(v), color=col, ls=":", lw=1)
    ax.set_xlabel("Unitary EPSP (mV)"); ax.set_ylabel("Pairs")
    ax.set_title(f"a  L5TTPC→L5TTPC, n={len(df)}  (-- mean, ·· median)", loc="left", fontsize=9)
    ax.legend(fontsize=6, frameon=False)
    m = hi
    bx.plot([0, m], [0, m], color="0.5", lw=1, ls="--")
    bx.scatter(df.epsp_og, df.epsp_antic, s=10, color="k", alpha=0.7)
    bx.set_xlim(0, m); bx.set_ylim(0, m); bx.set_aspect("equal")
    bx.set_xlabel("EPSP og-delta (mV)"); bx.set_ylabel("EPSP antic-delta (mV)")
    bx.set_title(f"b  median ratio {ratio.median():.2f}, Wilcoxon p={p:.1e}", loc="left", fontsize=9)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, "figs", f"epsp_hist_og_vs_antic.{ext}"), dpi=200)


if __name__ == "__main__":
    main()
