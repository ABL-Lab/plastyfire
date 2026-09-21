"""
Validate the EPSP basis by checking the assumptions compute_epsp_from_basis makes.

The basis stores, per pair: EPSP with all synapses depressed ("0,0,..,0"), EPSP with
exactly one synapse potentiated (each singleton), and EPSP with all potentiated
("1,1,..,1"). compute_epsp_from_basis then predicts ANY rho vector by linear
superposition:

    EPSP(rho) = e0 + sum_i rho_i * (e_i - e0)

The all-ones row is the one configuration that is BOTH measured and predictable, so
it is the only internal check on that superposition assumption. This script uses it.

Checks per pair
---------------
  superposition : predicted all-ones vs measured all-ones (the headline test)
  monotonic     : every singleton EPSP > all-zeros EPSP (potentiating one synapse
                  must not lower the EPSP)
  bracketed     : measured all-ones >= every singleton
  positive      : all means > 0
  trials        : same trial count on every row
  noise         : per-config coefficient of variation (std/mean)

Usage:
    python validate_basis.py --basis-dir basis_results_edges_ion_channels
    python validate_basis.py --basis-dir basis_results_edges_ion_channels --sample 10 --seed 0
"""

import argparse
import glob
import os
import random

import numpy as np
import pandas as pd


def check_pair(path):
    d = pd.read_csv(path)
    cfgs = {r.config: r for r in d.itertuples()}
    n = len(d.config.iloc[0].split(","))
    zeros = ",".join(["0"] * n)
    ones = ",".join(["1"] * n)
    out = {"file": os.path.basename(path), "n_syn": n, "n_rows": len(d)}

    if zeros not in cfgs or ones not in cfgs:
        out["status"] = "MISSING all-zeros or all-ones"
        return out

    e0 = cfgs[zeros].mean
    e1 = cfgs[ones].mean
    out["e0"] = e0
    out["e_allones_measured"] = e1

    singles, missing = [], 0
    for i in range(n):
        s = ["0"] * n
        s[i] = "1"
        k = ",".join(s)
        if k not in cfgs:
            missing += 1
            continue
        singles.append(cfgs[k].mean)
    out["n_singletons"] = len(singles)
    out["missing_singletons"] = missing
    if missing:
        out["status"] = "MISSING %d singleton(s)" % missing
        return out

    pred = e0 + sum(s - e0 for s in singles)
    out["e_allones_predicted"] = pred
    out["superposition_err_pct"] = 100.0 * (pred - e1) / e1 if e1 else np.nan

    out["monotonic"] = all(s > e0 for s in singles)
    out["n_nonmonotonic"] = sum(1 for s in singles if s <= e0)
    out["bracketed"] = e1 >= max(singles) - 1e-12
    out["positive"] = bool((d["mean"] > 0).all())
    out["trials_consistent"] = d["count"].nunique() == 1
    out["trials"] = int(d["count"].iloc[0])
    cv = (d["std"] / d["mean"]).replace([np.inf, -np.inf], np.nan)
    out["cv_median"] = float(cv.median())
    out["cv_max"] = float(cv.max())
    out["status"] = "ok"
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--basis-dir", required=True)
    ap.add_argument("--sample", type=int, default=0,
                    help="validate only N random pairs (0 = all)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--csv", default=None, help="write per-pair results here")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.basis_dir, "basis_*.csv")))
    if not files:
        raise SystemExit("no basis CSVs in %s" % args.basis_dir)
    if args.sample:
        random.seed(args.seed)
        files = sorted(random.sample(files, min(args.sample, len(files))))

    print("basis dir : %s" % args.basis_dir)
    print("pairs     : %d%s\n" % (len(files),
          "  (random sample, seed=%d)" % args.seed if args.sample else ""))

    rows = [check_pair(f) for f in files]
    df = pd.DataFrame(rows)

    bad = df[df.status != "ok"]
    if len(bad):
        print("STRUCTURAL PROBLEMS (%d):" % len(bad))
        for r in bad.itertuples():
            print("  %s: %s" % (r.file, r.status))
        print()
    ok = df[df.status == "ok"].copy()
    if not len(ok):
        raise SystemExit("no structurally valid pairs")

    print("%-30s %5s %8s %9s %9s %8s" % (
        "pair", "nsyn", "e0", "ones_meas", "ones_pred", "err%"))
    print("-" * 74)
    for r in ok.itertuples():
        print("%-30s %5d %8.4f %9.4f %9.4f %+8.1f" % (
            r.file.replace("basis_", "").replace(".csv", ""),
            r.n_syn, r.e0, r.e_allones_measured, r.e_allones_predicted,
            r.superposition_err_pct))

    e = ok.superposition_err_pct
    print("\n--- superposition error (predicted vs measured all-ones) ---")
    print("  median %+.1f%%   mean %+.1f%%   min %+.1f%%   max %+.1f%%" % (
        e.median(), e.mean(), e.min(), e.max()))
    print("  within +/-10%%: %d/%d      within +/-20%%: %d/%d" % (
        (e.abs() <= 10).sum(), len(e), (e.abs() <= 20).sum(), len(e)))
    print("  sign: %d over-predict, %d under-predict" % ((e > 0).sum(), (e < 0).sum()))
    if (e > 0).all():
        print("  NOTE: all over-predict -> systematic, consistent with EPSP saturation")
        print("        (synapses do not sum linearly; the model has no saturation term)")

    print("\n--- structural checks ---")
    print("  monotonic (every singleton > all-zeros) : %d/%d" % (
        ok.monotonic.sum(), len(ok)))
    if (~ok.monotonic).any():
        for r in ok[~ok.monotonic].itertuples():
            print("      %s: %d singleton(s) <= e0" % (r.file, r.n_nonmonotonic))
    print("  bracketed (all-ones >= max singleton)   : %d/%d" % (
        ok.bracketed.sum(), len(ok)))
    print("  all means positive                      : %d/%d" % (
        ok.positive.sum(), len(ok)))
    print("  trial count consistent within pair      : %d/%d  (trials=%s)" % (
        ok.trials_consistent.sum(), len(ok), sorted(ok.trials.unique())))
    print("  per-config CV: median %.3f  worst %.3f" % (
        ok.cv_median.median(), ok.cv_max.max()))

    print("\n--- error vs synapse count (saturation should grow with n) ---")
    for lo, hi in [(0, 4), (5, 8), (9, 12), (13, 99)]:
        sub = ok[(ok.n_syn >= lo) & (ok.n_syn <= hi)]
        if len(sub):
            print("  n_syn %2d-%2d : n=%2d  median err %+.1f%%" % (
                lo, hi, len(sub), sub.superposition_err_pct.median()))

    if args.csv:
        df.to_csv(args.csv, index=False)
        print("\nwrote %s" % args.csv)


if __name__ == "__main__":
    main()
