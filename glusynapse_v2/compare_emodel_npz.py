"""Plumbing check for run_new_emodel_chain.sh: compare the extraction npz (and EPSP basis csv) of a re-run with the
existing ones for the same pair/protocol.

    python glusynapse_v2/compare_emodel_npz.py --pair NEWDIR=REFDIR [--pair ...] [--basis NEWDIR=REFDIR ...]

Run on the og-delta emodel under a test name the numbers must be identical (same workdirs, same seeds). Per npz file it
prints max |d| of rho0, rho_obs, c_pre, c_post, the per-record ratio mean(rho_obs)/mean(rho0) (new, ref), the max
relative difference of effcai / shaft_cai / vdcc, and IDENTICAL / DIFF / NOREF. Exit 1 only when a NEWDIR has no npz at all
(a step did not run); a DIFF is reported, not an error (e.g. the Markram 10Hz workdirs of og-delta were calibrated on
modified_emodels_hoc, so a re-calibration on delta is expected to differ there).
"""
import argparse
import glob
import os
import sys

import numpy as np

TOL = 1e-9


def _maxabs(a, b):
    a, b = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    if a.shape != b.shape:
        return float("inf")
    if a.size == 0:
        return 0.0
    both_nan = np.isnan(a) & np.isnan(b)
    d = np.abs(a - b)
    d[both_nan] = 0.0
    return float(np.nanmax(d)) if d.size else 0.0


def _maxrel(a, b):
    a, b = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    if a.shape != b.shape:
        return float("inf")
    if a.size == 0:
        return 0.0
    scale = max(float(np.nanmax(np.abs(b))), 1e-30)
    return _maxabs(a, b) / scale


def _ratio(z):
    r0 = float(np.mean(z["rho0"]))
    return float(np.mean(z["rho_obs"])) / r0 if r0 > 0 else float("nan")


def compare_dir(new, ref):
    files = sorted(glob.glob(os.path.join(new, "*.npz")))
    print(f"== {new}\n   vs {ref}\n   {len(files)} new npz", flush=True)
    n_same = n_diff = n_noref = 0
    for f in files:
        name = os.path.basename(f)
        g = os.path.join(ref, name)
        if not os.path.isfile(g):
            n_noref += 1
            print(f"  NOREF     {name}")
            continue
        a, b = np.load(f), np.load(g)
        same_syn = a["syn"].shape == b["syn"].shape and bool(np.array_equal(a["syn"], b["syn"]))
        if not same_syn:
            n_diff += 1
            print(f"  DIFF      {name}  synapse sets differ ({a['syn'].size} vs {b['syn'].size})")
            continue
        d = {k: _maxabs(a[k], b[k]) for k in ("rho0", "rho_obs", "c_pre", "c_post", "prespikes", "postspikes")}
        sig = {k: _maxrel(a[k], b[k]) for k in ("effcai", "shaft_cai", "vdcc") if k in a.files and k in b.files}
        ok = max(d.values()) <= TOL and (not sig or max(sig.values()) <= 1e-6)
        n_same += ok
        n_diff += (not ok)
        print(f"  {'IDENTICAL' if ok else 'DIFF     '} {name}  n_syn {a['syn'].size}  "
              f"d_rho0 {d['rho0']:.2e} d_rho_obs {d['rho_obs']:.2e} d_cpre {d['c_pre']:.2e} d_cpost {d['c_post']:.2e} "
              f"d_spikes {max(d['prespikes'], d['postspikes']):.2e}  "
              f"ratio new {_ratio(a):.4f} ref {_ratio(b):.4f}  "
              + " ".join(f"rel_{k} {v:.1e}" for k, v in sig.items()), flush=True)
    print(f"   SUMMARY identical {n_same}, diff {n_diff}, no-reference {n_noref}", flush=True)
    return len(files), n_same, n_diff, n_noref


def compare_basis(new, ref):
    import pandas as pd
    files = sorted(glob.glob(os.path.join(new, "basis_*.csv")))
    print(f"== basis {new}\n   vs {ref}\n   {len(files)} new csv", flush=True)
    n_same = n_diff = n_noref = 0
    for f in files:
        name = os.path.basename(f)
        g = os.path.join(ref, name)
        if not os.path.isfile(g):
            n_noref += 1
            print(f"  NOREF     {name}")
            continue
        a, b = pd.read_csv(f), pd.read_csv(g)
        if "config" in a.columns and "config" in b.columns:   # row order differs between script versions: align on the config key
            a, b = a.sort_values("config").reset_index(drop=True), b.sort_values("config").reset_index(drop=True)
            if not a["config"].equals(b["config"]):
                n_diff += 1
                print(f"  DIFF      {name}  config sets differ")
                continue
        num = [c for c in a.columns if c in b.columns and np.issubdtype(a[c].dtype, np.number)]
        if a.shape[0] != b.shape[0] or not num:
            n_diff += 1
            print(f"  DIFF      {name}  rows {a.shape[0]} vs {b.shape[0]}")
            continue
        if a.shape[1] != b.shape[1]:   # older reference csv has fewer columns (script gained columns since): compare the common ones
            print(f"   note: columns only in new {sorted(set(a.columns)-set(b.columns))}, only in ref {sorted(set(b.columns)-set(a.columns))}")
        m = max(_maxabs(a[c].to_numpy(), b[c].to_numpy()) for c in num)
        rel = max(_maxrel(a[c].to_numpy(), b[c].to_numpy()) for c in num)
        ok = m <= 1e-6
        close = rel <= 1e-3
        n_same += ok
        n_diff += (not ok and not close)
        print(f"  {'IDENTICAL' if ok else ('CLOSE    ' if close else 'DIFF     ')} {name}  max|d| over {len(num)} common numeric columns {m:.2e} (max rel {rel:.1e})", flush=True)
    print(f"   SUMMARY identical {n_same}, diff {n_diff}, no-reference {n_noref}", flush=True)
    return len(files), n_same, n_diff, n_noref


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pair", action="append", default=[], help="NEWDIR=REFDIR of extraction npz (repeatable)")
    ap.add_argument("--basis", action="append", default=[], help="NEWDIR=REFDIR of basis csv (repeatable)")
    args = ap.parse_args()
    empty = []
    tot = [0, 0, 0, 0]
    for s in args.pair:
        new, ref = s.split("=", 1)
        r = compare_dir(new, ref)
        tot = [x + y for x, y in zip(tot, r)]
        if r[0] == 0:
            empty.append(new)
    for s in args.basis:
        new, ref = s.split("=", 1)
        r = compare_basis(new, ref)
        if r[0] == 0:
            empty.append(new)
    print(f"RESULT npz: identical {tot[1]}, diff {tot[2]}, no-reference {tot[3]} of {tot[0]} new files", flush=True)
    if empty:
        print("MISSING OUTPUT (step did not run): " + ", ".join(empty), flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
