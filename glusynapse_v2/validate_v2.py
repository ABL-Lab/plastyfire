"""Data checks that must pass before fitting on a new extraction (run on every new npz set).

  1. files: every npz has shaft_cai, t_ms monotonic, prespikes inside the grid, rho0 in {0,1}
  2. rho fast vs full (mean-field vs ODE) per protocol, at the given a-params: max |ratio diff|
  3. offline v1 rho (full) vs BCL rho_obs where the traces came from a run with the same a-params

    python glusynapse_v2/validate_v2.py --dirs glusynapse_v2/extracted/ebner_delta-prefire \
        --a fit_results/delta-cooker.json --pairs 10
"""
import argparse, glob, json, os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from batch_v2 import BatchV2   # noqa: E402


def check_files(dirs):
    bad = []
    for d in dirs:
        for f in sorted(glob.glob(os.path.join(d, "*.npz"))):
            z = np.load(f)
            t = z["t_ms"] if "t_ms" in z.files else None
            msg = []
            if "shaft_cai" not in z.files: msg.append("no shaft_cai")
            if t is not None and np.any(np.diff(t) <= 0): msg.append("t_ms not increasing")
            if t is not None and len(z["prespikes"]) and (z["prespikes"].max() > t[-1] or z["prespikes"].min() < t[0]):
                msg.append("prespikes outside grid")
            if not np.all(np.isin(z["rho0"], [0.0, 1.0])): msg.append("rho0 not binary")
            if "shaft_cai" in z.files and z["shaft_cai"].shape != z["effcai"].shape: msg.append("shape mismatch")
            if msg: bad.append((os.path.basename(f), msg))
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dirs", required=True)
    ap.add_argument("--a", required=True, help="json with a-params (fit_results/*.json)")
    ap.add_argument("--pairs", type=int, default=10, help="first N pairs for the fast/full check")
    ap.add_argument("--same-run", action="store_true", help="traces were generated with these a-params: compare rho to rho_obs")
    args = ap.parse_args()
    dirs = args.dirs.split(",")

    bad = check_files(dirs)
    print(f"[1] files: {len(bad)} bad" + "".join(f"\n    {b}" for b in bad[:20]), flush=True)

    a = json.load(open(args.a)); a = a.get("params", a.get("a", a))
    pairs = sorted({os.path.basename(f).split("__")[0] for d in dirs for f in glob.glob(os.path.join(d, "*.npz"))})
    B = BatchV2(dirs, pairs=set(pairs[: args.pairs]))
    rows = []
    for r in B.recs:
        td, tp = B.thetas(r, a)
        rf, rF = B.rho_fast(r, td, tp), B.rho_full(r, td, tp)
        b = B.basis(r); z = np.zeros(len(rf))
        rows.append(dict(proto=r["proto"], pair=r["pair"], ratio_fast=b.ratio(r["rho0"], rf, z),
                         ratio_full=b.ratio(r["rho0"], rF, z), rho_fast=rf.mean(), rho_full=rF.mean(),
                         rho_obs=np.mean(r["rho_obs"]), bin_flip=int(np.sum((rf >= .5) != (rF >= .5)))))
    df = pd.DataFrame(rows)
    df["d_ratio"] = (df.ratio_fast - df.ratio_full).abs()
    g = df.groupby("proto").agg(mean_full=("ratio_full", "mean"), mean_fast=("ratio_fast", "mean"),
                                max_d=("d_ratio", "max"), flips=("bin_flip", "sum"), n=("pair", "count"))
    g["d_mean"] = (g.mean_fast - g.mean_full).abs()
    print("[2] rho fast vs full (protocol mean ratio):\n" + g.round(4).to_string(), flush=True)
    worst = g["d_mean"].max()
    print(f"    worst |mean fast - mean full| = {worst:.4f}  -> {'fast OK' if worst < 0.02 else 'USE --rho full'}")
    if args.same_run:
        d = (df.rho_full - df.rho_obs).abs()
        print(f"[3] offline rho_full vs BCL rho_obs: mean |diff| {d.mean():.4f}, max {d.max():.4f}")
    out = os.path.join(HERE, "results", "validate_" + os.path.basename(dirs[0].rstrip("/")) + ".csv")
    os.makedirs(os.path.dirname(out), exist_ok=True); df.to_csv(out, index=False)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
