"""STDP curve from full-protocol bluecellulab runs (no basis).

Reads simulation_full_<param_hash>.pkl written by
    run_de_fit2_pool.py --full ...
and takes each run's measured EPSP ratio, mean(C02 EPSP) / mean(C01 EPSP) over
all test pulses. Aggregates per dt (mean +/- SEM over pairs), compares with
Markram 1997 and, with --basis-dir, with what the basis would have predicted
from the same run's rho0 -> snapped rho (a direct check of the basis shortcut).

    python curate_stdp_full.py --results-dir <.../simulations> --param-hash fixed_beam \
        --basis-dir basis_results_edges_sabrina_n120 --out stdp_full_fixed_beam.png
"""
import argparse, glob, os, pickle, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
INVITRO = {-10: (0.7922, 0.0259), 5: (1.2038, 0.0644), 10: (1.2013, 0.0626)}


def collect(results_dir, param_hash, basis_dir=None):
    if basis_dir:
        os.environ["ANALYTICAL_BASIS_DIR"] = os.path.abspath(basis_dir)
        sys.path.insert(0, os.path.join(HERE, "analytical_method"))
        import model as M
    rows, bad = [], []
    for f in sorted(glob.glob(os.path.join(results_dir, "*", "*",
                                           f"simulation_full_{param_hash}.pkl"))):
        proto = os.path.basename(os.path.dirname(f))
        pair = os.path.basename(os.path.dirname(os.path.dirname(f)))
        dt = int(proto.split("_")[1].replace("ms", ""))
        r = pickle.load(open(f, "rb"))
        if r.get("epsp_ratio") is None:
            bad.append((pair, proto, r.get("epsp_error", "no ratio")))
            continue
        row = dict(pair=pair, protocol=proto, dt=dt, ratio=r["epsp_ratio"],
                   epsp_before=float(np.mean(r["epsp_c01"])),
                   epsp_after=float(np.mean(r["epsp_c02"])),
                   n_c01=len(r["epsp_c01"]), n_c02=len(r["epsp_c02"]),
                   guardrail_forced=r.get("guardrail_forced", False))
        if basis_dir and "snapped_rho" in r:
            pre, post = pair.split("-")
            try:
                row["basis_ratio"] = M.epsp_ratio(M.load_basis(pre, post),
                                                  np.asarray(r["initial_rho"]),
                                                  np.asarray(r["snapped_rho"]))
            except Exception:
                row["basis_ratio"] = np.nan
        rows.append(row)
    return pd.DataFrame(rows), bad


def summarise(df, col):
    g = df.groupby("dt")[col]
    return pd.DataFrame({"mean": g.mean(), "sem": g.sem(), "n": g.count()}).reset_index()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--param-hash", required=True)
    ap.add_argument("--basis-dir", default=None,
                    help="also report the basis prediction from the same runs' rho")
    ap.add_argument("--label", default=None)
    ap.add_argument("--out", default=None, help="png path (CSV written next to it)")
    args = ap.parse_args()

    df, bad = collect(args.results_dir, args.param_hash, args.basis_dir)
    print(f"param_hash : {args.param_hash}")
    print(f"runs       : {len(df)} with a measured ratio, {len(bad)} unusable")
    for b in bad[:20]:
        print(f"   skip {b[0]}/{b[1]}: {b[2]}")
    if df.empty:
        sys.exit("no simulation_full pkls with a ratio")

    s = summarise(df, "ratio")
    sb = summarise(df, "basis_ratio") if "basis_ratio" in df else None
    print(f"\n  {'dt':>5} {'measured':>16} {'n':>4} {'basis(same rho)':>16} {'in vitro':>16}")
    for _, r in s.iterrows():
        dt = int(r["dt"])
        b = (f"{sb.loc[sb['dt'] == dt, 'mean'].values[0]:.4f}" if sb is not None else "-")
        iv = f"{INVITRO[dt][0]:.4f} ± {INVITRO[dt][1]:.4f}" if dt in INVITRO else "-"
        print(f"  {dt:>5} {r['mean']:>8.4f} ± {r['sem']:.4f} {int(r.n):>4} {b:>16} {iv:>16}")

    out = args.out or f"stdp_full_{args.param_hash}.png"
    csv = os.path.splitext(out)[0] + ".csv"
    df.to_csv(csv, index=False)
    print(f"\nper-run table -> {csv}")

    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    ax.errorbar(s["dt"], s["mean"], yerr=s["sem"], fmt="-o", capsize=3, color="#3a78c2",
                label=args.label or f"{args.param_hash} (full protocol, BCL)")
    if sb is not None:
        ax.errorbar(sb["dt"], sb["mean"], yerr=sb["sem"], fmt=":s", capsize=2, ms=4,
                    color="#999999", alpha=0.8, label="basis prediction, same rho")
    iv = sorted(INVITRO)
    ax.errorbar(iv, [INVITRO[d][0] for d in iv], yerr=[INVITRO[d][1] for d in iv],
                fmt="--^", capsize=3, color="#f08a24", label="In vitro (Markram 1997)")
    ax.axhline(1.0, ls="--", c="gray", lw=0.8); ax.axvline(0, ls="--", c="gray", lw=0.8)
    ax.set_xlabel(r"$\Delta t$ (ms)"); ax.set_ylabel("EPSP ratio")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(out, dpi=150)
    print(f"plot -> {out}")


if __name__ == "__main__":
    main()
