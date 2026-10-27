"""Neurodamus vs bluecellulab prefire agreement for one param set.

Per workdir: neurodamus <workdir>/out_<tag>/rho.h5 against bluecellulab
<workdir>/simulation_edges_<tag>.pkl. Synapses are matched by edge id (the ND
report's element_ids, sorted, against the pkl's global_ids, sorted), rho is
binarised at 0.5, and the per-synapse agreement plus the fraction of
potentiated / depressed synapses per dt are reported for both simulators.
With --basis-dir, each run's rho0 -> final rho also goes through the EPSP basis,
giving an ND and a BCL STDP curve (mean +/- SEM per dt) and a plot.

    python compare_nd_bcl_prefire.py --tag delta-cooker --basis-dir basis_results_edges_sabrina_n120_delta
"""
import argparse, glob, os, pickle, sys
import numpy as np, pandas as pd, h5py

HERE = os.path.dirname(os.path.abspath(__file__))
SIMS = os.path.join(HERE, "refitting_results/fitting/n120/seed20262009/"
                          "Sabrina_L5TTPC_L5TTPC_STDP/simulations")


def read_nd(path):
    with h5py.File(path, "r") as f:
        (pop,) = list(f["report"])
        g = f[f"report/{pop}"]
        data = g["data"][()]
        eids = g["mapping/element_ids"][()]
    o = np.argsort(eids)
    return eids[o], data[0, o], data[-1, o]


INVITRO = {-10: (0.7922, 0.0259), 5: (1.2038, 0.0644), 10: (1.2013, 0.0626)}
M = None   # analytical_method.model, imported once --basis-dir sets its env


def compare(workdir, tag):
    nd = os.path.join(workdir, f"out_{tag}", "rho.h5")
    bc = os.path.join(workdir, f"simulation_edges_{tag}.pkl")
    if not (os.path.isfile(nd) and os.path.isfile(bc)):
        return None
    r = pickle.load(open(bc, "rb"))
    gid = np.asarray(r["global_ids"]); o = np.argsort(gid)
    b0, b1 = np.asarray(r["initial_rho"])[o] >= 0.5, np.asarray(r["final_rho"])[o] >= 0.5
    eids, n0, n1 = read_nd(nd)
    n0, n1 = n0 >= 0.5, n1 >= 0.5
    if len(eids) != len(gid):
        return dict(error=f"n_syn nd={len(eids)} bcl={len(gid)}")
    out = dict(n_syn=len(gid), same_init=int((n0 == b0).sum()), same_final=int((n1 == b1).sum()),
               nd_pot=int((~n0 & n1).sum()), nd_dep=int((n0 & ~n1).sum()),
               bcl_pot=int((~b0 & b1).sum()), bcl_dep=int((b0 & ~b1).sum()))
    if M is not None:
        pre, post = os.path.basename(os.path.dirname(workdir)).split("-")
        try:   # basis rows are in ascending synapse-id order, like both sorted vectors here
            basis = M.load_basis(int(pre), int(post))
            out["nd_ratio"] = M.epsp_ratio(basis, n0.astype(float), n1.astype(float))
            out["bcl_ratio"] = M.epsp_ratio(basis, b0.astype(float), b1.astype(float))
        except Exception as e:
            out["basis_error"] = str(e)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--results-dir", default=SIMS)
    ap.add_argument("--csv", default=None)
    ap.add_argument("--basis-dir", default=None, help="also compute basis STDP curves for ND and BCL")
    ap.add_argument("--out", default=None, help="png for the curves (default compare_nd_bcl_<tag>.png)")
    a = ap.parse_args()
    global M
    if a.basis_dir:
        os.environ["ANALYTICAL_BASIS_DIR"] = os.path.abspath(a.basis_dir)
        sys.path.insert(0, os.path.join(HERE, "analytical_method"))
        import model
        M = model
    rows = []
    for w in sorted(glob.glob(os.path.join(a.results_dir, "*-*", "10Hz_*ms"))):
        c = compare(w, a.tag)
        if c is None:
            continue
        rows.append(dict(pair=os.path.basename(os.path.dirname(w)),
                         dt=int(os.path.basename(w).split("_")[1][:-2]), **c))
    df = pd.DataFrame(rows)
    if df.empty:
        raise SystemExit("no workdir has both outputs")
    bad = df[df.get("error").notna()] if "error" in df else df.iloc[:0]
    ok = df.drop(bad.index)
    print(f"tag {a.tag}: {len(ok)} workdirs compared, {len(bad)} with mismatched synapse counts")
    for _, b in bad.head(10).iterrows():
        print(f"   {b.pair} dt={b.dt}: {b.error}")
    n = ok.n_syn.sum()
    print(f"synapses: {n}  initial agree {ok.same_init.sum()/n:.4f}  final agree {ok.same_final.sum()/n:.4f}")
    g = ok.groupby("dt")[["n_syn", "same_final", "nd_pot", "bcl_pot", "nd_dep", "bcl_dep"]].sum()
    g["final_agree"] = g.same_final / g.n_syn
    print(g.to_string())
    out = a.csv or os.path.join(HERE, f"compare_nd_bcl_{a.tag}.csv")
    df.to_csv(out, index=False)
    print(f"per-workdir table -> {out}")
    if M is None or "nd_ratio" not in ok:
        return
    c = ok.dropna(subset=["nd_ratio", "bcl_ratio"])
    s = c.groupby("dt")[["nd_ratio", "bcl_ratio"]].agg(["mean", "sem", "count"])
    print("\nbasis STDP curve (ND vs BCL, same pairs):")
    print(f"  {'dt':>5} {'neurodamus':>16} {'bluecellulab':>16} {'n':>4} {'in vitro':>16}")
    for dt, r in s.iterrows():
        iv = f"{INVITRO[dt][0]:.4f} ± {INVITRO[dt][1]:.4f}" if dt in INVITRO else "-"
        print(f"  {dt:>5} {r[('nd_ratio','mean')]:>8.4f} ± {r[('nd_ratio','sem')]:.4f} "
              f"{r[('bcl_ratio','mean')]:>8.4f} ± {r[('bcl_ratio','sem')]:.4f} "
              f"{int(r[('nd_ratio','count')]):>4} {iv:>16}")
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    for col, lab, colr in (("nd_ratio", "neurodamus", "#c23a3a"), ("bcl_ratio", "bluecellulab", "#3a78c2")):
        ax.errorbar(s.index, s[(col, "mean")], yerr=s[(col, "sem")], fmt="-o", capsize=3,
                    color=colr, label=f"{lab} (prefire + basis)")
    iv = sorted(INVITRO)
    ax.errorbar(iv, [INVITRO[d][0] for d in iv], yerr=[INVITRO[d][1] for d in iv],
                fmt="--^", capsize=3, color="#f08a24", label="In vitro (Markram 1997)")
    ax.axhline(1.0, ls="--", c="gray", lw=0.8); ax.axvline(0, ls="--", c="gray", lw=0.8)
    ax.set_xlabel(r"$\Delta t$ (ms)"); ax.set_ylabel("EPSP ratio"); ax.set_title(a.tag)
    ax.legend(frameon=False, fontsize=8); fig.tight_layout()
    png = a.out or os.path.join(HERE, f"compare_nd_bcl_{a.tag}.png")
    fig.savefig(png, dpi=150)
    print(f"plot -> {png}")


if __name__ == "__main__":
    main()
