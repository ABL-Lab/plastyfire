"""Per-target plot of a fit table (fit_v2 / eval_v2 csv): experiment mean ± sem vs model ± sem across pairs.

    python glusynapse_v2/plot_fit.py glusynapse_v2/results/reduced_gpu_all_ckpt_best_allpairs.csv [--ref other.csv]

Panels: Markram 10 Hz dt curve, Sjostrom frequency (+10 / -10), Nevian 3AP 50 Hz dt curve, Nevian bursts and
single pairings, pharmacology (mGluR block, post NMDAR block), Letzkus, and predicted vs target (z coloured).
"""
import argparse, json, os, re
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def dt_of(name):
    m = re.search(r"dt([+-]\d+)ms", name) or re.search(r"^10Hz_(-?\d+)ms", name)
    return float(m.group(1)) if m else np.nan


LAB = {"model": "model", "ref": "reference"}


def errbar(ax, x, df, ref, dx=0.0, label=True):
    ax.errorbar(x - dx, df.target_mean, df.target_sem, fmt="o", color="k", ms=5, capsize=2, label="experiment" if label else None)
    ax.errorbar(x + dx, df.pred, df.get("pred_sem", 0 * df.pred), fmt="s", color="C3", ms=5, capsize=2,
                label=LAB["model"] if label else None)
    if ref is not None:
        r = ref.set_index(["target", "condition"]).reindex(list(zip(df.target, df.condition)))
        ax.errorbar(x + 2 * dx, r.pred.values, r.get("pred_sem", 0 * r.pred).values, fmt="^", color="C0", ms=5,
                    capsize=2, alpha=0.8, label=LAB["ref"] if label else None)
    ax.axhline(1.0, color="0.7", lw=0.8, zorder=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--ref", default=None, help="second table drawn as blue triangles (e.g. an earlier fit)")
    ap.add_argument("--title", default=None)
    ap.add_argument("--label", default="model")
    ap.add_argument("--ref-label", default="reference")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    LAB.update(model=args.label, ref=args.ref_label)
    T = pd.read_csv(args.csv); ref = pd.read_csv(args.ref) if args.ref else None
    T["dt"] = T.target.map(dt_of)
    chi2 = float(((T.pred - T.target_mean) / T.target_sem).pow(2).sum())
    js = args.csv.replace(".csv", ".json")
    npairs = json.load(open(js)).get("n_pairs", "") if os.path.exists(js) else ""
    ctl = T[T.condition == "control"]

    fig, axs = plt.subplots(2, 4, figsize=(18, 8.5))
    # Markram 10 Hz
    ax = axs[0, 0]; d = ctl[ctl.target.str.match(r"^10Hz_")].sort_values("dt")
    errbar(ax, d.dt.values, d, ref, dx=0.6); ax.set(title="Markram 1997, 10 Hz pairing", xlabel="dt (ms)", ylabel="EPSP ratio")
    # Sjostrom frequency
    ax = axs[0, 1]
    for sgn, c in (("+", "C2"), ("-", "C4")):
        d = ctl[ctl.target.str.match(rf"^sjostrom_.*dt\{sgn}10ms")].copy()
        d["f"] = d.target.str.extract(r"sjostrom_([\d.]+)hz")[0].astype(float); d = d.sort_values("f")
        x = np.arange(len(d)) + (0.15 if sgn == "+" else -0.15)
        errbar(ax, x, d, ref, dx=0.05, label=sgn == "+")
        ax.plot(x, d.target_mean, "-", color=c, lw=1, alpha=0.6); ax.plot(x, d.pred, "--", color=c, lw=1, alpha=0.8)
        ax.set_xticks(np.arange(len(d))); ax.set_xticklabels([f"{f:g}" for f in d.f])
    ax.set(title="Sjöström 2001: +10 ms (green), −10 ms (purple)", xlabel="pairing frequency (Hz)")
    # Nevian 3AP 50 Hz dt curve
    ax = axs[0, 2]; d = ctl[ctl.target.str.match(r"^nevian_3ap_50hz_dt")].sort_values("dt")
    errbar(ax, d.dt.values, d, ref, dx=2); ax.set(title="Nevian 2006, 3AP 50 Hz", xlabel="dt (ms)")
    # Nevian other
    ax = axs[0, 3]
    d = ctl[ctl.target.str.startswith("nevian") & ~ctl.target.str.match(r"^nevian_3ap_50hz_dt")]
    errbar(ax, np.arange(len(d)), d, ref, dx=0.12)
    ax.set_xticks(np.arange(len(d))); ax.set_xticklabels([t.replace("nevian_", "") for t in d.target], rotation=60, ha="right", fontsize=7)
    ax.set(title="Nevian: bursts, 1AP, controls")
    # pharmacology
    ax = axs[1, 0]; d = T[T.condition != "control"]
    errbar(ax, np.arange(len(d)), d, ref, dx=0.12)
    ax.set_xticks(np.arange(len(d)))
    ax.set_xticklabels([f"{t.replace('nevian_', '')}\n{c}" for t, c in zip(d.target, d.condition)], rotation=60, ha="right", fontsize=7)
    ax.set(title="Pharmacology (Nevian)", ylabel="EPSP ratio")
    # Letzkus
    ax = axs[1, 1]; d = ctl[ctl.target.str.startswith("letzkus")].sort_values("dt")
    errbar(ax, d.dt.values, d, ref, dx=0.8); ax.set(title="Letzkus 2006, 3AP 200 Hz (proximal)", xlabel="dt (ms)")
    # predicted vs target
    ax = axs[1, 2]; z = (T.pred - T.target_mean) / T.target_sem
    sc = ax.scatter(T.target_mean, T.pred, c=z.clip(-4, 4), cmap="coolwarm", vmin=-4, vmax=4, s=30, edgecolor="k", lw=0.4)
    ax.errorbar(T.target_mean, T.pred, xerr=T.target_sem, fmt="none", ecolor="0.7", lw=0.8, zorder=0)
    lo, hi = 0.4, 2.5; ax.plot([lo, hi], [lo, hi], "k-", lw=0.8); ax.set(xlim=(lo, hi), ylim=(lo, hi),
                                                                        xlabel="experiment", ylabel="model", title="all targets")
    plt.colorbar(sc, ax=ax, label="z")
    for _, r in T[z.abs() > 2.5].iterrows():
        ax.annotate(r.target.replace("nevian_", "nv_").replace("sjostrom_", "sj_").replace("letzkus_3ap_200hz_", "lz_")[:22],
                    (r.target_mean, r.pred), fontsize=6)
    # z bars
    ax = axs[1, 3]; o = z.abs().sort_values(ascending=False).index[:14]
    ax.barh(np.arange(len(o)), z[o], color=["C3" if v > 0 else "C0" for v in z[o]])
    ax.set_yticks(np.arange(len(o))); ax.set_yticklabels([f"{T.target[i][:28]} {T.condition[i][:4]}" for i in o], fontsize=7)
    ax.invert_yaxis(); ax.axvline(0, color="k", lw=0.8); ax.set(title="largest misses (z)", xlabel="z = (model − exp) / sem")
    fig.suptitle(args.title or f"{os.path.basename(args.csv)}: χ² {chi2:.1f} over {len(T)} targets"
                 + (f", {npairs} pairs" if npairs else ""), fontsize=13)
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.955), fontsize=10, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = args.out or args.csv.replace(".csv", ".png")
    fig.savefig(out, dpi=110); print(out, f"chi2 {chi2:.2f}")


if __name__ == "__main__":
    main()
