"""Compare several fit_v2 results against the same targets (experiment mean +- SEM vs each fit).

usage: python plot_fit_compare.py label=results/<save> label=results/<save> ... [--out figs/x.pdf]
"""
import argparse, json, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_fit_vs_data import STYLE, group, label

HERE = os.path.dirname(os.path.abspath(__file__))
COLORS = ["#68a8e0", "#f9a24b", "#9ad175", "#ff4747"]
MARKERS = ["D", "s", "^", "v"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fits", nargs="+", help="label=save_prefix")
    ap.add_argument("--out", default=os.path.join(HERE, "figs", "fit_compare.pdf"))
    args = ap.parse_args()
    fits = [f.split("=", 1) for f in args.fits]
    base = None
    for lab, save in fits:
        d = pd.read_csv(save + ".csv")[["target", "condition", "target_mean", "target_sem", "pred", "z"]]
        d = d.rename(columns={"pred": f"pred_{lab}", "z": f"z_{lab}"})
        base = d if base is None else base.merge(d, on=["target", "condition", "target_mean", "target_sem"])
    base["g"], base["gname"] = zip(*[group(t, c) for t, c in zip(base["target"], base["condition"])])
    base["freq"] = base["target"].str.extract(r"_(\d+\.?\d*)hz")[0].astype(float).fillna(10.0)
    base["dt"] = base["target"].str.extract(r"(-?\+?\d+)ms$")[0].astype(float)
    base["ci"] = base["condition"].map(["control", "mglu_block", "nmdar_block", "post_nmdar"].index)
    df = base.sort_values(["g", "ci", "freq", "dt"]).reset_index(drop=True)

    plt.style.use(STYLE)
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(7.0, 5.2), sharex=True, gridspec_kw=dict(height_ratios=[2.2, 1]))
    x = np.arange(len(df))
    w = 0.6 / len(fits)
    ax.axhline(1.0, color="0.6", lw=0.5, ls="--", zorder=0)
    ax.errorbar(x - 0.35, df["target_mean"], yerr=df["target_sem"], fmt="o", color="k", mfc="white", capsize=1.5,
                label="experiment (mean ± SEM)")
    for i, (lab, _) in enumerate(fits):
        chi2 = float((df[f"z_{lab}"] ** 2).sum())
        ax.plot(x - 0.15 + i * w, df[f"pred_{lab}"], MARKERS[i], color=COLORS[i], ms=3.2, label=f"{lab} (χ² {chi2:.1f})")
        bx.bar(x - 0.15 + i * w, df[f"z_{lab}"], width=w, color=COLORS[i])
    for axx in (ax, bx):
        for b in df.index[df["g"].diff().fillna(0) != 0]:
            axx.axvline(b - 0.5, color="0.8", lw=0.5)
    for g, sub in df.groupby("g"):
        ax.text(sub.index.to_numpy().mean(), 2.12, sub["gname"].iloc[0], ha="center", va="bottom", fontsize=6)
    ax.set_ylim(0.3, 2.3); ax.set_ylabel("EPSP ratio (after / before)")
    ax.legend(loc="upper left", ncol=2, fontsize=6, bbox_to_anchor=(0.0, 0.9))
    for yy in (-2, 2):
        bx.axhline(yy, color="#ff4747", lw=0.5, ls=":")
    bx.axhline(0, color="0.5", lw=0.5)
    bx.set_ylabel("z = (fit − exp) / SEM")
    lbl = [label(t) + ("" if c == "control" else {"mglu_block": " AM251", "nmdar_block": " APV"}.get(c, " " + c))
           for t, c in zip(df["target"], df["condition"])]
    bx.set_xticks(x); bx.set_xticklabels(lbl, rotation=70, ha="right")
    bx.set_xlim(-0.7, len(df) - 0.3)
    ax.set_title("a", loc="left", fontweight="bold"); bx.set_title("b", loc="left", fontweight="bold")
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out); fig.savefig(args.out.replace(".pdf", ".png"), dpi=300)
    print(args.out)
    cols = ["target", "condition", "target_mean"] + [f"pred_{l}" for l, _ in fits] + [f"z_{l}" for l, _ in fits]
    print(df[cols].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
