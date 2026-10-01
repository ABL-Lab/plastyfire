"""Fitted vs experimental EPSP ratio for one fit_v2 result (per-target csv + json).

usage: python plot_fit_vs_data.py results/<save>  [--out figs/<name>.pdf]
"""
import argparse, json, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
STYLE = os.path.join(HERE, "..", "onerule.mplstyle")
COND_COLOR = {"control": "#68a8e0", "mglu_block": "#f9a24b", "nmdar_block": "#9ad175", "post_nmdar": "#ff4747"}
COND_LABEL = {"control": "control", "mglu_block": "AM251 / ifenprodil (T block)",
              "nmdar_block": "APV (NMDAR block)", "post_nmdar": "post NMDAR block"}


def group(t, cond):
    if t.startswith("10Hz"):
        return 0, "Markram 1997"
    if cond != "control":
        return 3, "Sjöström 2003\ndrugs"
    if "burst" in t or t.split("_dt")[1] not in ("+10ms", "-10ms"):
        return 2, "Sjöström 2003\ntiming"
    return 1, "Sjöström 2001"


def label(t):
    if t.startswith("10Hz_"):
        return "10 Hz " + t[5:].replace("ms", " ms").replace("_", " ")
    s = t.replace("sjostrom_", "").replace("_dt", " ").replace("ms", " ms")
    return s.replace("burst5x20hz", "burst 5x20Hz").replace("hz", " Hz")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("save", help="fit_v2 --save prefix (reads .csv and .json)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    df = pd.read_csv(args.save + ".csv")
    meta = json.load(open(args.save + ".json"))
    df["g"], df["gname"] = zip(*[group(t, c) for t, c in zip(df["target"], df["condition"])])
    df["freq"] = df["target"].str.extract(r"_(\d+\.?\d*)hz")[0].astype(float).fillna(10.0)
    df["dt"] = df["target"].str.extract(r"(-?\+?\d+)ms$")[0].astype(float)
    df["ci"] = df["condition"].map(list(COND_COLOR).index)
    df = df.sort_values(["g", "ci", "freq", "dt"]).reset_index(drop=True)

    plt.style.use(STYLE)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.0, 3.9), gridspec_kw=dict(width_ratios=[2.6, 1]))

    # a: per target, data mean +- SEM (black) and fitted value (colour = condition)
    x = np.arange(len(df))
    ax.axhline(1.0, color="0.6", lw=0.5, ls="--", zorder=0)
    ax.errorbar(x - 0.15, df["target_mean"], yerr=df["target_sem"], fmt="o", color="k", mfc="white",
                capsize=1.5, label="experiment (mean ± SEM)")
    for c, sub in df.groupby("condition"):
        ax.plot(sub.index + 0.15, sub["pred"], "D", color=COND_COLOR[c], ms=3.5, label=f"fit: {COND_LABEL[c]}")
    for b in df.index[df["g"].diff().fillna(0) != 0]:
        ax.axvline(b - 0.5, color="0.8", lw=0.5)
    for g, sub in df.groupby("g"):
        ax.text(sub.index.to_numpy().mean(), 2.12, sub["gname"].iloc[0], ha="center", va="bottom", fontsize=6)
    ax.set_xticks(x)
    ax.set_xticklabels([label(t) for t in df["target"]], rotation=70, ha="right")
    ax.set_xlim(-0.7, len(df) - 0.3)
    ax.set_ylim(0.3, 2.3)
    ax.set_ylabel("EPSP ratio (after / before)")
    ax.legend(loc="upper center", ncol=4, bbox_to_anchor=(0.5, -0.42), fontsize=6)
    ax.set_title("a", loc="left", fontweight="bold")

    # b: fitted vs experiment, with |z| > 2 marked
    lim = (0.3, 2.2)
    bx.plot(lim, lim, color="0.6", lw=0.5, ls="--")
    for c, sub in df.groupby("condition"):
        bx.errorbar(sub["target_mean"], sub["pred"], xerr=sub["target_sem"], fmt="o", color=COND_COLOR[c],
                    capsize=1.2, elinewidth=0.4)
    bad = df[df["z"].abs() > 2]
    bx.plot(bad["target_mean"], bad["pred"], "o", mfc="none", mec="#ff4747", ms=7, mew=0.7)
    for (m, pr), sub in bad.groupby(["target_mean", "pred"]):
        txt = " / ".join(label(t) for t in sub["target"]).replace("burst 5x20Hz -200 ms / burst 5x20Hz -120 ms",
                                                                   "burst 5x20Hz -120/-200 ms")
        bx.annotate(txt, (m, pr), xytext=(6, -3), textcoords="offset points", fontsize=5, color="#ff4747")
    bx.set_xlim(lim); bx.set_ylim(lim); bx.set_aspect("equal")
    bx.set_xlabel("experiment"); bx.set_ylabel("fit")
    chi2 = float((df["z"] ** 2).sum())
    bx.text(0.35, 2.12, f"χ² = {chi2:.1f} / {len(df)} targets\n"
            f"|z| > 2: {len(bad)}   median |z| = {df['z'].abs().median():.2f}", va="top", fontsize=6)
    bx.set_title("b", loc="left", fontweight="bold")

    fig.suptitle(f"GluSynapse_v2 fit to paired L5→L5 data ({os.path.basename(args.save)})", fontsize=7)
    out = args.out or os.path.join(HERE, "figs", os.path.basename(args.save) + "_fit_vs_data.pdf")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out); fig.savefig(out.replace(".pdf", ".png"), dpi=300)
    print(out)
    print(df[["target", "condition", "target_mean", "target_sem", "pred", "z"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
