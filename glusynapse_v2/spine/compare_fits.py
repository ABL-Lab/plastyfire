"""og-delta vs delta-sv reduced v2 fits on subset24: chi2 per protocol group and per target, best seed of each.

    python glusynapse_v2/spine/compare_fits.py [--og reduced_gpu_subset_s1,..] [--sv reduced_gpu_subset_sv24_s1,..]
"""
import argparse, json, os
import numpy as np, pandas as pd

R = "glusynapse_v2/results"


def group(t):
    return "markram" if t[:1].isdigit() else t.split("_")[0]


def load(names):
    fits = []
    for n in names.split(","):
        if os.path.exists(f"{R}/{n}.json"):
            fits.append((json.load(open(f"{R}/{n}.json"))["chi2"], n))
    chi, best = min(fits)
    df = pd.read_csv(f"{R}/{best}.csv")
    return best, [round(c, 2) for c, _ in sorted(fits, key=lambda x: x[1])], df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--og", default="reduced_gpu_subset_s1,reduced_gpu_subset_s2,reduced_gpu_subset_s3")
    ap.add_argument("--sv", default="reduced_gpu_subset_sv24_s1,reduced_gpu_subset_sv24_s2,reduced_gpu_subset_sv24_s3")
    ap.add_argument("--out", default=f"{R}/spine_sv24_vs_og24.csv")
    a = ap.parse_args()
    out = {}
    for k in ("og", "sv"):
        best, chis, df = load(getattr(a, k))
        print(f"{k}: chi2 per seed {chis}  best {best}")
        df["group"] = df.target.map(group)
        out[k] = df.set_index(["target", "condition"])
    m = out["og"][["group", "target_mean", "target_sem", "pred", "z"]].join(out["sv"][["pred", "z"]], rsuffix="_sv")
    m = m.rename(columns={"pred": "pred_og", "z": "z_og"})
    m["z2_og"], m["z2_sv"] = m.z_og ** 2, m.z_sv ** 2
    print("\nchi2 by group (best seed)")
    print(m.groupby("group")[["z2_og", "z2_sv"]].sum().round(2).to_string())
    pd.set_option("display.width", 200, "display.max_rows", 100)
    print("\nper target")
    print(m[["group", "target_mean", "target_sem", "pred_og", "pred_sv", "z_og", "z_sv"]].round(3).to_string())
    m.to_csv(a.out)


if __name__ == "__main__":
    main()
