"""Score delta-bd screen results (screen_bd.py csvs) against the L5 basal targets; og-delta measured the same way.

Per cell, then median over cells (hinge penalties, 0 = target met):
  (all Ca ratios are per somatic AP: the soma drops pulses at high rates)
  kampa_R100   distal Ca(5 APs @100 Hz) / Ca(@66 Hz) = 5 +- 1.2          (Kampa & Stuart 2006, fig 1g)
  kampa_Rhi    same at 150 and 200 Hz >= 3.8 (supralinear above the critical frequency)
  R3_dist      Ca(3 APs @200 Hz) / Ca(1 AP), distal >= 130 um: >= 4.5   (supralinear distally)
  R3_prox      same, proximal < 75 um: <= 3.5                           (not proximally)
  bap200       single bAP dV at >= 200 um <= 30 mV                      (Nevian 2007; Kampa voltage imaging)
  fourap       4-AP (basal Ka off) raises distal 1-AP Ca more than proximal: ratio(dist)/ratio(prox) >= 1.5
  nickel       Ni (basal LVA off) removes distal supralinearity: R3_dist <= 3.3
  follow       fraction of cells firing 3/3 for the protocol pulses at 50, 100 (Nevian) and 200 Hz (Letzkus): 1
  fI           somatic spike counts within 15% of og-delta on the same cell (3 steps)

    python glusynapse_v2/emodel_bd/score_bd.py 'glusynapse_v2/emodel_bd/results/s1_*.csv' [--top 15]
"""
import argparse, glob
import numpy as np, pandas as pd

# kampa_R100 and nickel carry no weight (s1: no basal-channel candidate moves per-AP R_f100 off 1.1-1.4, and Ca ratios with
# LVA off rise because the 1-AP denominator is mostly LVA Ca; both still printed)
W = dict(kampa_R100=0.0, kampa_Rhi=1.0, R3_dist=1.0, R3_prox=1.0, bap200=1.0, fourap=0.5, nickel=0.0, follow=2.0, fI=2.0)


def per_cell(df):
    og = df[df.cand == "og"].set_index("pair")
    d = df.copy()
    # per-AP normalisation: the soma drops pulses at >= 150 Hz (and some at 66 Hz), Kampa compares equal AP counts
    nz = lambda c: d[c].where(d[c] > 0)
    for f in ("f100", "f150", "f200"):
        d["R_" + f] = (d[f"ca_{f}_dist"] / nz(f"nsp_{f}")) / (d["ca_f66_dist"] / nz("nsp_f66"))
    for g in ("prox", "mid", "dist"):
        d["R3_" + g] = 3 * d[f"ca_r3_{g}"] / (nz("nsp_r3") * d[f"ca_ap1_{g}"])
        if f"lva_off_ca_r3_{g}" in d:
            d["lva_off_R3_" + g] = 3 * d[f"lva_off_ca_r3_{g}"] / (nz("lva_off_nsp_r3") * d[f"lva_off_ca_ap1_{g}"])
    d["p_kampa_R100"] = ((d.R_f100 - 5.0) / 1.2).abs().clip(lower=1) - 1               # 0 inside +-1.2
    d["p_kampa_Rhi"] = np.clip((3.8 - d[["R_f150", "R_f200"]].min(axis=1)) / 1.2, 0, None)
    d["p_R3_dist"] = np.clip((4.5 - d.R3_dist) / 1.0, 0, None)
    d["p_R3_prox"] = np.clip((d.R3_prox - 3.5) / 1.0, 0, None)
    d["p_bap200"] = np.clip((d.bap_dv_d200 - 30.0) / 5.0, 0, None)
    fa = (d.ka_off_ca_ap1_dist / d.ca_ap1_dist) / (d.ka_off_ca_ap1_prox / d.ca_ap1_prox)
    d["fourap_ratio"] = fa
    d["p_fourap"] = np.clip((1.5 - fa) / 0.5, 0, None)
    d["p_nickel"] = np.clip((d.lva_off_R3_dist - 3.3) / 1.0, 0, None)
    d["n_follow"] = (d.nsp_n50 >= 3).astype(int) + (d.nsp_n100 >= 3) + (d.nsp_l200 >= 3)
    d["p_follow"] = (3 - d.n_follow) / 3.0
    fi = [c for c in d.columns if c.startswith("fi_")]
    ref = og.loc[d.pair, fi].to_numpy()
    rel = np.abs(d[fi].to_numpy() - ref) / np.maximum(ref, 3)
    d["p_fI"] = np.clip(rel.max(1) - 0.15, 0, None) / 0.15
    return d


def score(d):
    keys = ["p_" + k for k in W]
    g = d.groupby("cand")
    s = g[keys].median()
    s["p_follow"] = 1 - g.n_follow.apply(lambda x: (x == 3).mean())                        # fraction of cells, not median
    s["score"] = sum(W[k] * s["p_" + k] for k in W)
    params = [c for c in ("na", "ka0", "lam", "m", "dh") if c in d]
    med = g[["R_f100", "R_f200", "R3_prox", "R3_dist", "bap_dv_d200", "bap_dv_d50_100", "fourap_ratio",
             "lva_off_R3_dist", "nsp_f66", "nsp_f100", "nsp_f200", "nsp_r3"]].median().add_prefix("med_")
    return pd.concat([s, med, g[params].first(), g.pair.nunique().rename("ncell")], axis=1).sort_values("score")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("glob"); ap.add_argument("--top", type=int, default=15); ap.add_argument("--out")
    a = ap.parse_args()
    df = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(a.glob))], ignore_index=True)
    S = score(per_cell(df))
    pd.set_option("display.width", 260, "display.max_columns", 40)
    print(S.head(a.top).round(3).to_string())
    print("\nreferences:"); print(S.loc[[c for c in ("og", "antic") if c in S.index]].round(3).to_string())
    if a.out:
        S.to_csv(a.out)


if __name__ == "__main__":
    main()
