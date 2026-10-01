"""Overnight ion-channel summary figure -> glusynapse_v2/results/ion_channels_overnight.png
(a) synapse-level burst gain per variant (delta-bd basal candidates, spine-VDCC variants), og-delta pairs, spike-matched
(b) subset24 reduced-fit chi2 per protocol group, og-delta vs delta-sv (best of 3 seeds each)
(c) per-synapse peak effcai / c_pre for the protocols that conflict, og-delta vs delta-sv (24 pairs, 191 synapses)"""
import os
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

V2 = "glusynapse_v2"
b = pd.read_csv(f"{V2}/emodel_bd/results/syn1_compare.csv"); s = pd.read_csv(f"{V2}/emodel_bd/results/sp1_compare.csv")
fig, ax = plt.subplots(1, 3, figsize=(18, 4.8))
P = "3ap_50hz@+10"
rows = []
for df, kind in ((b, "basal (delta-bd)"), (s, "spine VDCC")):
    m = df[df.full & (df.proto == P)].groupby("variant").gain.median()
    for v, g in m.items():
        if kind == "spine VDCC" and v == "og":
            continue
        rows.append((kind, v, g))
r = pd.DataFrame(rows, columns=["kind", "v", "g"])
col = r.kind.map({"basal (delta-bd)": "C0", "spine VDCC": "C1"})
ax[0].bar(range(len(r)), r.g, color=col)
ax[0].set_xticks(range(len(r)), r.v, rotation=45, fontsize=8)
ax[0].axhline(2.01 / 1.04, color="r", ls="--"); ax[0].text(0, 1.95, "Nevian target 1.93", color="r", fontsize=8)
ax[0].set_ylim(1, 2.1); ax[0].set_ylabel("effcai peak, 3AP50+10 / 1AP+10 (median over pairs)")
ax[0].set_title("(a) burst gain at real synapses\nblue: basal Na/Ka/LVA, orange: spine VDCC (s10g3 = delta-sv)", fontsize=9)

c = pd.read_csv(f"{V2}/results/spine_sv24_vs_og24.csv").groupby("group")[["z2_og", "z2_sv"]].sum()
x = np.arange(len(c))
ax[1].bar(x - .2, c.z2_og, .4, label=f"og-delta (total {c.z2_og.sum():.1f})")
ax[1].bar(x + .2, c.z2_sv, .4, label=f"delta-sv (total {c.z2_sv.sum():.1f})")
ax[1].set_xticks(x, c.index); ax[1].set_ylabel("chi2"); ax[1].legend(fontsize=8)
ax[1].set_title("(b) reduced v2 refit, subset24, best of 3 seeds", fontsize=9)

pairs = open(f"{V2}/subset24_pairs.txt").read().strip().split(",")
PR = [("nevian_1ap_dt+10ms", "Nevian 1AP+10\n(1.04)"), ("nevian_3ap_50hz_dt+10ms", "Nevian 3AP50+10\n(LTP 2.01)"),
      ("sjostrom_20hz_dt-10ms", "Sjostrom 20Hz -10\n(LTD 0.65)"), ("nevian_3ap_50hz_dt-50ms", "Nevian 3AP50 -50\n(LTD 0.68)")]
for j, (d, lab, cc) in enumerate((("ebner_delta-prefire", "og-delta", "C0"), ("ebner_delta-sv", "delta-sv", "C1"))):
    data = []
    for p, _ in PR:
        v = []
        for q in pairs:
            f = f"{V2}/extracted/{d}/{q}__{p}.npz"
            if os.path.exists(f):
                z = np.load(f, allow_pickle=True); v += list(z["effcai"].max(1) / z["c_pre"])
        data.append(np.log10(v))
    bp = ax[2].boxplot(data, positions=np.arange(len(PR)) + (j - .5) * .35, widths=.3, showfliers=False,
                       patch_artist=True)
    for bx in bp["boxes"]:
        bx.set_facecolor(cc); bx.set_alpha(.6)
    ax[2].plot([], [], color=cc, lw=6, alpha=.6, label=lab)
ax[2].set_xticks(range(len(PR)), [l for _, l in PR], fontsize=8)
ax[2].set_ylabel("log10 peak effcai / c_pre (per synapse)"); ax[2].legend(fontsize=8)
ax[2].set_title("(c) per-synapse Ca drive: Sjostrom LTD protocol > Nevian LTP protocol", fontsize=9)
fig.suptitle("Overnight ion-channel test: basal channels do not reach the synapses; the spine VDCC raises the burst gain "
             "but the fit does not improve (Ca-amplitude conflict between protocols)", fontsize=10)
fig.tight_layout(); fig.savefig(f"{V2}/results/ion_channels_overnight.png", dpi=110)
print("saved")
