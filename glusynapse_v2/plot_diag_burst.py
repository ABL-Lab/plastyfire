"""og-delta vs antic-delta burst dissection figure, from results/diag_burst_antic_matched.csv (runs where both models fire
the full somatic AP count; written by the matching step in DECISIONS 2026-09-28) -> results/diag_burst_antic.png"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

df = pd.read_csv("glusynapse_v2/results/diag_burst_antic_matched.csv")
V, C, L = ["delta", "antic_delta"], ["C0", "C1"], ["og-delta", "antic-delta"]
P = ["3ap_50hz@+10", "3ap_50hz@-10", "3ap_100hz@+10", "3ap_100hz@-10"]
fig, ax = plt.subplots(1, 4, figsize=(17, 4.2))
one = df.drop_duplicates(["pair", "v"])
for j, v in enumerate(V):                                              # bAP peak at the synapse, single AP
    y = one[one.v == v].vpk_1ap.dropna()
    ax[0].bar(j, y.median(), color=C[j], alpha=.5); ax[0].scatter(np.full(len(y), j), y, s=16, color="k")
ax[0].set_xticks(range(2), L); ax[0].set_title("bAP peak at synapse, 1 AP (mV)", fontsize=10)
for i, (k, t) in enumerate([("vshare", "VDCC share of synaptic Ca in the pairing"),
                            ("gain", "effcai peak: 3AP pairing / 1AP pairing (same dt)")], start=1):
    a = ax[i]
    for j, v in enumerate(V):
        m = df[df.v == v].groupby("proto")[k].median().reindex(P)
        a.bar(np.arange(4) + (j - .5) * .38, m, .38, color=C[j], alpha=.7, label=L[j])
        d = df[df.v == v]
        a.scatter([P.index(p) + (j - .5) * .38 for p in d.proto], d[k], s=10, color="k", zorder=3)
    a.set_xticks(range(4), [p.replace("3ap_", "3AP ").replace("hz@", " Hz ") + " ms" for p in P], rotation=20, fontsize=8)
    a.set_title(t, fontsize=10)
ax[1].legend(fontsize=8)
ax[2].axhline(2.01 / 1.04, color="r", ls="--"); ax[2].text(-.4, 1.96, "Nevian +10: LTP 2.01 vs 1.04", color="r", fontsize=8)
ax[2].axhline(1, color=".6", lw=.8); ax[2].set_ylim(.9, 2.2)
for j, v in enumerate(V):                                              # the gain vs synapse-to-synapse spread
    d = df[(df.v == v) & (df.proto == "3ap_50hz@+10")]
    ax[3].scatter(d.spread, d.gain, color=C[j], label=L[j])
ax[3].plot([1, 6], [1, 6], color=".6", lw=.8); ax[3].set_xlim(1, 6); ax[3].set_ylim(1, 6)
ax[3].set_xlabel("within-pair spread of 1AP+10 effcai (max/min over synapses)"); ax[3].set_ylabel("burst gain, 3AP50+10")
ax[3].set_title("burst gain vs synapse-to-synapse spread", fontsize=10); ax[3].legend(fontsize=8)
fig.suptitle("Nevian burst pairing in the real post cells: og-delta vs antic-delta "
             f"(basal Na 0.003->0.0143, Ka 0.002->0.025; {df.pair.nunique()} pairs, spike-matched runs only)")
fig.tight_layout(); fig.savefig("glusynapse_v2/results/diag_burst_antic.png", dpi=110)
print(df.groupby(["proto", "v"])[["gain", "vshare"]].median().round(3))
