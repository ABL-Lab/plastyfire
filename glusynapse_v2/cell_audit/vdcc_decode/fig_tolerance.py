"""Figure: t_drive 2 rows vs a per-synapse scale lambda of the per-bAP drive (te 0.85, tauT 40) against the spread of
the per-bAP drive across the 191 synapses.  Login node -> figs/tolerance.png"""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import replay as RP
A = RP.A
recs = A.load("og"); t = recs[0]["t"]
lams = np.exp(np.linspace(np.log(0.5), np.log(2), 301))
L, N = [], []
for lam in lams:
    T = {}
    for p in ("ap1", "burst"):
        imp = np.zeros(len(t[p]))
        for sp in recs[0]["spk"][p]:
            imp[np.searchsorted(t[p], sp)] = lam
        T[p] = A.lp(np.maximum(A.imp_decay(imp[None], A.TAU_E1) - 0.85, 0), 40.0)[0]
    at = lambda p, x: T[p][np.searchsorted(t[p], x + 0.1)]
    L.append([at(p, x) for p, x in A.LTD]); N.append([at(p, x) for p, x in A.NOLTD])
L, N = np.array(L), np.array(N)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
for j, (p, x) in enumerate(A.LTD):
    ax[0].semilogy(lams, L[:, j], label=f"LTD {p}@{x}")
for j, (p, x) in enumerate(A.NOLTD):
    ax[0].semilogy(lams, N[:, j], "--", label=f"no-LTD {p}@{x}")
th = (0.122 + 0.2125) / 2
ax[0].axhline(th, c="k", lw=0.8); ax[0].axhline(2 * th, c="k", lw=0.8, ls=":")
ax[0].axvspan(0.965, 1.04, color="g", alpha=0.2, label="pass band, theta_Tg fixed")
ax[0].set_xscale("log"); ax[0].set_xlabel("lambda (per-bAP drive / fitted unit)"); ax[0].set_ylabel("T at pre arrival")
ax[0].set_ylim(1e-3, 100); ax[0].legend(fontsize=6)
gk = np.load(os.path.join(RP.HERE, "gK.npz")); g, K = gk["g"], gk["K"]
v = RP.stack(recs, "ap1", "v")[0]
for lab, key in (("bAP amplitude (V)", None), ("spine VDCC Ca, R as used", 0.0), ("spine VDCC Ca, R -25 mV", 25.0)):
    x = (v - v[:, :1]).max(1) if key is None else RP.spine(v, np.zeros_like(v), g, K, ljp=key)[1].max(1) - 7e-5
    x = np.sort(x / np.median(x))
    ax[1].semilogx(x, np.linspace(0, 1, len(x)), label=lab)
ax[1].axvspan(0.965, 1.04, color="g", alpha=0.2); ax[1].set_xlabel("per-bAP drive / median (191 synapses)"); ax[1].set_ylabel("cumulative fraction")
ax[1].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(RP.HERE, "figs", "tolerance.png"), dpi=110)
