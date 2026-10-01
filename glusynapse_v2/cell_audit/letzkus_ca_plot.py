"""Tables + figure for letzkus_ca.py outputs (trunk-path sites and the real L23->L5 synapses)."""
import glob, json, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = os.path.dirname(os.path.abspath(__file__))
plt.style.use("/project/rrg-emuller/dhuruva/plastyfire/onerule.mplstyle")
VN = {"delta": "og-delta", "antic_delta": "antic-delta"}
P3, P2, P1 = "3ap_200hz", "2ap_50hz", "1ap"
data = {}
for f in sorted(glob.glob(os.path.join(D, "letzkus_ca_*-*.json"))):
    d = json.load(open(f))
    for r in d["results"]:
        data[(d["pair"], r["variant"], r["proto"])] = r
pairs = sorted({k[0] for k in data})
L5L5 = {"181015-184976", "182339-200396"}


def sites(pair, v, p, real):
    return [s for s in data[(pair, v, p)]["sites"] if s["sec"].startswith("syn:") == real]


print("soma spikes per protocol:")
for pair in pairs:
    print(pair, {p: len(data[(pair, "delta", p)]["soma_spikes"]) for p in (P1, P2, P3, "3ap_100hz", "3ap_50hz")})

# trunk path: median across cells at each target distance
TG = [50, 150, 250, 350, 450, 550, 650, 750, 900]
for v in VN:
    print(f"\n### trunk path, {VN[v]} (median over {len(pairs)} cells)")
    print("| dist | bAP 1AP mV | 3AP200 spikes 1/2/3 mV | 3AP/1AP | 3AP200/2AP50 | 3AP200/3AP50 | V-int 3AP/1AP |")
    for i, tg in enumerate(TG):
        row = []
        for pair in pairs:
            a = {p: sites(pair, v, p, False) for p in (P1, P2, P3, "3ap_50hz")}
            if i >= len(a[P1]): continue
            row.append([a[P1][i]["per_spike"][0]["vpk"], *[e["vpk"] for e in a[P3][i]["per_spike"]],
                        a[P3][i]["vdcc_int"] / a[P1][i]["vdcc_int"], a[P3][i]["vdcc_int"] / a[P2][i]["vdcc_int"],
                        a[P3][i]["vdcc_int"] / a["3ap_50hz"][i]["vdcc_int"], a[P3][i]["v_int"] / a[P1][i]["v_int"]])
        m = np.median(np.array(row), 0)
        print(f"| {tg} | {m[0]:.0f} | {m[1]:.0f}/{m[2]:.0f}/{m[3]:.0f} | {m[4]:.2f} | {m[5]:.2f} | {m[6]:.2f} | {m[7]:.2f} |")

# real synapses of the L23->L5 pairs, binned by path distance
BINS = [(0, 200, "<200"), (200, 450, "200-450"), (450, 2000, ">450")]
for v in VN:
    print(f"\n### real L23->L5 synapses, {VN[v]}")
    print("| bin | n syn | bAP 1AP mV | VDCC 1AP | 3AP/1AP | 3AP200/2AP50 |")
    for lo, hi, lab in BINS:
        rr = []
        for pair in pairs:
            if pair in L5L5: continue
            a = {p: sites(pair, v, p, True) for p in (P1, P2, P3)}
            for j, s in enumerate(a[P1]):
                if lo <= s["dist"] < hi:
                    rr.append([s["per_spike"][0]["vpk"], s["vdcc_int"], a[P3][j]["vdcc_int"] / s["vdcc_int"],
                               a[P3][j]["vdcc_int"] / a[P2][j]["vdcc_int"]])
        if rr:
            m = np.median(np.array(rr), 0)
            print(f"| {lab} | {len(rr)} | {m[0]:.0f} | {m[1]:.2e} | {m[2]:.2f} | {m[3]:.2f} |")

# figure: pooled over cells, real synapses (dots) + trunk path (lines, median)
fig, ax = plt.subplots(1, 3, figsize=(12, 3.6))
for v, ls in (("delta", "-"), ("antic_delta", "--")):
    for p, c in ((P1, "C0"), (P2, "C1"), (P3, "C3")):
        tr = np.array([[s["vdcc_int"] for s in sites(pair, v, p, False)] for pair in pairs])
        ax[0].semilogy(TG, np.median(tr, 0), ls, color=c, label=f"{p} {VN[v]}")
        if v == "delta":
            xs = [s["dist"] for pair in pairs if pair not in L5L5 for s in sites(pair, v, p, True)]
            ys = [s["vdcc_int"] for pair in pairs if pair not in L5L5 for s in sites(pair, v, p, True)]
            ax[0].semilogy(xs, ys, ".", color=c, alpha=0.4, ms=3)
    t1 = np.array([[s["vdcc_int"] for s in sites(pair, v, P1, False)] for pair in pairs])
    t2 = np.array([[s["vdcc_int"] for s in sites(pair, v, P2, False)] for pair in pairs])
    t3 = np.array([[s["vdcc_int"] for s in sites(pair, v, P3, False)] for pair in pairs])
    ax[2].plot(TG, np.median(t3 / t1, 0), ls, color="C3", marker="o", ms=3, label=f"3AP200/1AP {VN[v]}")
    ax[2].plot(TG, np.median(t3 / t2, 0), ls, color="C2", marker="o", ms=3, label=f"3AP200/2AP50 {VN[v]}")
for k in range(3):
    y = np.array([[s["per_spike"][k]["vpk"] if k < len(s["per_spike"]) else np.nan for s in sites(pair, "delta", P3, False)]
                  for pair in pairs])
    ax[1].plot(TG, np.median(y, 0), color=f"C{k+4}", marker="o", ms=3, label=f"3AP200 spike {k+1}")
y = np.array([[s["per_spike"][0]["vpk"] for s in sites(pair, "delta", P1, False)] for pair in pairs])
ax[1].plot(TG, np.median(y, 0), "k-", marker="o", ms=3, label="1 AP")
ax[1].axvspan(600, 700, color="0.9", zorder=0)
ax[1].text(470, 60, "Letzkus Fig 6A: burst\nCa spikes at 660 um", fontsize=7)
ax[0].set(xlabel="path distance (um)", ylabel="GluSynapse ica_VDCC integral", title="spine VDCC charge (lines: trunk, dots: L23 syns)")
ax[1].set(xlabel="path distance (um)", ylabel="local peak - rest (mV)", title="bAP amplitude, og-delta")
ax[2].set(xlabel="path distance (um)", ylabel="ratio", title="burst / single VDCC charge (trunk)")
ax[2].axhline(1, color="0.5", lw=0.5)
for a in ax:
    a.axvline(450, color="0.6", lw=0.5, ls="--"); a.legend(fontsize=6)
fig.suptitle(f"L5TTPC post cells (n={len(pairs)}), somatic pulses only, median over cells", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(D, "figs", "letzkus_ca.png"), dpi=150)
