"""Per-synapse decoding of the bAP Ca signal (local_t/out, og and s5).  Login node, numpy only.
    python cell_audit/vdcc_decode/decode.py  -> figs/decode.png, decode.json, gK.npz (per-synapse g, K)"""
import json, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import replay as RP
A, C = RP.A, RP.C
HERE = RP.HERE
PROTOS = ["ap1", "burst", "epsp"] + A.TRAINS


def win(t, x, a, b):
    m = (t >= a) & (t < b)
    return x[:, m]


def metrics(recs, p, a, b, extra=None):
    t = recs[0]["t"][p]
    base = lambda x: x[:, t < -2].mean(1, keepdims=True)
    out = {}
    for key in ("v", "nvdcc", "cacr", "cai", "ica_nmda"):
        x, _ = RP.stack(recs, p, key)
        d = win(t, x - base(x), a, b)
        if key == "ica_nmda":
            d = -d
        out[key + "_pk"] = d.max(1); out[key + "_int"] = d.clip(0).sum(1) * A.DT
    out["vrest"] = base(RP.stack(recs, p, "v")[0])[:, 0]
    out["vabs_pk"] = out["vrest"] + out["v_pk"]
    return out


def main():
    res = {}
    recs = A.load("og"); recs5 = A.load("s5")
    dist = np.concatenate([r["dist"] for r in recs]); kind = np.concatenate([r["kind"] for r in recs])
    # --- 1. spine replay validation: fit g, K on og, predict og and s5 (ljp 5) ---
    V = {p: RP.stack(recs, p, "v")[0] for p in PROTOS}
    get = lambda rr, p, k: RP.stack(rr, p, k)[0]
    g, K = RP.fit_gK([V[p] for p in ("ap1", "epsp")], [get(recs, p, "ica_vdcc") for p in ("ap1", "epsp")],
                     [get(recs, p, "ica_nmda") for p in ("ap1", "epsp")], [get(recs, p, "cacr") for p in ("ap1", "epsp")])
    np.savez(os.path.join(HERE, "gK.npz"), g=g, K=K, dist=dist, kind=kind)
    val = {}
    for name, rr, ljp in (("og", recs, 0.0), ("s5", recs5, 5.0)):
        for p in ("ap1", "burst", "sj20@-10"):
            vv = get(rr, p, "v")
            I, Cc = RP.spine(vv, get(rr, p, "ica_nmda"), g, K, ljp=ljp)
            ir, cr = get(rr, p, "ica_vdcc"), get(rr, p, "cacr")
            qi = (I.sum(1) + 1e-30) / (ir.sum(1) + 1e-30); qc = (Cc - 70e-6).max(1) / np.maximum((cr - 70e-6).max(1), 1e-12)
            val[f"{name}|{p}"] = dict(vdcc_charge_ratio_median=float(np.median(qi)), vdcc_ratio_IQR=np.percentile(qi, [5, 95]).round(3).tolist(),
                                      cacr_peak_ratio_median=float(np.median(qc)), cacr_ratio_5_95=np.percentile(qc, [5, 95]).round(3).tolist())
    res["replay_validation"] = val
    # volume consistency: K ~ 1/vol, g ~ vol^(2/3)  -> g * K^(2/3) should be ~constant
    gk = g * K ** (2 / 3)
    res["g_K23_cv"] = float(np.std(gk) / np.mean(gk))
    # --- 2. per-synapse bAP and EPSP signals ---
    bap = metrics(recs, "ap1", -2, 60); ep = metrics(recs, "epsp", 0, 150)
    bap5 = metrics(recs5, "ap1", -2, 60)
    bins = [(0, 60), (60, 150), (150, 250), (250, 1e4)]
    tab = []
    for lo, hi in bins:
        m = (dist >= lo) & (dist < hi)
        row = dict(bin=f"{lo}-{hi if hi < 1e4 else ''}", n=int(m.sum()), n_basal=int((m & (kind == "basal")).sum()))
        for k in ("v_pk", "vabs_pk", "nvdcc_int", "cacr_pk", "cacr_int", "cai_pk", "cai_int"):
            row["bap_" + k] = float(np.median(bap[k][m]))
            row["epsp_" + k] = float(np.median(ep[k][m]))
            if k != "vabs_pk":
                row["ratio_bap_epsp_" + k] = float(np.median(bap[k][m] / np.maximum(ep[k][m], 1e-30)))
        row["s5_bap_nvdcc_int"] = float(np.median(bap5["nvdcc_int"][m]))
        tab.append(row)
    res["bins"] = tab
    # log-log slopes vs bAP amplitude and e-fold per mV of peak V
    sl = {}
    ok = bap["v_pk"] > 1
    for k in ("nvdcc_int", "cacr_pk", "cai_pk"):
        y = np.log(np.maximum(bap[k], 1e-30))
        m = ok & (bap[k] > 0)
        sl[k] = dict(n_per_log_amp=float(np.polyfit(np.log(bap["v_pk"][m]), y[m], 1)[0]),
                     mV_per_efold=float(1 / np.polyfit(bap["vabs_pk"][m], y[m], 1)[0]),
                     range_p5_p95_decades=float(np.log10(np.percentile(bap[k][m], 95) / np.percentile(bap[k][m], 5))))
    sl["v_pk"] = dict(range_p5_p95_decades=float(np.log10(np.percentile(bap["v_pk"], 95) / np.percentile(bap["v_pk"], 5))))
    # EPSP-dominated synapses: bAP signal below the own EPSP signal
    for k in ("nvdcc_pk", "cacr_pk", "cai_pk", "v_pk"):
        sl["frac_epsp_gt_bap_" + k] = float(np.mean(ep[k] > bap[k]))
    # what feeds cacr: fraction of spine Ca charge from VDCC in ap1 and epsp
    res["slopes"] = sl
    # shaft cai: locally driven? replay HVA2 + LVA on v, compare integral with cai rise
    ih, il = RP.dend(V["ap1"])
    t = recs[0]["t"]["ap1"]; m = (t > -2) & (t < 60)
    qh, ql = -(ih[:, m]).clip(max=0).sum(1), -(il[:, m]).clip(max=0).sum(1)
    q = qh + ql
    cc = np.corrcoef(np.log(q + 1e-30), np.log(np.maximum(bap["cai_pk"], 1e-30)))[0, 1]
    res["shaft"] = dict(lva_share_median=float(np.median(ql / np.maximum(q, 1e-30))),
                        corr_log_replay_charge_vs_cai_pk=float(cc),
                        cai_pk_over_cacr_pk_median=float(np.median(bap["cai_pk"] / np.maximum(bap["cacr_pk"], 1e-30))))
    for lo, hi in bins:
        mm = (dist >= lo) & (dist < hi)
        res["shaft"][f"lva_share_{lo}"] = float(np.median((ql / np.maximum(q, 1e-30))[mm]))
    json.dump(res, open(os.path.join(HERE, "decode.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))
    # --- figure ---
    fig, ax = plt.subplots(2, 3, figsize=(14, 8))
    cb = dict(c=dist, cmap="viridis", s=12)
    a = ax[0, 0]; sc = a.scatter(dist, bap["v_pk"], **cb); a.set_yscale("log"); a.set_xlabel("path dist (um)"); a.set_ylabel("bAP amplitude (mV)")
    a.scatter(dist, ep["v_pk"], c="r", s=6, marker="x", label="own EPSP"); a.legend(fontsize=7)
    for a, k, lab in ((ax[0, 1], "nvdcc_int", "spine VDCC charge (nA ms)"), (ax[0, 2], "cacr_pk", "spine d[Ca] peak (mM)"),
                      (ax[1, 0], "cai_pk", "shaft d[Ca] peak (mM)")):
        a.scatter(bap["vabs_pk"], np.maximum(bap[k], 1e-12), **cb); a.scatter(ep["vabs_pk"], np.maximum(ep[k], 1e-12), c="r", s=6, marker="x", label="own EPSP")
        a.set_yscale("log"); a.set_xlabel("peak local V (mV)"); a.set_ylabel(lab + " (bAP dots, EPSP x)"); a.legend(fontsize=7)
    a = ax[1, 1]
    a.scatter(dist, bap["cacr_pk"] / np.maximum(ep["cacr_pk"], 1e-12), **cb); a.set_yscale("log"); a.axhline(1, c="k", lw=0.5)
    a.set_xlabel("path dist (um)"); a.set_ylabel("spine Ca peak bAP / own EPSP")
    a = ax[1, 2]
    a.scatter(q + 1e-12, np.maximum(bap["cai_pk"], 1e-12), **cb); a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel("replayed HVA2+LVA charge density at the segment"); a.set_ylabel("recorded shaft d[Ca] peak")
    fig.colorbar(sc, ax=ax[0, 0], label="dist")
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "figs", "decode.png"), dpi=110)


if __name__ == "__main__":
    main()
