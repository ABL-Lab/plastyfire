"""Low-threshold (T-type) spine VDCC counterfactual + window tolerance to per-synapse drive spread.  Login node.
Spine channel variants (replayed on recorded local V, spine Ca from VDCC only, NMDAR share excluded):
  R0   : GluSynapse R-type as used            R25 : R-type, 110 mM Ba -> 2 mM Ca surface-charge shift (-25 mV)
  Tu   : Ca_LVAst kinetics as used (-10 mV)   Ta  : Ca_LVAst at Avery & Johnston 1996 values (no -10 mV)
  R0+Tu: half of the proximal bAP Ca each
Every variant's density is set so the median bAP spine Ca at < 60 um equals R0's (the Chindemi/Sabatini-calibrated
proximal value), i.e. calibrated to Ca data, not to plasticity.
    python cell_audit/vdcc_decode/ttype.py -> ttype.json, figs/ttype.png"""
import itertools, json, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import replay as RP
import counterfactual as CF
A, C = RP.A, RP.C
HERE = RP.HERE
PROTOS = CF.PROTOS
A.MARGIN = 1.0
BINS = ((0, 60), (60, 150), (150, 250), (250, 1e9))


def spine_multi(v, g, K, chans, cmin=70e-6, tau=12.0):
    """VDCC-only spine Ca for a sum of m^2 h channels [(fn, gscale)], Eca from the spine Ca."""
    st = []
    for fn, _ in chans:
        mi, mt, hi, ht = fn(v[:, 0]); st.append([mi.copy(), hi.copy()])
    c = np.full(v.shape[0], cmin)
    Cc, I = np.empty_like(v), np.empty_like(v)
    dt = A.DT / RP.NSUB
    for k in range(v.shape[1]):
        for s in range(RP.NSUB if k else 1):
            vv = v[:, k - 1] + (v[:, k] - v[:, k - 1]) * (s + 1) / RP.NSUB if k else v[:, 0]
            i = 0.0
            for (fn, gs), mh in zip(chans, st):
                mi, mt, hi, ht = fn(vv)
                if k:
                    mh[0] += (1 - np.exp(-dt / mt)) * (mi - mh[0]); mh[1] += (1 - np.exp(-dt / ht)) * (hi - mh[1])
                i = i + gs * g * mh[0] ** 2 * mh[1] * (vv - RP.nernst(c))
            if k:
                c += dt * (-K * i - (c - cmin) / tau)
        Cc[:, k], I[:, k] = c, i
    return I, Cc


def main():
    recs = A.load("og"); t = recs[0]["t"]
    dist = np.concatenate([r["dist"] for r in recs])
    pid = np.concatenate([np.full(len(r["dist"]), i) for i, r in enumerate(recs)])
    V = {p: RP.stack(recs, p, "v")[0] for p in PROTOS}
    gk = np.load(os.path.join(HERE, "gK.npz")); g, K = gk["g"], gk["K"]
    prox = dist < 60
    R0 = lambda x: C.rtype(x, 0.0); R25 = lambda x: C.rtype(x, 25.0)
    Tu = lambda x: C.lva(x, 10.0); Ta = lambda x: C.lva(x, 0.0)
    variants = {"R0": [(R0, 1.0)], "R25": [(R25, 1.0)], "Tu": [(Tu, 1.0)], "Ta": [(Ta, 1.0)], "R0+Tu": [(R0, 0.5), (Tu, 1.0)]}
    bap_amp = (V["ap1"] - V["ap1"][:, :1]).max(1); vpk = V["ap1"].max(1)
    out = {}
    target = None
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    for name, ch in variants.items():
        # calibrate: scale so median proximal bAP Ca equals R0's (T part only for the mix: T gets half)
        for it in range(3):
            I, Cc = spine_multi(V["ap1"], g, K, ch)
            pk = Cc.max(1) - 70e-6
            if name == "R0":
                target = np.median(pk[prox]); break
            if name == "R0+Tu":
                Ir, Cr = spine_multi(V["ap1"], g, K, [ch[0]])
                want = target - np.median((Cr.max(1) - 70e-6)[prox]); _, Ct = spine_multi(V["ap1"], g, K, [ch[1]])
                ch = [ch[0], (ch[1][0], ch[1][1] * want / np.median((Ct.max(1) - 70e-6)[prox]))]
            else:
                ch = [(ch[0][0], ch[0][1] * target / np.median(pk[prox]))]
        I, Cc = spine_multi(V["ap1"], g, K, ch)
        pk = Cc.max(1) - 70e-6
        Ie, Ce = spine_multi(V["epsp"], g, K, ch)
        pe = Ce.max(1) - 70e-6
        row = dict(gscale=[c[1] for c in ch],
                   bap_uM_by_bin=[float(np.median(pk[(dist >= a) & (dist < b)]) * 1e3) for a, b in BINS],
                   epsp_uM_by_bin=[float(np.median(pe[(dist >= a) & (dist < b)]) * 1e3) for a, b in BINS],
                   bap_uM_amp40_60=float(np.median(pk[(bap_amp >= 40) & (bap_amp < 60)]) * 1e3),
                   bap_uM_amp20_40=float(np.median(pk[(bap_amp >= 20) & (bap_amp < 40)]) * 1e3),
                   decades_p5_p95=float(np.log10(np.percentile(pk, 95) / max(np.percentile(pk, 5), 1e-15))),
                   mV_per_efold=float(1 / np.polyfit(vpk[pk > 0], np.log(pk[pk > 0]), 1)[0]),
                   frac_epsp_gt_bap=float(np.mean(pe > pk)))
        # window scores for Ca drives built on this channel
        S = {}
        for p in PROTOS:
            I, Cc = spine_multi(V[p], g, K, ch)
            S[p] = {"c": Cc, "i": -I}
        sc = {}
        for key in ("c", "i"):
            for fn, p in CF.grids(S, key):
                if not p.get("w") or (key == "i" and p["kind"] != "ev"):
                    continue
                L, N, E = CF.rows(S, key, fn, p, t)
                k = f"{key}:{p['kind']}+v"
                fs = A.score(L, N, E, strict=True)
                if k not in sc or fs[1] > sc[k]["f_syn"]:
                    th = fs[2]
                    ok = (L.min(1) >= 2 * th) & (np.maximum(N.max(1), E) <= th) if np.isfinite(th) else np.zeros(len(L), bool)
                    sc[k] = dict(f_syn=fs[1], by_bin=[float(ok[(dist >= a) & (dist < b)].mean()) for a, b in BINS],
                                 f_conn=A.score(CF.conn(L, pid), CF.conn(N, pid), CF.conn(E, pid))[1])
                for (a, b) in ((0, 60), (0, 150)):                 # fit restricted to proximal synapses
                    m = (dist >= a) & (dist < b)
                    kk = f"{k}|fit<{b}"
                    f = A.score(L[m], N[m], E[m], strict=True)[1]
                    if kk not in sc or f > sc[kk]:
                        sc[kk] = f
        row["scores"] = sc
        out[name] = row
        print(name, json.dumps(row), flush=True)
        ax[0].scatter(vpk, np.maximum(pk * 1e3, 1e-6), s=8, label=name)
        ax[1].scatter(dist, np.maximum(pk * 1e3, 1e-6), s=8, label=name)
    # tolerance of the t_drive 2 window to a per-synapse scale lambda of the per-bAP drive (soma impulses x lambda)
    recs1 = recs[:1]
    tol = {}
    for te, tt in ((0.85, 40.0), (0.6, 40.0)):
        base = None; res = []
        for lam in np.logspace(-1, 1, 81):
            T = {}
            for p in ("ap1", "burst", "epsp"):
                imp = np.zeros(len(t[p]))
                for sp in recs1[0]["spk"][p]:
                    imp[np.searchsorted(t[p], sp)] = lam
                Sx = A.imp_decay(imp[None], A.TAU_E1); T[p] = A.lp(np.maximum(Sx - te, 0), tt)[0]
            at = lambda p, x: T[p][np.searchsorted(t[p], x + 0.1)]
            L = np.array([at(p, x) for p, x in A.LTD]); N = np.array([at(p, x) for p, x in A.NOLTD])
            if base is None and abs(lam - 1) < 1e-9:
                pass
            res.append((lam, L.min(), N.max()))
        res = np.array(res)
        i1 = np.argmin(abs(res[:, 0] - 1)); th = res[i1, 1] / 2           # theta_Tg at lambda = 1 (weakest LTD row = 2 theta)
        ok = (res[:, 1] >= 2 * th) & (res[:, 2] <= th)
        # widest theta that passes at lambda=1: max(N) <= theta <= min(L)/2 ; best tolerance: pick theta in that interval
        best = (0, None)
        for th in np.linspace(res[i1, 2], res[i1, 1] / 2, 50):
            okk = (res[:, 1] >= 2 * th) & (res[:, 2] <= th)
            lo_i = i1; hi_i = i1
            while lo_i > 0 and okk[lo_i - 1]: lo_i -= 1
            while hi_i < len(okk) - 1 and okk[hi_i + 1]: hi_i += 1
            if okk[i1] and res[hi_i, 0] / res[lo_i, 0] > best[0]:
                best = (res[hi_i, 0] / res[lo_i, 0], (float(res[lo_i, 0]), float(res[hi_i, 0])))
        tol[f"te{te}_tau{tt}"] = dict(fold_range=best[0], lam_range=best[1])
    out["tolerance_t_drive2"] = tol
    print(json.dumps(tol))
    json.dump(out, open(os.path.join(HERE, "ttype.json"), "w"), indent=1)
    for a in ax[:2]:
        a.set_yscale("log"); a.set_ylabel("bAP spine Ca from VDCC (uM), proximal-calibrated"); a.legend(fontsize=7)
    ax[0].set_xlabel("peak local V (mV)"); ax[1].set_xlabel("path dist (um)")
    names = list(variants)
    for j, lab in enumerate(("<60", "60-150", "150-250", ">250")):
        ax[2].bar(np.arange(len(names)) + 0.2 * j - 0.3, [out[n]["scores"].get("i:ev+v", {}).get("by_bin", [0] * 4)[j] for n in names], 0.2, label=lab)
    ax[2].set_xticks(range(len(names))); ax[2].set_xticklabels(names); ax[2].set_ylabel("f_win, VDCC-current event + veto"); ax[2].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "figs", "ttype.png"), dpi=110)


if __name__ == "__main__":
    main()
