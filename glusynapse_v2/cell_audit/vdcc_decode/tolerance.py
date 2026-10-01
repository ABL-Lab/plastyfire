"""(1) How much per-synapse spread of the per-bAP drive does the t_drive 2 window tolerate?  S += lambda per post AP,
T = lp(pos(S - theta_Te), tau_T); find the widest lambda range that one theta_Tg passes (LTD rows >= 2 theta_Tg,
no-LTD rows <= theta_Tg).  (2) Integrated Ca drives fitted inside distance bins (bin-normalised) and with an oracle
per-synapse normalisation (u / own bAP area), to separate 'spread across synapses' from 'shape of the drive'.
Login node.  python cell_audit/vdcc_decode/tolerance.py -> tolerance.json"""
import itertools, json, os
import numpy as np
import replay as RP
import counterfactual as CF
A = RP.A
HERE = RP.HERE
A.MARGIN = 1.0
PROTOS = CF.PROTOS


def lam_tolerance(recs, t):
    spk = {p: recs[0]["spk"][p] for p in ("ap1", "burst")}
    lams = np.exp(np.linspace(np.log(0.25), np.log(4), 801))
    out = {}
    for te, tt in itertools.product((0.3, 0.6, 0.85, 1.0, 1.5), (20.0, 40.0, 80.0)):
        LN = []
        for lam in lams:
            T = {}
            for p in spk:
                imp = np.zeros(len(t[p]))
                for sp in spk[p]:
                    imp[np.searchsorted(t[p], sp)] = lam
                T[p] = A.lp(np.maximum(A.imp_decay(imp[None], A.TAU_E1) - te, 0), tt)[0]
            at = lambda p, x: T[p][np.searchsorted(t[p], x + 0.1)]
            LN.append((min(at(p, x) for p, x in A.LTD), max(at(p, x) for p, x in A.NOLTD)))
        LN = np.array(LN)
        best = (1.0, None)
        for th in np.unique(np.concatenate([LN[:, 1], LN[:, 0] / 2])):
            if th <= 0:
                continue
            ok = (LN[:, 0] >= 2 * th) & (LN[:, 1] <= th)
            if not ok.any():
                continue
            # longest contiguous run
            idx = np.where(ok)[0]; runs = np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)
            r = max(runs, key=len)
            f = lams[r[-1]] / lams[r[0]]
            if f > best[0]:
                best = (float(f), (float(lams[r[0]]), float(lams[r[-1]])))
        out[f"te{te}_tau{tt:.0f}"] = dict(fold=best[0], lam=best[1])
    return out


def fit(S, key, t, mask, oracle=False):
    x0 = S["ap1"][key][mask]; rest = x0[:, :1]
    pk = (x0 - rest).max(1); area = np.maximum(x0 - rest, 0).sum(1) * A.DT
    a = area[:, None] if oracle else np.median(area)
    best = (0.0, None)
    for q, te, tt, w in itertools.product((0.0, 0.03, 0.1, 0.3), (0.3, 0.5, 0.7, 0.85, 1.0, 1.3), (20.0, 40.0, 80.0), (0.0, 15.0)):
        th = (pk[:, None] if oracle else np.median(pk)) * q
        T = {}
        for p in PROTOS:
            if p not in ("ap1", "burst", "epsp"):
                continue
            x = S[p][key][mask]
            u = np.maximum(x - rest - th, 0) / a
            if w:
                u = u * A.veto_mask(t[p], A.pre_times(p) + 0.1, w)[None, :]
            T[p] = A.lp(np.maximum(A.lp(u, A.TAU_E1) - te, 0), tt)
        at = lambda p, xx: T[p][:, np.searchsorted(t[p], xx + 0.1)]
        L = np.stack([at(p, xx) for p, xx in A.LTD], 1); N = np.stack([at(p, xx) for p, xx in A.NOLTD], 1)
        f = A.score(L, N, at("epsp", 50.0), strict=True)[1]
        if f > best[0]:
            best = (float(f), dict(q=q, te=te, tauT=tt, w=w))
    return best


def main():
    recs = A.load("og"); t = recs[0]["t"]
    dist = np.concatenate([r["dist"] for r in recs])
    res = dict(lambda_tolerance=lam_tolerance(recs, t))
    print(json.dumps(res["lambda_tolerance"]))
    gk = np.load(os.path.join(HERE, "gK.npz")); g, K = gk["g"], gk["K"]
    S = {}
    for p in ("ap1", "burst", "epsp"):
        S[p] = {k: RP.stack(recs, p, k)[0] for k in ("v", "cai", "cacr")}
        z = np.zeros_like(S[p]["v"])
        for ljp in (0, 25):
            S[p][f"cvd{ljp}"] = RP.spine(S[p]["v"], z, g, K, ljp=float(ljp))[1]
        S[p]["vrel"] = S[p]["v"]
    bins = {"<60": dist < 60, "60-150": (dist >= 60) & (dist < 150), ">150": dist >= 150, "all": dist >= 0}
    res["fits"] = {}
    for key in ("cvd0", "cvd25", "cai", "cacr", "vrel"):
        for b, m in bins.items():
            for orc in (False, True):
                f = fit(S, key, t, m, orc)
                res["fits"][f"{key}|{b}|{'oracle' if orc else 'bin-norm'}"] = f
                print(key, b, "oracle" if orc else "bin-norm", f, flush=True)
    json.dump(res, open(os.path.join(HERE, "tolerance.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
