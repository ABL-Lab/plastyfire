"""Follow-ups (CA_DRIVE_DESIGN.md): Hill onset on influx with a wider K/n grid, the same with K scaled by the synapse's own
bAP influx peak (candidates 1+2), failure classes by distance bin, and own-history (NRM) robustness.  Login node.
    CA_CACHE=<scratch> python cell_audit/ca_drive/extra.py -> out/extra.json"""
import itertools, json, os
import numpy as np
import common as C
import cands as K
A = C.A


def hon(x, Kv, n):
    h = K.hill(x, Kv, n); d = np.zeros_like(h); d[:, 1:] = np.maximum(np.diff(h, axis=1), 0); return d


def classify(S, t, dist, pid, incf, te, tt, veto):
    inc = {p: incf(p) * C.veto_weight(veto, t[p], C.arrivals(p))[None, :] for p in C.PROTOS}
    L, N, E, TR = C.rows({p: C.rule(inc[p], te, tt) for p in C.PROTOS}, t)
    th = A.score(L, N, E)[2]
    ok = (L.min(1) >= 2 * th) & (np.maximum(N.max(1), E) <= th)
    silent = ~ok & (L.max(1) < 2 * th) & (np.maximum(N.max(1), E) <= th)      # no LTD anywhere: benign loss
    early = ~ok & ~silent & (L[:, 0] < 2 * th) & (L[:, 1:].min(1) >= 2 * th) & (np.maximum(N.max(1), E) <= th)  # only -10 missing
    wrong = ~ok & ~silent & ~early
    return {f"{a}-{b:.0f}": dict(n=int(m.sum()), ok=round(float(ok[m].mean()), 2), silent=round(float(silent[m].mean()), 2),
                                 only_minus10=round(float(early[m].mean()), 2), wrong=round(float(wrong[m].mean()), 2))
            for a, b in C.BINS for m in [(dist >= a) & (dist < b)]}


def main():
    S, t, dist, pid = C.build()
    out = {}
    for lj in (0, 25):
        s = S[lj]; qmed = np.median(s["ap1"]["q"].max(1)); qpk = s["ap1"]["q"].max(1)
        for lab, Kf in (("HOq global K", lambda k: k * qmed), ("HOq own K (K_i = k x own bAP peak)", lambda k: k * qpk[:, None])):
            best = None
            for k, n in itertools.product(10.0 ** np.arange(-5.0, -0.49, 0.5), (4.0, 8.0, 16.0)):
                r = K.evaluate(s, t, dist, pid, lambda p, Kv=Kf(k), n=n: hon(s[p]["q"], Kv, n), veto="nmda")
                if best is None or (r["f_syn"], r["conn_meanT"]) > (best["f_syn"], best["conn_meanT"]):
                    best = dict(r, k=float(k), n=n)
            best["classes"] = classify(s, t, dist, pid, lambda p, Kv=Kf(best["k"]), n=best["n"]: hon(s[p]["q"], Kv, n),
                                       best["te"], best["tauT"], "nmda")
            out[f"{lab}|ljp{lj}"] = best
            print(lab, lj, {kk: best[kk] for kk in ("f_syn", "bins", "conn_major", "conn_meanT", "k", "n", "te", "tauT", "trains")},
                  "\n   ", best["classes"], flush=True)
        # EV reference failure classes
        th = 0.003 * qmed
        out[f"EV classes|ljp{lj}"] = classify(s, t, dist, pid, lambda p: A.crossings(s[p]["q"], th), 0.85, 40.0, "nmda")
        print("EV classes", lj, out[f"EV classes|ljp{lj}"])
        # NRM: E row (own EPSP read 50 ms later) under pre-only history, by floor, for each veto
        Ab = K.charge(s["ap1"]["q"])
        for veto in ("none", "nmda", "hand"):
            wE = C.veto_weight(veto, t["epsp"], C.arrivals("epsp"))[None, :]
            Ae = K.charge(s["epsp"]["q"] * wE)
            for q0 in (0, 10, 25, 50):
                M0 = np.percentile(Ab, q0) if q0 else 0.0
                lamE = Ae / np.maximum(Ae, M0 if q0 else 1e-300)        # EPSP drive per event under pre-only history
                lamB = Ab / np.maximum(Ab, M0) if q0 else np.ones_like(Ab)  # bAP drive after a bAP history
                out[f"NRM pre-only|ljp{lj}|{veto}|floor p{q0}"] = dict(
                    frac_EPSP_unit_drive_ge_0p5=round(float(np.mean(lamE >= 0.5)), 2),
                    frac_bAP_drive_below_0p96=round(float(np.mean(lamB < 0.96)), 2))
                print("NRM", lj, veto, q0, out[f"NRM pre-only|ljp{lj}|{veto}|floor p{q0}"])
    json.dump(out, open(os.path.join(C.HERE, "out", "extra.json"), "w"), indent=0)


if __name__ == "__main__":
    main()
