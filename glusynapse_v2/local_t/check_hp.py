"""D1 variants with a rest-free high-pass V detector (v - vbar, vbar' = (v - vbar)/tau_b). Login node."""
import sys
import numpy as np
import analyze as A
variant = sys.argv[1] if len(sys.argv) > 1 else "og"
recs = A.load(variant)
for r in recs:                     # x - rest before the filter = the mod's INITIAL vbar = v
    for p in r["sig"]:
        x = r["sig"][p]["v"] - r["sig"]["ap1"]["v"][:, :1]
        for tb in (10.0, 30.0):
            r["sig"][p][f"vhp{tb:.0f}"] = x - A.lp(x, tb) / tb
        r["sig"][p]["vrel"] = x
rows = []
for key in ("vhp10", "vhp30", "vrel"):
    for th in (3, 4, 6):
        for w in (0.0, 15.0):
            for te, tt in ((0.85, 40.0), (0.6, 40.0), (0.9, 20.0), (1.2, 40.0)):
                p = dict(th=th, te=te, tauT=tt, w=w)
                L, N, E, TR, D = A.evaluate(recs, key, A.c_event, p)
                A.MARGIN = 1.0
                fo, fw, th_g = A.score(L, N, E, strict=True)
                fl = A.score(L, N, E, strict=False)[1]
                far = D >= 150
                ok = (L.min(1) >= 2 * th_g) & (np.maximum(N.max(1), E) <= th_g)
                gate = {k: round(float(np.median(np.tanh(np.maximum(v - th_g, 0) / th_g).mean(1))), 2) for k, v in TR.items()}
                rows.append((fw, key, th, w, te, tt, fl, ok[far].mean(), gate,
                             np.round(np.median(L, 0), 2).tolist(), np.round(np.median(N, 0), 3).tolist()))
for r in sorted(rows, key=lambda r: -r[0])[:14]:
    print(f"f_win {r[0]:.3f} (loose {r[6]:.3f}, >150um {r[7]:.2f})  {r[1]} th {r[2]} w {r[3]} te {r[4]} tauT {r[5]}\n"
          f"    medT LTD {r[9]} noLTD {r[10]}  trains {r[8]}")
