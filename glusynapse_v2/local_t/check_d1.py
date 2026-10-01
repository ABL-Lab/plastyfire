"""Per-row T and robustness of the recommended local drives (LOCAL_T_DRIVE.md). Login node."""
import sys
import numpy as np
import analyze as A
variant = sys.argv[1] if len(sys.argv) > 1 else "og"
recs = A.load(variant)
sets = [("A0 ref t_drive 2", "vsoma", A.c_soma, dict(te=0.85, tauT=40.0)),
        ("A1h th4", "v", A.c_event, dict(th=4, te=0.85, tauT=40.0, rel=1, hyst=2.0)),
        ("D1 th4 w10", "v", A.c_event, dict(th=4, te=0.85, tauT=40.0, rel=1, hyst=2.0, w=10.0)),
        ("D1 th4 w15", "v", A.c_event, dict(th=4, te=0.85, tauT=40.0, rel=1, hyst=2.0, w=15.0)),
        ("D1 th6 w10", "v", A.c_event, dict(th=6, te=0.85, tauT=40.0, rel=1, hyst=2.0, w=10.0)),
        ("D1 th8 w10", "v", A.c_event, dict(th=8, te=0.85, tauT=40.0, rel=1, hyst=2.0, w=10.0)),
        ("D1 th4 w10 te.9 tau20", "v", A.c_event, dict(th=4, te=0.9, tauT=20.0, rel=1, hyst=2.0, w=10.0)),
        ("C1v VDCC event + veto", "nvdcc", A.c_event, None)]
med = np.median(np.concatenate([(-r["sig"]["ap1"]["ica_vdcc"]).max(1) for r in recs]))
sets[-1] = (sets[-1][0], "nvdcc", A.c_event, dict(th=med * 0.003, te=0.85, tauT=40.0, rel=1, w=10.0))
print("rows LTD", A.LTD, "noLTD", A.NOLTD)
for name, key, fn, p in sets:
    L, N, E, TR, D = A.evaluate(recs, key, fn, p)
    out = []
    for m in (0.25, 1.0):
        A.MARGIN = m
        out.append(A.score(L, N, E, strict=True)[1:])
    th = out[1][1]
    gate = {k: round(float(np.median(np.tanh(np.maximum(v - th, 0) / th).mean(1))), 2) for k, v in TR.items()} \
        if np.isfinite(th) else {}
    far = D >= 150
    print(f"{name:24s} f_win(m.25) {out[0][0]:.2f} f_win(m1) {out[1][0]:.2f} [>150um {np.mean((L[far].min(1) >= 2 * th) & (np.maximum(N[far].max(1), E[far]) <= th)):.2f}]"
          f"  medT LTD {np.round(np.median(L, 0), 2).tolist()} noLTD {np.round(np.median(N, 0), 3).tolist()} "
          f"epsp {np.median(E):.3f}\n{'':24s} train gate tanh((T-th)/th) {gate}")
