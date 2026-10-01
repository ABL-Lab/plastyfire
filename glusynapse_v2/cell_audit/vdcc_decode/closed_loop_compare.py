"""Compare closed-loop ljp_VDCC = 25 sims (out/*_s25.npz, job 22088814) with the open-loop replay on og V.  Login node."""
import json, os
import numpy as np
import replay as RP
A = RP.A
res = {}
gk = np.load(os.path.join(RP.HERE, "gK.npz"))
recs = A.load("og")
off = np.cumsum([0] + [len(r["dist"]) for r in recs])
for i, r in enumerate(recs):
    f = os.path.join(RP.HERE, "out", f"{r['pair']}_s25.npz")
    if not os.path.exists(f):
        continue
    z = np.load(f); g, K = gk["g"][off[i]:off[i + 1]], gk["K"][off[i]:off[i + 1]]
    for p in ("ap1", "burst", "epsp", "sj20@-10"):
        v_og = r["sig"][p]["v"]; v25 = z[f"s25|{p}|v"].astype(float)
        I, Cc = RP.spine(v_og, r["sig"][p]["ica_nmda"], g, K, ljp=25.0)
        i25 = z[f"s25|{p}|ica_vdcc"].astype(float); c25 = z[f"s25|{p}|cacr"].astype(float)
        res[f"{r['pair']}|{p}"] = dict(max_abs_dV_mV=float(np.abs(v25 - v_og).max()),
                                     vdcc_charge_ratio_open_over_closed=np.round(np.percentile(I.sum(1) / i25.sum(1), [0, 50, 100]), 3).tolist(),
                                     cacr_peak_ratio=np.round(np.percentile((Cc.max(1) - 7e-5) / (c25.max(1) - 7e-5), [0, 50, 100]), 3).tolist())
print(json.dumps(res, indent=1))
json.dump(res, open(os.path.join(RP.HERE, "closed_loop_compare.json"), "w"), indent=1)
