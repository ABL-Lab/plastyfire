"""LTP ceiling per protocol: EPSP ratio if every synapse ends potentiated (rho=1), with dpre 0 or dmax."""
import sys, numpy as np
sys.path.insert(0, "glusynapse_v2")
from batch_v2 import BatchV2
B = BatchV2(["glusynapse_v2/extracted/ebner_delta-prefire"], protocols=["nevian_3ap_50hz_dt+10ms"],
            keep_traces=False, fast=False, signals=())
out = {"rho1": [], "rho1_dmax": [], "rho1_d0.3": []}
for r in B.recs:
    b = B.basis(r); n = len(r["syn"]); one = np.ones(n)
    out["rho1"].append(b.ratio(r["rho0"], one, np.zeros(n)))
    out["rho1_dmax"].append(b.ratio(r["rho0"], one, np.ones(n)))
    out["rho1_d0.3"].append(b.ratio(r["rho0"], one, np.full(n, 0.3)))
    out.setdefault("rho0_frac_pot", []).append(np.mean(r["rho0"] >= 0.5))
for k, v in out.items():
    v = np.array(v); print(f"{k:14s} mean {np.nanmean(v):.3f}  median {np.nanmedian(v):.3f}  p10 {np.nanpercentile(v,10):.3f} p90 {np.nanpercentile(v,90):.3f}  n={len(v)}")
