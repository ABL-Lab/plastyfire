"""(b) protocol encoding on L2/3 letzkus_1ap_dt+10ms: nreps truncation and a 0.1 Hz emulation (9 s of extra
rho relaxation per pairing, inserted 600 ms after each pairing's post AP). dpre is 0 for this protocol."""
import json, os, sys, ast
import numpy as np, pandas as pd
G = "/project/rrg-emuller/dhuruva/plastyfire/glusynapse_v2"; ROOT = os.path.dirname(G)
sys.path.insert(0, G)
import batch_v2 as BV
FIT = "reduced_gpu_subset_pl5r50_td2_s2"
fit = json.load(open(f"{G}/results/{FIT}.json")); a = fit["a"]
g = pd.read_csv(f"{ROOT}/ebner/pair_geometry_L23PC_L5TTPC.csv"); g["pair"] = g.pregid.astype(str) + "-" + g.postgid.astype(str)
L23 = sorted(g[g.all_protocols].pair)
BV.BASIS_DIR = os.path.join(ROOT, "basis_results_edges_ebner_l23l5_delta")
rows = []; bad = 0
for c0 in range(0, len(L23), 20):
    B = BV.BatchV2([f"{G}/extracted/ebner_l23l5_delta-prefire"], protocols=["letzkus_1ap_dt+10ms"],
                   pairs=set(L23[c0:c0 + 20]), fast=False, signals=(), verbose=False)
    th = [B.thetas(r, a) for r in B.recs]
    td = np.concatenate([t[0] for t in th]); tp = np.concatenate([t[1] for t in th])
    E, H = B._stack(); H0 = H.copy(); L0 = B._len.copy()
    res = {}
    for nrep in (100, 50, 15):
        for slow in (False, True):
            H[:] = H0; L = L0.copy()
            for j, r in enumerate(B.recs):
                ps = np.sort(r["postspikes"])[:nrep]
                t = r["t"]; kcut = np.searchsorted(t, ps[-1] + 600.0)
                L[B._rec == j] = min(kcut, len(t) - 1)
                if slow:
                    ks = np.searchsorted(t, ps + 600.0)
                    ks = ks[ks < len(t) - 1]
                    H[j, ks] += 9.0 / BV.TAU_IND_GB     # 9 s more between pairings (0.1 Hz)
                    sl = r["sl"]; e = r["effcai"][:, ks]
                    bad += int((e > np.minimum(td[sl], tp[sl])[:, None]).sum())
            rho = BV._rho_kernel(E, H, B._rec, L, td, tp, B._rho0, float(BV.GAMMA_D_GB), float(BV.GAMMA_P_GB),
                                 float(BV.RHO_STAR_GB))
            for j, r in enumerate(B.recs):
                b = B.basis(r)
                rows.append(dict(pair=r["pair"], nrep=nrep, slow=slow,
                                 R=b.ratio(r["rho0"], rho[r["sl"]], np.zeros(len(r["syn"])))))
    del B
df = pd.DataFrame(rows)
print("steps above a threshold at inserted gaps:", bad)
print(df.groupby(["nrep", "slow"]).R.agg(["mean", "sem", "count"]).round(3))
