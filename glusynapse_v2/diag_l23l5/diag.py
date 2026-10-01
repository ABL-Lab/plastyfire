"""L2/3->L5 vs L5->L5 transfer diagnosis: per-synapse post-rule quantities + EPSP mapping, fit td2_s2 (r50).
Login node: chunks of pairs, numpy + the numba rho kernel only, no jax."""
import json, os, sys, ast
import numpy as np, pandas as pd
G = "/project/rrg-emuller/dhuruva/plastyfire/glusynapse_v2"; ROOT = os.path.dirname(G)
sys.path.insert(0, G)
import batch_v2 as BV
import model_v2 as MV

FIT = sys.argv[1] if len(sys.argv) > 1 else "reduced_gpu_subset_pl5r50_td2_s2"
fit = json.load(open(f"{G}/results/{FIT}.json")); fa = ast.literal_eval(fit["args"]) if isinstance(fit["args"], str) else fit["args"]
a = fit["a"]; P = {**MV.DEFAULTS, **json.loads(fa["filters"]), **fit["pre"]}
g = pd.read_csv(f"{ROOT}/ebner/pair_geometry_L23PC_L5TTPC.csv"); g["pair"] = g.pregid.astype(str) + "-" + g.postgid.astype(str)
L23 = sorted(g[g.all_protocols].pair)
S24 = open(f"{G}/subset24_pairs.txt").read().strip().split(",")
l5all = sorted({f.split("__")[0] for f in os.listdir(f"{G}/extracted/ebner_delta-prefire") if "sjostrom_50hz_dt+10ms" in f})
JOBS = [("L23", "ebner_l23l5_delta-prefire", "basis_results_edges_ebner_l23l5_delta", L23,
         ["sjostrom_50hz_dt+10ms", "letzkus_1ap_dt+10ms"]),
        ("L5", "ebner_delta-prefire", "basis_results_edges_sabrina_n120_delta", l5all,
         ["sjostrom_50hz_dt+10ms", "sjostrom_0.1hz_dt+10ms"])]
syn_rows, rec_rows = [], []
for path, d, bdir, pairs, protos in JOBS:
    BV.BASIS_DIR = os.path.join(ROOT, bdir)
    for c0 in range(0, len(pairs), 15):
        B = BV.BatchV2([f"{G}/extracted/{d}"], protocols=protos, pairs=set(pairs[c0:c0 + 15]), fast=False,
                       signals=("vdcc",), verbose=False)
        if not B.recs:
            continue
        th = [B.thetas(r, a) for r in B.recs]
        rho_all = B.rho_full_all(np.concatenate([t[0] for t in th]), np.concatenate([t[1] for t in th]))
        feats = B.features(P)
        for i, r in enumerate(B.recs):
            rho = rho_all[r["sl"]]; td, tp = th[i]; tT, K = feats[i]
            dp = MV.dpre_final(tT, K, P["A_mglu"], P["A_NO"], P["dpre_min"], P["dpre_max"], P["dpre0"])
            b = B.basis(r); z = np.zeros(len(rho))
            E = r["effcai"]; h = np.append(np.diff(r["t"]), 0.0) / 1000.0
            b_m, b_s = b.epsp(r["rho0"], z)
            rec_rows.append(dict(path=path, pair=r["pair"], proto=r["proto"], nsyn=len(rho), e0=b.e0, b_m=b_m,
                                 cv2=(b_s / b_m) ** 2, R=b.ratio(r["rho0"], rho, dp), R_post=b.ratio(r["rho0"], rho, z),
                                 R_pre=b.ratio(r["rho0"], r["rho0"], dp), R_frozen=b.ratio(r["rho0"], r["rho0"], z),
                                 R_allpot=b.ratio(r["rho0"], np.ones(len(rho)), z),
                                 n_post=len(r["postspikes"]), n_pre=len(r["prespikes"])))
            for j in range(len(rho)):
                syn_rows.append(dict(path=path, pair=r["pair"], proto=r["proto"], syn=int(r["syn"][j]),
                                     rho0=r["rho0"][j], rho_f=rho[j], dpre=dp[j], c_pre=r["c_pre"][j],
                                     c_post=r["c_post"][j], peak=float(E[j].max()), td=td[j], tp=tp[j],
                                     t_above_p=float(h[E[j] > tp[j]].sum()), t_above_d=float(h[E[j] > td[j]].sum()),
                                     gain=b.delta[j] / b_m, use_d=b.Use_d[j], use_p=b.Use_p[j]))
        print(path, c0, len(B.recs), flush=True)
        del B
out = os.path.dirname(os.path.abspath(__file__))
pd.DataFrame(syn_rows).to_csv(f"{out}/syn_{FIT}.csv", index=False)
pd.DataFrame(rec_rows).to_csv(f"{out}/rec_{FIT}.csv", index=False)
