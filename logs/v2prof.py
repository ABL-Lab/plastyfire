import sys, time, json, cProfile, pstats
sys.path.insert(0, "/project/rrg-emuller/dhuruva/plastyfire/glusynapse_v2")
import numpy as np
from batch_v2 import BatchV2
from targets import load_targets
import model_v2 as MV
V2 = "/project/rrg-emuller/dhuruva/plastyfire/glusynapse_v2"
fit = json.load(open(f"{V2}/results/preview_fit_v22.json"))
T = load_targets(("markram", "nevian", "ebner"))
B = BatchV2([f"{V2}/extracted/ebner_preview", f"{V2}/extracted/markram_delta-cooker"],
            protocols=sorted({k[0].split("@")[0] for k in T}), pairs=set(fit["args"]["pairs"].split(",")))
have = {r["proto"] for r in B.recs}; T = {k: v for k, v in T.items() if k[0].split("@")[0] in have}
P = {**json.loads(fit["args"]["filters"]), **fit["pre"]}
t0 = time.time(); B._stack(); print("stack", time.time() - t0)
t0 = time.time(); B.features(P); print("one feature pass (T+K, all recs)", time.time() - t0)
th = [B.thetas(r, fit["a"]) for r in B.recs]
td = np.concatenate([t[0] for t in th]); tp = np.concatenate([t[1] for t in th])
B.rho_full_all(td, tp); t0 = time.time(); B.rho_full_all(td, tp); print("rho numba serial", time.time() - t0)
B.rho_full_all(td, tp, parallel=True); t0 = time.time(); B.rho_full_all(td, tp, parallel=True); print("rho numba parallel", time.time() - t0)
E, H = B._E.T, B._H[B._rec].T; t0 = time.time(); rho = B._rho0.copy()
for k in range(E.shape[0]):
    e = E[k]; dep = e > td; pot = e > tp
    rho += H[k] * (-rho * (1 - rho) * (0.5 - rho) + pot * 101.5 * (1 - rho) - dep * (1 - pot) * 199.77 * rho)
    np.clip(rho, 0.0, 1.0, out=rho)
print("rho old numpy loop", time.time() - t0, "N", len(td), "T", E.shape[0])
chi2, _ = B.objective(fit["a"], P, T, mode="full"); print("chi2", chi2)
pr = cProfile.Profile(); pr.enable(); B.objective(fit["a"], P, T, mode="full"); pr.disable()
pstats.Stats(pr).sort_stats("cumulative").print_stats(15)
