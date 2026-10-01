"""Debug the v2.4 gate case: mod N_GB and cai_CR vs offline N from cacr recovered from effcai."""
import os, sys, numpy as np
sys.path.insert(0, "glusynapse_v2/tests"); sys.path.insert(0, "glusynapse_v2")
import test_v2_single_synapse as S, test_offline_gate as G, model_v2 as MV
from batch_v2 import cacr_from_effcai
S._snapshot_defaults()
dt_pp, v, par = G.CASES_VDCC["cacr NO thN"]
pre, steps = S.pairing(dt_pp, v=v, dur=2.0)
r = S.run(pre, steps, par)
P = {k[:-3]: val for k, val in {**S.DEFAULTS, **par}.items()}
n = len(r["v2.effcai_GB"]); t = MV.grid(n, 0.025)
u_true = r["v2.cai_CR"] - 70e-6
e = r["v2.effcai_GB"].astype(np.float64); h = np.diff(t); a = np.exp(-h / 200.0); b = 200.0 * (1 - a)
u64 = np.zeros(n); u64[:-1] = (e[1:] - a * e[:-1]) / b
u32 = cacr_from_effcai(r["v2.effcai_GB"][None].astype(np.float32), t)[0]
for name, u in (("f64", u64), ("f32", u32)):
    for lag in (0, 1):
        ut = np.roll(u_true, -lag)
        print(f"{name} lag {lag}: max|u - (cai_CR - min)| {np.abs(u[:-2] - ut[:-2]).max():.2e}  (peak {u_true.max():.2e})")
arr = MV.arrival_index(pre, [0.1], t)
tT = MV.feature_T((-r["v2.ica_VDCC"])[None], t, arr, P)
u_lib = cacr_from_effcai(r["v2.effcai_GB"][None].astype(np.float64), t)[0]
print("lib dtype", u_lib.dtype, "effcai dtype", r["v2.effcai_GB"].dtype, "max|lib - u64|", np.abs(u_lib - u64).max())
for name, sig in (("lib f64 in", u_lib), ("true cai_CR", u_true), ("true lag1", np.roll(u_true, -1)), ("recovered f64", u64), ("recovered f32", u32)):
    N = MV._lowpass_grid(MV.drive(sig[None], P, "NO"), t, P["tau_NO"])[0]
    K = MV.feature_K(sig[None], t, arr, P)
    d = MV.dpre_final(tT, K, P["A_mglu"], P["A_NO"], P["dpre_min"], P["dpre_max"])[0]
    print(f"{name:14s} N max {N.max():.4f} vs mod {r['v2.N_GB'].max():.4f}; dpre {d:+.5f} vs mod {r['v2.dpre_GB'][-1]:+.5f}")
print("P used:", {k: P[k] for k in ("no_drive", "theta_NOc", "c_scale", "theta_N", "tau_NO", "tau_Z", "theta_Z", "A_NO")})
