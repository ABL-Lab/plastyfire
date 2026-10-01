"""Seconds-long GPU check: gpu_v3 kernel vs jax_v2's own scan (GPUModel._make_scan) on synthetic data, for every supported
(t_drive, no_drive, pre_drive). 1 record, 3 synapses, 3 candidates; v2 lanes = 3 control + 3 A_mglu 0 + 3 A_NO 0.
    sbatch glusynapse_v2/run_smoke_gpu_v3.sh
"""
import os, sys
from types import SimpleNamespace
import numpy as np
V2 = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, V2)
import model_v2 as MV
import gpu_v3
from batch_v2 import TAU_IND_GB

rng = np.random.default_rng(0)
T, n, P = 3000, 3, 3
h = np.where(np.arange(T) < 1500, 0.025, 0.1) / 1000.0 / TAU_IND_GB; h[-1] = 0.0      # two step sizes, last 0
E = (rng.random((n, T)) ** 8).astype(np.float32)                                           # effcai, spiky
S = (rng.random((n, T)) ** 6 * 5e-5).astype(np.float32)                                    # vdcc drive
C = np.zeros((n, T), np.int8); C[:, rng.integers(0, T, 40)] = 1; C[0, 100] = 2; C[1, T - 1] = 1
PC0 = np.zeros((n, T), np.float32); PC0[:, rng.integers(0, T, 30)] = rng.random(30).astype(np.float32) * 2; PC0[2, 7] = 0.4
rho0 = np.array([0.2, 0.7, 1.0]); cpost = np.array([0.3, 0.5, 0.8])
td = rng.uniform(0.2, 0.4, (P, n)); tp = td + rng.uniform(0.0, 0.3, (P, n))
worst = 0.0
MODES = ((4, 0, 1), (3, 0, 1), (2, 0, 1), (0, 0, 1), (0, 0, 0), (1, 1, 1), (1, 0, 0))
CASES = {}
# phase 1: v3 only, before jax is imported (checks numba-cuda alone)
for td_mode, no_mode, pre in MODES:
    PC = np.tile(PC0[0], (n, 1)) if td_mode == 2 else PC0          # t_drive 2: one impulse train per record
    P0 = {**MV.DEFAULTS, "pre_drive": pre, "i_scale": 1e-5, "t_drive": td_mode, "no_drive": no_mode, "tau_E1": 100.0}
    Ps = [dict(A_mglu=10 ** rng.uniform(-1, 0), A_NO=10 ** rng.uniform(0, 2), theta_Ti=1e-6, theta_NOi=2e-6,
               theta_Te=rng.uniform(0.1, 1.0), tau_T=rng.uniform(3, 50), tau_NO=rng.uniform(3, 20), tau_Z=rng.uniform(25, 60),
               theta_Z=rng.uniform(0, 0.8), theta_Tg=rng.uniform(0, 0.5), dpre_min=-0.5, theta_N=rng.uniform(0, 0.5),
               theta_NOe=0.01, e_scale=0.01) for _ in range(P)]
    # v3 kernel
    g = SimpleNamespace(P0=P0, t_mode=td_mode, no_mode=no_mode)
    prm = gpu_v3.GPUModelV3._params(g, Ps)
    use_s = td_mode == 0 or no_mode == 0
    kern = gpu_v3._make_kernel(td_mode, no_mode, use_s)
    a_st = [np.flatnonzero(C[i]) for i in range(n)]; p_st = [np.flatnonzero(PC[i]) for i in range(n)]
    aptr = np.r_[0, np.cumsum([len(x) for x in a_st])]; pptr = np.r_[0, np.cumsum([len(x) for x in p_st])]
    orr = gpu_v3.cuda.device_array((P, n)); od = gpu_v3.cuda.device_array((P, n, 3))
    imp = td_mode in (2, 3, 4)
    kern[(1, 1), (32, 4)](E.ravel(), S.ravel(), h, np.arange(n + 1, dtype=np.int64) * T, np.full(n, T, np.int64),
                          np.zeros(n, np.int64), aptr.astype(np.int64), np.concatenate(a_st).astype(np.int64),
                          np.concatenate([C[i, a_st[i]] for i in range(n)]).astype(np.float64), pptr.astype(np.int64),
                          np.concatenate(p_st).astype(np.int64),
                          np.concatenate([PC[i, p_st[i]] for i in range(n)]).astype(np.float64),
                          1.0 / cpost, rho0, np.arange(n, dtype=np.int64), td, tp, prm, orr, od)
    r3, d3 = orr.copy_to_host(), od.copy_to_host()
    CASES[(td_mode, no_mode, pre)] = (P0, Ps, PC, r3, d3)
    print(f'v3 ok t_drive {td_mode} no_drive {no_mode} pre_drive {pre}', flush=True)

# phase 2: jax in the same process (the equality test needs both)
import jax, jax.numpy as jnp
import jax_v2 as JV
print("jax devices", jax.devices(), flush=True)
for (td_mode, no_mode, pre), (P0, Ps, PC, r3, d3) in CASES.items():
    # v2 scan on 3n lanes
    obj = SimpleNamespace(P0=P0)
    run = jax.jit(JV.GPUModel._make_scan(obj, 8))
    Tb = -(-T // 8) * 8; pad = lambda x: np.concatenate([x, np.zeros((Tb - T,) + x.shape[1:], x.dtype)])
    lane = lambda x: np.tile(x, 3)
    xs = (pad(np.tile(E.T, 3)), pad(np.tile(S.T, 3)), pad(np.tile(C.T, 3)), pad(h[:, None]), pad(np.tile(PC.T, 3)))
    mg = np.r_[np.ones(n), np.zeros(n), np.ones(n)]; no = np.r_[np.ones(n), np.ones(n), np.zeros(n)]
    keys = ("A_mglu", "A_NO", "theta_Ti", "theta_NOi", "i_scale", "ca_sh_rest", "theta_T", "theta_NO", "ca_scale", "tau_T",
            "tau_NO", "tau_Z", "theta_Z", "dpre_min", "dpre_max", "dpre0", "theta_NOe", "e_scale", "theta_NOc", "c_scale",
            "theta_N", "theta_Te", "Te_scale", "tau_E1", "theta_Tg")
    pr = {k: jnp.asarray(np.array([[float({**P0, **p}[k])] for p in Ps])) for k in keys}
    r2, d2 = run(jnp.asarray(np.tile(td, 3)), jnp.asarray(np.tile(tp, 3)), pr, tuple(jnp.asarray(x) for x in xs),
                 jnp.zeros(3 * n, jnp.int32), jnp.asarray(mg), jnp.asarray(no), jnp.asarray(lane(rho0)),
                 jnp.asarray(lane(1.0 / cpost)))
    r2, d2 = np.asarray(r2), np.asarray(d2)
    er = np.abs(r3 - r2[:, :n]).max()
    ed = max(np.abs(d3[:, :, v] - d2[:, v * n:(v + 1) * n]).max() for v in range(3))
    worst = max(worst, er, ed)
    print(f"t_drive {td_mode} no_drive {no_mode} pre_drive {pre}: max|rho| {er:.1e} max|dpre| {ed:.1e}  "
          f"rho {np.round(r3[0], 4)} dpre {np.round(d3[0, :, 0], 4)} (mglu0 {np.round(d3[0, :, 1], 4)}, NO0 {np.round(d3[0, :, 2], 4)})",
          flush=True)
assert worst < 1e-10, worst
print(f"PASS smoke gpu_v3 == jax_v2 scan (max abs {worst:.1e})")
