"""Seconds-long GPU check of the v4 post-rule prototype (synthetic data, as tests/smoke_gpu_v3.py):
  1. thetas with rho_gamma 1 == batch_v2.BatchV2.thetas bit for bit; gamma 0 == a00 c_pre + a01 C_REF.
  2. v4 kernel (fast / vgate / both) == CPU reference rho_v4d.rho_loop_v4d (rho) and == gpu_v3 (dpre, 3 variants);
  3. gate open at every step == gpu_v3 bit for bit; CPU with every option off == gpu_v3.
  4. rates kernel (gamma_d, gamma_p per candidate) at the GB values == constant-rate kernel (<1e-12, 4 variants);
     other rates (gp only, and gd+gp at the bounds) change rho and == CPU rho_v4d.rho_loop_v4d (<1e-10).
  (22108540 passed the first version: fast kernel only, tau_fast = tau_effca == v3 to 0.0.)
    sbatch glusynapse_v2/rho_redesign/run_smoke_gpu_v4.sh
"""
import os, sys
from types import SimpleNamespace
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); V2 = os.path.dirname(HERE)
sys.path.insert(0, V2); sys.path.insert(0, HERE)
import model_v2 as MV
import gpu_v3, gpu_v4_rho, rho_v4
from batch_v2 import BatchV2, TAU_IND_GB, TAU_EFFCA, GAMMA_D_GB, GAMMA_P_GB, RHO_STAR_GB

rng = np.random.default_rng(0)
T, n, P = 4000, 4, 3
dt = np.where(np.arange(T) < 2000, 0.025, np.where(np.arange(T) < 3000, 0.25, 5.0)); dt[-1] = 0.0   # ms
h = dt / 1000.0 / TAU_IND_GB
# effcai from a spiky spine Ca through the exact tau_effca filter (as extract.effcai_from_cai)
ca = (rng.random((n, T)) ** 30) * 2e-2
E = np.zeros((n, T))
for k in range(T - 1):
    if dt[k] > 0:
        a = np.exp(-dt[k] / TAU_EFFCA); E[:, k + 1] = a * E[:, k] + TAU_EFFCA * (1 - a) * ca[:, k]
E = E.astype(np.float32)
S = (rng.random((n, T)) ** 6 * 5e-5).astype(np.float32)
C = np.zeros((n, T), np.int8); C[:, rng.integers(0, T - 1, 40)] = 1
PC = np.zeros((n, T), np.float32); PC[:, rng.integers(0, T - 1, 30)] = 1.0
rho0 = np.array([0.0, 1.0, 0.0, 1.0]); cpost = np.array([3e-3, 1e-4, 2e-2, 8e-3]); cpre = np.array([0.05, 0.03, 0.1, 0.06])
m = float(np.median(E.max(1)))
td = rng.uniform(0.2, 0.5, (P, n)) * m; tp = td + rng.uniform(0.1, 0.5, (P, n)) * m
worst = {}

# 1. thetas
a = dict(a00=1.04, a01=0.65, a10=1.16, a11=4.8)
r = dict(c_pre=cpre, c_post=cpost, syn=np.arange(n))
t3 = BatchV2.thetas(r, {**a, "a20": 1, "a21": 1, "a30": 1, "a31": 1})
t4 = rho_v4.thetas(cpre, cpost, a, 1.0)
assert all(np.array_equal(x, y) for x, y in zip(t3, t4)), "gamma 1 != v3 thetas"
t0 = rho_v4.thetas(cpre, cpost, a, 0.0)
assert np.allclose(t0[1], 1.16 * cpre + 4.8 * rho_v4.C_REF, rtol=1e-15)
print("1. thetas: gamma 1 bitwise v3, gamma 0 ok", flush=True)

# kernels
P0 = {**MV.DEFAULTS, "pre_drive": 1, "i_scale": 1e-5, "t_drive": 4, "no_drive": 0, "tau_E1": 100.0}
Ps = [dict(A_mglu=10 ** rng.uniform(-2, -1), A_NO=10 ** rng.uniform(0, 2), theta_NOi=2e-6, theta_Te=0.6,
           tau_T=rng.uniform(20, 60), tau_NO=10.0, tau_Z=25.0, theta_Z=0.4, theta_Tg=0.5, dpre_min=-0.3,
           tau_fast=tf) for tf in (10.0, 40.0, float(TAU_EFFCA))]
g = SimpleNamespace(P0=P0, t_mode=4, no_mode=0, fast=True)
prm3 = gpu_v3.GPUModelV3._params(g, Ps)
prm4 = np.ascontiguousarray(np.concatenate([prm3, np.array([[p["tau_fast"]] for p in Ps])], axis=1))
a_st = [np.flatnonzero(C[i]) for i in range(n)]; p_st = [np.flatnonzero(PC[i]) for i in range(n)]
aptr = np.r_[0, np.cumsum([len(x) for x in a_st])].astype(np.int64); pptr = np.r_[0, np.cumsum([len(x) for x in p_st])].astype(np.int64)
args = (E.ravel(), S.ravel(), h, np.arange(n + 1, dtype=np.int64) * T, np.full(n, T, np.int64), np.zeros(n, np.int64), aptr,
        np.concatenate(a_st).astype(np.int64), np.concatenate([C[i, a_st[i]] for i in range(n)]).astype(np.float64), pptr,
        np.concatenate(p_st).astype(np.int64), np.concatenate([PC[i, p_st[i]] for i in range(n)]).astype(np.float64),
        1.0 / cpost, rho0, np.arange(n, dtype=np.int64), td, tp)


def run(kern, prm):
    orr = gpu_v3.cuda.device_array((P, n)); od = gpu_v3.cuda.device_array((P, n, 3))
    kern[(1, 1), (32, 4)](*args, prm, orr, od)
    return orr.copy_to_host(), od.copy_to_host()


# ---- final kernel signature (fast / vgate flags), vs gpu_v3 and vs the CPU reference rho_v4d
import rho_v4d
gl = [np.unique(rng.integers(0, T - 1, 12)) for _ in range(n)]
gptr = np.r_[0, np.cumsum([len(x) for x in gl])].astype(np.int64); gstep = np.concatenate(gl).astype(np.int64)
gcnt = np.ones(len(gstep)); gcnt[::5] = 2.0
nog = (np.zeros(n + 1, np.int64), np.zeros(1, np.int64), np.zeros(1))


def run2(kern, prm, gev):
    orr = gpu_v3.cuda.device_array((P, n)); od = gpu_v3.cuda.device_array((P, n, 3))
    kern[(1, 1), (32, 4)](*args, prm, orr, od, *gev)
    return orr.copy_to_host(), od.copy_to_host()


r3, d3 = run(gpu_v3._make_kernel(4, 0, True), prm3)
H = h[None, :].astype(np.float64); rec = np.zeros(n, np.int64); L = np.full(n, T - 1, np.int64)
cpu = lambda k, tf, vg, gev: rho_v4d.rho_loop_v4d(E, H, rec, L, td[k], tp[k], rho0, float(GAMMA_D_GB), float(GAMMA_P_GB),
                                                  float(RHO_STAR_GB), tf, float(TAU_EFFCA), 1e3 * TAU_IND_GB, vg, 100.0, *gev)
for fast, vg in ((True, False), (False, True), (True, True)):
    kern = gpu_v4_rho._make_kernel_fast(4, 0, True, fast, vg)
    gev = (gptr, gstep, gcnt) if vg else nog
    r4, d4 = run2(kern, prm4 if fast else prm3, gev)
    assert np.array_equal(d3, d4), "dpre changed"
    for k, p in enumerate(Ps):
        worst[f"f{int(fast)}g{int(vg)}_{k}"] = float(np.abs(cpu(k, p["tau_fast"] if fast else 0.0, vg, gev) - r4[k]).max())
    if vg:
        print(f"   fast {fast} vgate: rho {np.round(r4, 4).tolist()}  (v3 {np.round(r3, 4).tolist()})", flush=True)
for k in range(P):
    worst[f"cpu_off_v3_{k}"] = float(np.abs(cpu(k, 0.0, False, nog) - r3[k]).max())
print("GPU vs CPU reference:", {k: f"{v:.1e}" for k, v in worst.items()}, flush=True)
# gate with an event at every step and tau_fast off == v3
allg = (np.arange(n + 1, dtype=np.int64) * T, np.tile(np.arange(T, dtype=np.int64), n), np.ones(n * T))
r4, _ = run2(gpu_v4_rho._make_kernel_fast(4, 0, True, False, True), prm3, allg)
e0 = float(np.abs(r4 - r3).max()); print(f"gate always open == v3: {e0:.1e}")
assert max(worst.values()) < 1e-10 and e0 == 0.0, (worst, e0)
print("PASS smoke v4 gate", flush=True)

# 4. free rho rates (fit_v4 --fit-gamma): gamma_d, gamma_p as prm columns after tau_fast
GB = np.array([[float(GAMMA_D_GB), float(GAMMA_P_GB)]] * P)
RT = np.array([[float(GAMMA_D_GB), float(GAMMA_P_GB)], [50.0, 300.0], [200.0, 150.0]])   # per candidate, bounds
args_short = args[:4] + (np.full(n, 2400, np.int64),) + args[5:]
cat = lambda *xs: np.ascontiguousarray(np.concatenate(xs, axis=1))
wr = {}
for fast, vg in ((False, False), (True, False), (False, True), (True, True)):
    gev = (gptr, gstep, gcnt) if vg else nog
    base = prm4 if fast else prm3
    ref = r3 if not (fast or vg) else run2(gpu_v4_rho._make_kernel_fast(4, 0, True, fast, vg), base, gev)[0]
    kr = gpu_v4_rho._make_kernel_fast(4, 0, True, fast, vg, True)
    rg, dg = run2(kr, cat(base, GB), gev)
    tag = f"f{int(fast)}g{int(vg)}"
    wr[f"GB_{tag}"] = float(np.abs(rg - ref).max())
    assert np.array_equal(dg, d3), "dpre changed (rates)"
    rr, _ = run2(kr, cat(base, RT), gev)
    for k, p in enumerate(Ps):
        c = rho_v4d.rho_loop_v4d(E, H, rec, L, td[k], tp[k], rho0, RT[k, 0], RT[k, 1], float(RHO_STAR_GB),
                                 p["tau_fast"] if fast else 0.0, float(TAU_EFFCA), 1e3 * TAU_IND_GB, vg, 100.0, *gev)
        wr[f"cpu_{tag}_{k}"] = float(np.abs(c - rr[k]).max())
    # a changed gp / gd must change rho (full run, or mid-run when the full run settles at the same 0/1 states)
    orr = gpu_v3.cuda.device_array((P, n)); od = gpu_v3.cuda.device_array((P, n, 3))
    kr[(1, 1), (32, 4)](*args_short, cat(base, GB), orr, od, *gev); s0 = orr.copy_to_host()
    kr[(1, 1), (32, 4)](*args_short, cat(base, RT), orr, od, *gev); s1 = orr.copy_to_host()
    ch = max(float(np.abs(rr[1:] - rg[1:]).max()), float(np.abs(s1[1:] - s0[1:]).max()))
    wr[f"change_{tag}"] = ch
    print(f"   rates {tag}: rho GB {np.round(rg, 4).tolist()}  changed {np.round(rr, 4).tolist()}  "
          f"(mid-run change {float(np.abs(s1[1:] - s0[1:]).max()):.2e})", flush=True)
    assert ch > 1e-6, f"gamma change has no effect ({tag})"
# gp only (gd at GB): must change rho too, and match the CPU reference
RP = np.array([[float(GAMMA_D_GB), float(GAMMA_P_GB)], [float(GAMMA_D_GB), 150.0], [float(GAMMA_D_GB), 300.0]])
kr = gpu_v4_rho._make_kernel_fast(4, 0, True, False, False, True)
rp, _ = run2(kr, cat(prm3, RP), nog)
orr = gpu_v3.cuda.device_array((P, n)); od = gpu_v3.cuda.device_array((P, n, 3))
kr[(1, 1), (32, 4)](*args_short, cat(prm3, GB), orr, od, *nog); s0 = orr.copy_to_host()
kr[(1, 1), (32, 4)](*args_short, cat(prm3, RP), orr, od, *nog); s1 = orr.copy_to_host()
wr["gp_only_change"] = max(float(np.abs(rp[1:] - r3[1:]).max()), float(np.abs(s1[1:] - s0[1:]).max()))
for k in range(P):
    c = rho_v4d.rho_loop_v4d(E, H, rec, L, td[k], tp[k], rho0, RP[k, 0], RP[k, 1], float(RHO_STAR_GB), 0.0,
                             float(TAU_EFFCA), 1e3 * TAU_IND_GB, False, 100.0, *nog)
    wr[f"gp_only_cpu_{k}"] = float(np.abs(c - rp[k]).max())
print("rates:", {k: f"{v:.1e}" for k, v in wr.items()}, flush=True)
assert all(wr[k] < 1e-12 for k in wr if k.startswith("GB_")), wr
assert all(wr[k] < 1e-10 for k in wr if "cpu" in k), wr
assert wr["gp_only_change"] > 1e-6, wr
print("PASS smoke v4 rates")

# 5. ROUND2 C1 / C2 VDCC-amplitude gate (vamp kernel; prm columns theta_V, i_scale after the rates columns) vs the CPU
#    ground truth scan_vgate_amp.rho_rec, on the synthetic traces (S = the VDCC current). theta_V 0 == gate off bitwise.
import scan_vgate_amp as SV
ISC, TAUE = 1e-5, 100.0
aV = np.exp(-dt[:-1] / TAUE); bV = TAUE * (1.0 - aV) / ISC
VD = S.astype(np.float64)


def cpu_vamp(k, mode, th, gd, gp):
    out = np.empty((2, 1, n)); vm = np.empty(n)
    SV.rho_rec(E, VD, h[:-1], aV, bV, td[k].astype(float), tp[k].astype(float), rho0, gd, gp, float(RHO_STAR_GB),
               np.array([th]), out, vm)
    return out[mode - 1, 0], vm


vm0 = np.median(cpu_vamp(0, 1, 0.0, float(GAMMA_D_GB), float(GAMMA_P_GB))[1])
THV = np.array([0.0, 0.3 * vm0, 1.0 * vm0])          # per candidate: off, partly closed, mostly closed
print(f"5. vamp: median V at crossings {vm0:.4g}, theta_V per candidate {THV.tolist()}", flush=True)
wv = {}
r_rt, d_rt = run2(gpu_v4_rho._make_kernel_fast(4, 0, True, False, False, True), cat(prm3, RT), nog)
r_nr, d_nr = run2(gpu_v4_rho._make_kernel_fast(4, 0, True, False, False, False), prm3, nog)   # v4 kernel, all off
wv["nr_vs_v3"] = float(np.abs(r_nr - r3).max())
for mode in (1, 2):
    for rates_ in (False, True):
        kv = gpu_v4_rho._make_kernel_fast(4, 0, True, False, False, rates_, mode)
        base = cat(prm3, RT) if rates_ else prm3
        tag = f"C{mode}{'_rates' if rates_ else ''}"
        # gate off (theta_V 0 for every candidate) == the kernel without vamp, bit for bit
        r0, d0 = run2(kv, cat(base, np.array([[0.0, ISC]] * P)), nog)
        ref = r_rt if rates_ else r_nr
        wv[f"off_{tag}"] = float(np.abs(r0 - ref).max())
        assert np.array_equal(r0, ref) and np.array_equal(d0, d3), f"vamp {tag} theta_V 0 != no-vamp kernel"
        rv, dv = run2(kv, cat(base, np.c_[THV, np.full(P, ISC)]), nog)
        assert np.array_equal(dv, d3), "dpre changed (vamp)"
        for k in range(P):
            gd, gp = (RT[k] if rates_ else (float(GAMMA_D_GB), float(GAMMA_P_GB)))
            wv[f"cpu_{tag}_{k}"] = float(np.abs(cpu_vamp(k, mode, THV[k], float(gd), float(gp))[0] - rv[k]).max())
        wv[f"change_{tag}"] = float(np.abs(rv[1:] - ref[1:]).max())
        print(f"   {tag}: rho {np.round(rv, 4).tolist()}  (off {np.round(ref, 4).tolist()})", flush=True)
print("vamp:", {k: f"{v:.1e}" for k, v in wv.items()}, flush=True)
assert all(wv[k] < 1e-10 for k in wv if k.startswith("cpu_")), wv
assert any(wv[k] > 1e-6 for k in wv if k.startswith("change_")), "theta_V > 0 has no effect"
print("PASS smoke v4 vamp", flush=True)

# 6. --real: C1 / C2 at theta_V 2, 5, 10 on the A0g_s3 parameters, L5 (29 targets): GPU fit_v4 objective == CPU
#    scan_vgate_amp.run (computed here) to 1e-9; vamp off (gpu_v4 rates kernel) == C1 / C2 theta_V 0 == json chi2.
if "--real" in sys.argv:
    import json, gc
    import batch_v2, fit_v4
    from targets import load_targets
    FIT = os.path.join(HERE, "results", "v4_A0g_s3.json")
    fit = json.load(open(FIT)); fa = fit["args"]; batch_v2.BAP_GATE = fit.get("bap_gate")
    fil = {**json.loads(fa["filters"]), **json.loads(fa.get("set", "{}"))}
    Pc = {**MV.DEFAULTS, **fil, **fit["pre"]}
    conds = set(fa["conditions"].split(","))
    pairs = set(fa["pairs"].split(","))
    TH = [0.0, 2.0, 5.0, 10.0]
    # CPU ground truth first (freed before the GPU model is built)
    T5 = {k: v for k, v in load_targets(tuple(fa["groups"].split(","))).items() if k[1] in conds}
    batch_v2.BASIS_DIR = fit_v4.L5_BASIS
    SV.THETAS = np.array(TH)
    Bc = BatchV2(fa["dirs"].split(","), protocols=sorted({k[0].split("@")[0] for k in T5}), pairs=pairs, fast=False,
                 signals=("vdcc",))
    for rr in Bc.recs:
        Bc.basis(rr)
    R5, _ = SV.run(Bc, fit, Pc, T5, sorted({c for _, c in T5}), "L5", {}, {})
    cpu5 = {(c, float(t)): float((x.z ** 2).sum()) for (c, t), x in R5.groupby(["cand", "theta_V"])}
    ncpu = {(c, float(t)): len(x) for (c, t), x in R5.groupby(["cand", "theta_V"])}
    del Bc, R5; gc.collect()
    # GPU: one model (traces resident), kernels swapped per mode
    G, Tg, _ = fit_v4.build(fa["dirs"], fa["groups"], conds, {**fil, "vamp_mode": 1}, pairs, fit_v4.L5_BASIS, False,
                            rates=True)
    a = fit["a"]; Pg = {**fil, **fit["pre"]}
    G.vamp = 0; G._kern = gpu_v4_rho._make_kernel_fast(G.t_mode, G.no_mode, G.use_s, False, False, True)
    c_off = float(G.evaluate([a], [Pg])[0])
    print(f"6. real L5 ({len(Tg)} targets GPU): vamp off {c_off:.9f}  json chi2 {fit['chi2']:.9f}", flush=True)
    res = {}
    for mode in (1, 2):
        G.vamp = mode; G._kern = gpu_v4_rho._make_kernel_fast(G.t_mode, G.no_mode, G.use_s, False, False, True, mode)
        cs = G.evaluate([a] * len(TH), [{**Pg, "vamp_mode": mode, "theta_V": t} for t in TH])
        for t, c in zip(TH, cs):
            k = (f"C{mode}", t); res[k] = float(c)
            print(f"   C{mode} theta_V {t:5.1f}: GPU {c:.9f}  CPU {cpu5[k]:.9f} ({ncpu[k]} targets)  diff {abs(c - cpu5[k]):.1e}",
                  flush=True)
    assert res[("C1", 0.0)] == c_off and res[("C2", 0.0)] == c_off, "theta_V 0 != vamp off"
    assert abs(c_off - fit["chi2"]) < 1e-9, (c_off, fit["chi2"])
    assert all(abs(res[k] - cpu5[k]) < 1e-9 for k in res), "GPU != CPU scan"
    print("PASS smoke v4 vamp real", flush=True)
