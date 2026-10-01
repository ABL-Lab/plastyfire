"""jax_v2.GPUModel vs the CPU path (batch_v2 + model_v2) on the 8 preview pairs. GPU node only.

  per synapse: rho vs rho_full_all, dpre vs dpre_final(features) for the v2.2 preview fit and random candidates
  chi2: GPU readout vs batch_v2.objective (preview v2.2 must give 84.16)
  timing: one population of 80 candidates

    sbatch glusynapse_v2/run_jax_test.sh
"""
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); V2 = os.path.join(HERE, "..")
sys.path.insert(0, V2)
from batch_v2 import BatchV2
from targets import load_targets
import model_v2 as MV
import jax_v2 as JV

PAIRS = "180351-198084,181015-184976,182339-200396,182381-186941,182900-195902,183240-189246,183769-207759,186261-208490"


def cpu_dpre(B, P):
    out = []
    for (tT, K), r in zip(B.features(P), B.recs):
        out.append(MV.dpre_final(tT, K, P["A_mglu"], P["A_NO"], P["dpre_min"], P["dpre_max"], P["dpre0"]))
    return np.concatenate(out)


def main():
    import jax; print("devices", jax.devices(), flush=True)
    fit = json.load(open(os.path.join(V2, "results/preview_fit_v22.json")))
    T = load_targets(("markram", "nevian", "ebner"))
    B = BatchV2([os.path.join(V2, "extracted/ebner_preview"), os.path.join(V2, "extracted/markram_delta-cooker")],
                protocols=sorted({k[0].split("@")[0] for k in T}), pairs=set(PAIRS.split(",")),
                signals=("vdcc", "cacr"))
    have = {r["proto"] for r in B.recs}; T = {k: v for k, v in T.items() if k[0].split("@")[0] in have}
    filters = json.loads(fit["args"]["filters"])
    t0 = time.time(); G = JV.GPUModel(B, T, filters); print(f"GPUModel build {time.time()-t0:.1f}s, lanes {G.L}, buckets {[(b['T'], len(b['idx'])) for b in G.buckets]}, "
          f"padded work / real {G.pad_work:.2f}", flush=True)

    rng = np.random.default_rng(0)
    cands = [(fit["a"], {**filters, **fit["pre"]})]
    for _ in range(3):
        a = {k: float(v) for k, v in zip(("a00", "a01", "a10", "a11"), rng.uniform(0.5, 3, 4))}
        a.update(a20=a["a00"], a21=a["a01"], a30=a["a10"], a31=a["a11"])
        p = dict(filters, A_mglu=10 ** rng.uniform(-3, -1), A_NO=10 ** rng.uniform(0, 3),
                 theta_Ti=10 ** rng.uniform(-7, -5), theta_NOi=10 ** rng.uniform(-7, -5),
                 tau_T=rng.uniform(2, 50), tau_NO=rng.uniform(2, 30), tau_Z=rng.uniform(30, 100), theta_Z=rng.uniform(0, 3))
        cands.append((a, p))
    A = [c[0] for c in cands]; Ps = [c[1] for c in cands]
    t0 = time.time(); rho, d = G.rho_dpre(A, Ps); print(f"first call (compile + run) {time.time()-t0:.1f}s", flush=True)
    worst_r = worst_d = 0.0
    nc = G.n_ctrl
    for i, (a, p) in enumerate(cands):
        th = [B.thetas(r, a) for r in B.recs]
        rc = B.rho_full_all(np.concatenate([t[0] for t in th]), np.concatenate([t[1] for t in th]))
        P = {**MV.DEFAULTS, **p}
        dc = cpu_dpre(B, P)
        er, ed = np.abs(rho[i, :nc] - rc).max(), np.abs(d[i, :nc] - dc).max()
        flips = int(np.sum((rho[i, :nc] >= 0.5) != (rc >= 0.5)))
        print(f"cand {i}: max|rho| {er:.2e}  max|dpre| {ed:.2e}  binary flips {flips}  mean dpre {dc.mean():+.3f}", flush=True)
        worst_r, worst_d = max(worst_r, er), max(worst_d, ed)
        # mglu_block lanes == CPU dpre with A_mglu = 0
        if G.dup:
            j = next(iter(G.dup)); r = B.recs[j]
            (tT, K) = B.features(P)[j]
            ref = MV.dpre_final(tT, K, 0.0, P["A_NO"], P["dpre_min"], P["dpre_max"], P["dpre0"])
            assert np.abs(d[i, G.dup[j]] - ref).max() < 1e-9, "mglu_block lanes"
    chi_g, _ = G.chi2(rho, d, return_pred=True)
    chi_c = [B.objective(a, p, T, mode="full")[0] for a, p in cands]
    print("chi2 gpu", np.round(chi_g, 4), "\nchi2 cpu", np.round(chi_c, 4), flush=True)
    assert worst_r < 1e-9 and worst_d < 1e-9, (worst_r, worst_d)
    assert np.allclose(chi_g, chi_c, rtol=1e-8, atol=1e-8)
    assert abs(chi_g[0] - 84.1625) < 1e-3
    # stream mode (chunk_gb, all 120 pairs): small chunks force many record-aligned, lane-padded chunks
    Gs = JV.GPUModel(B, T, filters, chunk_gb=0.5)
    rs_, ds_ = Gs.rho_dpre(A, Ps)
    es = max(np.abs(rs_ - rho).max(), np.abs(ds_ - d).max())
    print(f"stream mode: {len(Gs.buckets)} chunks vs {len(G.buckets)} buckets, max|diff| {es:.1e}", flush=True)
    assert es < 1e-12 and np.allclose(Gs.chi2(rs_, ds_), chi_g, rtol=1e-12)
    del Gs
    # v2.3: NO driven by effcai (no_drive 1), post_nmdar lanes with A_NO = 0
    f3 = dict(filters, no_drive=1)
    G3 = JV.GPUModel(B, T, f3)
    c3 = []
    for k in range(3):
        a, p = cands[k + 1]
        c3.append((a, dict(p, no_drive=1, theta_NOe=10 ** rng.uniform(-2, -0.7), e_scale=10 ** rng.uniform(-3, -1))))
    rho3, d3 = G3.rho_dpre([c[0] for c in c3], [c[1] for c in c3])
    for i, (a, p) in enumerate(c3):
        P = {**MV.DEFAULTS, **p}
        dc = cpu_dpre(B, P); ed = np.abs(d3[i, :nc] - dc).max()
        print(f"v2.3 cand {i}: max|dpre| {ed:.2e}  mean dpre {dc.mean():+.3f}  post_nmdar lanes {len(G3.dup_nm)}", flush=True)
        worst_d = max(worst_d, ed)
        if G3.dup_nm:
            j = next(iter(G3.dup_nm)); (tT, K) = B.features(P)[j]
            ref = MV.dpre_final(tT, K, P["A_mglu"], 0.0, P["dpre_min"], P["dpre_max"], P["dpre0"])
            assert np.abs(d3[i, G3.dup_nm[j]] - ref).max() < 1e-9, "post_nmdar lanes"
    chi3 = G3.chi2(rho3, d3)
    chi3c = [B.objective(a, p, T, mode="full")[0] for a, p in c3]
    print("v2.3 chi2 gpu", np.round(chi3, 4), "\nv2.3 chi2 cpu", np.round(chi3c, 4), flush=True)
    assert worst_d < 1e-9 and np.allclose(chi3, chi3c, rtol=1e-8, atol=1e-8)
    del G3
    # v2.4: NO driven by spine Ca (cacr), gate tanh(pos(N - theta_N))
    G4 = JV.GPUModel(B, T, dict(filters, no_drive=2))
    c4 = []
    for k in range(3):
        a, p = cands[k + 1]
        c4.append((a, dict(p, no_drive=2, theta_NOc=10 ** rng.uniform(-5, -2.5), c_scale=1e-3,
                           theta_N=rng.uniform(0, 20), tau_NO=rng.uniform(5, 60), tau_Z=rng.uniform(60, 150),
                           theta_Z=0.0)))
    rho4, d4 = G4.rho_dpre([c[0] for c in c4], [c[1] for c in c4])
    w4 = 0.0
    for i, (a, p) in enumerate(c4):
        P = {**MV.DEFAULTS, **p}
        dc = cpu_dpre(B, P); ed = np.abs(d4[i, :nc] - dc).max()
        print(f"v2.4 cand {i}: max|dpre| {ed:.2e}  mean dpre {dc.mean():+.3f}", flush=True)
        w4 = max(w4, ed)
    chi4 = G4.chi2(rho4, d4); chi4c = [B.objective(a, p, T, mode="full")[0] for a, p in c4]
    print("v2.4 chi2 gpu", np.round(chi4, 4), "\nv2.4 chi2 cpu", np.round(chi4c, 4), flush=True)
    # the GPU recovers spine Ca inside the scan (XLA exp vs numpy exp, then the float32 cast): a rare 1-ulp
    # float32 difference in u, so this path is checked to 1e-6 instead of bit level
    assert w4 < 1e-6 and np.allclose(chi4, chi4c, rtol=1e-6, atol=1e-6), w4
    del G4
    # T25: T driven by effcai / c_post (t_drive 1), free dpre_min; paired_l5 targets (mglu_block, nmdar_block)
    T5 = {k: v for k, v in load_targets(("paired_l5",)).items() if k[0].split("@")[0] in have}
    G5 = JV.GPUModel(B, T5, dict(filters, t_drive=1))
    c5 = []
    for k in range(3):
        a, p = cands[k + 1]
        c5.append((a, dict(p, t_drive=1, theta_Te=10 ** rng.uniform(-1, 0.5), tau_T=10 ** rng.uniform(0.5, 2.5),
                           dpre_min=rng.uniform(-0.8, -0.1), A_mglu=10 ** rng.uniform(-2, -0.5))))
    rho5, d5 = G5.rho_dpre([c[0] for c in c5], [c[1] for c in c5])
    w5 = 0.0
    for i, (a, p) in enumerate(c5):
        P = {**MV.DEFAULTS, **p}
        dc = cpu_dpre(B, P); ed = np.abs(d5[i, :nc] - dc).max()
        print(f"T25 cand {i}: max|dpre| {ed:.2e}  mean dpre {dc.mean():+.3f}", flush=True)
        w5 = max(w5, ed)
    chi5 = G5.chi2(rho5, d5); chi5c = [B.objective(a, p, T5, mode="full")[0] for a, p in c5]
    print("T25 chi2 gpu", np.round(chi5, 4), "\nT25 chi2 cpu", np.round(chi5c, 4), flush=True)
    assert w5 < 1e-9 and np.allclose(chi5, chi5c, rtol=1e-8, atol=1e-8), w5
    del G5
    # T25: cell-level eCB from post APs (t_drive 2), gate threshold theta_Tg
    G6 = JV.GPUModel(B, T5, dict(filters, t_drive=2))
    c6 = []
    for k in range(3):
        a, p = cands[k + 1]
        c6.append((a, dict(p, t_drive=2, theta_Te=rng.uniform(0.3, 1.5), tau_T=10 ** rng.uniform(1, 2.3),
                           theta_Tg=rng.uniform(0, 3), dpre_min=rng.uniform(-0.8, -0.1), A_mglu=10 ** rng.uniform(-1.5, 0))))
    rho6, d6 = G6.rho_dpre([c[0] for c in c6], [c[1] for c in c6])
    w6 = 0.0
    for i, (a, p) in enumerate(c6):
        P = {**MV.DEFAULTS, **p}
        dc = cpu_dpre(B, P); ed = np.abs(d6[i, :nc] - dc).max()
        print(f"T25 t_drive 2 cand {i}: max|dpre| {ed:.2e}  mean dpre {dc.mean():+.3f}", flush=True)
        w6 = max(w6, ed)
    chi6 = G6.chi2(rho6, d6); chi6c = [B.objective(a, p, T5, mode="full")[0] for a, p in c6]
    print("t_drive 2 chi2 gpu", np.round(chi6, 4), "\nt_drive 2 chi2 cpu", np.round(chi6c, 4), flush=True)
    assert w6 < 1e-9 and np.allclose(chi6, chi6c, rtol=1e-8, atol=1e-8), w6
    del G6
    A80 = [A[i % 4] for i in range(80)]; P80 = [Ps[i % 4] for i in range(80)]
    G.evaluate(A80, P80)
    t0 = time.time(); G.evaluate(A80, P80); t80 = time.time() - t0
    t0 = time.time(); rho, d = G.rho_dpre(A80, P80); tg = time.time() - t0
    print(f"PASS: rho/dpre match to {max(worst_r, worst_d):.1e}; chi2 equal; 80 candidates: {t80:.1f}s "
          f"(GPU scan {tg:.1f}s, readout {t80 - tg:.1f}s)")


if __name__ == "__main__":
    main()
