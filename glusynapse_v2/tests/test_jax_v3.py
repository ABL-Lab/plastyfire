"""gpu_v3.GPUModelV3 vs jax_v2.GPUModel (as the running fits: CHUNK 4) on the pl5sj07 t_drive 4 setup (5 -vca dirs,
paired_l5 + sjostrom07, 5 conditions). Also profiles v2 (compile vs run, cost vs number of candidates P) and
counts admissible candidates in the initial and checkpoint populations. GPU node only.

    sbatch glusynapse_v2/run_test_jax_v3.sh          (TDRIVE=3 for the t_drive 3 setup)
Pass: chi2 relative difference < 1e-6 and per-target prediction difference < 1e-8 at every candidate.
"""
import os, sys, json, time, glob, shutil, tempfile
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); V2 = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, V2)
from batch_v2 import BatchV2
from targets import load_targets
import fit_v2 as F
import fit_v3
import gpu_v3

X = os.path.join(V2, "extracted")
DIRS = [os.path.join(X, d) for d in ("ebner_delta-prefire-vca", "markram_delta-prefire-vca", "sj03_delta-prefire-vca",
                                     "sj03r50_delta-prefire-vca", "sj07_delta-prefire-vca")]
FREE = "theta_Te,tau_T,theta_Tg,dpre_min,theta_NOi,tau_NO,tau_Z,theta_Z".split(",")
TD = int(os.environ.get("TDRIVE", 4))
FIL = {"pre_drive": 1, "i_scale": 1e-5, "t_drive": TD, "tau_E1": 100.0}
RES = os.path.join(V2, "results")


def load_ckpt(f):
    tmp = tempfile.mktemp(suffix=".npz")
    for _ in range(5):                    # running fits rewrite these every generation
        try:
            shutil.copy(f, tmp); z = np.load(tmp); return {k: z[k] for k in z.files}
        except Exception as e:  # noqa: BLE001
            print("retry", f, e, flush=True); time.sleep(2)
    raise RuntimeError(f)


def main():
    t0 = time.time()
    conds = {"control", "mglu_block", "post_nmdar", "nmdar_block", "no_block"}
    T = {k: v for k, v in load_targets(("paired_l5", "sjostrom07")).items() if k[1] in conds}
    B = BatchV2(DIRS, protocols=sorted({k[0].split("@")[0] for k in T}), pairs=set(open(os.path.join(V2, "subset24_pairs.txt")).read().strip().split(",")),
                fast=False, signals=("vdcc",))
    have = {r["proto"] for r in B.recs}; T = {k: v for k, v in T.items() if k[0].split("@")[0] in have}
    peak = np.concatenate([r["peak"] for r in B.recs]); B._stack()
    cfg = dict(model="v2", rho="full", filters=FIL, free_filters=FREE, filter_steps=8, loc=None, peak=peak,
               tie={"theta_NOi": "theta_Ti"}, continuous=True)
    F._B, F._T, F._CFG = B, T, cfg
    print(f"{len(B.recs)} records, {len(T)} targets, loaded in {time.time() - t0:.0f} s", flush=True)

    # ---- admissible candidates per population (why the v2 call grows with the generation)
    D = 6 + len(FREE); bounds = [(0.0, 5.0)] * 4 + [(-4.0, 0.0), (-1.0, 4.0)] + [(0.0, 1.0)] * len(FREE)
    lo, hi = np.array(bounds).T
    rng = np.random.default_rng(1); pop0 = lo + rng.random((8 * D, D)) * (hi - lo)
    a0 = json.load(open(os.path.join(V2, "..", "fit_results", "delta-cooker.json"))); a0 = a0.get("params", a0)
    pop0[0, :4] = [a0[k] for k in F.NAMES_A]
    pops = {"init (seed 1)": pop0}
    for f in sorted(glob.glob(os.path.join(RES, "reduced_gpu_subset_pl5sj07_td*_s*_ckpt.npz"))):
        z = load_ckpt(f); pops[f"{os.path.basename(f)} nit {int(z['nit'])}"] = z["x"]

    def admissible(pop):
        return [i for i, x in enumerate(pop) if not F.penalty(*F.unpack(x, cfg))]
    for k, pop in pops.items():
        print(f"admissible {len(admissible(pop)):3d} / {len(pop)}  {k}", flush=True)

    # ---- candidates: seed fit, best of each ckpt, then the admissible td4_s1 population
    fixed = [F.pack(json.load(open(os.path.join(RES, "reduced_gpu_subset_pl5r50_td3_s5.json")))["a"],
                    json.load(open(os.path.join(RES, "reduced_gpu_subset_pl5r50_td3_s5.json")))["pre"], cfg)]
    for f in sorted(glob.glob(os.path.join(RES, "reduced_gpu_subset_pl5sj07_td*_s*_ckpt.npz"))):
        fixed.append(load_ckpt(f)["best"])
    fixed = np.clip(np.array(fixed), lo, hi)
    popk = [k for k in pops if "td%d_s1" % TD in k][0]
    popx = pops[popk][admissible(pops[popk])]
    ok = [not F.penalty(*F.unpack(x, cfg)) for x in fixed]
    fixed = fixed[ok]
    print(f"{len(fixed)} fixed candidates (admissible), population sample {len(popx)} ({popk})", flush=True)
    AP = lambda xs: ([F.unpack(x, cfg)[0] for x in xs], [{**FIL, **F.unpack(x, cfg)[1]} for x in xs])

    # ---- v3
    t1 = time.time(); G3 = gpu_v3.GPUModelV3(B, T, FIL)
    res3 = {}
    for name, xs in (("fixed", fixed), ("pop", popx)):
        A, Ps = AP(xs)
        for rep in range(2):
            t = time.time(); rho, d = G3.rho_dpre(A, Ps); tg = time.time() - t
            t = time.time(); c, pr = G3.chi2(rho, d, return_pred=True); tc = time.time() - t
            print(f"v3 {name} P={len(xs)} rep {rep}: rho_dpre {tg:.2f} s (incl. first-call compile on rep 0), "
                  f"readout {tc:.2f} s", flush=True)
        res3[name] = (rho, d, c, pr)
    # full generation as fit_v3.objective_gpu (penalties + GPU + readout) on the whole ckpt population
    fit_v3._G, fit_v3._CFG = G3, cfg
    t = time.time(); o3 = fit_v3.objective_gpu(pops[popk].T); tgen = time.time() - t
    F._G = None
    print(f"v3 fit_v3.objective_gpu on {len(pops[popk])} members: {tgen:.2f} s", flush=True)

    # ---- v2 as the running fits (stream, CHUNK 4), compile vs run
    import jax, jax_v2
    jax.config.update("jax_log_compiles", True)
    G2 = jax_v2.GPUModel(B, T, FIL, chunk_gb=4)
    print("v2 buckets", [(b["T"], b["Lb"]) for b in G2.buckets], flush=True)
    res2 = {}
    for name, xs in (("fixed", fixed), ("pop", popx)):
        A, Ps = AP(xs)
        for rep in range(2):
            t = time.time(); rho, d = G2.rho_dpre(A, Ps); tg = time.time() - t
            t = time.time(); c, pr = G2.chi2(rho, d, return_pred=True); tc = time.time() - t
            print(f"v2 {name} P={len(xs)} rep {rep}: rho_dpre {tg:.2f} s (rep 0 compiles for a new P), "
                  f"readout {tc:.2f} s", flush=True)
        res2[name] = (rho, d, c, pr)
    A, Ps = AP(popx[:len(popx) - 1])                    # one candidate fewer: new shape
    t = time.time(); G2.rho_dpre(A, Ps); print(f"v2 P={len(popx) - 1} (new shape): {time.time() - t:.2f} s", flush=True)
    t = time.time(); th = [np.concatenate([B.thetas(r, a)[0] for r in B.recs]) for a in A]
    print(f"v2 host td/tp loop for P={len(A)}: {2 * (time.time() - t):.2f} s", flush=True)
    F._CFG = cfg; F._G = G2
    t = time.time(); o2 = F.objective_gpu(pops[popk].T); tgen2 = time.time() - t
    print(f"v2 fit_v2.objective_gpu on {len(pops[popk])} members: {tgen2:.2f} s", flush=True)

    # ---- equality
    worst = 0.0
    for name in ("fixed", "pop"):
        r2, d2, c2, p2 = res2[name]; r3, d3, c3, p3 = res3[name]
        n = G2.n_ctrl
        rel = np.abs(c3 - c2) / np.maximum(np.abs(c2), 1e-12)
        dp = max(np.nanmax(np.abs(p3[k] - p2[k])) for k in p2)
        print(f"{name}: chi2 v2 {np.round(c2[:6], 4)}\n{name}: chi2 v3 {np.round(c3[:6], 4)}", flush=True)
        print(f"{name}: max rel chi2 diff {rel.max():.2e}, max |pred diff| {dp:.2e}, "
              f"max |rho diff| {np.abs(r3 - r2).max():.2e}, max |dpre diff| {np.abs(d3 - d2).max():.2e} "
              f"(control lanes {np.abs(d3[:, :n] - d2[:, :n]).max():.2e}), "
              f"rho binary flips {int(((r3 >= 0.5) != (r2 >= 0.5)).sum())}", flush=True)
        worst = max(worst, rel.max())
    print(f"objective_gpu v2 vs v3 on the full population: max |diff| {np.abs(o3 - o2).max():.2e}", flush=True)
    assert worst < 1e-6, worst
    print(f"PASS v3 == v2 (chi2 rel {worst:.1e}); total {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
