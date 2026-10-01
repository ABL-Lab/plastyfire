"""jax_v2.GPUModel vs the CPU path for t_drive 4 (-ica_VDCC events weighted 1 - b_NMDA) on the pilot pairs (vca dirs), plus a
no_block (A_NO = 0) target to check its duplicate lanes. GPU node only.

    sbatch glusynapse_v2/run_jax_tdrive4_test.sh
"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); V2 = os.path.join(HERE, "..")
sys.path.insert(0, V2)
from batch_v2 import BatchV2
from targets import load_targets
import model_v2 as MV
import jax_v2 as JV

PAIRS = "180351-198084,181455-195199"          # pilot 22073044
DIRS = [os.path.join(V2, "extracted", d) for d in ("sj03_delta-prefire-vca", "sj03r50_delta-prefire-vca")]


def main():
    import jax; print("devices", jax.devices(), flush=True)
    fit = json.load(open(os.path.join(V2, "results/reduced_gpu_subset_pl5r50_v22_s2.json")))
    filters = json.loads(fit["args"]["filters"])
    sig = ("vdcc",) if filters.get("pre_drive") else ("shaft_cai",)
    T = load_targets(("paired_l5",))
    B = BatchV2(DIRS, protocols=sorted({k[0].split("@")[0] for k in T}), pairs=set(PAIRS.split(",")), signals=sig)
    have = {r["proto"] for r in B.recs}; T = {k: v for k, v in T.items() if k[0].split("@")[0] in have}
    T[(sorted(have)[0], "no_block")] = (1.0, 0.1, 1, "T30 no_block lanes (A_NO = 0), test only")
    assert all(r.get("cev") is not None for r in B.recs), "records without cev"
    print(f"{len(B.recs)} records, {len(T)} targets", flush=True)
    G = JV.GPUModel(B, T, dict(filters, t_drive=4))
    rng = np.random.default_rng(3)
    a0, p0 = fit["a"], {**filters, **fit["pre"]}
    cands = []
    for k in range(3):
        cands.append((a0, dict(p0, t_drive=4, theta_Te=rng.uniform(0.3, 1.5), tau_T=10 ** rng.uniform(1, 2.3),
                               theta_Tg=rng.uniform(0, 3), dpre_min=rng.uniform(-0.8, -0.1),
                               A_mglu=10 ** rng.uniform(-1.5, 0))))
    rho, d = G.rho_dpre([c[0] for c in cands], [c[1] for c in cands])
    worst = 0.0
    for i, (a, p) in enumerate(cands):
        P = {**MV.DEFAULTS, **p}
        dc = np.concatenate([MV.dpre_final(tT, K, P["A_mglu"], P["A_NO"], P["dpre_min"], P["dpre_max"], P["dpre0"])
                             for (tT, K), r in zip(B.features(P), B.recs)])
        ed = np.abs(d[i, :len(dc)] - dc).max(); worst = max(worst, ed)
        print(f"t_drive 4 cand {i}: max|dpre| {ed:.2e}  mean dpre {dc.mean():+.3f}  min {dc.min():+.3f}", flush=True)
    chi = G.chi2(rho, d); chic = [B.objective(a, p, T, mode="full")[0] for a, p in cands]
    print("t_drive 4 chi2 gpu", np.round(chi, 4), "\nt_drive 4 chi2 cpu", np.round(chic, 4), flush=True)
    assert worst < 1e-9 and np.allclose(chi, chic, rtol=1e-8, atol=1e-8), worst
    print("PASS t_drive 4 GPU == CPU")


if __name__ == "__main__":
    main()
