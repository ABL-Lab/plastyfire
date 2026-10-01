"""CPU scoring (numba, no jax, no GPU) of a fit json (v3 or v4) with the v4 post rule: the L2/3 -> L5 transfer test
(eval_l23l5.py selections) or the L5 targets. The v4 options and the rho rates gamma_d, gamma_p (fit_v4 --fit-gamma)
come from the json (filters + pre); absent = v3 / GAMMA_D_GB, GAMMA_P_GB. ROUND2 C1 / C2 (vamp_mode 1 / 2, theta_V) use
the CPU ground truth scan_vgate_amp.rho_rec.

    ANALYTICAL_BASIS_DIR=basis_results_edges_ebner_l23l5_delta_rs python glusynapse_v2/rho_redesign/eval_v4.py --l23 \
        --fit <json> --dirs glusynapse_v2/extracted/ebner_l23l5_delta-prefire-vseg-rs --save <prefix>
"""
import argparse, json, os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); V2 = os.path.dirname(HERE); ROOT = os.path.dirname(V2)
sys.path.insert(0, V2); sys.path.insert(0, HERE)
import batch_v2                          # noqa: E402
from batch_v2 import BatchV2            # noqa: E402
from targets import load_targets        # noqa: E402
import model_v2 as MV                   # noqa: E402
import rho_v4                           # noqa: E402
import rho_v4d                          # noqa: E402

GEOM = os.path.join(ROOT, "ebner", "pair_geometry_L23PC_L5TTPC.csv")


def rho_vamp(B, a, P, mode, theta_V):
    """ROUND2 C1 / C2 rho_f per record (scan_vgate_amp.rho_rec, the CPU ground truth; tau_fast 0, no vgate)."""
    import scan_vgate_amp as SV
    gamma, tau_fast, _ = rho_v4.opts(P)
    assert tau_fast == 0.0 and not P.get("vgate"), "vamp: tau_fast 0, no vgate"
    gd, gp = rho_v4.rates(P); tau = float(P["tau_E1"]); isc = float(P["i_scale"])
    out_ = []
    for r in B.recs:
        td, tp = rho_v4.thetas(r["c_pre"], r["c_post"], a, gamma)
        h = np.diff(r["t"]); aV = np.exp(-h / tau); bV = tau * (1.0 - aV) / isc
        n = len(r["syn"]); o = np.empty((2, 1, n)); vm = np.empty(n)
        SV.rho_rec(np.ascontiguousarray(r["effcai"]), np.ascontiguousarray(r["vdcc"], dtype=np.float64),
                   h / 1000.0 / batch_v2.TAU_IND_GB, aV, bV, np.asarray(td, float), np.asarray(tp, float),
                   r["rho0"].astype(float), gd, gp, float(batch_v2.RHO_STAR_GB), np.array([theta_V]), o, vm)
        out_.append(o[mode - 1, 0])
    return out_


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True); ap.add_argument("--dirs", required=True); ap.add_argument("--save", required=True)
    ap.add_argument("--l23", action="store_true"); ap.add_argument("--groups", default="paired_l5,sjostrom07")
    a_ = ap.parse_args()
    fit = json.load(open(a_.fit)); fa = fit["args"]
    batch_v2.BAP_GATE = fit.get("bap_gate")
    fil = {**json.loads(fa["filters"]), **json.loads(fa.get("set", "{}"))}
    P = {**MV.DEFAULTS, **fil, **fit["pre"]}
    gamma, tau_fast, sigma = rho_v4.opts(P)
    gd, gp = rho_v4.rates(P)
    T = load_targets(("paired_l23l5",) if a_.l23 else tuple(a_.groups.split(",")))
    T = {k: v for k, v in T.items() if k[1] in set(fa["conditions"].split(","))}
    pairs, dist, sh = None, {}, {}
    if a_.l23:
        g = pd.read_csv(GEOM); g["pair"] = g.pregid.astype(str) + "-" + g.postgid.astype(str); g = g[g.all_protocols]
        pairs = set(g.pair); dist = dict(zip(g.pair, g.letzkus_distal)); sh = dict(zip(g.pair, g.sh_distal))
    B = BatchV2(a_.dirs.split(","), protocols=sorted({k[0].split("@")[0] for k in T}), pairs=pairs, fast=False,
                signals=("vdcc",) if fil.get("pre_drive") else ("shaft_cai",))
    vmode, theta_V = rho_v4.vamp(P)
    if vmode:
        rhos = rho_vamp(B, fit["a"], P, vmode, theta_V)
    else:
        rho = rho_v4d.rho_all(B, fit["a"], P)
        rhos = [rho[r["sl"]] for r in B.recs]
    feats = B.features(P)
    rows = []
    conds = sorted({c for _, c in T})
    for i, r in enumerate(B.recs):
        b = B.basis(r); rf = rhos[i]; tT, K = feats[i]
        for c in conds:
            rc = r["rho0"] if c in ("post_nmdar", "nmdar_block") else rf
            if c == "nmdar_block":
                d = np.zeros(len(rc))
            else:
                d = MV.dpre_final(tT, K, 0.0 if c == "mglu_block" else P["A_mglu"],
                                  0.0 if (c == "no_block" or (c == "post_nmdar" and P["no_drive"])) else P["A_NO"],
                                  P["dpre_min"], P["dpre_max"], P["dpre0"])
            rows.append(dict(pair=r["pair"], proto=r["proto"], condition=c,
                             ratio=rho_v4.ratio(b, r["rho0"], rc, d, 0.0 if c in ("post_nmdar", "nmdar_block") else sigma),
                             dpre=float(np.mean(d)), rho=float(np.mean(rc))))
    df = pd.DataFrame(rows); df.to_csv(a_.save + "_records.csv", index=False)
    out = []
    for (pid, cond), (m, s, n, src) in sorted(T.items()):
        proto, _, where = pid.partition("@")
        sel = df[(df.proto == proto) & (df.condition == cond)].dropna(subset=["ratio"])
        if a_.l23 and where and proto.startswith("letzkus"):
            sel = sel[sel.pair.map(dist) == (where == "distal")]
        elif a_.l23 and where == "distal":
            sel = sel[sel.pair.map(sh).fillna(False).astype(bool)]
        elif where == "distal":
            continue
        if sel.empty:
            continue
        pred = sel.ratio.mean()
        out.append(dict(target=pid, condition=cond, target_mean=m, target_sem=s, pred=pred, pred_sem=sel.ratio.sem(),
                        n_pairs=len(sel), z=(pred - m) / s, dpre=sel.dpre.mean(), rho=sel.rho.mean()))
    o = pd.DataFrame(out); o.to_csv(a_.save + ".csv", index=False)
    chi2 = float((o.z ** 2).sum())
    json.dump(dict(fit=a_.fit, chi2=chi2, n_targets=len(o), v4=dict(rho_gamma=gamma, tau_fast=tau_fast, rho_sigma=sigma, vgate=int(P.get("vgate", 0)), vamp_mode=vmode, theta_V=theta_V,
                                                                gamma_d=gd, gamma_p=gp),
                   l23=a_.l23), open(a_.save + ".json", "w"), indent=1)
    print(o.round(3).to_string(index=False))
    print(f"chi2 {chi2:.2f} over {len(o)} targets (v4 gamma {gamma} tau_fast {tau_fast} sigma {sigma}; "
          f"gamma_d {gd:.3f} gamma_p {gp:.3f})")


if __name__ == "__main__":
    main()
