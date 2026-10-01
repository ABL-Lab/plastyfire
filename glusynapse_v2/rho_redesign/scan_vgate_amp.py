"""ROUND2.md C1 / C2 theta_V scan: potentiation gated on the synapse's own filtered spine VDCC current. CPU/numba, no jax,
no refit (A0g_s3 parameters). Scores the 29 L5 targets (paired_l5 + sjostrom07) and the 9 L2/3->L5 targets.

    V' = -V / tau_E1 + (-ica_VDCC) / i_scale          (tau_E1 100 ms, i_scale 1e-5: the t_drive 4 constants)
    pot = [effcai > theta_p] * [V > theta_V]           (theta_V <= 0: gate off = A0)
    C1: dep = [effcai > theta_d]                       (a crossing of theta_p without VDCC Ca drives depression)
    C2: dep = [effcai > theta_d] * [effcai <= theta_p] (the gated-out crossing is neutral)
    rho' = ... + pot gp (1 - rho) - dep (1 - pot) gd rho   (rho_v4.rho_loop_v4 with tau_fast 0, same skip rule)
Unlike candidate D (rho_v4d: unweighted K_ca event count G > 1/2 in a 69 ms window, i.e. a binary "any VDCC event"),
C1 thresholds the amplitude-weighted VDCC charge, so it can sit between the L2/3 and L5 levels.

Also: vmax = max V over the steps where effcai > theta_p (A0 crossings), per synapse, for the between-levels check.

    python glusynapse_v2/rho_redesign/scan_vgate_amp.py --fit <json> --save <prefix> --fig <png>
"""
import argparse, json, os, sys
import numpy as np, pandas as pd
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__)); V2 = os.path.dirname(HERE); ROOT = os.path.dirname(V2)
sys.path.insert(0, V2); sys.path.insert(0, HERE)
import batch_v2                          # noqa: E402
from batch_v2 import BatchV2, RHO_STAR_GB, TAU_IND_GB   # noqa: E402
from targets import load_targets        # noqa: E402
import model_v2 as MV                   # noqa: E402
import rho_v4                           # noqa: E402

THETAS = np.array([0.0, 2.0, 5.0, 10.0, 20.0, 35.0, 60.0, 100.0, 160.0, 250.0, 400.0, 700.0])
CANDS = ("C1", "C2")
L23_DIR = os.path.join(V2, "extracted", "ebner_l23l5_delta-prefire-vseg-rs")
L23_BASIS = os.path.join(ROOT, "basis_results_edges_ebner_l23l5_delta_rs")
L5_BASIS = os.path.join(ROOT, "basis_results_edges_sabrina_n120_delta")
GEOM = os.path.join(ROOT, "ebner", "pair_geometry_L23PC_L5TTPC.csv")
KEY = ["letzkus_1ap_dt+10ms|control", "letzkus_3ap_200hz_dt+10ms@distal|control", "letzkus_3ap_200hz_dt+10ms@proximal|control",
       "sjostrom_50hz_dt+10ms@distal|control", "sjostrom_50hz_dt+10ms|control", "10Hz_10ms|control", "10Hz_5ms|control",
       "10Hz_-10ms|control", "sjostrom_10hz_dt+10ms|control", "sjostrom_20hz_dt+10ms|control", "sjostrom_50hz_dt+10ms|control",
       "sjostrom_0.1hz_dt+10ms|control", "sjostrom_0.1hz_dt-10ms|control", "sjostrom07_step200ms_pair|control"]
L5_LTP = ("10Hz_10ms", "10Hz_5ms", "sjostrom_10hz_dt+10ms", "sjostrom_20hz_dt+10ms", "sjostrom_40hz_dt+10ms",
          "sjostrom_50hz_dt+10ms", "sjostrom07_step200ms_pair")


@njit(cache=True)
def rho_rec(E, VD, hs, aV, bV, td, tp, rho0, gd, gp, rs, ths, out, vmax):
    """out (2, n_theta, n_syn): rho_f for C1 / C2 at each theta_V; vmax (n_syn,)."""
    n, T = E.shape; nt = ths.shape[0]
    R = np.empty((2, nt))
    for i in range(n):
        for c in range(2):
            for q in range(nt):
                R[c, q] = rho0[i]
        a = td[i]; b = tp[i]; V = 0.0; vm = 0.0
        for k in range(T - 1):
            x = E[i, k]
            ab = x > b; aa = x > a
            if ab and V > vm:
                vm = V
            h = hs[k]
            for c in range(2):
                for q in range(nt):
                    r = R[c, q]
                    if not (aa or ab) and (r == 0.0 or r == 1.0):
                        continue
                    pot = 1.0 if (ab and (ths[q] <= 0.0 or V > ths[q])) else 0.0
                    if c == 0:
                        dep = 1.0 if aa else 0.0
                    else:
                        dep = 1.0 if (aa and not ab) else 0.0
                    r += h * (-r * (1 - r) * (rs - r) + pot * gp * (1 - r) - dep * (1 - pot) * gd * r)
                    R[c, q] = min(max(r, 0.0), 1.0)
            V = aV[k] * V + bV[k] * VD[i, k]
        for c in range(2):
            for q in range(nt):
                out[c, q, i] = R[c, q]
        vmax[i] = vm


def run(B, fit, P, T, conds, path, loc, sh):
    gamma = rho_v4.opts(P)[0]; gd, gp = rho_v4.rates(P); tau = float(P["tau_E1"]); isc = float(P["i_scale"])
    feats = B.features(P)
    for r in B.recs:
        r.pop("_c4", None)
    rows, syn = [], []
    for j, r in enumerate(B.recs):
        td, tp = rho_v4.thetas(r["c_pre"], r["c_post"], fit["a"], gamma)
        t = r["t"]; h = np.diff(t)
        hs = h / 1000.0 / TAU_IND_GB; aV = np.exp(-h / tau); bV = tau * (1.0 - aV) / isc
        n = len(r["syn"]); out = np.empty((2, len(THETAS), n)); vmax = np.empty(n)
        rho_rec(np.ascontiguousarray(r["effcai"]), np.ascontiguousarray(r["vdcc"], dtype=np.float64), hs, aV, bV,
                np.asarray(td, float), np.asarray(tp, float), r["rho0"].astype(float), gd, gp, float(RHO_STAR_GB), THETAS,
                out, vmax)
        b = B.basis(r); tT, K = feats[j]
        for c in conds:
            if c == "nmdar_block":
                d = np.zeros(n)
            else:
                d = MV.dpre_final(tT, K, 0.0 if c == "mglu_block" else P["A_mglu"],
                                  0.0 if (c == "no_block" or (c == "post_nmdar" and P["no_drive"])) else P["A_NO"],
                                  P["dpre_min"], P["dpre_max"], P["dpre0"])
            for ci, cn in enumerate(CANDS):
                for q, th in enumerate(THETAS):
                    rc = r["rho0"] if c in ("post_nmdar", "nmdar_block") else out[ci, q]
                    rows.append(dict(path=path, pair=r["pair"], proto=r["proto"], condition=c, cand=cn, theta_V=th,
                                     ratio=rho_v4.ratio(b, r["rho0"], rc, d)))
        up0 = (r["rho0"] < 0.5) & (out[0, 0] >= 0.5)
        for i in range(n):
            syn.append(dict(path=path, pair=r["pair"], proto=r["proto"], syn=int(r["syn"][i]), rho0=float(r["rho0"][i]),
                            rho_A0=float(out[0, 0, i]), up=bool(up0[i]), vmax=float(vmax[i]),
                            letzkus_distal=loc.get(r["pair"]), sh_distal=sh.get(r["pair"])))
    df = pd.DataFrame(rows)
    res = []
    for (pid, cond), (m, s, _n, _src) in sorted(T.items()):
        proto, _, where = pid.partition("@")
        sel = df[(df.proto == proto) & (df.condition == cond)].dropna(subset=["ratio"])
        if path == "L23" and where and proto.startswith("letzkus"):
            sel = sel[sel.pair.map(loc) == (where == "distal")]
        elif path == "L23" and where == "distal":
            sel = sel[sel.pair.map(sh).fillna(False).astype(bool)]
        elif where == "distal":
            continue
        if sel.empty:
            continue
        for (cn, th), x in sel.groupby(["cand", "theta_V"]):
            pred = x.ratio.mean()
            res.append(dict(path=path, target=f"{pid}|{cond}", cand=cn, theta_V=th, data=m, sem=s, pred=pred, z=(pred - m) / s))
    return pd.DataFrame(res), pd.DataFrame(syn)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True); ap.add_argument("--save", required=True); ap.add_argument("--fig", required=True)
    a_ = ap.parse_args()
    fit = json.load(open(a_.fit)); fa = fit["args"]
    batch_v2.BAP_GATE = fit.get("bap_gate")
    fil = {**json.loads(fa["filters"]), **json.loads(fa.get("set", "{}"))}
    P = {**MV.DEFAULTS, **fil, **fit["pre"]}
    for k in ("gamma_d", "gamma_p"):
        if k in fit:
            P[k] = fit[k]
    assert rho_v4.opts(P)[1] == 0.0 and rho_v4.opts(P)[2] == 0.0, "scan assumes tau_fast 0, rho_sigma 0"
    conds_ok = set(fa["conditions"].split(","))
    print("rule: rho_gamma", rho_v4.opts(P)[0], "rates", rho_v4.rates(P), "tau_E1", P["tau_E1"], "i_scale", P["i_scale"],
          "theta_V grid", THETAS.tolist(), flush=True)
    # L5
    T5 = {k: v for k, v in load_targets(tuple(fa["groups"].split(","))).items() if k[1] in conds_ok}
    batch_v2.BASIS_DIR = L5_BASIS
    B = BatchV2(fa["dirs"].split(","), protocols=sorted({k[0].split("@")[0] for k in T5}), pairs=set(fa["pairs"].split(",")),
                fast=False, signals=("vdcc",))
    for r in B.recs:
        B.basis(r)
    R5, S5 = run(B, fit, P, T5, sorted({c for _, c in T5}), "L5", {}, {})
    del B
    # L2/3 -> L5
    T23 = {k: v for k, v in load_targets(("paired_l23l5",)).items() if k[1] in conds_ok}
    g = pd.read_csv(GEOM); g["pair"] = g.pregid.astype(str) + "-" + g.postgid.astype(str); g = g[g.all_protocols]
    loc = dict(zip(g.pair, g.letzkus_distal.astype(bool))); sh = dict(zip(g.pair, g.sh_distal.astype(bool)))
    batch_v2.BASIS_DIR = L23_BASIS
    B = BatchV2([L23_DIR], protocols=sorted({k[0].split("@")[0] for k in T23}), pairs=set(g.pair), fast=False,
                signals=("vdcc",))
    for r in B.recs:
        B.basis(r)
    R23, S23 = run(B, fit, P, T23, sorted({c for _, c in T23}), "L23", loc, sh)
    del B
    R = pd.concat([R5, R23], ignore_index=True); S = pd.concat([S5, S23], ignore_index=True)
    R.to_csv(a_.save + "_targets.csv", index=False); S.to_csv(a_.save + "_syn.csv", index=False)
    # summary: chi2 per path and key targets per (cand, theta_V)
    summ = []
    for (cn, th), x in R.groupby(["cand", "theta_V"]):
        q = dict(cand=cn, theta_V=th, chi2_L5=float((x[x.path == "L5"].z ** 2).sum()), n_L5=int((x.path == "L5").sum()),
                 chi2_L23=float((x[x.path == "L23"].z ** 2).sum()), n_L23=int((x.path == "L23").sum()))
        for k in KEY:
            for pth in ("L5", "L23"):
                y = x[(x.target == k) & (x.path == pth)]
                if len(y):
                    q[f"{pth}:{k.split('|')[0]}"] = float(y.pred.iloc[0])
        summ.append(q)
    M = pd.DataFrame(summ); M.to_csv(a_.save + ".csv", index=False)
    pd.set_option("display.width", 300); pd.set_option("display.max_columns", 60)
    print("=== chi2 per candidate and theta_V (theta_V 0 must reproduce A0g_s3: L5 34.16, L2/3 220.4)")
    print(M[["cand", "theta_V", "chi2_L5", "n_L5", "chi2_L23", "n_L23"]].to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    print("=== key targets (pred)")
    print(M.drop(columns=["chi2_L5", "n_L5", "chi2_L23", "n_L23"]).to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    # between-levels check: vmax (V at A0 theta_p crossings) of up-flipping synapses
    grp = {"L5 LTP protos, up": S[(S.path == "L5") & S.proto.isin(L5_LTP) & S.up],
           "L5 all, crossing": S[(S.path == "L5") & (S.vmax > 0)],
           "L2/3 1AP+10, up": S[(S.path == "L23") & (S.proto == "letzkus_1ap_dt+10ms") & S.up],
           "L2/3 3AP+10 distal, up": S[(S.path == "L23") & (S.proto == "letzkus_3ap_200hz_dt+10ms") & (S.letzkus_distal == True) & S.up],  # noqa: E712
           "L2/3 3AP+10 prox, up": S[(S.path == "L23") & (S.proto == "letzkus_3ap_200hz_dt+10ms") & (S.letzkus_distal == False) & S.up],  # noqa: E712
           "L2/3 SH50 distal, up": S[(S.path == "L23") & (S.proto == "sjostrom_50hz_dt+10ms") & (S.sh_distal == True) & S.up]}  # noqa: E712
    qs = [0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95]
    Qv = pd.DataFrame([dict(group=k, n=len(v), **{f"q{int(100 * p)}": float(v.vmax.quantile(p)) if len(v) else np.nan for p in qs})
                       for k, v in grp.items()])
    Qv.to_csv(a_.save + "_vmax.csv", index=False)
    print("=== V at A0 theta_p crossings (vmax) quantiles: can one theta_V sit above L2/3 and below L5?")
    print(Qv.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    fig_(M, grp, a_.fig)
    print(f"saved {a_.save}.csv, _targets.csv, _syn.csv, _vmax.csv, {a_.fig}")


def fig_(M, grp, path):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(17, 5))
    xs = lambda th: np.where(th == 0, 1.0, th)
    for cn, ls in (("C1", "-"), ("C2", "--")):
        x = M[M.cand == cn].sort_values("theta_V")
        ax[0].plot(xs(x.theta_V), x.chi2_L5, ls, color="#2e86c1", marker="o", ms=3, label=f"{cn} L5 (29)")
        ax[0].plot(xs(x.theta_V), x.chi2_L23, ls, color="#c0392b", marker="o", ms=3, label=f"{cn} L2/3 (9)")
    ax[0].set(xscale="log", yscale="log", xlabel=r"$\theta_{\mathrm{VDCC}}$, spine VDCC Ca gate (0 plotted at 1)", ylabel=r"$\chi^2$", title=r"a  $\chi^2$ vs $\theta_{\mathrm{VDCC}}$ (no refit)")
    ax[0].legend(fontsize=7)
    cols = {"L23:letzkus_1ap_dt+10ms": ("#c0392b", 0.72), "L23:letzkus_3ap_200hz_dt+10ms@distal": ("#e67e22", 0.79),
            "L23:letzkus_3ap_200hz_dt+10ms@proximal": ("#8e44ad", 1.30), "L23:sjostrom_50hz_dt+10ms@distal": ("#d35400", 0.86),
            "L5:10Hz_10ms": ("#2e86c1", 1.2013), "L5:sjostrom_10hz_dt+10ms": ("#1abc9c", 1.16),
            "L5:sjostrom_20hz_dt+10ms": ("#16a085", 1.31), "L5:sjostrom_50hz_dt+10ms": ("#34495e", 1.57)}
    x = M[M.cand == "C1"].sort_values("theta_V")
    for k, (c, d) in cols.items():
        if k in x:
            ax[1].plot(xs(x.theta_V), x[k], "-", color=c, marker="o", ms=3, label=k.replace("L23:", "L2/3 ").replace("L5:", "L5 "))
            ax[1].axhline(d, color=c, lw=0.6, ls=":")
    ax[1].set(xscale="log", xlabel=r"$\theta_{\mathrm{VDCC}}$", ylabel="ratio", title="b  C1 key targets (dotted = data)"); ax[1].legend(fontsize=6)
    for k, v in grp.items():
        if len(v):
            vv = np.sort(np.maximum(v.vmax.values, 1e-2))
            ax[2].plot(vv, np.linspace(0, 1, len(vv)), label=f"{k} ({len(v)})", ls="-" if k.startswith("L5") else "--")
    for th in M.theta_V.unique():
        if th > 0:
            ax[2].axvline(th, color="k", lw=0.3)
    ax[2].set(xscale="log", xlabel=r"spine VDCC Ca $c_{\mathrm{VDCC}}$ at $\theta_p$ crossings", ylabel="CDF", title="c  own spine VDCC Ca at A0 crossings")
    ax[2].legend(fontsize=6)
    fig.tight_layout(); fig.savefig(path, dpi=130)


if __name__ == "__main__":
    main()
