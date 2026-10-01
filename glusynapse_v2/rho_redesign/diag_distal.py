"""Distal Letzkus diagnosis (DISTAL_LTP.md step 1). CPU/numba, no jax. Fixed fit (A0 joint s4), no refit.

Question: which synapse-local signal separates distal 3AP 200 Hz dt -10 (data LTP 1.42) from distal 3AP +10 (0.79) and
1AP +10 (0.72)? Proximal is the control (3AP +10 LTP 1.30, -10 LTD 0.89): a usable signal should favour -10 at distal
and +10 at proximal synapses.

Per synapse and pairing (window [arr - 30, arr + 150] ms around the synapse's own pre arrival arr = prespike + delay),
averaged over pairings:
  eff_pk, m_p, m_d     peak effcai, peak / theta_p, peak / theta_d
  ca_pk, ca_int        spine Ca above rest (cacr, recovered from effcai) peak and integral (mM ms)
  ca_bef, ca_aft       cacr integral in [arr - 30, arr) and [arr, arr + 50)    (pre-locked spine Ca)
  ca_late              cacr integral in [arr + 10, arr + 60)                     (NMDA-window spine Ca)
  sh_pk, sh_int        shaft Ca above the window-start value: peak, integral
  vd_pk, vd_int        -ica_VDCC (spine VDCC current, a local-V proxy) peak, integral
  vd_aft               -ica_VDCC integral in [arr, arr + 50)
  n_cev, t_cev         VDCC Ca events (cev, K_ca crossings) per pairing, mean time rel. arr
  n_cev_aft            events in [arr, arr + 50)
  b_cev, bsum_cev      own glutamate-bound NMDA state b at the events (model_v2.glu_weight: w = 1 - b), mean and sum
  n_vev, t_vev         local depolarisation events (vev, t_drive 3 crossings of v_seg - vbar > theta_V) per pairing
Raw local V is not stored in the npz (extract_v2 keeps only vev); vdcc and vev are the local-V proxies.
Also rho0, rho_f (A0 rule), theta_d, theta_p, c_pre, c_post, and the record ratio (sanity check vs the fit's _l23.csv).

    ANALYTICAL_BASIS_DIR=basis_results_edges_ebner_l23l5_delta_rs python glusynapse_v2/rho_redesign/diag_distal.py \
        --fit <json> --dirs <dir> --save <prefix> --fig <png>
"""
import argparse, json, os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); V2 = os.path.dirname(HERE); ROOT = os.path.dirname(V2)
sys.path.insert(0, V2); sys.path.insert(0, HERE)
import batch_v2                                   # noqa: E402
from batch_v2 import BatchV2, cacr_from_effcai    # noqa: E402
import model_v2 as MV                             # noqa: E402
import rho_v4, rho_v4d                            # noqa: E402

PROTOS = ("letzkus_3ap_200hz_dt-10ms", "letzkus_3ap_200hz_dt+10ms", "letzkus_1ap_dt+10ms")
SHORT = {PROTOS[0]: "3AP-10", PROTOS[1]: "3AP+10", PROTOS[2]: "1AP+10"}
W0, W1 = -30.0, 150.0
FEATS = ("eff_pk", "m_p", "m_d", "ca_pk", "ca_int", "ca_bef", "ca_aft", "ca_late", "sh_pk", "sh_int", "vd_pk",
         "vd_int", "vd_aft", "n_cev", "t_cev", "n_cev_aft", "b_cev", "bsum_cev", "n_vev", "t_vev")


def integ(t, y, lo, hi):
    m = (t >= lo) & (t < hi)
    if m.sum() < 2:
        return 0.0
    tt, yy = t[m], y[m]
    return float(np.sum(yy[:-1] * np.diff(tt)))


def syn_features(r, i, arr, cacr, td, tp, tau_d):
    t = r["t"]; E = r["effcai"][i]; sh = r["shaft_cai"][i] if r.get("shaft_cai") is not None else None
    vd = r["vdcc"][i] if r.get("vdcc") is not None else None
    cev = r["cev"][i] if r.get("cev") is not None else np.zeros(0)
    cev = np.asarray(cev, float); cev = cev[np.isfinite(cev)]
    vev = r["vev"][i] if r.get("vev") is not None else np.zeros(0)
    vev = np.asarray(vev, float); vev = vev[np.isfinite(vev)]
    ev_all, w_all = MV.glu_weight(cev, arr, tau_d) if len(cev) else (np.zeros(0), np.zeros(0))
    rows = []
    for a in arr:
        m = (t >= a + W0) & (t < a + W1)
        if m.sum() < 2:
            continue
        f = dict(eff_pk=float(E[m].max()))
        f["m_p"] = f["eff_pk"] / tp if np.isfinite(tp) and tp > 0 else np.nan
        f["m_d"] = f["eff_pk"] / td if np.isfinite(td) and td > 0 else np.nan
        c = cacr[i]
        f["ca_pk"] = float(c[m].max()); f["ca_int"] = integ(t, c, a + W0, a + W1)
        f["ca_bef"] = integ(t, c, a + W0, a); f["ca_aft"] = integ(t, c, a, a + 50); f["ca_late"] = integ(t, c, a + 10, a + 60)
        if sh is not None:
            s0 = float(sh[m][0]); f["sh_pk"] = float(sh[m].max() - s0); f["sh_int"] = integ(t, sh - s0, a + W0, a + W1)
        if vd is not None:
            f["vd_pk"] = float(vd[m].max()); f["vd_int"] = integ(t, vd, a + W0, a + W1); f["vd_aft"] = integ(t, vd, a, a + 50)
        k = (ev_all >= a + W0) & (ev_all < a + W1)
        f["n_cev"] = int(k.sum()); f["t_cev"] = float(np.mean(ev_all[k] - a)) if k.any() else np.nan
        f["n_cev_aft"] = int(((ev_all >= a) & (ev_all < a + 50)).sum())
        b = 1.0 - w_all[k]
        f["b_cev"] = float(b.mean()) if k.any() else np.nan; f["bsum_cev"] = float(b.sum())
        kv = (vev >= a + W0) & (vev < a + W1)
        f["n_vev"] = int(kv.sum()); f["t_vev"] = float(np.mean(vev[kv] - a)) if kv.any() else np.nan
        rows.append(f)
    if not rows:
        return {}
    D = pd.DataFrame(rows)
    return {k: float(D[k].mean()) for k in D.columns}


def auc(x, y):
    """P(x > y) + 0.5 P(x == y), NaN-safe."""
    x = np.asarray(x, float); y = np.asarray(y, float); x = x[np.isfinite(x)]; y = y[np.isfinite(y)]
    if not len(x) or not len(y):
        return np.nan
    return float(((x[:, None] > y[None, :]).mean() + 0.5 * (x[:, None] == y[None, :]).mean()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True); ap.add_argument("--dirs", required=True)
    ap.add_argument("--save", required=True); ap.add_argument("--fig", required=True)
    a_ = ap.parse_args()
    fit = json.load(open(a_.fit)); fa = fit["args"]
    batch_v2.BAP_GATE = fit.get("bap_gate")
    fil = {**json.loads(fa["filters"]), **json.loads(fa.get("set", "{}"))}
    P = {**MV.DEFAULTS, **fil, **fit["pre"]}
    for k in ("gamma_d", "gamma_p"):
        if k in fit:
            P[k] = fit[k]
    gamma = rho_v4.opts(P)[0]; tau_d = float(P["tau_d_NMDA"])
    print("rule: rho_gamma", gamma, "rates", rho_v4.rates(P), "tau_d_NMDA", tau_d, flush=True)
    g = pd.read_csv(os.path.join(ROOT, "ebner", "pair_geometry_L23PC_L5TTPC.csv"))
    g["pair"] = g.pregid.astype(str) + "-" + g.postgid.astype(str); g = g[g.all_protocols]
    dist = dict(zip(g.pair, g.letzkus_distal.astype(bool)))
    B = BatchV2(a_.dirs.split(","), protocols=list(PROTOS), pairs=set(g.pair), fast=False, signals=("vdcc", "shaft_cai"))
    # npz keys of one record (no login-node python)
    z = np.load(os.path.join(a_.dirs.split(",")[0], f"{B.recs[0]['pair']}__{B.recs[0]['proto']}.npz"))
    print("npz keys:", {k: z[k].shape for k in z.files}, flush=True)
    rho = rho_v4d.rho_all(B, fit["a"], P)
    feats = B.features(P)
    rows, recs, ex = [], [], {}
    for j, r in enumerate(B.recs):
        td, tp = rho_v4.thetas(r["c_pre"], r["c_post"], fit["a"], gamma)
        rf = rho[r["sl"]]; cacr = cacr_from_effcai(r["effcai"], r["t"])
        delay = MV.edge_params(r["syn"])["delay"]
        pre = np.asarray(r["prespikes"], float); post = np.asarray(r["postspikes"], float)
        # sign check: pre time minus the first post spike of the same pairing (nearest post onset)
        dts = []
        for p in pre:
            if len(post):
                q = post[np.argmin(np.abs(post - p))]; on = post[(post <= q) & (post > q - 20)].min()
                dts.append(p - on)
        tT, K = feats[j]
        d = MV.dpre_final(tT, K, P["A_mglu"], P["A_NO"], P["dpre_min"], P["dpre_max"], P["dpre0"])
        b = B.basis(r)
        recs.append(dict(pair=r["pair"], proto=r["proto"], distal=dist.get(r["pair"]), n_syn=len(rf),
                         ratio=rho_v4.ratio(b, r["rho0"], rf, d), dt_pre_minus_post=float(np.median(dts)) if dts else np.nan,
                         n_pairings=len(pre)))
        for i in range(len(rf)):
            arr = pre + delay[i]
            f = syn_features(r, i, arr, cacr, td[i], tp[i], tau_d)
            rows.append(dict(pair=r["pair"], proto=r["proto"], P=SHORT[r["proto"]], distal=dist.get(r["pair"]),
                             syn=int(r["syn"][i]), c_pre=float(r["c_pre"][i]), c_post=float(r["c_post"][i]),
                             td=float(td[i]), tp=float(tp[i]), rho0=float(r["rho0"][i]), rho_f=float(rf[i]), **f))
        # example traces: first distal pair that has all protocols, one pairing in the middle of the train
        if dist.get(r["pair"]) and len(pre) and len(ex) < 6 and not any(len(v) == 3 for v in ex.values()):
            ex.setdefault(r["pair"], {})[SHORT[r["proto"]]] = dict(
                t=r["t"], cacr=cacr, E=r["effcai"], vd=r.get("vdcc"), tp=tp, td=td, pre=pre, post=post, delay=delay,
                cev=r.get("cev"), syn=r["syn"])
    S = pd.DataFrame(rows); R = pd.DataFrame(recs)
    S.to_csv(a_.save + "_syn.csv", index=False); R.to_csv(a_.save + "_rec.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    print("=== records: mean ratio (check vs v4_A0_joint_s4_l23.csv), sign check dt = pre - post onset (ms)")
    print(R.groupby(["proto", "distal"]).agg(n=("ratio", "size"), ratio=("ratio", "mean"),
                                             dt=("dt_pre_minus_post", "median"), npair=("n_pairings", "median")).round(3))
    S["floor"] = S.c_post < 1e-3
    S["up"] = (S.rho0 < 0.5) & (S.rho_f >= 0.5); S["down"] = (S.rho0 >= 0.5) & (S.rho_f < 0.5)
    print("=== per protocol x location x floor: medians")
    cols = ["m_p", "m_d", "ca_int", "ca_aft", "ca_late", "vd_int", "vd_aft", "n_cev", "n_cev_aft", "b_cev", "bsum_cev",
            "n_vev", "t_cev"]
    gq = S.groupby(["distal", "floor", "P"])
    T = gq[cols].median()
    T.insert(0, "n", gq.size()); T["P_up"] = gq.apply(lambda x: x.up.sum() / max((x.rho0 < 0.5).sum(), 1))
    T["P_down"] = gq.apply(lambda x: x.down.sum() / max((x.rho0 >= 0.5).sum(), 1))
    print(T.to_string(float_format=lambda v: f"{v:.3g}"))
    # separation: paired (same pair, synapse) log ratio and AUC; want distal 3AP-10 > 3AP+10 and > 1AP+10,
    # proximal 3AP+10 > 3AP-10
    sep = []
    piv = {p: S[S.P == p].set_index(["pair", "syn"]) for p in SHORT.values()}
    for loc in (True, False):
        for fl in ("all", True, False):
            for f in FEATS:
                if f not in S:
                    continue
                def get(p):
                    x = piv[p]; x = x[x.distal == loc]
                    return x if fl == "all" else x[x.floor == fl]
                a, b, c = get("3AP-10"), get("3AP+10"), get("1AP+10")
                ab = a[[f]].join(b[[f]], rsuffix="_b", how="inner").dropna()
                lr = np.log((ab[f].abs() + 1e-12) / (ab[f + "_b"].abs() + 1e-12)) if len(ab) else np.array([])
                sep.append(dict(distal=loc, floor=fl, feat=f, n=len(ab), med_m10=float(a[f].median()),
                                med_p10=float(b[f].median()), med_1ap=float(c[f].median()),
                                frac_m10_gt_p10=float((ab[f] > ab[f + "_b"]).mean()) if len(ab) else np.nan,
                                med_logratio=float(np.median(lr)) if len(lr) else np.nan,
                                auc_m10_vs_p10=auc(a[f], b[f]), auc_m10_vs_1ap=auc(a[f], c[f])))
    Q = pd.DataFrame(sep); Q.to_csv(a_.save + "_sep.csv", index=False)
    # a feature separates if distal favours -10 (AUC > 0.5 vs both +10) and proximal favours +10 (AUC < 0.5)
    qa = Q[Q.floor == "all"].pivot(index="feat", columns="distal", values=["auc_m10_vs_p10", "auc_m10_vs_1ap", "frac_m10_gt_p10"])
    print("=== separation (all synapses): AUC(-10 > +10), AUC(-10 > 1AP) for distal True / proximal False")
    print(qa.round(3).to_string())
    print("=== separation (distal only) by floor (c_post < 1e-3)")
    print(Q[Q.distal == True].pivot(index="feat", columns="floor", values="auc_m10_vs_p10").round(3).to_string())  # noqa: E712
    fig_(S, Q, ex, a_.fig)
    print(f"saved {a_.save}_syn.csv / _rec.csv / _sep.csv, {a_.fig}  ({len(S)} synapse rows, {len(R)} records)")


def fig_(S, Q, ex, path):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    col = {"3AP-10": "#c0392b", "3AP+10": "#2e86c1", "1AP+10": "#7f8c8d"}
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))
    # (a, b) example traces: distal pair with all three protocols, synapse with the largest -10 cacr
    pair = next((p for p, d in ex.items() if len(d) == 3), next(iter(ex), None))
    if pair:
        d = ex[pair]; i = int(np.argmax([np.max(d["3AP-10"]["cacr"][k]) for k in range(len(d["3AP-10"]["syn"]))])) \
            if "3AP-10" in d else 0
        for p, e in d.items():
            k = len(e["pre"]) // 2; a = e["pre"][k] + e["delay"][i]; t = e["t"]; m = (t >= a - 40) & (t < a + 120)
            ax[0, 0].plot(t[m] - a, 1e3 * e["cacr"][i][m], color=col[p], label=p)
            ax[0, 1].plot(t[m] - a, e["E"][i][m] * 1e3, color=col[p], label=p)
            ax[0, 1].axhline(e["tp"][i] * 1e3, color=col[p], ls="--", lw=0.8)
            if e["vd"] is not None:
                ax[0, 2].plot(t[m] - a, e["vd"][i][m], color=col[p], label=p)
            for q in e["post"][(e["post"] > a - 40) & (e["post"] < a + 120)]:
                ax[0, 0].axvline(q - a, color=col[p], lw=0.6, alpha=0.6)
        ax[0, 0].set(title=f"distal {pair} syn {i}: spine Ca above rest (uM)", xlabel="ms from own pre arrival")
        ax[0, 1].set(title="effcai (uM), dashed = theta_p", xlabel="ms from own pre arrival")
        ax[0, 2].set(title="-ica_VDCC (spine VDCC current)", xlabel="ms from own pre arrival")
        ax[0, 0].legend(fontsize=8)
    # (c) distal vs proximal AUC(-10 > +10) per feature
    q = Q[Q.floor == "all"].pivot(index="feat", columns="distal", values="auc_m10_vs_p10").dropna()
    y = np.arange(len(q))
    ax[1, 0].barh(y - 0.2, q.get(True, np.nan), 0.4, color="#c0392b", label="distal (want > 0.5)")
    ax[1, 0].barh(y + 0.2, q.get(False, np.nan), 0.4, color="#2e86c1", label="proximal (want < 0.5)")
    ax[1, 0].axvline(0.5, color="k", lw=0.8); ax[1, 0].set_yticks(y, q.index, fontsize=8)
    ax[1, 0].set(title="AUC(3AP -10 > 3AP +10)", xlim=(0, 1)); ax[1, 0].legend(fontsize=7, loc="lower right")
    # (d, e) distal per-protocol distributions of two key features
    D = S[S.distal == True]  # noqa: E712
    for axx, f, lab in ((ax[1, 1], "ca_late", "spine Ca integral [arr+10, arr+60) (mM ms)"),
                        (ax[1, 2], "m_p", "pairing effcai peak / theta_p")):
        data = [D[D.P == p][f].dropna().values for p in col]
        axx.boxplot(data, showfliers=False); axx.set_xticks([1, 2, 3], list(col))
        for k, p in enumerate(col):
            v = D[D.P == p][f].dropna().values
            axx.scatter(np.full(len(v), k + 1) + np.random.uniform(-0.15, 0.15, len(v)), v, s=4, color=col[p], alpha=0.4)
        axx.set(title=f"distal: {lab}")
    ax[1, 2].axhline(1.0, color="k", lw=0.8)
    fig.tight_layout(); fig.savefig(path, dpi=130)


if __name__ == "__main__":
    main()
