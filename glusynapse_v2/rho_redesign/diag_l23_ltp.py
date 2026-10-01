"""Round-2 diagnosis (ROUND2.md): why A0g_s3 over-potentiates L2/3->L5 pre-before-post (+10) while L5->L5 +10 must
stay LTP. CPU/numba, no jax. Fixed fit, no refit.

Targets (data):  L2/3->L5  Letzkus 1AP +10 all 0.72, 3AP 200 Hz +10 distal 0.79 (proximal 1.30 = LTP control),
                           S&H 50 Hz +10 distal 0.86, all 1.06
                 L5->L5    Markram 10 Hz +10 1.20, Sjostrom 10/20/50 Hz +10 1.16/1.31/1.57, 0.1 Hz +10 0.97
1. Ratio decomposition per target (record means, as eval_v4):
     full; rho_only (dpre 0); pre_only (rho = rho0); no_ecb (A_mglu 0); no_no (A_NO 0); pot_off (rho_f <= rho0);
     dep_off (rho_f >= rho0); ecb_unw (t_drive 4 events unweighted, i.e. 1 - b -> 1: is eCB-LTD missing because of the
     glutamate weight?).
2. Per-synapse local signals (mean over pairings; pairing = pre spikes closer than 300 ms; window [first pre - 30,
   last pre + 150] ms):
     m_p, m_d      pairing effcai peak / theta_p, / theta_d (A0: thresholds = a c_pre + const, gamma 0)
     ca_pk, ca_int spine Ca above rest (cacr) peak and integral (mM ms)
     sh_pk, vd_int shaft Ca peak above window start; -ica_VDCC integral
     n_cev, b_cev, w_cev   VDCC Ca events (cev, K_ca) per pairing, mean NMDA-bound b at them, sum of 1 - b
     n_pre         own pre spikes per pairing; n_pair pairings; f_pair pairing rate (Hz)
     t_p, t_d      s per pairing above theta_p, in the depression band (theta_d < effcai <= theta_p)
     pot_frac      gamma_p t_p / (gamma_p t_p + gamma_d t_d)  (local drive balance)
     tT_max, tTu_max   max eCB gate tanh(pos(T - theta_Tg)) at own arrivals, weighted / unweighted events
     c_pre, c_post, theta_d, theta_p, rho0, rho_f, sec_type (syn_section_type.npz, else edges afferent_section_type; here L5 = 2, L2/3 = 3),
     dist (L5: local_t path distance um), path_mean (L2/3: pair mean path distance um)
   AUC(L2/3 should-depress > L5 should-potentiate) per feature, all synapses and per rho0 state.

    python glusynapse_v2/rho_redesign/diag_l23_ltp.py --fit <json> --save <prefix> --fig <png>
"""
import argparse, glob, json, os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); V2 = os.path.dirname(HERE); ROOT = os.path.dirname(V2)
sys.path.insert(0, V2); sys.path.insert(0, HERE)
import batch_v2                                   # noqa: E402
from batch_v2 import BatchV2, cacr_from_effcai    # noqa: E402
import model_v2 as MV                             # noqa: E402
import rho_v4, rho_v4d                            # noqa: E402

EXT = os.path.join(V2, "extracted")
L23_DIR = os.path.join(EXT, "ebner_l23l5_delta-prefire-vseg-rs")
L23_BASIS = os.path.join(ROOT, "basis_results_edges_ebner_l23l5_delta_rs")
L5_BASIS = os.path.join(ROOT, "basis_results_edges_sabrina_n120_delta")
L23_PROTOS = ("letzkus_1ap_dt+10ms", "letzkus_3ap_200hz_dt+10ms", "sjostrom_50hz_dt+10ms")
L5_DIRS = (os.path.join(EXT, "markram_delta-prefire-vca"), os.path.join(EXT, "ebner_delta-prefire-vca"))
L5_PROTOS = ("10Hz_10ms", "sjostrom_0.1hz_dt+10ms", "sjostrom_10hz_dt+10ms", "sjostrom_20hz_dt+10ms",
             "sjostrom_50hz_dt+10ms")
# (label, path, proto, selection, data, role) role: dep = L2/3 should depress, pot = should potentiate, ctl = other
TARGETS = [("L23 1AP+10 all", "L23", "letzkus_1ap_dt+10ms", None, 0.72, "dep"),
           ("L23 3AP+10 distal", "L23", "letzkus_3ap_200hz_dt+10ms", "ldist", 0.79, "dep"),
           ("L23 3AP+10 prox", "L23", "letzkus_3ap_200hz_dt+10ms", "lprox", 1.30, "ctl"),
           ("L23 SH50+10 distal", "L23", "sjostrom_50hz_dt+10ms", "sh", 0.86, "dep"),
           ("L23 SH50+10 all", "L23", "sjostrom_50hz_dt+10ms", None, 1.06, "dep"),
           ("L5 Mk10Hz+10", "L5", "10Hz_10ms", None, 1.2013, "pot"),
           ("L5 Sj0.1Hz+10", "L5", "sjostrom_0.1hz_dt+10ms", None, 0.97, "ctl"),
           ("L5 Sj10Hz+10", "L5", "sjostrom_10hz_dt+10ms", None, 1.16, "pot"),
           ("L5 Sj20Hz+10", "L5", "sjostrom_20hz_dt+10ms", None, 1.31, "pot"),
           ("L5 Sj50Hz+10", "L5", "sjostrom_50hz_dt+10ms", None, 1.57, "pot")]
FEATS = ("m_p", "m_d", "ca_pk", "ca_int", "sh_pk", "vd_int", "n_cev", "b_cev", "w_cev", "n_pre", "n_pair", "f_pair",
         "t_p", "t_d", "pot_frac", "tT_max", "tTu_max", "c_pre", "c_post", "td", "tp", "sec_type", "dist", "path_mean")


def auc(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float); x = x[np.isfinite(x)]; y = y[np.isfinite(y)]
    if not len(x) or not len(y):
        return np.nan
    return float(((x[:, None] > y[None, :]).mean() + 0.5 * (x[:, None] == y[None, :]).mean()))


def pairings(pre, gap=300.0):
    pre = np.sort(np.asarray(pre, float))
    if not len(pre):
        return []
    cut = np.where(np.diff(pre) > gap)[0] + 1
    return np.split(pre, cut)


def load(dirs, protos, basis_dir, pairs):
    batch_v2.BASIS_DIR = basis_dir                 # PairBasis and the loader read the module global
    B = BatchV2(list(dirs), protocols=list(protos), pairs=pairs, fast=False, signals=("vdcc", "shaft_cai"))
    for r in B.recs:
        B.basis(r)                                 # cache now, under this BASIS_DIR
    return B


def unweighted_tT(r, P):
    """feature_T with the t_drive 4 impulses unweighted (every VDCC event counts 1, no 1 - b)."""
    t = r["t"]; imp = np.zeros((len(r["syn"]), len(t)))
    for i in range(len(r["syn"])):
        ev = np.asarray(r["cev"][i], float); ev = ev[np.isfinite(ev)]
        if len(ev):
            np.add.at(imp[i], np.clip(np.searchsorted(t, ev - 1e-9), 0, len(t) - 1), 1.0)
    return MV.feature_T(imp, t, r["arr"], P)


def run_path(name, B, fit, P, loc, styp, dist5):
    gamma = rho_v4.opts(P)[0]; gd, gp = rho_v4.rates(P); tau_d = float(P["tau_d_NMDA"])
    rho = rho_v4d.rho_all(B, fit["a"], P)
    feats = B.features(P)
    recs, rows, ex = [], [], {}
    for j, r in enumerate(B.recs):
        td, tp = rho_v4.thetas(r["c_pre"], r["c_post"], fit["a"], gamma)
        rf = rho[r["sl"]]; r0 = r["rho0"]; b = B.basis(r); tT, K = feats[j]
        dp = lambda tt, am, ano: MV.dpre_final(tt, K, am, ano, P["dpre_min"], P["dpre_max"], P["dpre0"])
        d = dp(tT, P["A_mglu"], P["A_NO"]); z = np.zeros(len(rf))
        tTu = unweighted_tT(r, P); du = dp(tTu, P["A_mglu"], P["A_NO"])
        R = lambda rr, dd: rho_v4.ratio(b, r0, rr, dd)
        L = loc.get(r["pair"], {})
        recs.append(dict(path=name, pair=r["pair"], proto=r["proto"], n_syn=len(rf), **L,
                         full=R(rf, d), rho_only=R(rf, z), pre_only=R(r0, d), no_ecb=R(rf, dp(tT, 0.0, P["A_NO"])),
                         no_no=R(rf, dp(tT, P["A_mglu"], 0.0)), pot_off=R(np.minimum(rf, r0), d),
                         dep_off=R(np.maximum(rf, r0), d), ecb_unw=R(rf, du), ecb_unw_pre=R(r0, du),
                         dpre=float(d.mean()), dpre_unw=float(du.mean()), rho0=float(r0.mean()), rho_f=float(rf.mean())))
        # per-synapse local signals
        t = r["t"]; E = r["effcai"]; cacr = cacr_from_effcai(E, t); h = np.append(np.diff(t), 0.0) / 1e3
        sh = r.get("shaft_cai"); vd = r.get("vdcc")
        delay = MV.edge_params(r["syn"])["delay"]
        pre = np.asarray(r["prespikes"], float); prs = pairings(pre)
        f_pair = 1e3 / np.median(np.diff([p[0] for p in prs])) if len(prs) > 1 else np.nan
        ab_p = (E > tp[:, None]); band = (E > td[:, None]) & ~ab_p
        t_p = (ab_p * h).sum(1); t_d = (band * h).sum(1)
        npair = max(len(prs), 1)
        acc = {k: np.zeros(len(rf)) for k in ("eff_pk", "ca_pk", "ca_int", "sh_pk", "vd_int", "n_cev", "w_cev", "bsum")}
        evw = []
        for i in range(len(rf)):
            ev, w = MV.glu_weight(r["cev"][i], pre + delay[i], tau_d) if r.get("cev") is not None else (np.zeros(0),) * 2
            evw.append((ev, w))
        for p in prs:
            lo, hi = np.searchsorted(t, p[0] - 30.0), np.searchsorted(t, p[-1] + 150.0)
            if hi - lo < 2:
                continue
            sl = slice(lo, hi); hh = h[sl]
            acc["eff_pk"] += E[:, sl].max(1); acc["ca_pk"] += cacr[:, sl].max(1)
            acc["ca_int"] += (cacr[:, sl] * hh * 1e3).sum(1)
            if sh is not None:
                acc["sh_pk"] += sh[:, sl].max(1) - sh[:, lo]
            if vd is not None:
                acc["vd_int"] += (vd[:, sl] * hh * 1e3).sum(1)
            for i, (ev, w) in enumerate(evw):
                k = (ev >= p[0] - 30.0) & (ev < p[-1] + 150.0)
                acc["n_cev"][i] += k.sum(); acc["w_cev"][i] += w[k].sum(); acc["bsum"][i] += (1 - w[k]).sum()
        for k in acc:
            acc[k] /= npair
        for i in range(len(rf)):
            s = int(r["syn"][i])
            rows.append(dict(path=name, pair=r["pair"], proto=r["proto"], syn=s, **L,
                             c_pre=float(r["c_pre"][i]), c_post=float(r["c_post"][i]), td=float(td[i]), tp=float(tp[i]),
                             rho0=float(r0[i]), rho_f=float(rf[i]), dpre=float(d[i]), dpre_unw=float(du[i]),
                             m_p=acc["eff_pk"][i] / tp[i], m_d=acc["eff_pk"][i] / td[i], ca_pk=acc["ca_pk"][i],
                             ca_int=acc["ca_int"][i], sh_pk=acc["sh_pk"][i], vd_int=acc["vd_int"][i],
                             n_cev=acc["n_cev"][i], w_cev=acc["w_cev"][i],
                             b_cev=acc["bsum"][i] / acc["n_cev"][i] if acc["n_cev"][i] > 0 else np.nan,
                             n_pre=len(pre) / npair, n_pair=len(prs), f_pair=f_pair,
                             t_p=t_p[i] / npair, t_d=t_d[i] / npair,
                             pot_frac=gp * t_p[i] / (gp * t_p[i] + gd * t_d[i]) if t_p[i] + t_d[i] > 0 else np.nan,
                             tT_max=float(tT[i].max()) if tT.shape[1] else 0.0,
                             tTu_max=float(tTu[i].max()) if tTu.shape[1] else 0.0,
                             sec_type=styp.get(s, np.nan), dist=dist5.get(s, np.nan)))
        # example: first up-flipping synapse per proto (middle pairing)
        up = np.where((r0 < 0.5) & (rf >= 0.5))[0]
        if r["proto"] not in ex.get(name, {}) and len(up) and len(prs):
            i = int(up[0]); p = prs[len(prs) // 2]; m = (t >= p[0] - 40) & (t < p[-1] + 200)
            ex.setdefault(name, {})[r["proto"]] = dict(t=t[m] - p[0], E=E[i, m], tp=tp[i], td=td[i], pair=r["pair"])
    return pd.DataFrame(recs), pd.DataFrame(rows), ex


def select(df, path, proto, sel):
    x = df[(df.path == path) & (df.proto == proto)]
    if sel == "ldist":
        x = x[x.letzkus_distal == True]   # noqa: E712
    elif sel == "lprox":
        x = x[x.letzkus_distal == False]  # noqa: E712
    elif sel == "sh":
        x = x[x.sh_distal == True]        # noqa: E712
    return x


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
    print("rule: rho_gamma", rho_v4.opts(P)[0], "rates", rho_v4.rates(P), "t_drive", P["t_drive"],
          "theta_Te", P["theta_Te"], "theta_Tg", P["theta_Tg"], "A_mglu", P["A_mglu"], "A_NO", P["A_NO"], flush=True)
    g = pd.read_csv(os.path.join(ROOT, "ebner", "pair_geometry_L23PC_L5TTPC.csv"))
    g["pair"] = g.pregid.astype(str) + "-" + g.postgid.astype(str); g = g[g.all_protocols]
    loc23 = {p: dict(letzkus_distal=bool(a), sh_distal=bool(b), path_mean=float(c))
             for p, a, b, c in zip(g.pair, g.letzkus_distal, g.sh_distal, g.path_mean)}
    st = np.load(os.path.join(V2, "syn_section_type.npz")); styp = dict(zip(st["syn"].tolist(), st["section_type"].tolist()))
    dist5 = {}
    for f in glob.glob(os.path.join(V2, "local_t", "out", "*.npz")):
        zz = np.load(f)
        for s, dd in zip(zz["og|ap1|sid"], zz["og|ap1|dist"]):
            dist5[int(s)] = float(dd)
    B5 = load(L5_DIRS, L5_PROTOS, L5_BASIS, set(fa["pairs"].split(",")))
    R5, S5, ex5 = run_path("L5", B5, fit, P, {}, styp, dist5)
    del B5
    B23 = load([L23_DIR], L23_PROTOS, L23_BASIS, set(g.pair))
    R23, S23, ex23 = run_path("L23", B23, fit, P, loc23, styp, dist5)
    del B23
    R = pd.concat([R5, R23], ignore_index=True); S = pd.concat([S5, S23], ignore_index=True)
    # section type for synapses not in syn_section_type.npz: edges afferent_section_type
    miss = sorted(set(S.syn[S.sec_type.isna()].astype(int)))
    if miss:
        import h5py
        with h5py.File(MV.EDGES, "r") as f:
            v = f[f"edges/{MV.EDGE_POP}/0/afferent_section_type"][np.array(miss)]
        mm = dict(zip(miss, v.astype(float))); S["sec_type"] = S.sec_type.fillna(S.syn.map(mm))
    S["up"] = (S.rho0 < 0.5) & (S.rho_f >= 0.5); S["down"] = (S.rho0 >= 0.5) & (S.rho_f < 0.5)
    S.to_csv(a_.save + "_syn.csv", index=False); R.to_csv(a_.save + "_rec.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 50)
    cols = ["full", "rho_only", "pre_only", "no_ecb", "no_no", "pot_off", "dep_off", "ecb_unw", "ecb_unw_pre"]
    dec, grp = [], {}
    for lab, path, proto, sel, data, role in TARGETS:
        x = select(R, path, proto, sel); y = select(S, path, proto, sel)
        grp[lab] = (y, role)
        dec.append(dict(target=lab, role=role, data=data, n_rec=len(x), n_syn=len(y), **x[cols].mean().to_dict(),
                        dpre=x.dpre.mean(), dpre_unw=x.dpre_unw.mean(), frac_rho0_1=(y.rho0 >= 0.5).mean(),
                        P_up=y.up.sum() / max((y.rho0 < 0.5).sum(), 1), P_down=y.down.sum() / max((y.rho0 >= 0.5).sum(), 1)))
    D = pd.DataFrame(dec); D.to_csv(a_.save + "_decomp.csv", index=False)
    print("=== ratio decomposition (record means; full should match v4_A0g_s3.csv / l23_v4_A0g_s3.csv)")
    print(D.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    med = []
    for lab, (y, role) in grp.items():
        q = {"target": lab, "role": role}
        for f in FEATS:
            if f in y:
                q[f] = y[f].median()
        q["up_m_p"] = y[y.up].m_p.median() if y.up.any() else np.nan
        med.append(q)
    M = pd.DataFrame(med); M.to_csv(a_.save + "_med.csv", index=False)
    print("=== per-target medians of local signals")
    print(M.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    dep = pd.concat([y for y, ro in grp.values() if ro == "dep"]).drop_duplicates(["path", "pair", "proto", "syn"])
    pot = pd.concat([y for y, ro in grp.values() if ro == "pot"])
    sep = []
    for sub, fx in (("all", lambda x: x), ("rho0=0", lambda x: x[x.rho0 < 0.5]), ("rho0=1", lambda x: x[x.rho0 >= 0.5]),
                    ("up", lambda x: x[x.up])):
        a, b = fx(dep), fx(pot)
        for f in FEATS:
            if f in a:
                sep.append(dict(subset=sub, feat=f, n_dep=len(a), n_pot=len(b), med_dep=a[f].median(), med_pot=b[f].median(),
                                auc_dep_gt_pot=auc(a[f], b[f])))
    Q = pd.DataFrame(sep); Q.to_csv(a_.save + "_sep.csv", index=False)
    Q["sepn"] = (Q.auc_dep_gt_pot - 0.5).abs()
    print("=== separation L2/3 should-depress vs L5 should-potentiate: AUC(dep > pot), sorted by |AUC - 0.5|")
    for sub in Q.subset.unique():
        print(f"--- {sub}"); print(Q[Q.subset == sub].sort_values("sepn", ascending=False).drop(columns="sepn")
                                    .to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    # within L2/3: proximal 3AP+10 (LTP 1.30) vs distal 3AP+10 (LTD 0.79), same protocol
    a, b = grp["L23 3AP+10 distal"][0], grp["L23 3AP+10 prox"][0]
    print("=== L2/3 3AP+10 distal vs proximal: AUC(distal > prox)")
    print(pd.DataFrame([dict(feat=f, med_dist=a[f].median(), med_prox=b[f].median(), auc=auc(a[f], b[f]))
                        for f in FEATS if f in a]).to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    fig_(D, S, grp, Q, {**{("L5", k): v for k, v in ex5.get("L5", {}).items()},
                        **{("L23", k): v for k, v in ex23.get("L23", {}).items()}}, a_.fig)
    print(f"saved {a_.save}_{{rec,syn,decomp,med,sep}}.csv, {a_.fig}  ({len(S)} synapse rows, {len(R)} records)")


def fig_(D, S, grp, Q, ex, path):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 3, figsize=(17, 9.5))
    x = np.arange(len(D)); w = 0.16
    for k, (c, lab, colr) in enumerate((("data", "data", "k"), ("full", "model", "#c0392b"), ("rho_only", "post only (dpre 0)", "#e59866"),
                                        ("pre_only", "pre only (rho0)", "#2e86c1"), ("ecb_unw", "eCB unweighted", "#27ae60"))):
        ax[0, 0].bar(x + (k - 2) * w, D[c], w, color=colr, label=lab)
    ax[0, 0].axhline(1, color="k", lw=0.6); ax[0, 0].set_xticks(x, D.target, rotation=45, ha="right", fontsize=7)
    ax[0, 0].set(title="a  EPSP ratio decomposition", ylabel="ratio"); ax[0, 0].legend(fontsize=7)
    ax[0, 1].bar(x - 0.2, D.P_up, 0.4, color="#c0392b", label="P(up | rho0 = 0)")
    ax[0, 1].bar(x + 0.2, D.P_down, 0.4, color="#2e86c1", label="P(down | rho0 = 1)")
    ax[0, 1].plot(x, D.frac_rho0_1, "ko", ms=4, label="frac rho0 = 1")
    ax[0, 1].set_xticks(x, D.target, rotation=45, ha="right", fontsize=7); ax[0, 1].set(title="b  rho flips"); ax[0, 1].legend(fontsize=7)
    q = Q[Q.subset == "rho0=0"].set_index("feat").auc_dep_gt_pot.dropna(); q2 = Q[Q.subset == "all"].set_index("feat").auc_dep_gt_pot
    y = np.arange(len(q))
    ax[0, 2].barh(y - 0.2, q.values, 0.4, color="#c0392b", label="rho0 = 0")
    ax[0, 2].barh(y + 0.2, q2.reindex(q.index).values, 0.4, color="#7f8c8d", label="all")
    ax[0, 2].axvline(0.5, color="k", lw=0.8); ax[0, 2].set_yticks(y, q.index, fontsize=7)
    ax[0, 2].set(title="c  AUC(L2/3 should-depress > L5 +10 LTP)", xlim=(0, 1)); ax[0, 2].legend(fontsize=7)
    colr = {"dep": "#c0392b", "pot": "#2e86c1", "ctl": "#7f8c8d"}
    for lab, (yy, role) in grp.items():
        z = yy[yy.rho0 < 0.5]
        ax[1, 0].scatter(z.m_p, z.pot_frac, s=4, alpha=0.35, color=colr[role], label=f"{lab} ({role})" if role != "ctl" else None)
    ax[1, 0].set(xscale="log", title="d  rho0 = 0 synapses: pairing peak / theta_p vs pot_frac", xlabel="m_p", ylabel="pot_frac")
    ax[1, 0].axvline(1, color="k", lw=0.6); ax[1, 0].legend(fontsize=6, markerscale=3)
    data = [yy.b_cev.dropna().values for yy, _ in grp.values()]
    ax[1, 1].boxplot(data, showfliers=False); ax[1, 1].set_xticks(np.arange(1, len(grp) + 1), list(grp), rotation=45, ha="right", fontsize=7)
    ax[1, 1].set(title="e  NMDA-bound b at own VDCC events (eCB weight = 1 - b)")
    for (pth, pr), e in ex.items():
        ax[1, 2].plot(e["t"], e["E"] / e["tp"], lw=0.9, label=f"{pth} {pr}", color="#2e86c1" if pth == "L5" else "#c0392b",
                      alpha=0.8, ls="-" if "50" in pr or "10Hz" in pr else "--")
    ax[1, 2].axhline(1, color="k", lw=0.6)
    ax[1, 2].set(title="f  effcai / theta_p, first up-flipping synapse (mid pairing)", xlabel="ms from pairing onset")
    ax[1, 2].legend(fontsize=6)
    fig.tight_layout(); fig.savefig(path, dpi=130)


if __name__ == "__main__":
    main()
