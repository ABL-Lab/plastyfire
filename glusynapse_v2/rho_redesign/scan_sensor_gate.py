"""CA_DECODE.md section 3: replace the C1 / C2 c_VDCC gate by a gate read from ONE spine Ca pool. CPU/numba, no jax,
no refit (A0g_s3 parameters; the A0 rows must reproduce L5 34.16, L2/3 220.4). Scores the 29 L5 targets
(paired_l5 + sjostrom07) and the 9 L2/3->L5 targets, as scan_vgate_amp.py.

Input: the synapse's own spine Ca, Ca = cacr (cai_CR - min, recovered exactly from effcai, as diag_ca_decode.py).
    sensor   S' = kon Ca (1 - S) - koff S       n = 1, tau_off = 1/koff = 3 ms, Kd = koff/kon in {0.5, 1} x ca_ref
    pool S   G' = -G / tau_g + S                (G in ms: leaky integral of sensor occupancy)
    pool T   G' = -G / tau_g + H(Ca - Kd)       (G in ms: leaky time above Kd)
    tau_g in {0.1, 0.5, 2, 10} s
    pot = [effcai > theta_p] * [G > theta_g]
    C1: dep = [effcai > theta_d]                        (a crossing without enough sensor drive depresses)
    C2: dep = [effcai > theta_d] * [effcai <= theta_p]  (the gated-out crossing is neutral)
    rho' as scan_vgate_amp.rho_rec (rho_v4 with tau_fast 0, same skip rule).
theta_g grid per pool: quantiles {2,5,10,20,30,40,50,60,75,90}% of G_max at theta_p crossings over all crossing
synapses of both paths (pre-pass), plus the geometric mid of the L2/3 should-depress / L5 should-potentiate up-flipper
medians. The grid is one set of numbers per pool, applied to both paths (uniform).

Window check (pre-pass): the diag_ca_decode `_int` feature integrates S over [first pre - 30, last pre + 150] ms, a
protocol-dependent length. Recomputed here per pairing as (var) that integral, (mean) that integral / window length,
(fix200) the integral over [first pre, first pre + 200] ms; for S and for time above Kd; plus window length and pre
spikes per pairing alone. AUCs as diag_ca_decode: auc_up (L2/3 dep vs L5 pot up-flippers), auc_3ap (3AP+10 distal >
proximal, all synapses).

    python glusynapse_v2/rho_redesign/scan_sensor_gate.py --fit <json> --save <prefix> --fig <png> [--md CA_DECODE.md]
        [--kd 0.5,1] [--taug 100,500,2000,10000] [--types S,T]
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import scan_vgate_amp as SV                       # noqa: E402  (sets sys.path; constants, targets, loaders)
from scan_vgate_amp import (BatchV2, RHO_STAR_GB, TAU_IND_GB, load_targets, MV, rho_v4, KEY,  # noqa: E402
                            L23_DIR, L23_BASIS, L5_BASIS, GEOM)
import batch_v2                                   # noqa: E402
from batch_v2 import cacr_from_effcai             # noqa: E402
from diag_l23_ltp import auc, pairings, select, TARGETS   # noqa: E402
from diag_ca_decode import _sensor, spearman      # noqa: E402

TAU_OFF = 3.0
QS = (0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 0.9)
FIX_WIN = 200.0
BUDGETS = (40.0, 55.0, 90.0, 140.0)
MARK = "<!-- GATE_SCAN -->"
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f} s] {msg}", flush=True)


@njit(cache=False)
def gate_rec(E, CA, h, hs, td, tp, rho0, gd, gp, rs, kds, koff, taus, AG, BG, pk, pt, pty, TH, do_rho,
             out, out0, gmax, xing):
    """One pass per synapse. out (2, NP, nth, n): rho_f for C1 / C2 per pool and theta; out0 (n,): A0 rho_f;
    gmax (NP, n): max G over steps with effcai > theta_p; xing (n,): any theta_p crossing."""
    n, T = E.shape; nk = kds.shape[0]; NP = pk.shape[0]; nth = TH.shape[1]
    kon = koff / kds
    S = np.empty(nk); G = np.empty(NP); R = np.empty((2, NP, nth))
    for i in range(n):
        c0 = max(CA[i, 0], 0.0)
        for j in range(nk):
            x = kon[j] * c0
            S[j] = x / (x + koff)
        for p in range(NP):
            u = S[pk[p]] if pty[p] == 0 else (1.0 if c0 > kds[pk[p]] else 0.0)
            G[p] = taus[pt[p]] * u                     # steady state of the initial input
        for c in range(2):
            for p in range(NP):
                for q in range(nth):
                    R[c, p, q] = rho0[i]
        r0 = rho0[i]; free = (r0 != 0.0) and (r0 != 1.0)
        a = td[i]; b = tp[i]; crossed = False
        for k in range(T - 1):
            x = E[i, k]
            ab = x > b; aa = x > a
            if ab:
                crossed = True
                for p in range(NP):
                    if G[p] > gmax[p, i]:
                        gmax[p, i] = G[p]
            hk = hs[k]
            if aa or ab or (r0 != 0.0 and r0 != 1.0):
                pot = 1.0 if ab else 0.0
                dep = 1.0 if aa else 0.0
                r0 += hk * (-r0 * (1 - r0) * (rs - r0) + pot * gp * (1 - r0) - dep * (1 - pot) * gd * r0)
                r0 = min(max(r0, 0.0), 1.0)
            if do_rho and (aa or ab or free):
                nf = False
                for p in range(NP):
                    g = G[p]
                    for q in range(nth):
                        po = 1.0 if (ab and g > TH[p, q]) else 0.0
                        for c in range(2):
                            r = R[c, p, q]
                            if not (aa or ab) and (r == 0.0 or r == 1.0):
                                continue
                            if c == 0:
                                dep = 1.0 if aa else 0.0
                            else:
                                dep = 1.0 if (aa and not ab) else 0.0
                            r += hk * (-r * (1 - r) * (rs - r) + po * gp * (1 - r) - dep * (1 - po) * gd * r)
                            r = min(max(r, 0.0), 1.0)
                            R[c, p, q] = r
                            if r != 0.0 and r != 1.0:
                                nf = True
                free = nf
            # advance sensor (exact for step-constant Ca) and pools (input = step-mean S or H(Ca - Kd))
            ck = max(CA[i, k], 0.0); hm = h[k]
            for j in range(nk):
                x = kon[j] * ck; rr = x + koff; si = x / rr; e = np.exp(-rr * hm)
                sav = si + (S[j] - si) * (1.0 - e) / (rr * hm) if hm > 0.0 else S[j]
                for p in range(NP):
                    if pk[p] == j:
                        u = sav if pty[p] == 0 else (1.0 if ck > kds[j] else 0.0)
                        G[p] = AG[pt[p], k] * G[p] + BG[pt[p], k] * u
                S[j] = si + (S[j] - si) * e
        for c in range(2):
            for p in range(NP):
                for q in range(nth):
                    out[c, p, q, i] = R[c, p, q]
        out0[i] = r0; xing[i] = crossed


class Pools:
    def __init__(self, types, kdm, taug, ca_ref):
        self.kdm = np.asarray(kdm, float); self.kds = self.kdm * ca_ref; self.taus = np.asarray(taug, float)
        self.list = [(ty, ik, it) for ty in types for ik in range(len(kdm)) for it in range(len(taug))]
        self.pk = np.array([x[1] for x in self.list], np.int64); self.pt = np.array([x[2] for x in self.list], np.int64)
        self.pty = np.array([0 if x[0] == "S" else 1 for x in self.list], np.int64)
        self.names = [f"{ty}_kd{self.kdm[ik]:g}_tg{self.taus[it] / 1000:g}s" for ty, ik, it in self.list]

    def meta(self, p):
        ty, ik, it = self.list[p]
        return dict(pool=self.names[p], ptype=ty, kd_mult=float(self.kdm[ik]), tau_g_s=float(self.taus[it] / 1000))


def window_feats(r, tp, E, cacr, PO):
    """Per-synapse mean over pairings: S integral / mean / fixed-window integral, time above Kd likewise."""
    t = r["t"]; m = len(r["syn"]); koff = 1.0 / TAU_OFF
    prs = pairings(np.asarray(r["prespikes"], float)); npair = max(len(prs), 1)
    acc = {}

    def add(k, v):
        acc[k] = acc.get(k, 0.0) + v

    for p in prs:
        for wn, (t_lo, t_hi) in (("var", (p[0] - 30.0, p[-1] + 150.0)), ("fix200", (p[0], p[0] + FIX_WIN))):
            lo, hi = np.searchsorted(t, t_lo), np.searchsorted(t, t_hi)
            if hi - lo < 3:
                continue
            hw = np.ascontiguousarray(np.diff(t[lo:hi]), dtype=np.float64)
            ca = np.ascontiguousarray(cacr[:, lo:hi], dtype=np.float64); Ew = np.ascontiguousarray(E[:, lo:hi])
            wl = float(t[hi - 1] - t[lo])
            if wn == "var":
                add("win_len", np.full(m, wl)); add("n_pre", np.full(m, float(len(p))))
            for ik, kd in enumerate(PO.kds):
                o = _sensor(ca, hw, koff / kd, koff, 1, 100.0, 0.5, Ew, tp)
                ta = ((ca[:, :-1] > kd) * hw[None, :]).sum(1)
                tag = f"kd{PO.kdm[ik]:g}"
                add(f"S_{tag}_{wn}", o[:, 1]); add(f"T_{tag}_{wn}", ta)
                if wn == "var":
                    add(f"S_{tag}_mean", o[:, 1] / wl); add(f"T_{tag}_mean", ta / wl)
    return {k: np.asarray(v, float) / npair for k, v in acc.items()}


def run(B, fit, P, T, path, loc, sh, PO, TH, do_rho, win_protos=()):
    """do_rho False: pre-pass (A0 rho, G_max at crossings, window features) -> synapse table.
    do_rho True: scored replay -> long table of target predictions per (pool, cand, theta)."""
    gamma = rho_v4.opts(P)[0]; gd, gp = rho_v4.rates(P)
    NP = len(PO.list); nth = TH.shape[1]
    need = {}
    for pid, cond in T:
        need.setdefault(pid.split("@")[0], set()).add(cond)
    feats = B.features(P) if do_rho else None
    for r in B.recs:
        r.pop("_c4", None)
    rows, syn = [], []
    for j, r in enumerate(B.recs):
        td, tp = rho_v4.thetas(r["c_pre"], r["c_post"], fit["a"], gamma)
        td = np.asarray(td, np.float64); tp = np.asarray(tp, np.float64)
        t = r["t"]; h = np.diff(t); hs = h / 1000.0 / TAU_IND_GB
        E = np.ascontiguousarray(r["effcai"]); CA = np.ascontiguousarray(cacr_from_effcai(E, t))
        AG = np.exp(-h[None, :] / PO.taus[:, None]); BG = PO.taus[:, None] * (1.0 - AG)
        n = len(r["syn"]); out = np.empty((2, NP, nth, n)); out0 = np.empty(n)
        gmax = np.zeros((NP, n)); xing = np.zeros(n, np.bool_)
        gate_rec(E, CA, h, hs, td, tp, r["rho0"].astype(float), gd, gp, float(RHO_STAR_GB), PO.kds, 1.0 / TAU_OFF,
                 PO.taus, AG, BG, PO.pk, PO.pt, PO.pty, TH, do_rho, out, out0, gmax, xing)
        if do_rho:
            b = B.basis(r); tT, K = feats[j]
            for c in sorted(need.get(r["proto"], ())):
                if c == "nmdar_block":
                    d = np.zeros(n)
                else:
                    d = MV.dpre_final(tT, K, 0.0 if c == "mglu_block" else P["A_mglu"],
                                      0.0 if (c == "no_block" or (c == "post_nmdar" and P["no_drive"])) else P["A_NO"],
                                      P["dpre_min"], P["dpre_max"], P["dpre0"])
                fixed = c in ("post_nmdar", "nmdar_block")
                v0 = rho_v4.ratio(b, r["rho0"], r["rho0"] if fixed else out0, d)
                rows.append((path, r["pair"], r["proto"], c, "A0", "A0", -1, 0.0, v0))
                for p in range(NP):
                    for q in range(nth):
                        for ci, cn in enumerate(("C1", "C2")):
                            v = v0 if fixed else rho_v4.ratio(b, r["rho0"], out[ci, p, q], d)
                            rows.append((path, r["pair"], r["proto"], c, PO.names[p], cn, q, float(TH[p, q]), v))
        else:
            wf = window_feats(r, tp, E, CA, PO) if r["proto"] in win_protos else {}
            up = (r["rho0"] < 0.5) & (out0 >= 0.5)
            for i in range(n):
                q = dict(path=path, pair=r["pair"], proto=r["proto"], syn=int(r["syn"][i]), rho0=float(r["rho0"][i]),
                         rho_A0=float(out0[i]), up=bool(up[i]), xing=bool(xing[i]),
                         letzkus_distal=loc.get(r["pair"]), sh_distal=sh.get(r["pair"]))
                for p in range(NP):
                    q["g_" + PO.names[p]] = float(gmax[p, i]) if xing[i] else np.nan
                for k, v in wf.items():
                    q[k] = float(v[i])
                syn.append(q)
    if not do_rho:
        return pd.DataFrame(syn)
    df = pd.DataFrame(rows, columns=["path", "pair", "proto", "condition", "pool", "cand", "iq", "theta", "ratio"])
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
        for (pn, cn, iq), x in sel.groupby(["pool", "cand", "iq"]):
            pred = x.ratio.mean()
            res.append(dict(path=path, target=f"{pid}|{cond}", pool=pn, cand=cn, iq=iq, theta=float(x.theta.iloc[0]),
                            data=m, sem=s, pred=pred, z=(pred - m) / s))
    return pd.DataFrame(res)


def groups(S):
    grp = {lab: (select(S, path, proto, sel), role) for lab, path, proto, sel, data, role in TARGETS}
    dep = pd.concat([y for y, ro in grp.values() if ro == "dep"]).drop_duplicates(["path", "pair", "proto", "syn"])
    pot = pd.concat([y for y, ro in grp.values() if ro == "pot"])
    return dep, pot, grp["L23 3AP+10 distal"][0], grp["L23 3AP+10 prox"][0]


def theta_grid(S, PO):
    dep, pot, _, _ = groups(S)
    TH = np.empty((len(PO.list), len(QS) + 1)); info = []
    for p, nm in enumerate(PO.names):
        g = S["g_" + nm].dropna().values
        qv = np.quantile(g, QS) if len(g) else np.zeros(len(QS))
        md, mp = dep[dep.up]["g_" + nm].median(), pot[pot.up]["g_" + nm].median()
        mid = np.sqrt(md * mp) if (md > 0 and mp > 0) else 0.5 * (np.nan_to_num(md) + np.nan_to_num(mp))
        TH[p] = np.sort(np.append(qv, mid))
        for q in range(TH.shape[1]):
            info.append(dict(table="grid", pool=nm, iq=q, theta=float(TH[p, q])))
    return TH, pd.DataFrame(info)


def checks(S, PO):
    dep, pot, dist, prox = groups(S)
    du, pu = dep[dep.up], pot[pot.up]
    log(f"n up-flippers dep {len(du)} pot {len(pu)} (diag_ca_decode: 235 / 112); 3AP distal {len(dist)} prox {len(prox)}")
    out = []
    wfe = [c for c in S.columns if c.startswith(("S_kd", "T_kd")) or c in ("win_len", "n_pre")]
    for f in wfe:
        out.append(dict(table="window", feat=f, n_dep_up=len(du), n_pot_up=len(pu), auc_up=auc(du[f], pu[f]),
                        auc_all=auc(dep[f], pot[f]), auc_3ap=auc(dist[f], prox[f]),
                        rho_winlen_up=spearman(pd.concat([du, pu])[f], pd.concat([du, pu]).win_len),
                        med_dep_up=du[f].median(), med_pot_up=pu[f].median()))
    for p, nm in enumerate(PO.names):
        f = "g_" + nm
        out.append(dict(table="gmax", feat=f, **PO.meta(p), n_dep_up=int(du[f].notna().sum()),
                        n_pot_up=int(pu[f].notna().sum()), auc_up=auc(du[f], pu[f]),
                        auc_3ap=auc(dist[f].dropna(), prox[f].dropna()), med_dep_up=du[f].median(),
                        med_pot_up=pu[f].median()))
    return pd.DataFrame(out), (du, pu, dist, prox)


def summarize(R, PO):
    meta = {nm: PO.meta(p) for p, nm in enumerate(PO.names)}
    meta["A0"] = dict(pool="A0", ptype="A0", kd_mult=np.nan, tau_g_s=np.nan)
    summ = []
    for (pn, cn, iq), x in R.groupby(["pool", "cand", "iq"]):
        q = dict(**meta[pn], cand=cn, iq=int(iq), theta=float(x.theta.iloc[0]),
                 chi2_L5=float((x[x.path == "L5"].z ** 2).sum()), n_L5=int((x.path == "L5").sum()),
                 chi2_L23=float((x[x.path == "L23"].z ** 2).sum()), n_L23=int((x.path == "L23").sum()))
        q["chi2_tot"] = q["chi2_L5"] + q["chi2_L23"]
        for k in KEY:
            for pth in ("L5", "L23"):
                y = x[(x.target == k) & (x.path == pth)]
                if len(y):
                    q[f"{pth}:{k.split('|')[0]}"] = float(y.pred.iloc[0])
        summ.append(q)
    return pd.DataFrame(summ)


def budget_table(M, ref):
    sets = [("c_VDCC C1", ref[ref.cand == "C1"], "ref"), ("c_VDCC C2", ref[ref.cand == "C2"], "ref")]
    for ty, lab in (("S", "sensor pool S"), ("T", "time>Kd pool T")):
        if (M.ptype == ty).any():
            sets.append((lab, M[M.ptype == ty], "new"))
    out = []
    for bd in BUDGETS:
        for lab, X, kind in sets:
            y = X[X.chi2_L5 <= bd]
            if not len(y):
                out.append(dict(budget=bd, gate=lab, chi2_L5=np.nan, chi2_L23=np.nan, at="none")); continue
            z = y.loc[y.chi2_L23.idxmin()]
            at = (f"{z.cand} theta_V {z.theta_V:g}" if kind == "ref" else
                  f"{z.cand} Kd {z.kd_mult:g}xca_ref tau_g {z.tau_g_s:g} s theta {z.theta:.4g} ms (grid {int(z.iq)})")
            out.append(dict(budget=bd, gate=lab, chi2_L5=float(z.chi2_L5), chi2_L23=float(z.chi2_L23), at=at))
    return pd.DataFrame(out)


def md_tab(df, cols, fmt="{:.3g}"):
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, q in df.iterrows():
        lines.append("| " + " | ".join(q[c] if isinstance(q[c], str) else
                                       (fmt.format(q[c]) if pd.notna(q[c]) else "-") for c in cols) + " |")
    return "\n".join(lines)


def verdict(M, ref, BT, C):
    sens = BT[BT.gate.isin(["sensor pool S", "time>Kd pool T"])]; vd = BT[BT.gate.str.startswith("c_VDCC")]
    cmp_ = []
    for bd in (55.0, 90.0):
        s = sens[sens.budget == bd].chi2_L23.min(); r = vd[vd.budget == bd].chi2_L23.min()
        cmp_.append((bd, s, r, (s / r if np.isfinite(s) else np.inf) if r > 0 else np.nan))
    ok = [c[3] <= 1.1 for c in cmp_ if not np.isnan(c[3])]
    new = M[M.ptype != "A0"]; b = new.loc[new.chi2_tot.idxmin()]; rb = ref.loc[(ref.chi2_L5 + ref.chi2_L23).idxmin()]
    if ok and all(ok):
        v = "YES: the single-pool sensor gate matches or beats c_VDCC"
    elif any(ok):
        v = "PARTLY: the single-pool sensor gate matches c_VDCC at one L5 budget only"
    else:
        v = "NO: the single-pool sensor gate does not reach the c_VDCC trade-off"
    txt = [f"**Verdict (auto, no refit).** {v} (min L2/3 chi2 at L5 chi2 <= budget; sensor / c_VDCC: "
           + "; ".join(f"<= {bd:g}: {s:.1f} / {r:.1f} = {q:.2f}" for bd, s, r, q in cmp_) + ").",
           f"Best sensor point by L5 + L2/3 chi2: {b.cand}, pool {b.ptype}, Kd {b.kd_mult:g} x ca_ref, tau_g {b.tau_g_s:g} s, "
           f"theta {b.theta:.4g} ms: L5 {b.chi2_L5:.1f}, L2/3 {b.chi2_L23:.1f} (c_VDCC best: {rb.cand} theta_V "
           f"{rb.theta_V:g}: L5 {rb.chi2_L5:.1f}, L2/3 {rb.chi2_L23:.1f}; A0 34.16 / 220.4)."]
    w = C[C.table == "window"].set_index("feat")
    if "S_kd0.5_fix200" in w.index:
        a_v, a_m, a_f = (w.loc[f"S_kd0.5_{s}", "auc_up"] for s in ("var", "mean", "fix200"))
        a_wl = w.loc["win_len", "auc_up"] if "win_len" in w.index else np.nan
        conf = ("survives a fixed window: synapse-local" if a_f <= 0.27 and a_m <= 0.3 else
                "is largely a window-length / spike-count effect" if min(a_f, a_m) >= 0.35 else "is partly window-driven")
        txt.append(f"Window check (Kd 0.5 x ca_ref): auc_up var {a_v:.3f}, per-ms mean {a_m:.3f}, fixed 200 ms {a_f:.3f}; "
                   f"window length alone {a_wl:.3f}. The single-pool separation {conf}.")
    return "\n".join(txt)


def write_md(path, M, ref, BT, C, TH_info, save):
    txt = open(path).read() if os.path.isfile(path) else ""
    head = txt.split(MARK)[0].rstrip() + "\n\n"
    W = C[C.table == "window"].copy(); G = C[C.table == "gmax"].copy()
    top = M[M.ptype != "A0"].sort_values("chi2_tot").head(12)
    keyc = [c for c in ("L23:letzkus_1ap_dt+10ms", "L23:letzkus_3ap_200hz_dt+10ms@distal",
                        "L23:letzkus_3ap_200hz_dt+10ms@proximal", "L23:sjostrom_50hz_dt+10ms@distal", "L5:10Hz_10ms",
                        "L5:sjostrom_10hz_dt+10ms", "L5:sjostrom_50hz_dt+10ms") if c in M]
    a0 = M[M.ptype == "A0"].iloc[0]
    body = [MARK, "## 3. Gate scan: single-pool sensor gate in place of c_VDCC (no refit)", "",
            f"`scan_sensor_gate.py` (job {os.environ.get('SLURM_JOB_ID', '?')}), fit A0g_s3, no refit; outputs "
            f"`{os.path.relpath(save, HERE)}{{,_targets,_check}}.csv`, `figs/fig10_sensor_gate.png`.",
            "Rule: own spine Ca (cacr) drives one sensor S' = kon Ca (1 - S) - koff S (n 1, tau_off 3 ms, Kd = {0.5, 1} x "
            "ca_ref). Pool S: G' = -G/tau_g + S; pool T: G' = -G/tau_g + H(Ca - Kd) (G in ms). pot = H(c* - theta_p) "
            "H(G - theta_g); C1 converts a gated-out crossing to depression, C2 leaves it neutral. theta_g grid per pool: "
            "quantiles 2-90% of G_max at theta_p crossings (both paths) plus the dep/pot up-flipper mid. No membrane "
            "voltage, no source label, one parameter set for both pathways.",
            f"A0 check: L5 {a0.chi2_L5:.2f}, L2/3 {a0.chi2_L23:.2f} (must be 34.16 / 220.4).", "",
            "### 3a. Window check (is the `_int` separation window length?)", "",
            md_tab(W, ["feat", "auc_up", "auc_all", "auc_3ap", "rho_winlen_up", "med_dep_up", "med_pot_up"]), "",
            "kd = multiple of ca_ref; S = sensor integral, T = time above Kd; var = [first pre - 30, last pre + 150] ms (diag_ca_decode), mean = var / "
            "window length, fix200 = [first pre, first pre + 200] ms. Reference vd_int: auc_up 0.20, auc_3ap 0.205.", "",
            "### 3b. Gate pool at theta_p crossings (G_max), per pool", "",
            md_tab(G, ["pool", "auc_up", "auc_3ap", "med_dep_up", "med_pot_up", "n_dep_up", "n_pot_up"]), "",
            "### 3c. chi2: best L2/3 at a given L5 budget, sensor gate vs c_VDCC (results/vgate_amp_A0g_s3.csv)", "",
            md_tab(BT, ["budget", "gate", "chi2_L5", "chi2_L23", "at"], "{:.4g}"), "",
            "### 3d. Best 12 sensor points by L5 + L2/3 chi2", "",
            md_tab(top, ["cand", "ptype", "kd_mult", "tau_g_s", "theta", "chi2_L5", "chi2_L23"] + keyc), "",
            "Data: L2/3 1AP 0.72, 3AP distal 0.79, 3AP proximal 1.30, S&H distal 0.86; L5 Mk10Hz 1.20, Sj10Hz 1.16, "
            "Sj50Hz 1.57.", "", verdict(M, ref, BT, C), ""]
    open(path, "w").write(head + "\n".join(body))


def fig_(M, ref, C, Sall, path):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(15, 11))
    a = ax[0, 0]
    for cn, ls in (("C1", "-"), ("C2", "--")):
        x = ref[ref.cand == cn].sort_values("theta_V")
        a.plot(x.chi2_L5, x.chi2_L23, ls, color="k", marker="o", ms=3, label=f"c_VDCC {cn}")
    taus = sorted(M.tau_g_s.dropna().unique()); cm = plt.get_cmap("viridis")
    for (ty, cn, kd, tg), x in M[M.ptype != "A0"].groupby(["ptype", "cand", "kd_mult", "tau_g_s"]):
        col = cm(taus.index(tg) / max(len(taus) - 1, 1))
        a.scatter(x.chi2_L5, x.chi2_L23, s=18 if kd == 0.5 else 34, marker="o" if ty == "S" else "^",
                  facecolors=col if cn == "C1" else "none", edgecolors=col, linewidths=0.8)
    for tg in taus:
        a.scatter([], [], color=cm(taus.index(tg) / max(len(taus) - 1, 1)), label=f"tau_g {tg:g} s")
    a.scatter([], [], marker="o", color="grey", label="pool S (o) / T (^); filled C1, open C2; small Kd 0.5")
    a.scatter([34.16], [220.4], marker="*", s=150, color="#c0392b", label="A0", zorder=5)
    a.set(xscale="log", yscale="log", xlabel=r"$\chi^2$ L5 (29)", ylabel=r"$\chi^2$ L2/3$\to$L5 (9)",
          title="a  trade-off, no refit: sensor gate vs c_VDCC"); a.legend(fontsize=6)
    new = M[M.ptype != "A0"]; b = new.loc[new.chi2_tot.idxmin()]
    x = new[(new.pool == b.pool) & (new.cand == b.cand)].sort_values("theta")
    cols = {"L23:letzkus_1ap_dt+10ms": ("#c0392b", 0.72), "L23:letzkus_3ap_200hz_dt+10ms@distal": ("#e67e22", 0.79),
            "L23:letzkus_3ap_200hz_dt+10ms@proximal": ("#8e44ad", 1.30), "L23:sjostrom_50hz_dt+10ms@distal": ("#d35400", 0.86),
            "L5:10Hz_10ms": ("#2e86c1", 1.2013), "L5:sjostrom_10hz_dt+10ms": ("#1abc9c", 1.16),
            "L5:sjostrom_20hz_dt+10ms": ("#16a085", 1.31), "L5:sjostrom_50hz_dt+10ms": ("#34495e", 1.57)}
    for k, (c, d) in cols.items():
        if k in x:
            ax[0, 1].plot(np.maximum(x.theta, 1e-6), x[k], "-", color=c, marker="o", ms=3,
                          label=k.replace("L23:", "L2/3 ").replace("L5:", "L5 "))
            ax[0, 1].axhline(d, color=c, lw=0.6, ls=":")
    ax[0, 1].set(xscale="log", xlabel=r"$\theta_g$ (ms)", ylabel="ratio",
                 title=f"b  key targets, best point's pool: {b.cand} {b.pool} (dotted = data)"); ax[0, 1].legend(fontsize=6)
    W = C[C.table == "window"].reset_index(drop=True); y = np.arange(len(W))
    ax[1, 0].barh(y - 0.2, W.auc_up, 0.4, color="#2e86c1", label="auc_up (L2/3 dep vs L5 pot up-flippers)")
    ax[1, 0].barh(y + 0.2, W.auc_3ap, 0.4, color="#e67e22", label="auc_3ap (distal > prox)")
    ax[1, 0].set_yticks(y, W.feat, fontsize=7); ax[1, 0].axvline(0.5, color="k", lw=0.8)
    ax[1, 0].axvline(0.20, color="#2e86c1", lw=0.8, ls="--"); ax[1, 0].axvline(0.205, color="#e67e22", lw=0.8, ls=":")
    ax[1, 0].set(xlim=(0, 1), title="c  window check (dashed: vd_int auc_up 0.20)"); ax[1, 0].legend(fontsize=7)
    f = "g_" + b.pool
    dep, pot, dist, prox = groups(Sall)
    for lab, z, ls in (("L2/3 dep, up", dep[dep.up], "--"), ("L5 pot, up", pot[pot.up], "-"),
                       ("3AP+10 distal, cross", dist, ":"), ("3AP+10 prox, cross", prox, "-.")):
        v = np.sort(z[f].dropna().values)
        if len(v):
            ax[1, 1].plot(np.maximum(v, 1e-6), np.linspace(0, 1, len(v)), ls, label=f"{lab} ({len(v)})")
    for th in x.theta.unique():
        ax[1, 1].axvline(max(th, 1e-6), color="k", lw=0.3)
    ax[1, 1].axvline(max(b.theta, 1e-6), color="#c0392b", lw=1.2)
    ax[1, 1].set(xscale="log", xlabel=f"G_max at theta_p crossings, {b.pool} (ms)", ylabel="CDF",
                 title="d  gate pool at crossings (red = best theta)"); ax[1, 1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(path, dpi=120)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True); ap.add_argument("--save", required=True); ap.add_argument("--fig", required=True)
    ap.add_argument("--md", default=None)
    ap.add_argument("--kd", default="0.5,1"); ap.add_argument("--taug", default="100,500,2000,10000")
    ap.add_argument("--types", default="S,T")
    ap.add_argument("--ref-csv", default=os.path.join(HERE, "results", "vgate_amp_A0g_s3.csv"))
    ap.add_argument("--ref-syn", default=os.path.join(HERE, "results", "l23ltp_A0g_s3_syn.csv"))
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
    ca_ref = 3.5e-3
    if os.path.isfile(a_.ref_syn):
        z = pd.read_csv(a_.ref_syn, usecols=["ca_pk", "up"]); ca_ref = float(z.ca_pk[z.up.astype(bool)].median())
    PO = Pools(a_.types.split(","), [float(v) for v in a_.kd.split(",")], [float(v) for v in a_.taug.split(",")], ca_ref)
    log(f"ca_ref {ca_ref:.4g} mM; pools {PO.names}; tau_off {TAU_OFF} ms; {len(QS) + 1} thetas per pool")
    T5 = {k: v for k, v in load_targets(tuple(fa["groups"].split(","))).items() if k[1] in conds_ok}
    T23 = {k: v for k, v in load_targets(("paired_l23l5",)).items() if k[1] in conds_ok}
    g = pd.read_csv(GEOM); g["pair"] = g.pregid.astype(str) + "-" + g.postgid.astype(str); g = g[g.all_protocols]
    loc = dict(zip(g.pair, g.letzkus_distal.astype(bool))); sh = dict(zip(g.pair, g.sh_distal.astype(bool)))
    win_protos = {proto for _, _, proto, _, _, _ in TARGETS}
    p5 = sorted({k[0].split("@")[0] for k in T5}); p23 = sorted({k[0].split("@")[0] for k in T23})
    TH0 = np.zeros((len(PO.list), 0))

    def load23():
        batch_v2.BASIS_DIR = L23_BASIS
        return BatchV2([L23_DIR], protocols=p23, pairs=set(g.pair), fast=False, signals=("vdcc",))

    # 1. L2/3 pre-pass (A0 rho, G_max, window features)
    B = load23(); log("L2/3 loaded")
    S23 = run(B, fit, P, T23, "L23", loc, sh, PO, TH0, False, win_protos); del B; log("L2/3 pre-pass done")
    # 2. L5 pre-pass, theta grid, L5 replay
    batch_v2.BASIS_DIR = L5_BASIS
    B = BatchV2(fa["dirs"].split(","), protocols=p5, pairs=set(fa["pairs"].split(",")), fast=False, signals=("vdcc",))
    for r in B.recs:
        B.basis(r)
    log("L5 loaded")
    S5 = run(B, fit, P, T5, "L5", {}, {}, PO, TH0, False, win_protos); log("L5 pre-pass done")
    Sall = pd.concat([S5, S23], ignore_index=True)
    TH, grid = theta_grid(Sall, PO)
    C, _ = checks(Sall, PO)
    pd.set_option("display.width", 300); pd.set_option("display.max_columns", 60); pd.set_option("display.max_rows", 400)
    print("=== window check + G_max at crossings (AUCs)")
    print(C.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print("=== theta grid (ms)"); print(pd.DataFrame(TH, index=PO.names).to_string(float_format=lambda v: f"{v:.4g}"))
    R5 = run(B, fit, P, T5, "L5", {}, {}, PO, TH, True); del B; log("L5 replay done")
    # 3. L2/3 reload and replay
    B = load23()
    for r in B.recs:
        B.basis(r)
    log("L2/3 reloaded")
    R23 = run(B, fit, P, T23, "L23", loc, sh, PO, TH, True); del B; log("L2/3 replay done")
    R = pd.concat([R5, R23], ignore_index=True)
    R.to_csv(a_.save + "_targets.csv", index=False)
    M = summarize(R, PO); M.to_csv(a_.save + ".csv", index=False)
    pd.concat([C, grid], ignore_index=True).to_csv(a_.save + "_check.csv", index=False)
    a0 = M[M.ptype == "A0"].iloc[0]
    print(f"=== A0 check: L5 {a0.chi2_L5:.2f} (n {a0.n_L5}), L2/3 {a0.chi2_L23:.2f} (n {a0.n_L23}); expect 34.16 / 220.4")
    print("=== chi2 per pool / cand / theta")
    print(M[["pool", "cand", "iq", "theta", "chi2_L5", "chi2_L23", "chi2_tot"]].to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    ref = pd.read_csv(a_.ref_csv); ref = ref[ref.theta_V > 0]
    BT = budget_table(M, ref)
    print("=== best L2/3 chi2 at L5 chi2 budget (sensor vs c_VDCC)"); print(BT.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print("=== best 12 sensor points (key targets)")
    print(M[M.ptype != "A0"].sort_values("chi2_tot").head(12).to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print(verdict(M, ref, BT, C))
    fig_(M, ref, C, Sall, a_.fig)
    if a_.md:
        write_md(a_.md, M, ref, BT, C, grid, a_.save)
    log(f"saved {a_.save}{{,_targets,_check}}.csv, {a_.fig}" + (f", {a_.md}" if a_.md else ""))


if __name__ == "__main__":
    main()
