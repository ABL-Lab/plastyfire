"""Can the spine decode the Ca source from its single pool? (CA_DECODE.md). CPU/numba, no jax, fixed fit.

The C1 gate reads c_VDCC (own spine VDCC current, low-passed 100 ms); among up-flipping synapses vd_int separates
L2/3 should-depress from L5 should-potentiate (AUC 0.20) while total-pool amplitude/integral do not (ca_pk 0.56,
ca_int 0.35). Here every feature is computed from the synapse's own spine Ca trace only (cacr, recovered from effcai),
per pairing window [first pre - 30, last pre + 150] ms, mean over pairings (as diag_l23_ltp.py):
  dca_*      dCa/dt peak; integral of pos(dCa/dt - k), k = 0 and ca_ref / {1, 3, 10} ms
  s_*        sensor dS/dt = kon Ca^n (1 - S) - koff S, n in {1, 2}, Kd = koff/kon (half saturation, Kd^n = koff/kon)
             in {0.5, 1, 2} x ca_ref, tau_off = 1/koff in {1, 3, 10} ms; readouts pk, int, tab (time S > 0.5),
             R (peak of 100 ms low-pass of S, as tauE1) and Ratp (that low-pass where effcai > theta_p)
  hp*        Ca - lowpass(Ca, tau), tau in {10, 30, 100} ms: peak and positive integral
  sharp      pairing peak / pairing integral
  refs       vd_int, ca_pk, ca_int, sh_pk, m_p (diag_l23_ltp code path); cv_pk, cv_atp (c_VDCC-like low-pass of own
             VDCC current: peak and value while effcai > theta_p)
ca_ref = median ca_pk of up-flippers in results/l23ltp_A0g_s3_syn.csv (fallback 3.5e-3 mM).
Per feature: auc_up (dep vs pot up-flippers), auc_all, auc_3ap (3AP+10 distal > proximal, all synapses),
Spearman with vd_int (all rows; up-flippers in dep + pot).

    python glusynapse_v2/rho_redesign/diag_ca_decode.py --fit <json> --save <prefix> --fig <png> --md <CA_DECODE.md>
"""
import argparse, glob, json, os, sys
import numpy as np, pandas as pd
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import diag_l23_ltp as D                          # noqa: E402  (sets sys.path for batch_v2/model_v2/rho_v4)
from diag_l23_ltp import auc, pairings, select, TARGETS, V2, ROOT  # noqa: E402
from batch_v2 import cacr_from_effcai             # noqa: E402
import model_v2 as MV                             # noqa: E402
import rho_v4                                     # noqa: E402

TAU_E1 = 100.0                                     # ms, readout low-pass (as c_VDCC)
S_THR = 0.5
KD_MULT = (("lo", 0.5), ("mid", 1.0), ("hi", 2.0))
TAU_OFF = (1.0, 3.0, 10.0)
NS = (1, 2)
HP_TAU = (10.0, 30.0, 100.0)
DCA_TAU = (1.0, 3.0, 10.0)
KEYS = ["path", "pair", "proto", "syn"]
REFS = ("vd_int", "ca_pk", "ca_int", "sh_pk", "m_p", "cv_pk", "cv_atp")


@njit(cache=False)
def _sensor(ca, h, kon, koff, n, tau, thr, E, tp):
    """Exact step update for piecewise-constant Ca (cacr is a step mean). Returns (m, 5):
    peak S, integral S (ms), time S > thr (ms), peak R, peak R while E > tp; R = low-pass(S, tau)."""
    m, L = ca.shape
    out = np.zeros((m, 5))
    for i in range(m):
        c = max(ca[i, 0], 0.0); x = kon * c ** n
        S = x / (x + koff); R = S
        spk = S; sint = 0.0; stab = 0.0; rpk = R; ratp = 0.0
        for k in range(L - 1):
            hk = h[k]
            c = max(ca[i, k], 0.0); x = kon * c ** n; r = x + koff; si = x / r
            sint += S * hk
            if S > thr:
                stab += hk
            S = si + (S - si) * np.exp(-r * hk)
            R = S + (R - S) * np.exp(-hk / tau)
            if S > spk:
                spk = S
            if R > rpk:
                rpk = R
            if E[i, k + 1] > tp[i] and R > ratp:
                ratp = R
        out[i, 0] = spk; out[i, 1] = sint; out[i, 2] = stab; out[i, 3] = rpk; out[i, 4] = ratp
    return out


@njit(cache=False)
def _lowhigh(x, h, tau, E, tp):
    """lp = low-pass(x, tau) (exact for step-constant x, starts at steady state x[0]). Returns (m, 4):
    peak(x - lp), integral pos(x - lp) (ms), peak lp, peak lp while E > tp."""
    m, L = x.shape
    out = np.zeros((m, 4))
    for i in range(m):
        lp = x[i, 0]; hpk = 0.0; hpi = 0.0; lpk = lp; latp = 0.0
        for k in range(L - 1):
            c = x[i, k]; d = c - lp
            if d > hpk:
                hpk = d
            if d > 0.0:
                hpi += d * h[k]
            lp = c + (lp - c) * np.exp(-h[k] / tau)
            if lp > lpk:
                lpk = lp
            if E[i, k + 1] > tp[i] and lp > latp:
                latp = lp
        out[i, 0] = hpk; out[i, 1] = hpi; out[i, 2] = lpk; out[i, 3] = latp
    return out


def sensor_names():
    out = []
    for n in NS:
        for kl, _ in KD_MULT:
            for to in TAU_OFF:
                out.append(f"s_n{n}_kd{kl}_to{int(to)}")
    return out


def decode(name, B, fit, P, ca_ref):
    gamma = rho_v4.opts(P)[0]
    rows = []
    r0 = B.recs[0]; dt = np.diff(r0["t"])
    print(f"{name}: grid min/median dt {dt.min():.4g}/{np.median(dt):.4g} ms, max {dt.max():.4g} ms "
          f"({len(r0['t'])} samples, {r0['pair']} {r0['proto']})", flush=True)
    sens = [(f"s_n{n}_kd{kl}_to{int(to)}", n, ca_ref * km, 1.0 / to) for n in NS for kl, km in KD_MULT for to in TAU_OFF]
    kk = [("k0", 0.0)] + [(f"r{int(tt)}ms", ca_ref / tt) for tt in DCA_TAU]
    for r in B.recs:
        td, tp = rho_v4.thetas(r["c_pre"], r["c_post"], fit["a"], gamma)
        tp = np.asarray(tp, np.float64)
        t = r["t"]; E = r["effcai"]; cacr = cacr_from_effcai(E, t); vd = r.get("vdcc")
        m = len(r["syn"])
        pre = np.asarray(r["prespikes"], float); prs = pairings(pre); npair = max(len(prs), 1)
        acc = {}

        def add(k, v):
            acc[k] = acc.get(k, 0.0) + v

        for p in prs:
            lo, hi = np.searchsorted(t, p[0] - 30.0), np.searchsorted(t, p[-1] + 150.0)
            if hi - lo < 3:
                continue
            hw = np.ascontiguousarray(np.diff(t[lo:hi]), dtype=np.float64)
            ca = np.ascontiguousarray(cacr[:, lo:hi], dtype=np.float64)
            Ew = np.ascontiguousarray(E[:, lo:hi])
            pk = ca.max(1); it = (ca[:, :-1] * hw).sum(1)
            add("ca_pk2", pk); add("ca_int2", it); add("sharp", np.where(it > 0, pk / np.maximum(it, 1e-30), np.nan))
            d = np.diff(ca, axis=1) / hw[None, :]
            add("dca_pk", d.max(1))
            for kl, kv in kk:
                add(f"dca_pos_{kl}", (np.maximum(d - kv, 0.0) * hw).sum(1))
            for nm, n, kd, koff in sens:
                kon = koff / kd ** n
                o = _sensor(ca, hw, kon, koff, n, TAU_E1, S_THR, Ew, tp)
                for j, sfx in enumerate(("pk", "int", "tab", "R", "Ratp")):
                    add(f"{nm}_{sfx}", o[:, j])
            for tau in HP_TAU:
                o = _lowhigh(ca, hw, tau, Ew, tp)
                add(f"hp{int(tau)}_pk", o[:, 0]); add(f"hp{int(tau)}_pint", o[:, 1])
            if vd is not None:
                o = _lowhigh(np.ascontiguousarray(vd[:, lo:hi], dtype=np.float64), hw, TAU_E1, Ew, tp)
                add("cv_pk", o[:, 2]); add("cv_atp", o[:, 3])
        for i in range(m):
            q = dict(path=name, pair=r["pair"], proto=r["proto"], syn=int(r["syn"][i]))
            for k, v in acc.items():
                q[k] = float(np.asarray(v)[i]) / npair
            rows.append(q)
    return pd.DataFrame(rows)


def spearman(x, y):
    z = pd.DataFrame({"x": np.asarray(x, float), "y": np.asarray(y, float)}).replace([np.inf, -np.inf], np.nan).dropna()
    if len(z) < 3 or z.x.nunique() < 2 or z.y.nunique() < 2:
        return np.nan
    return float(np.corrcoef(z.x.rank(), z.y.rank())[0, 1])


def family(f):
    if f in REFS or f in ("ca_pk2", "ca_int2"):
        return "ref"
    for p in ("dca", "s_n1", "s_n2", "hp", "sharp"):
        if f.startswith(p):
            return p
    return "other"


def md_table(T):
    cols = ["feat", "family", "auc_up", "auc_all", "auc_3ap", "rho_vd_all", "rho_vd_up", "med_dep_up", "med_pot_up"]
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, q in T.iterrows():
        cells = []
        for c in cols:
            v = q[c]
            cells.append(v if isinstance(v, str) else (f"{v:.3g}" if np.isfinite(v) else "nan"))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True); ap.add_argument("--save", required=True); ap.add_argument("--fig", required=True)
    ap.add_argument("--md", default=None); ap.add_argument("--ref-syn", default=os.path.join(HERE, "results", "l23ltp_A0g_s3_syn.csv"))
    a_ = ap.parse_args()
    fit = json.load(open(a_.fit)); fa = fit["args"]
    D.batch_v2.BAP_GATE = fit.get("bap_gate")
    fil = {**json.loads(fa["filters"]), **json.loads(fa.get("set", "{}"))}
    P = {**MV.DEFAULTS, **fil, **fit["pre"]}
    for k in ("gamma_d", "gamma_p"):
        if k in fit:
            P[k] = fit[k]
    ca_ref = 3.5e-3
    if os.path.isfile(a_.ref_syn):
        z = pd.read_csv(a_.ref_syn, usecols=["ca_pk", "up"])
        ca_ref = float(z.ca_pk[z.up.astype(bool)].median())
    print(f"ca_ref {ca_ref:.4g} mM -> Kd {[round(ca_ref * k, 6) for _, k in KD_MULT]} mM; "
          f"dCa/dt k {[round(ca_ref / t, 7) for t in DCA_TAU]} mM/ms", flush=True)
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
    B5 = D.load(D.L5_DIRS, D.L5_PROTOS, D.L5_BASIS, set(fa["pairs"].split(",")))
    _, S5, _ = D.run_path("L5", B5, fit, P, {}, styp, dist5)
    F5 = decode("L5", B5, fit, P, ca_ref)
    del B5
    B23 = D.load([D.L23_DIR], D.L23_PROTOS, D.L23_BASIS, set(g.pair))
    _, S23, _ = D.run_path("L23", B23, fit, P, loc23, styp, dist5)
    F23 = decode("L23", B23, fit, P, ca_ref)
    del B23
    S = pd.concat([S5, S23], ignore_index=True); F = pd.concat([F5, F23], ignore_index=True)
    S["up"] = (S.rho0 < 0.5) & (S.rho_f >= 0.5)
    keep = KEYS + ["rho0", "rho_f", "up", "letzkus_distal", "sh_distal"] + [c for c in ("vd_int", "ca_pk", "ca_int", "sh_pk", "m_p") if c in S]
    S = S[keep].merge(F, on=KEYS, how="left")
    feats = [c for c in S.columns if c not in KEYS + ["rho0", "rho_f", "up", "letzkus_distal", "sh_distal"]]
    S.to_csv(a_.save + "_syn.csv", index=False)
    grp = {lab: (select(S, path, proto, sel), role) for lab, path, proto, sel, data, role in TARGETS}
    dep = pd.concat([y for y, ro in grp.values() if ro == "dep"]).drop_duplicates(KEYS)
    pot = pd.concat([y for y, ro in grp.values() if ro == "pot"])
    du, pu = dep[dep.up], pot[pot.up]
    dist, prox = grp["L23 3AP+10 distal"][0], grp["L23 3AP+10 prox"][0]
    upall = pd.concat([du, pu])
    out = []
    for f in feats:
        out.append(dict(feat=f, family=family(f), n_dep_up=len(du), n_pot_up=len(pu),
                        auc_up=auc(du[f], pu[f]), auc_all=auc(dep[f], pot[f]), auc_3ap=auc(dist[f], prox[f]),
                        rho_vd_all=spearman(S[f], S.vd_int), rho_vd_up=spearman(upall[f], upall.vd_int),
                        med_dep_up=du[f].median(), med_pot_up=pu[f].median()))
    T = pd.DataFrame(out)
    T["sep_up"] = (T.auc_up - 0.5).abs()
    T["useful"] = (T.sep_up >= 0.25) & (((T.auc_up < 0.5) & (T.auc_3ap <= 0.3)) | ((T.auc_up > 0.5) & (T.auc_3ap >= 0.7)))
    T = T.sort_values("sep_up", ascending=False)
    T.to_csv(a_.save + "_auc.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 20); pd.set_option("display.max_rows", 300)
    print(f"=== n up-flippers dep {len(du)} pot {len(pu)} (expect 235 / 112); 3AP distal {len(dist)} prox {len(prox)}")
    print(T.drop(columns="sep_up").to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print("useful (|auc_up - .5| >= .25 and 3AP same side):", T[T.useful].feat.tolist())
    if a_.md:
        top = pd.concat([T[T.family != "ref"].head(25), T[T.family == "ref"]])
        txt = open(a_.md).read() if os.path.isfile(a_.md) else ""
        head = txt.split("<!-- TABLE -->")[0]
        body = (f"<!-- TABLE -->\nca_ref = {ca_ref:.4g} mM; n up-flippers dep {len(du)}, pot {len(pu)}; 3AP distal {len(dist)}, "
                f"prox {len(prox)}. Top 25 single-pool features by |auc_up - 0.5|, then references. Full table: "
                f"`{os.path.relpath(a_.save, HERE)}_auc.csv`.\n\n" + md_table(top) + "\n\n"
                f"Useful (|auc_up - 0.5| >= 0.25 and 3AP on the same side): {', '.join(T[T.useful].feat) or 'none'}\n")
        open(a_.md, "w").write(head + body)
    fig_(T, du, pu, dist, prox, a_.fig)
    print(f"saved {a_.save}_{{syn,auc}}.csv, {a_.fig}" + (f", {a_.md}" if a_.md else ""))


def fig_(T, du, pu, dist, prox, path):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    col = {"dca": "#e67e22", "s_n1": "#2e86c1", "s_n2": "#8e44ad", "hp": "#27ae60", "sharp": "#c0392b", "ref": "#7f8c8d",
           "other": "k"}
    fig, ax = plt.subplots(2, 2, figsize=(15, 11))
    vd = T.set_index("feat").loc["vd_int"] if "vd_int" in set(T.feat) else None
    top = T.head(35).iloc[::-1]; y = np.arange(len(top))
    ax[0, 0].barh(y, top.auc_up, color=[col[f] for f in top.family])
    ax[0, 0].set_yticks(y, top.feat, fontsize=6); ax[0, 0].axvline(0.5, color="k", lw=0.8)
    for v in (0.25, 0.75):
        ax[0, 0].axvline(v, color="k", lw=0.6, ls=":")
    if vd is not None:
        ax[0, 0].axvline(vd.auc_up, color="#c0392b", lw=1, ls="--", label=f"vd_int {vd.auc_up:.2f}")
    ax[0, 0].set(xlim=(0, 1), title="a  AUC(L2/3 dep > L5 pot), up-flippers, top 35 by |AUC - 0.5|")
    ax[0, 0].legend(fontsize=7)
    for fam, c in col.items():
        z = T[T.family == fam]
        if len(z):
            ax[0, 1].scatter(z.auc_up, z.auc_3ap, s=14, color=c, label=fam)
            ax[1, 0].scatter(z.auc_up, z.rho_vd_up, s=14, color=c, label=fam)
    if vd is not None:
        ax[0, 1].scatter([vd.auc_up], [vd.auc_3ap], s=80, marker="*", color="#c0392b", label="vd_int")
    ax[0, 1].axvline(0.5, color="k", lw=0.6); ax[0, 1].axhline(0.5, color="k", lw=0.6)
    ax[0, 1].axvspan(0, 0.25, color="#c0392b", alpha=0.07); ax[0, 1].axhspan(0, 0.3, color="#c0392b", alpha=0.07)
    ax[0, 1].set(xlim=(0, 1), ylim=(0, 1), xlabel="auc_up (dep vs pot up-flippers)", ylabel="auc_3ap (distal > prox)",
                 title="b  both separations"); ax[0, 1].legend(fontsize=7)
    ax[1, 0].axvline(0.5, color="k", lw=0.6)
    ax[1, 0].set(xlim=(0, 1), ylim=(-1, 1), xlabel="auc_up", ylabel="Spearman with vd_int (up-flippers)",
                 title="c  is the feature a proxy for own VDCC Ca?"); ax[1, 0].legend(fontsize=7)
    best = T[T.family != "ref"].iloc[0].feat
    data = [du[best].dropna(), pu[best].dropna(), dist[best].dropna(), prox[best].dropna()]
    ax[1, 1].boxplot(data, showfliers=False)
    ax[1, 1].set_xticks([1, 2, 3, 4], ["L2/3 dep up", "L5 pot up", "3AP+10 distal", "3AP+10 prox"], fontsize=8)
    if all((d > 0).all() and len(d) for d in data):
        ax[1, 1].set_yscale("log")
    ax[1, 1].set(title=f"d  best single-pool feature: {best}")
    fig.tight_layout(); fig.savefig(path, dpi=120)


if __name__ == "__main__":
    main()
