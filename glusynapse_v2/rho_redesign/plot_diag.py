"""Figures + summary tables for the rho diagnosis (diag_rho.py outputs in rho_redesign/out). CPU, matplotlib.
    python glusynapse_v2/rho_redesign/plot_diag.py [--tag td4_s2]
fig1_margins: per-synapse effcai peak / theta_p vs c_post, per protocol (L5 and L2/3), marker = transition.
fig2_decomp:  per target: full, rho only, pre only vs data (L5 and L2/3 with the loc selections).
fig3_traces:  effcai / theta_p around the first pairings, 50 Hz +-10 (L5) and Letzkus 3AP +-10 (L2/3).
"""
import argparse, os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); V2 = os.path.dirname(HERE)
sys.path.insert(0, V2)
from targets import load_targets   # noqa: E402

OUT = os.path.join(HERE, "out"); FIG = os.path.join(HERE, "figs")
C_UP, C_DN, C_NO = "#1f77b4", "#d62728", "#999999"


def sel(R, pid, l23):
    proto, _, where = pid.partition("@")
    s = R[R.proto == proto]
    if l23 and where and proto.startswith("letzkus"):
        s = s[s.letzkus_distal == (where == "distal")]
    elif l23 and where == "distal":
        s = s[s.sh_distal.astype(bool)]
    elif where == "distal":
        return s.iloc[:0]
    return s


def summarise(S, name):
    S = S.copy()
    S["floor"] = S.c_post < 2e-4; S["m_p"] = S.peak / S.tp
    S["up"] = (S.rho0 < 0.5) & (S.rho_f >= 0.5); S["down"] = (S.rho0 >= 0.5) & (S.rho_f < 0.5)
    S["sat"] = (S.rho_f < 0.02) | (S.rho_f > 0.98)
    grp = ["proto", "floor"] if "letzkus_distal" not in S else ["proto", "letzkus_distal", "floor"]
    rows = []
    for k, x in S.groupby(grp):
        n0 = (x.rho0 < 0.5).sum(); n1 = (x.rho0 >= 0.5).sum()
        rows.append(dict(zip(grp, k if isinstance(k, tuple) else (k,)), n=len(x), m_p_med=x.m_p.median(),
                         f_mp_0p8_1p25=((x.m_p > 0.8) & (x.m_p < 1.25)).mean(), f_above=(x.m_p > 1).mean(),
                         P_up=x.up.sum() / max(n0, 1), P_down=x.down.sum() / max(n1, 1), f_sat=x.sat.mean(),
                         t_tp=x.t_above_tp.median(), t_band=x.t_band.median(), dpre=x.dpre.mean()))
    T = pd.DataFrame(rows); T.to_csv(os.path.join(OUT, f"{name}_margins.csv"), index=False)
    print(f"=== {name}: margins / transitions"); print(T.round(3).to_string(index=False))
    if S.dist.notna().any():
        S["dbin"] = pd.cut(S.dist, [0, 60, 150, 250, 2000])
        D = S[S.dist.notna()].groupby(["proto", "dbin"], observed=True).agg(n=("m_p", "size"), m_p=("m_p", "median"),
                                                                              c_post=("c_post", "median"),
                                                                              up=("up", "mean"), down=("down", "mean"))
        print(f"=== {name}: by soma distance (local_t synapses, {S.dist.notna().sum()} rows)"); print(D.round(4).to_string())
    return S


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tag", default="td4_s2"); a_ = ap.parse_args()
    os.makedirs(FIG, exist_ok=True)
    S5 = summarise(pd.read_csv(os.path.join(OUT, f"l5_{a_.tag}_syn.csv")), f"l5_{a_.tag}")
    S23 = summarise(pd.read_csv(os.path.join(OUT, f"l23_{a_.tag}_syn.csv")), f"l23_{a_.tag}")
    R5 = pd.read_csv(os.path.join(OUT, f"l5_{a_.tag}_rec.csv")); R23 = pd.read_csv(os.path.join(OUT, f"l23_{a_.tag}_rec.csv"))

    # ---- fig 1: margins
    P5 = ["sjostrom_0.1hz_dt+10ms", "sjostrom_0.1hz_dt-10ms", "sjostrom_20hz_dt-10ms", "sjostrom_50hz_dt+10ms",
          "sjostrom_50hz_dt-10ms", "sjostrom07_step200ms_post_only"]
    P23 = ["letzkus_1ap_dt+10ms", "letzkus_3ap_200hz_dt+10ms", "letzkus_3ap_200hz_dt-10ms", "sjostrom_50hz_dt+10ms"]
    fig, ax = plt.subplots(2, 6, figsize=(22, 7.5), sharey=True)
    for row, (S, protos, lab) in enumerate(((S5, P5, "L5-L5"), (S23, P23, "L2/3-L5"))):
        for j in range(6):
            A = ax[row, j]
            if j >= len(protos):
                A.axis("off"); continue
            x = S[S.proto == protos[j]]
            cp = np.maximum(x.c_post, 5e-5)
            col = np.where(x.up, C_UP, np.where(x.down, C_DN, C_NO))
            mk = x.letzkus_distal.map({True: "^", False: "o"}) if "letzkus_distal" in x else None
            if mk is None:
                A.scatter(cp, x.m_p, s=7, c=col, alpha=0.6, lw=0)
            else:
                for m_ in ("o", "^"):
                    q = (mk == m_).values
                    A.scatter(cp[q], x.m_p[q], s=9, c=col[q], alpha=0.6, lw=0, marker=m_)
            A.axhline(1, color="k", lw=0.8); A.axhline((x.td / x.tp).median(), color="k", ls=":", lw=0.8)
            A.axvline(2e-4, color="0.6", ls="--", lw=0.6)
            A.set_xscale("log"); A.set_yscale("log"); A.set_ylim(0.2, 30)
            A.set_title(f"{lab}\n{protos[j]}", fontsize=8)
            A.set_xlabel("c_post (mM, bAP Ca at the synapse)", fontsize=8)
        ax[row, 0].set_ylabel("effcai peak / theta_p")
    ax[0, 0].text(0.02, 0.97, "blue 0->1, red 1->0, grey no flip; dotted theta_d/theta_p; dashed c_post floor\n"
                  "L2/3: triangles = distal pairs (rise > 2.7 ms)", transform=ax[0, 0].transAxes, fontsize=7, va="top")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig1_margins.png"), dpi=130); plt.close(fig)

    # ---- fig 2: decomposition per target
    fig, ax = plt.subplots(1, 2, figsize=(20, 5.5), gridspec_kw=dict(width_ratios=[3, 1.2]))
    for A, R, groups, l23 in ((ax[0], R5, ("paired_l5", "sjostrom07"), False), (ax[1], R23, ("paired_l23l5",), True)):
        T = {k: v for k, v in load_targets(groups).items() if k[1] == "control"}
        rows = []
        for (pid, _), (m, s, n, src) in sorted(T.items()):
            x = sel(R, pid, l23)
            if len(x):
                rows.append(dict(target=pid, data=m, sem=s, full=x.full.mean(), rho_only=x.rho_only.mean(),
                                 pre_only=x.pre_only.mean(), cont=x.cont.mean()))
        D = pd.DataFrame(rows); D.to_csv(os.path.join(OUT, f"decomp_{'l23' if l23 else 'l5'}_{a_.tag}.csv"), index=False)
        print(f"=== decomposition {'L2/3' if l23 else 'L5'}"); print(D.round(3).to_string(index=False))
        xx = np.arange(len(D)); w = 0.2
        A.bar(xx - 1.5 * w, D.full, w, label="model (full)", color="0.3")
        A.bar(xx - 0.5 * w, D.rho_only, w, label="post rho only (dpre 0)", color=C_UP)
        A.bar(xx + 0.5 * w, D.pre_only, w, label="pre only (rho frozen)", color="#ff7f0e")
        A.errorbar(xx + 1.5 * w, D["data"], D["sem"], fmt="ks", ms=4, label="data")
        A.axhline(1, color="k", lw=0.6); A.set_xticks(xx)
        A.set_xticklabels([t.replace("sjostrom_", "sj_").replace("letzkus_", "lz_") for t in D.target], rotation=70, fontsize=7)
        A.set_ylabel("EPSP ratio"); A.set_title("L5-L5 fit targets (td4_s2)" if not l23 else "L2/3-L5 transfer (no refit)")
    ax[0].legend(fontsize=7, ncol=4)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig2_decomp.png"), dpi=130); plt.close(fig)

    # ---- fig 3: traces
    fig, ax = plt.subplots(2, 2, figsize=(13, 7))
    for A, (who, proto) in zip(ax.ravel(), (("l5", "sjostrom_50hz_dt+10ms"), ("l5", "sjostrom_50hz_dt-10ms"),
                                            ("l23", "letzkus_3ap_200hz_dt+10ms"), ("l23", "letzkus_3ap_200hz_dt-10ms"))):
        z = np.load(os.path.join(OUT, f"{who}_{a_.tag}_traces.npz"))
        if f"{proto}|t" not in z.files:
            A.set_title(f"{proto}: no trace"); continue
        t = z[f"{proto}|t"]; E = z[f"{proto}|E"]; tp = z[f"{proto}|tp"]; td = z[f"{proto}|td"]; cpo = z[f"{proto}|c_post"]
        win = t < (1200 if who == "l5" else 400)
        for i in range(E.shape[0]):
            col = C_UP if z[f"{proto}|rho_f"][i] >= 0.5 else C_DN
            A.plot(t[win], E[i, win] / tp[i], color=col, lw=0.8, alpha=0.8,
                   ls="-" if cpo[i] >= 2e-4 else "--")
        A.axhline(1, color="k", lw=0.8); A.axhline(np.median(td / tp), color="k", ls=":", lw=0.8)
        for s in z[f"{proto}|pre"]:
            if s < t[win][-1]:
                A.axvline(s, color="g", lw=0.4, alpha=0.5)
        for s in z[f"{proto}|post"]:
            if s < t[win][-1]:
                A.axvline(s, color="m", lw=0.4, alpha=0.5)
        pr = [p.split("|")[1] for p in z["pairs"] if p.startswith(proto + "|")]
        A.set_title(f"{who} {proto} (pair {pr[0] if pr else '?'}); blue rho_f>=0.5, dashed = c_post floor", fontsize=8)
        A.set_xlabel("ms from first pre spike"); A.set_ylabel("effcai / theta_p")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig3_traces.png"), dpi=130); plt.close(fig)
    print("figs in", FIG)


if __name__ == "__main__":
    main()
