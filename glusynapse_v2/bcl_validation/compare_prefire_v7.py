"""Prefire BCL validation of a v7 fit: live GluSynapseV7 (live_v7.py) vs the offline v7 rule. compare_prefire_v5.py
is reused unchanged (imported); the offline reference is v7_rec, the CPU port of compare_prefire_v5.v5_rec plus the
two gpu_v7_rho options, with the kernel's arithmetic:
  veto: Q per own arrival = gpu_v7_rho.veto_q (sum bV s from the arrival's grid sample to the first sample past
        t_a + veto_T, bV = tau_E1 (1 - exp(-h/tau_E1)) / i_scale, s = the record's unweighted -ica_VDCC); a triggered
        arrival steps only if not (Q > theta_eCB,i);
  theta_eCB,i = theta_eCB uE_i (theta_eCB <= 0: theta_V); uE 1 (ecb_ref 0) or live_v7.ue_map (ecb_ref 2, the
        fit_v6 SCALE_E rule, the same values the live run sets as uE_GB).
Per record it adds the eCB step counts: necb_live / necb_off (sums over the common synapses), necb_syn_mismatch
(synapses whose counts differ), nveto_live / nveto_off, nvpend_max (live decisions left pending; must be 0).
Tables and verdict as compare_prefire_v5 (<save>_records.csv, <save>_targets.csv, figure when there are targets).

    python glusynapse_v2/bcl_validation/compare_prefire_v7.py --fit <json> --results <jsonl> --save <prefix>
"""
import argparse, json, os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare_prefire_v5 as C5                           # noqa: E402
from compare_prefire_v5 import C, MV, rho_v4, batch_v2, PATHS, V5_DEFAULTS   # noqa: E402
import live_v7                                            # noqa: E402
from numba import njit                                    # noqa: E402


@njit(cache=False)
def v7_rec(E, S, hs, cnt, td, tp, rho0, gd, gp, rs, kt, thV, isc, tauE1, dmin, Ae, dp0, thE, tauD, Q, veto,
           out_rho, out_d, out_necb, out_nveto):
    """compare_prefire_v5.v5_rec (weighted, v5_mode 2) + gpu_v7_rho: thE (n,) per synapse, Q (n, T) veto charge at
    the arrival samples, veto = veto_T > 0. out_necb / out_nveto count arrivals (c per sample), as necb_GB / nveto_GB."""
    n, T = E.shape
    for i in range(n):
        a_td = td[i]; a_tp = tp[i]
        r = rho0[i]; d0 = dp0
        V = 0.0; chV = -1.0; aV = 1.0; bV = 0.0
        W = 0.0; Bg = 0.0; chD = -1.0; aD = 1.0
        ne = 0.0; nv = 0.0
        thEi = thE[i]
        for k in range(T):
            e = E[i, k] * 1.0
            s = S[i, k] * 1.0
            h = hs[k]
            hm = h * kt
            c = cnt[i, k]
            pot = 1.0 if e > a_tp else 0.0
            dep = 1.0 if e > a_td else 0.0
            if pot > 0.0 and thV > 0.0 and not (V > thV):
                pot = 0.0
            r = r + h * (-r * (1 - r) * (rs - r) + pot * gp * (1 - r) - dep * (1 - pot) * gd * r)
            if r < 0.0:
                r = 0.0
            elif r > 1.0:
                r = 1.0
            if c > 0.0:
                Bg = Bg + c
            trig = W > thEi
            if c > 0.0 and Ae > 0.0 and trig:
                if veto and Q[i, k] > thEi:
                    nv += c
                else:
                    d0 = dmin + (d0 - dmin) * (1.0 - Ae) ** c
                    ne += c
            if hm != chV:
                chV = hm
                aV = np.exp(-hm / tauE1); bV = tauE1 * (1.0 - aV) / isc
            wb = 1.0 - (Bg if Bg < 1.0 else 1.0)
            W = aV * W + bV * (wb * s)
            if hm != chD:
                chD = hm
                aD = np.exp(-hm / tauD)
            Bg = Bg * aD
            V = aV * V + bV * s
        out_rho[i] = r
        out_d[i] = d0
        out_necb[i] = ne
        out_nveto[i] = nv


def veto_Q(r, cnt, Tv, tauE1, isc):
    """gpu_v7_rho.veto_q for one record, as an (n, T) array (values at the arrival samples, 0 elsewhere)."""
    t = np.asarray(r["t"], np.float64); n, T = r["effcai"].shape
    hm = np.append(np.diff(t), 0.0); bv = tauE1 * (1.0 - np.exp(-hm / tauE1)) / isc
    Cc = np.concatenate([[0.0], np.cumsum(hm)])[:T]
    S = np.asarray(r["vdcc"], np.float32).astype(np.float64)
    Q = np.zeros((n, T))
    for ii in range(n):
        nz = np.flatnonzero(cnt[ii])
        cq = np.concatenate([[0.0], np.cumsum(bv * S[ii])])
        je = np.searchsorted(Cc, Cc[nz] + Tv, side="right")
        Q[ii, nz] = cq[je] - cq[nz]
    return Q


def offline_v7(path, fit, P, pairs, protos, conds):
    batch_v2.BASIS_DIR = C.BASIS[path]
    B = batch_v2.BatchV2(PATHS[path]["dirs"], protocols=sorted(protos), pairs=set(pairs), fast=False, signals=("vdcc",))
    assert int(P["v5_mode"]) == 2, "v7 = v5_mode 2"
    gamma, _, sigma = rho_v4.opts(P); assert gamma == 1.0
    gd, gp = rho_v4.rates(P); thV = rho_v4.vamp(P)[1]
    kt = 1e3 * float(batch_v2.TAU_IND_GB); rs = float(batch_v2.RHO_STAR_GB)
    Tv = float(P.get("veto_T", 0.0)); tauE1 = float(P["tau_E1"]); isc = float(P["i_scale"])
    thE0 = float(P["theta_eCB"]); um = live_v7.ue_map(fit, path)
    out = {}
    for r in B.recs:
        n, T = r["effcai"].shape
        hs = np.zeros(T); hs[:T - 1] = np.diff(r["t"]) / 1000.0 / batch_v2.TAU_IND_GB
        cnt = np.zeros((n, T)); np.add.at(cnt, (np.broadcast_to(np.arange(n)[:, None], r["arr"].shape), r["arr"]), 1.0)
        td, tp = rho_v4.thetas(r["c_pre"], r["c_post"], fit["a"], 1.0)
        pre, post = map(int, r["pair"].split("-"))
        uE = np.ones(n) if um is None else np.array([um[(pre, post, int(s))] for s in r["syn"]])
        thE = uE * thE0 if thE0 > 0 else np.full(n, thV)
        Q = veto_Q(r, cnt, Tv, tauE1, isc) if Tv > 0 else np.zeros((n, T))
        rho_f = np.empty(n); d_f = np.empty(n); ne = np.empty(n); nv = np.empty(n)
        v7_rec(np.ascontiguousarray(r["effcai"], dtype=np.float32), np.ascontiguousarray(r["vdcc"], dtype=np.float32),
               hs, cnt, np.asarray(td, np.float64), np.asarray(tp, np.float64), r["rho0"].astype(np.float64), gd, gp, rs,
               kt, thV, isc, tauE1, float(P["dpre_min"]), float(P["A_eCB"]), 0.0 + float(P["dpre0"]),
               np.ascontiguousarray(thE, np.float64), float(P["tau_d_NMDA"]), Q, Tv > 0, rho_f, d_f, ne, nv)
        b = B.basis(r)
        ncev = np.isfinite(r["cev"]).sum(axis=1) if r.get("cev") is not None else None
        for c in conds:
            frozen = c in ("post_nmdar", "nmdar_block")
            blocked = c in ("mglu_block", "nmdar_block")           # A_eCB 0 live: no step, no veto decision
            rc = r["rho0"] if frozen else rho_f
            d = (np.zeros(n) if c == "nmdar_block" else np.full(n, 0.0 + float(P["dpre0"])) if c == "mglu_block"
                 else d_f)
            out[(r["pair"], r["proto"], c)] = dict(
                syn=np.asarray(r["syn"]), rho0=r["rho0"], rho_obs=np.asarray(r["rho_obs"], float), rho=np.asarray(rc),
                dpre=np.asarray(d), ncev=ncev, b=b, ratio=rho_v4.ratio(b, r["rho0"], rc, d, 0.0 if frozen else sigma),
                necb=np.zeros(n) if blocked else ne, nveto=np.zeros(n) if blocked else nv, uE=uE)
    return out


def ecb_cols(x, o):
    """eCB step / veto counts of one live result against its offline record (common synapses)."""
    if not x.get("ok") or o is None or "v7" not in x:
        return {}
    st = x["end"] if x["phase"] == "prefire" else x["snap"]
    live = {int(s): i for i, s in enumerate(x["syn"])}
    idx = np.array([live.get(int(s), -1) for s in o["syn"]]); ok = idx >= 0
    nl = np.array(st["necb_GB"])[idx[ok]]; vl = np.array(x["v7"]["nveto"])[idx[ok]]
    no = o["necb"][ok]; vo = o["nveto"][ok]
    return dict(necb_live=float(nl.sum()), necb_off=float(no.sum()), necb_syn_mismatch=int(np.sum(nl != no)),
                nveto_live=float(vl.sum()), nveto_off=float(vo.sum()), nveto_syn_mismatch=int(np.sum(vl != vo)),
                nvpend_max=float(np.max(x["v7"]["nvpend"])),
                uE_maxdiff=float(np.max(np.abs(np.array(x["v7"]["uE"])[idx[ok]] - o["uE"][ok]))))


def records(fit, res, show_syn=False):
    fa = fit["args"]
    batch_v2.BAP_GATE = fit.get("bap_gate")
    P = {**MV.DEFAULTS, **V5_DEFAULTS, **json.loads(fa["filters"]), **json.loads(fa.get("set", "{}")), **fit["pre"]}
    sigma = rho_v4.opts(P)[2]
    last = {}
    for x in res:
        last[tuple(x[k] for k in ("path", "pair", "proto", "cond", "phase", "express"))] = x
    res = [x for x in last.values() if x["phase"] == "prefire" and x["express"] == "cooker"]
    dfs = []
    for path in sorted({x["path"] for x in res}):
        rp = [x for x in res if x["path"] == path]
        off = offline_v7(path, fit, P, {x["pair"] for x in rp}, {x["proto"] for x in rp}, sorted({x["cond"] for x in rp}))
        d = C.compare(rp, off, sigma)
        ex = pd.DataFrame([ecb_cols(x, off.get((x["pair"], x["proto"], x["cond"]))) for x in rp], index=d.index)
        dfs.append(pd.concat([d, ex], axis=1))
        if show_syn:
            for x in rp:
                o = off.get((x["pair"], x["proto"], x["cond"]))
                if not x.get("ok") or o is None:
                    continue
                live = {int(s): i for i, s in enumerate(x["syn"])}
                print(f"\n{path} {x['pair']} {x['proto']} {x['cond']}: syn | necb live/off | nveto live/off | "
                      f"rho live/off | dpre live/off")
                for j, s in enumerate(o["syn"]):
                    i = live.get(int(s))
                    if i is None:
                        continue
                    print(f"  {int(s)} | {x['end']['necb_GB'][i]:.0f}/{o['necb'][j]:.0f} | {x['v7']['nveto'][i]:.0f}/"
                          f"{o['nveto'][j]:.0f} | {x['end']['rho_GB'][i]:.4f}/{o['rho'][j]:.4f} | "
                          f"{x['end']['dpre_GB'][i]:.4f}/{o['dpre'][j]:.4f}")
    df = pd.concat(dfs, ignore_index=True)
    for c in ("ratio_basis_live", "ratio_off", "rho_disagree", "n_common", "dpre_maxdiff", "necb_live", "necb_off",
              "necb_syn_mismatch", "nveto_live", "nveto_off", "nveto_syn_mismatch", "nvpend_max", "uE_maxdiff"):
        if c not in df:
            df[c] = np.nan
    df["d_basis"] = df.ratio_basis_live - df.ratio_off
    return df


def targets(df, fit):
    """compare_prefire_v5.targets for up to three pathways: L5 (fa groups), L23 (--joint), L23L23 (the fit's l23l23
    extra groups, fit csv <save>_l23l23.csv) + live_v7.VAL_GROUPS (validation_only, no fit csv row)."""
    from compare_v4 import load_targets, GEOM, ROOT
    fa = fit["args"]
    g = pd.read_csv(GEOM); g["pair"] = g.pregid.astype(str) + "-" + g.postgid.astype(str)
    dist = dict(zip(g.pair, g.letzkus_distal)); sh = dict(zip(g.pair, g.sh_distal))
    drop = set(t for t in (fa.get("drop_targets") or "").split(",") if t)
    conds = set(fa["conditions"].split(","))
    fitcsv = {}
    for suf, path in (("", "L5"), ("_l23", "L23"), ("_l23l23", "L23L23")):
        f = fa["save"] + suf + ".csv"
        f = f if os.path.isabs(f) else os.path.join(ROOT, f)
        if os.path.isfile(f):
            for _, r in pd.read_csv(f).iterrows():
                fitcsv[(path, r.target, r.condition)] = r.pred
    d0 = df[df.ok.fillna(False).astype(bool)].dropna(subset=["ratio_basis_live", "ratio_off"])
    groups = [("L5", tuple(fa["groups"].split(",")), False)] + ([("L23", ("paired_l23l5",), False)] if fa.get("joint") else [])
    gx, _ = live_v7.extra_groups(fit)
    if gx:
        groups += [("L23L23", gx, False), ("L23L23", live_v7.VAL_GROUPS["L23L23"], True)]
    out = []
    for path, grp, val in groups:
        T = {k: v for k, v in load_targets(grp).items() if k[1] in conds}
        d = d0[d0.path == path]
        for (pid, cond), (m, s, n, src) in sorted(T.items()):
            proto, _, where = pid.partition("@")
            sel = d[(d.proto == proto) & (d.cond == cond)]
            if path == "L23" and where and proto.startswith("letzkus"):
                sel = sel[sel.pair.map(dist) == (where == "distal")]
            elif path == "L23" and where == "distal":
                sel = sel[sel.pair.map(sh).fillna(False).astype(bool)]
            elif where:
                continue
            if sel.empty:
                continue
            live, offm = sel.ratio_basis_live.mean(), sel.ratio_off.mean()
            of = fitcsv.get((path, pid, cond), np.nan)
            tol = max(C5.TOL_ABS, C5.TOL_SEM * s)
            out.append(dict(path=path, target=pid, cond=cond, exp=m, SEM=s, offline_fit=of, offline=offm, BCL=live,
                            BCL_sem=sel.ratio_basis_live.sem(), d=live - offm, d_fit=live - of, d_SEM=(live - offm) / s,
                            tol=tol, PASS=bool(abs(live - offm) <= tol), rho_disagree=int(sel.rho_disagree.sum()),
                            n_syn=int(sel.n_common.sum()), dpre_maxdiff=float(sel.dpre_maxdiff.max()),
                            necb_live=float(sel.necb_live.sum()), necb_off=float(sel.necb_off.sum()),
                            z_fit2=((of - m) / s) ** 2, z_off2=((offm - m) / s) ** 2, z_live2=((live - m) / s) ** 2,
                            n_pairs=len(sel), validation_only=bool(val or f"{pid}|{cond}" in drop)))
    return pd.DataFrame(out)


def figure(o, name, out):
    """Per pathway: data (SEM), offline fit csv, offline on the same records, BCL live; last panel: live - offline per
    target with the tolerance whiskers."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    COL = {"L5": "#4c78a8", "L23": "#e45756", "L23L23": "#54a24b"}
    LAB = {"L5": "L5->L5", "L23": "L2/3->L5", "L23L23": "L2/3->L2/3"}
    paths = [p for p in COL if (o.path == p).any()]
    fig, axs = plt.subplots(len(paths) + 1, 1, figsize=(15, 3.8 * (len(paths) + 1)))
    for ax, p in zip(axs, paths):
        d = o[o.path == p].reset_index(drop=True); x = np.arange(len(d)); val = d.validation_only.values
        ax.errorbar(x, d.exp, yerr=d.SEM, fmt="_", color="k", ms=12, capsize=3, lw=1.5, label="data (SEM)")
        ax.plot(x[~val] - 0.15, d.offline_fit[~val], "o", color=COL[p], ms=6, label=f"offline fit {name}")
        ax.plot(x - 0.15, d.offline, "s", color="0.55", ms=3, label="offline, same records")
        if val.any():
            ax.plot(x[val] - 0.15, d.offline[val], "o", mfc="none", mec=COL[p], mew=1.5, ms=8, label="validation only")
        ax.errorbar(x + 0.15, d.BCL, yerr=d.BCL_sem, fmt="D", color="k", mfc="w", ms=5, lw=0.8, label="BCL live")
        ax.axhline(1, color="0.7", lw=0.8, zorder=0)
        ax.set_xticks(x); ax.set_xticklabels([C5.short(t + "|" + c) + (" [val]" if v else "")
                                              for t, c, v in zip(d.target, d.cond, val)], rotation=60, ha="right", fontsize=6)
        f = d[~d.validation_only]
        ax.set_ylabel("weight ratio"); ax.legend(fontsize=6, loc="best")
        ax.set_title(f"{LAB[p]}: chi2 BCL live {f.z_live2.sum():.2f} vs offline same records {f.z_off2.sum():.2f} over "
                     f"{len(f)} (fit csv {np.nansum(f.z_fit2):.2f}); PASS {int(d.PASS.sum())}/{len(d)}", fontsize=8)
    a = o.reset_index(drop=True); x = np.arange(len(a)); ax = axs[-1]
    ax.bar(x, a.d, color=[COL[p] for p in a.path], edgecolor=["k" if v else "none" for v in a.validation_only])
    ax.errorbar(x, np.zeros(len(a)), yerr=a.tol, fmt="none", ecolor="0.4", capsize=2, lw=0.8)
    ax.axhline(0, color="0.5", lw=0.6)
    ax.set_xticks(x); ax.set_xticklabels([f"{p} " + C5.short(t + "|" + c) for p, t, c in zip(a.path, a.target, a.cond)],
                                         rotation=60, ha="right", fontsize=5)
    ax.set_ylabel("BCL live - offline")
    ax.set_title(f"per-target difference, whiskers = max({C5.TOL_ABS}, {C5.TOL_SEM} SEM): {int(a.PASS.sum())}/{len(a)} "
                 "within", fontsize=8)
    fig.tight_layout(); fig.savefig(out + ".png", dpi=200); fig.savefig(out + ".pdf")
    print("wrote", out + ".png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True); ap.add_argument("--results", required=True)
    ap.add_argument("--save", required=True)
    ap.add_argument("--show-syn", action="store_true", help="per-synapse necb / nveto / rho table (smoke)")
    a = ap.parse_args()
    fit = json.load(open(a.fit)); name = os.path.basename(a.fit).replace(".json", "").replace("v7_", "")
    _, bx = live_v7.extra_groups(fit)
    C.BASIS["L23L23"] = os.environ.get("L23L23_BASIS_DIR") or bx or "/scratch/dhuruva/split1/basis_l23l23"
    print(f"v7 compare: {a.fit}, ecb_ref {live_v7.ecb_ref(fit)}", flush=True)
    res = [json.loads(l) for l in open(a.results)]
    df = records(fit, res, a.show_syn)
    os.makedirs(os.path.dirname(os.path.abspath(a.save)), exist_ok=True)
    df.to_csv(a.save + "_records.csv", index=False)
    ok = df.ok.fillna(False).astype(bool)
    print(f"\nrecords: {int(ok.sum())}/{len(df)} ok; failed: "
          f"{df[~ok][['path', 'pair', 'proto', 'cond', 'error']].to_string(index=False) if (~ok).any() else 'none'}")
    for path in sorted(df.path.unique()):
        d = df[ok & (df.path == path)].dropna(subset=["ratio_basis_live"])
        if len(d):
            print(f"{path}: {len(d)} records, rho disagreements {int(d.rho_disagree.sum())} / {int(d.n_common.sum())} syn, "
                  f"|d_basis| mean {d.d_basis.abs().mean():.4f} max {d.d_basis.abs().max():.4f}, "
                  f"dpre_maxdiff max {d.dpre_maxdiff.max():.4f}; eCB steps live {d.necb_live.sum():.0f} vs offline "
                  f"{d.necb_off.sum():.0f} ({int(d.necb_syn_mismatch.sum())} syn differ), vetoes live "
                  f"{d.nveto_live.sum():.0f} vs offline {d.nveto_off.sum():.0f} ({int(d.nveto_syn_mismatch.sum())} syn "
                  f"differ), nvpend max {d.nvpend_max.max():.0f}, uE maxdiff {d.uE_maxdiff.max():.2e}")
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        print(df[ok][["path", "pair", "proto", "cond", "ratio_basis_live", "ratio_off", "rho_disagree", "n_common",
                      "dpre_maxdiff", "necb_live", "necb_off", "nveto_live", "nveto_off", "nvpend_max"]]
              .round(4).to_string(index=False))
    o = targets(df, fit)
    if o.empty:
        print("no targets covered by these records"); return
    o.to_csv(a.save + "_targets.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(o.round(4).to_string(index=False))
    f = o[~o.validation_only]
    print(f"chi2 over {len(f)} fitted targets: BCL live {f.z_live2.sum():.2f}, offline same records {f.z_off2.sum():.2f}, "
          f"fit csv {np.nansum(f.z_fit2):.2f} (fit de_fun {fit.get('de_fun', float('nan')):.2f}); "
          f"PASS {int(o.PASS.sum())}/{len(o)}; max |offline - offline_fit| {np.nanmax(np.abs(o.offline - o.offline_fit)):.4f}")
    print("\n| target | data | offline fit | offline (same records) | BCL live | delta (live - offline) | verdict |")
    print("|---|---|---|---|---|---|---|")
    for _, r in o.iterrows():
        t = f"{'L2/3' if r.path == 'L23' else 'L5'} {r.target} {r.cond}" + (" [val]" if r.validation_only else "")
        print(f"| {t} | {r.exp:.2f} +/- {r.SEM:.2f} | {r.offline_fit:.3f} | {r.offline:.3f} | {r.BCL:.3f} | {r.d:+.3f} "
              f"| {'match' if r.PASS else 'MISMATCH'} (n {r.n_pairs}) |")
    figure(o, name, a.save)


if __name__ == "__main__":
    main()
