"""Live BCL (live_v4.py results jsonl) vs the offline v4 rule on the fit's extracted prefire records.

Offline per record and condition exactly as rho_redesign/eval_v4.py (C1 / C2: scan_vgate_amp.rho_rec; dpre:
model_v2.dpre_final on batch_v2 features; ratio: rho_v4.ratio on the pair's EPSP basis).

--mode equiv: one row per live run -> <save>.csv. Per synapse of the record: rho flips (live rho_GB at the end of a
  prefire run, or just before the snap of a full run, vs offline rho_f, binary at 0.5), max |dpre_live - dpre_off|,
  shadow v1 rule rhox vs the cooker's rho_obs (express cooker), t_drive 4 event count vs the record's cev (prefire),
  and three ratios: offline, basis_live (the live rho / dpre through the same basis readout) and bcl (measured
  C02 / C01 EPSPs of a full run). PASS: 0 flips and |basis_live - offline| <= 0.02 (and |bcl - offline| reported).
--mode full: target table (fit selections: L2/3 Letzkus @proximal/@distal = letzkus_distal, sjostrom @distal =
  sh_distal) -> <save>.csv (target, cond, exp, SEM, offline, BCL, BCL-offline, BCL z2), per-record <save>_records.csv,
  total BCL chi2 over the fitted targets vs the fit's de_fun, validation-only targets (the fit's drop_targets)
  separately, per-pair correlation and --fig.
"""
import argparse, json, os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); V2 = os.path.dirname(HERE); ROOT = os.path.dirname(V2)
sys.path.insert(0, V2); sys.path.insert(0, os.path.join(V2, "rho_redesign")); sys.path.insert(0, HERE)
import batch_v2                     # noqa: E402
from batch_v2 import BatchV2       # noqa: E402
import model_v2 as MV              # noqa: E402
import rho_v4                      # noqa: E402
import eval_v4                     # noqa: E402
from targets import load_targets   # noqa: E402
from live_v4 import PATHS          # noqa: E402

BASIS = {"L5": os.path.join(ROOT, "basis_results_edges_sabrina_n120_delta"),
         "L23": os.path.join(ROOT, "basis_results_edges_ebner_l23l5_delta_rs")}
GEOM = os.path.join(ROOT, "ebner", "pair_geometry_L23PC_L5TTPC.csv")


def offline(path, fit, P, pairs, protos, conds):
    """{(pair, proto, cond): dict(syn, rho0, rho_obs, rho, dpre, ratio, ncev)} and the BatchV2 (for its bases)."""
    batch_v2.BASIS_DIR = BASIS[path]
    B = BatchV2(PATHS[path]["dirs"], protocols=sorted(protos), pairs=set(pairs), fast=False, signals=("vdcc",))
    vmode, thV = rho_v4.vamp(P)
    assert vmode, "C1 / C2 fits only (eval_v4.rho_vamp)"
    rhos = eval_v4.rho_vamp(B, fit["a"], P, vmode, thV)
    feats = B.features(P); sigma = rho_v4.opts(P)[2]
    out = {}
    for i, r in enumerate(B.recs):
        b = B.basis(r); tT, K = feats[i]
        ncev = np.isfinite(r["cev"]).sum(axis=1) if r.get("cev") is not None else None
        for c in conds:
            rc = r["rho0"] if c in ("post_nmdar", "nmdar_block") else rhos[i]
            d = np.zeros(len(rc)) if c == "nmdar_block" else MV.dpre_final(
                tT, K, 0.0 if c == "mglu_block" else P["A_mglu"],
                0.0 if (c == "no_block" or (c == "post_nmdar" and P["no_drive"])) else P["A_NO"],
                P["dpre_min"], P["dpre_max"], P["dpre0"])
            out[(r["pair"], r["proto"], c)] = dict(
                syn=np.asarray(r["syn"]), rho0=r["rho0"], rho_obs=np.asarray(r["rho_obs"], float), rho=np.asarray(rc),
                dpre=np.asarray(d), ncev=ncev, b=b,
                ratio=rho_v4.ratio(b, r["rho0"], rc, d, 0.0 if c in ("post_nmdar", "nmdar_block") else sigma))
    return out


def compare(res, off, sigma):
    rows = []
    for x in res:
        row = {k: x.get(k) for k in ("path", "pair", "proto", "cond", "phase", "express", "ok", "wall_s", "error")}
        o = off.get((x["pair"], x["proto"], x["cond"]))
        if not x.get("ok") or o is None:
            row["note"] = "no live result" if not x.get("ok") else "no offline record"
            rows.append(row); continue
        st = x["end"] if x["phase"] == "prefire" else x["snap"]
        live = {int(s): i for i, s in enumerate(x["syn"])}
        idx = np.array([live.get(int(s), -1) for s in o["syn"]])
        ok = idx >= 0
        rl = np.array(st["rho_GB"])[idx[ok]]; dl = np.array(st["dpre_GB"])[idx[ok]]
        rx = np.array(st["rhox_GB"])[idx[ok]]; nv = np.array(st["nev_GB"])[idx[ok]]
        ro = o["rho"][ok]; do = o["dpre"][ok]
        rbl, rbo = rl >= 0.5, ro >= 0.5
        row.update(n_syn_live=len(x["syn"]), n_syn_rec=len(o["syn"]), n_common=int(ok.sum()),
                   n_missing_in_record=len(x.get("missing_in_record", [])), rho0_mismatch=x.get("rho0_mismatch"),
                   flips_live=int(np.sum(rbl != (o["rho0"][ok] >= 0.5))), flips_off=int(np.sum(rbo != (o["rho0"][ok] >= 0.5))),
                   rho_disagree=int(np.sum(rbl != rbo)), rho_maxdiff=float(np.max(np.abs(rl - ro))) if ok.any() else np.nan,
                   dpre_live=float(np.mean(dl)), dpre_off=float(np.mean(do)), dpre_maxdiff=float(np.max(np.abs(dl - do))),
                   ratio_off=float(o["ratio"]))
        if x["express"] == "cooker":
            row.update(rhox_disagree=int(np.sum((rx >= 0.5) != (o["rho_obs"][ok] >= 0.5))),
                       rhox_maxdiff=float(np.max(np.abs(rx - o["rho_obs"][ok]))))
        if x["phase"] == "prefire" and o["ncev"] is not None:
            row.update(nev_live=int(nv.sum()), nev_off=int(o["ncev"][ok].sum()))
        if ok.all():
            row["ratio_basis_live"] = float(rho_v4.ratio(o["b"], o["rho0"], rl, dl,
                                                         0.0 if x["cond"] in ("post_nmdar", "nmdar_block") else sigma))
        if x["phase"] == "full":
            se = x["end"]
            row.update(ratio_bcl=x.get("ratio"), ratio_bcl_last10=x.get("ratio_last"), n_post_test=x.get("n_post_test"),
                       n_post_ind=x.get("n_post_ind"),
                       rho_end_change=float(np.max(np.abs(np.array(se["rho_GB"]) -
                                                          (np.array(x["snap"]["rho_GB"]) >= 0.5)))),
                       dpre_end_change=float(np.max(np.abs(np.array(se["dpre_GB"]) - np.array(x["snap"]["dpre_GB"])))),
                       matches_written_full_config=x.get("matches_written_full_config"))
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True); ap.add_argument("--results", required=True)
    ap.add_argument("--save", required=True); ap.add_argument("--mode", choices=("equiv", "full"), default="equiv")
    ap.add_argument("--fig", default=None)
    a = ap.parse_args()
    fit = json.load(open(a.fit)); fa = fit["args"]
    batch_v2.BAP_GATE = fit.get("bap_gate")
    P = {**MV.DEFAULTS, **json.loads(fa["filters"]), **json.loads(fa.get("set", "{}")), **fit["pre"]}
    sigma = rho_v4.opts(P)[2]
    res = [json.loads(l) for l in open(a.results)]
    last = {}
    for x in res:                                   # last result per task
        last[tuple(x[k] for k in ("path", "pair", "proto", "cond", "phase", "express"))] = x
    res = list(last.values())
    dfs = []
    for path in sorted({x["path"] for x in res}):
        rp = [x for x in res if x["path"] == path]
        off = offline(path, fit, P, {x["pair"] for x in rp}, {x["proto"] for x in rp}, sorted({x["cond"] for x in rp}))
        dfs.append(compare(rp, off, sigma))
    df = pd.concat(dfs, ignore_index=True)
    os.makedirs(os.path.dirname(os.path.abspath(a.save)), exist_ok=True)
    if a.mode == "equiv":
        if "ratio_basis_live" in df:
            df["d_basis"] = df.ratio_basis_live - df.ratio_off
        if "ratio_bcl" in df:
            df["d_bcl"] = df.ratio_bcl - df.ratio_off
        nan = pd.Series(np.nan, index=df.index)
        df["PASS"] = (df.get("rho_disagree", nan) == 0) & (df.get("d_basis", nan).abs() <= 0.02)
        df.to_csv(a.save + ".csv", index=False)
        cols = [c for c in ("path", "pair", "proto", "cond", "phase", "express", "n_common", "flips_off", "flips_live",
                            "rho_disagree", "rho_maxdiff", "dpre_off", "dpre_live", "dpre_maxdiff", "rhox_disagree",
                            "nev_off", "nev_live", "ratio_off", "ratio_basis_live", "ratio_bcl", "d_basis", "d_bcl",
                            "PASS", "wall_s", "error") if c in df]
        with pd.option_context("display.width", 250, "display.max_columns", 40):
            print(df[cols].round(4).to_string(index=False))
        print(f"equivalence: {int(df.PASS.sum())}/{len(df)} runs pass (0 rho disagreements, |basis_live - offline| <= 0.02)")
        return
    full_mode(a, fit, fa, df)


def full_mode(a, fit, fa, df):
    """Target table, chi2, per-pair correlation and the figure (full validation)."""
    df = df[(df.phase == "full") & (df.express == "v4") & df.ok.fillna(False).astype(bool)].copy()
    df.to_csv(a.save + "_records.csv", index=False)
    g = pd.read_csv(GEOM); g["pair"] = g.pregid.astype(str) + "-" + g.postgid.astype(str)
    dist = dict(zip(g.pair, g.letzkus_distal)); sh = dict(zip(g.pair, g.sh_distal))
    drop = set(t for t in (fa.get("drop_targets") or "").split(",") if t)
    conds = set(fa["conditions"].split(","))
    fitcsv = {}
    for suf, path in (("", "L5"), ("_l23", "L23")):
        f = os.path.join(ROOT, fa["save"] + suf + ".csv") if not os.path.isabs(fa["save"]) else fa["save"] + suf + ".csv"
        if os.path.isfile(f):
            for _, r in pd.read_csv(f).iterrows():
                fitcsv[(path, r.target, r.condition)] = r.pred
    out = []
    for path, groups in (("L5", tuple(fa["groups"].split(","))), ("L23", ("paired_l23l5",))):
        T = {k: v for k, v in load_targets(groups).items() if k[1] in conds}
        d = df[df.path == path]
        for (pid, cond), (m, s, n, src) in sorted(T.items()):
            proto, _, where = pid.partition("@")
            sel = d[(d.proto == proto) & (d.cond == cond)].dropna(subset=["ratio_bcl"])
            if path == "L23" and where and proto.startswith("letzkus"):
                sel = sel[sel.pair.map(dist) == (where == "distal")]
            elif path == "L23" and where == "distal":
                sel = sel[sel.pair.map(sh).fillna(False).astype(bool)]
            elif where == "distal":
                continue
            if sel.empty:
                continue
            bcl, offm = sel.ratio_bcl.mean(), sel.ratio_off.mean()
            out.append(dict(path=path, target=pid, cond=cond, exp=m, SEM=s, offline=offm, offline_fit=fitcsv.get((path, pid, cond)),
                            BCL=bcl, BCL_sem=sel.ratio_bcl.sem(), basis_live=sel.ratio_basis_live.mean(),
                            **{"BCL-offline": bcl - offm}, z_off2=((offm - m) / s) ** 2, **{"BCL z2": ((bcl - m) / s) ** 2},
                            n_pairs=len(sel), validation_only=f"{pid}|{cond}" in drop))
    o = pd.DataFrame(out); o.to_csv(a.save + ".csv", index=False)
    fitted, val = o[~o.validation_only], o[o.validation_only]
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(o.round(4).to_string(index=False))
    print(f"BCL chi2 {fitted['BCL z2'].sum():.2f} over {len(fitted)} fitted targets (offline on the same pairs "
          f"{fitted.z_off2.sum():.2f}; fit de_fun {fit.get('de_fun', float('nan')):.2f})")
    if len(val):
        print("validation-only:\n" + val[["target", "cond", "exp", "SEM", "offline", "BCL", "BCL z2", "n_pairs"]].round(4).to_string(index=False))
    for path in ("L5", "L23"):
        d = df[df.path == path].dropna(subset=["ratio_bcl", "ratio_off"])
        if len(d) > 2:
            print(f"{path}: per-pair r(BCL, offline) = {np.corrcoef(d.ratio_bcl, d.ratio_off)[0, 1]:.3f}, "
                  f"r(basis_live, offline) = {np.corrcoef(d.ratio_basis_live, d.ratio_off)[0, 1]:.3f} (n {len(d)}), "
                  f"mean BCL-offline {np.mean(d.ratio_bcl - d.ratio_off):+.3f}, rho disagreements {int(d.rho_disagree.sum())}")
    if a.fig:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(17, 5.5), gridspec_kw=dict(width_ratios=[2.2, 1, 1]))
        x = np.arange(len(o))
        ax[0].errorbar(x, o.exp, o.SEM, fmt="o", color="k", label="experiment")
        ax[0].plot(x - 0.15, o.offline, "s", color="tab:blue", label="offline (same pairs)")
        ax[0].errorbar(x + 0.15, o.BCL, o.BCL_sem, fmt="D", color="tab:red", label="BCL live")
        for i in np.flatnonzero(o.validation_only.values):
            ax[0].axvspan(i - 0.5, i + 0.5, color="0.9", zorder=0)
        ax[0].set_xticks(x); ax[0].set_xticklabels([f"{p}:{t}|{c}" for p, t, c in zip(o.path, o.target, o.cond)],
                                                   rotation=90, fontsize=6)
        ax[0].axhline(1, color="0.6", lw=0.5); ax[0].set_ylabel("EPSP ratio"); ax[0].legend(fontsize=7)
        ax[0].set_title(f"BCL chi2 {fitted['BCL z2'].sum():.1f} vs fit {fit.get('de_fun', float('nan')):.1f} (grey: validation only)")
        for k, path in enumerate(("L5", "L23")):
            d = df[df.path == path].dropna(subset=["ratio_bcl", "ratio_off"])
            if d.empty:
                continue
            ax[k + 1].scatter(d.ratio_off, d.ratio_bcl, s=6, alpha=0.5)
            lo, hi = min(d.ratio_off.min(), d.ratio_bcl.min()), max(d.ratio_off.max(), d.ratio_bcl.max())
            ax[k + 1].plot([lo, hi], [lo, hi], "k--", lw=0.6)
            r = np.corrcoef(d.ratio_bcl, d.ratio_off)[0, 1] if len(d) > 2 else np.nan
            ax[k + 1].set_title(f"{path} per pair x protocol: r = {r:.3f} (n {len(d)})")
            ax[k + 1].set_xlabel("offline ratio"); ax[k + 1].set_ylabel("BCL ratio")
        fig.tight_layout(); fig.savefig(a.fig, dpi=150)
        print("wrote", a.fig)


if __name__ == "__main__":
    main()
