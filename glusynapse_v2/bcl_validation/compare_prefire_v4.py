"""Prefire BCL validation of a v4 fit: the live GluSynapseV4 rule (live_v4.py --from-fit: prefire runs, cooker
expression so the induction reproduces the fit's own traces, V4_BIN=0 = the continuous live rule) vs the offline fit.

Per record (compare_v4.compare): the live rho_GB / dpre_GB at the end of the prefire run through the pair's EPSP basis
(ratio_basis_live, the fit's own readout) vs the offline ratio (compare_v4.offline = eval_v4), rho disagreements.
Per target (the fit's selections, as compare_v4.full_mode): data mean / SEM, offline_fit (the fit's csv pred, over all
its pairs), offline (the same records as the live runs), BCL live (mean ratio_basis_live), d = live - offline,
d_fit = live - offline_fit. PASS: |d| <= max(TOL_ABS, TOL_SEM * SEM).
Writes <save>_targets.csv, <save>_records.csv, <save>.png / .pdf (fig8_joint_fit style) and prints the markdown
verdict table (target | data | offline fit | BCL live | delta).

    python glusynapse_v2/bcl_validation/compare_prefire_v4.py --fit <json> --results <jsonl> --save <prefix>
"""
import argparse, json, os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare_v4 as C                                   # noqa: E402
from compare_v4 import MV, rho_v4, batch_v2, load_targets, GEOM, ROOT   # noqa: E402

TOL_ABS, TOL_SEM = 0.02, 0.25


def records(fit, res):
    fa = fit["args"]
    batch_v2.BAP_GATE = fit.get("bap_gate")
    P = {**MV.DEFAULTS, **json.loads(fa["filters"]), **json.loads(fa.get("set", "{}")), **fit["pre"]}
    sigma = rho_v4.opts(P)[2]
    last = {}
    for x in res:                                         # last result per task
        last[tuple(x[k] for k in ("path", "pair", "proto", "cond", "phase", "express"))] = x
    res = [x for x in last.values() if x["phase"] == "prefire" and x["express"] == "cooker"]
    dfs = []
    for path in sorted({x["path"] for x in res}):
        rp = [x for x in res if x["path"] == path]
        off = C.offline(path, fit, P, {x["pair"] for x in rp}, {x["proto"] for x in rp}, sorted({x["cond"] for x in rp}))
        dfs.append(C.compare(rp, off, sigma))
    df = pd.concat(dfs, ignore_index=True)
    for c in ("ratio_basis_live", "ratio_off", "rho_disagree", "n_common", "dpre_maxdiff"):
        if c not in df:
            df[c] = np.nan
    df["d_basis"] = df.ratio_basis_live - df.ratio_off
    return df


def targets(df, fit):
    fa = fit["args"]
    g = pd.read_csv(GEOM); g["pair"] = g.pregid.astype(str) + "-" + g.postgid.astype(str)
    dist = dict(zip(g.pair, g.letzkus_distal)); sh = dict(zip(g.pair, g.sh_distal))
    drop = set(t for t in (fa.get("drop_targets") or "").split(",") if t)
    conds = set(fa["conditions"].split(","))
    fitcsv = {}
    for suf, path in (("", "L5"), ("_l23", "L23")):
        f = fa["save"] + suf + ".csv"
        f = f if os.path.isabs(f) else os.path.join(ROOT, f)
        if os.path.isfile(f):
            for _, r in pd.read_csv(f).iterrows():
                fitcsv[(path, r.target, r.condition)] = r.pred
    d0 = df[df.ok.fillna(False).astype(bool)].dropna(subset=["ratio_basis_live", "ratio_off"])
    out = []
    groups = [("L5", tuple(fa["groups"].split(",")))] + ([("L23", ("paired_l23l5",))] if fa.get("joint") else [])
    for path, grp in groups:
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
            tol = max(TOL_ABS, TOL_SEM * s)
            out.append(dict(path=path, target=pid, cond=cond, exp=m, SEM=s, offline_fit=of, offline=offm, BCL=live,
                            BCL_sem=sel.ratio_basis_live.sem(), d=live - offm, d_fit=live - of, d_SEM=(live - offm) / s,
                            tol=tol, PASS=bool(abs(live - offm) <= tol), rho_disagree=int(sel.rho_disagree.sum()),
                            n_syn=int(sel.n_common.sum()), dpre_maxdiff=float(sel.dpre_maxdiff.max()),
                            z_fit2=((of - m) / s) ** 2, z_off2=((offm - m) / s) ** 2, z_live2=((live - m) / s) ** 2,
                            n_pairs=len(sel), validation_only=f"{pid}|{cond}" in drop))
    return pd.DataFrame(out)


def short(k):
    return (k.replace("sjostrom07_step200ms_", "S07 ").replace("sjostrom_", "")
             .replace("|control", "").replace("|", " "))


def figure(o, df, fit, name, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    COL = {"L5": "#4c78a8", "L23": "#e45756"}; CLIVE, COFF = "k", "0.55"
    fig = plt.figure(figsize=(14, 13))
    gs = fig.add_gridspec(3, 6, height_ratios=[1.2, 1.0, 1.0], hspace=1.0, wspace=1.0)
    axes = {"L5": fig.add_subplot(gs[0, :]), "L23": fig.add_subplot(gs[1, 0:2])}
    sc = {"L5": fig.add_subplot(gs[1, 2:4]), "L23": fig.add_subplot(gs[1, 4:6])}
    axd = fig.add_subplot(gs[2, :])
    for path, ax in axes.items():
        d = o[o.path == path].reset_index(drop=True)
        if d.empty:
            ax.set_visible(False); continue
        x = np.arange(len(d)); keys = d.target + "|" + d.cond; val = d.validation_only.values
        ax.errorbar(x, d.exp, yerr=d.SEM, fmt="_", color="k", ms=12, capsize=3, lw=1.5, zorder=3, label="data (SEM)")
        ok = ~val
        ax.plot(x[ok] - 0.15, d.offline_fit[ok], "o", color=COL[path], ms=6, zorder=4, label=f"offline fit {name}")
        ax.plot(x - 0.15, d.offline, "s", color=COFF, ms=3, zorder=4, label="offline, same records")
        if val.any():
            ax.plot(x[val] - 0.15, d.offline[val], "o", mfc="none", mec=COL[path], mew=1.5, ms=8, zorder=4,
                    label="validation only (not fitted)")
        ax.errorbar(x + 0.15, d.BCL, yerr=d.BCL_sem, fmt="D", color=CLIVE, mfc="w", ms=5, lw=0.8, zorder=5,
                    label="BCL live (prefire)")
        ax.axhline(1, color="0.7", lw=0.8, zorder=0)
        ax.set_xticks(x); ax.set_xticklabels([short(k) + (" [val]" if v else "") for k, v in zip(keys, val)],
                                             rotation=60, ha="right", fontsize=6)
        f = d[~d.validation_only]
        ax.set_ylabel("weight ratio"); ax.legend(fontsize=6, loc="best")
        ax.set_title(f"({'A' if path == 'L5' else 'B'}) {'L5->L5' if path == 'L5' else 'L2/3->L5'}: chi2 BCL live "
                     f"{f.z_live2.sum():.2f} vs offline same records {f.z_off2.sum():.2f} over {len(f)} "
                     f"(fit csv {np.nansum(f.z_fit2):.2f})", fontsize=8)
    d0 = df[df.ok.fillna(False).astype(bool)].dropna(subset=["ratio_basis_live", "ratio_off"])
    for path, ax in sc.items():
        d = d0[d0.path == path]
        if d.empty:
            ax.set_visible(False); continue
        for c, mk in zip(sorted(d.cond.unique()), "osD^v"):
            e = d[d.cond == c]
            ax.scatter(e.ratio_off, e.ratio_basis_live, s=10, marker=mk, alpha=0.6, color=COL[path] if c == "control" else None,
                       label=c)
        lo = min(d.ratio_off.min(), d.ratio_basis_live.min()); hi = max(d.ratio_off.max(), d.ratio_basis_live.max())
        ax.plot([lo, hi], [lo, hi], "k--", lw=0.6)
        r = np.corrcoef(d.ratio_off, d.ratio_basis_live)[0, 1] if len(d) > 2 else np.nan
        ax.set_title(f"{'L5->L5' if path == 'L5' else 'L2/3->L5'} per record: r {r:.4f}, max |d| "
                     f"{d.d_basis.abs().max():.3f}, rho disagree {int(d.rho_disagree.sum())} (n {len(d)})", fontsize=7)
        ax.set_xlabel("offline ratio"); ax.set_ylabel("BCL live ratio"); ax.legend(fontsize=6)
    a = o.reset_index(drop=True); x = np.arange(len(a))
    axd.bar(x, a.d, color=[COL[p] for p in a.path], edgecolor=["k" if v else "none" for v in a.validation_only])
    axd.errorbar(x, np.zeros(len(a)), yerr=a.tol, fmt="none", ecolor="0.4", capsize=2, lw=0.8)
    axd.axhline(0, color="0.5", lw=0.6)
    axd.set_xticks(x); axd.set_xticklabels([("L23 " if p == "L23" else "") + short(t + "|" + c) + (" [val]" if v else "")
                                            for p, t, c, v in zip(a.path, a.target, a.cond, a.validation_only)],
                                           rotation=60, ha="right", fontsize=6)
    axd.set_ylabel("BCL live - offline (same records)")
    axd.set_title(f"(C) per-target difference, whiskers = tolerance max({TOL_ABS}, {TOL_SEM} SEM): "
                  f"{int(a.PASS.sum())}/{len(a)} within", fontsize=8)
    fig.subplots_adjust(top=0.96, bottom=0.08, left=0.06, right=0.98)
    fig.savefig(out + ".png", dpi=200); fig.savefig(out + ".pdf")
    print("wrote", out + ".png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True); ap.add_argument("--results", required=True)
    ap.add_argument("--save", required=True)
    a = ap.parse_args()
    fit = json.load(open(a.fit)); name = os.path.basename(a.fit).replace(".json", "").replace("v4_", "")
    res = [json.loads(l) for l in open(a.results)]
    df = records(fit, res)
    os.makedirs(os.path.dirname(os.path.abspath(a.save)), exist_ok=True)
    df.to_csv(a.save + "_records.csv", index=False)
    ok = df.ok.fillna(False).astype(bool)
    print(f"records: {int(ok.sum())}/{len(df)} ok; failed: {df[~ok][['path', 'pair', 'proto', 'cond', 'error']].to_string(index=False) if (~ok).any() else 'none'}")
    for path in sorted(df.path.unique()):
        d = df[ok & (df.path == path)].dropna(subset=["ratio_basis_live"])
        if len(d):
            print(f"{path}: {len(d)} records, rho disagreements {int(d.rho_disagree.sum())} / {int(d.n_common.sum())} syn, "
                  f"|d_basis| mean {d.d_basis.abs().mean():.4f} max {d.d_basis.abs().max():.4f}, "
                  f"dpre_maxdiff max {d.dpre_maxdiff.max():.4f}, mean d_basis {d.d_basis.mean():+.4f}")
    o = targets(df, fit)
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
    figure(o, df, fit, name, a.save)


if __name__ == "__main__":
    main()
