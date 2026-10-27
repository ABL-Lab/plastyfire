"""T27: collect full-protocol BCL (simulation_full_<hash>.pkl) EPSP ratios of the L23PC->L5TTPC Ebner workdirs, merge with
pair geometry and the offline prediction, and summarise per target (selections as in eval_cooker_l23l5.py).
    python ebner/eval_full_l23l5.py [--hash delta-cooker-full]"""
import argparse, csv, glob, os, pickle
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser()
ap.add_argument("--hash", default="delta-cooker-full")
ap.add_argument("--out", default=os.path.join(ROOT, "ebner", "full_delta-cooker_l23l5"))
args = ap.parse_args()
SIMS = os.path.join(ROOT, "refitting_results/fitting/n120/seed20262009/Ebner2019_L23PC_L5TTPC/simulations")
rows, bad = [], []
for f in sorted(glob.glob(f"{SIMS}/*/*/simulation_full_{args.hash}.pkl")):
    proto = os.path.basename(os.path.dirname(f)); pair = os.path.basename(os.path.dirname(os.path.dirname(f)))
    r = pickle.load(open(f, "rb"))
    if r.get("epsp_ratio") is None:
        bad.append((pair, proto, r.get("epsp_error"))); continue
    sr = r.get("snapped_rho")
    rows.append(dict(pair=pair, proto=proto, ratio_full=r["epsp_ratio"], epsp_c01_mV=float(np.mean(r["epsp_c01"])),
                     epsp_c02_mV=float(np.mean(r["epsp_c02"])), n_c01=len(r["epsp_c01"]), n_c02=len(r["epsp_c02"]),
                     rho_snapped_mean=float(np.mean(sr)) if sr is not None and len(sr) else np.nan,
                     guardrail_forced=r.get("guardrail_forced", False)))
full = pd.DataFrame(rows)
print(len(full), "with ratio;", len(bad), "unusable", bad[:5])
geo = pd.read_csv(os.path.join(ROOT, "ebner/pair_geometry_L23PC_L5TTPC.csv"))
geo["pair"] = geo.pregid.astype(str) + "-" + geo.postgid.astype(str)
off = pd.read_csv(os.path.join(ROOT, "ebner/cooker_l23l5_offline_records.csv"))[["pair", "proto", "ratio", "rho"]]\
    .rename(columns={"ratio": "ratio_offline", "rho": "rho_offline"})
df = full.merge(geo[["pair", "pregid", "postgid", "all_protocols", "rise_10_90_ms", "letzkus_distal", "sh_distal",
                     "epsp_mV", "path_mean"]], on="pair", how="left").merge(off, on=["pair", "proto"], how="left")
df.to_csv(args.out + "_records.csv", index=False)
print(df.groupby("proto")[["ratio_full", "ratio_offline"]].agg(["mean", "count"]).round(3))
for p, g in df.groupby("proto"):
    g = g.dropna(subset=["ratio_offline"])
    print(p, "per-pair corr full vs offline: %.2f (n=%d)" % (np.corrcoef(g.ratio_full, g.ratio_offline)[0, 1], len(g)))

off_sum = pd.read_csv(os.path.join(ROOT, "ebner/cooker_l23l5_offline.csv"))
out = []
def add(target, pid, sel, mask, data, sem, n_exp, loc):
    for scope in ("all_protocols", "any"):
        sub = df[(df.proto == pid)]
        sub = sub[mask(sub) & (sub.all_protocols.astype(bool) if scope == "all_protocols" else True)]
        o = off_sum[(off_sum.target == f"{pid}@{loc}") & (off_sum.selection == sel) & (off_sum.scope == scope)]
        out.append(dict(target=f"{pid}@{loc}", selection=sel, scope=scope, data=data, data_sem=sem, n_exp=n_exp,
                        offline=o.pred.iloc[0] if len(o) else np.nan,
                        offline_n=o.n_pairs.iloc[0] if len(o) else np.nan,
                        full=sub.ratio_full.mean() if len(sub) else np.nan,
                        full_sem=sub.ratio_full.sem() if len(sub) > 1 else np.nan, n_full=len(sub),
                        frac_ltp=(sub.ratio_full > 1.05).mean(), frac_ltd=(sub.ratio_full < 0.95).mean(),
                        offline_same_pairs=sub.ratio_offline.mean() if len(sub) else np.nan,
                        baseline_epsp_mV=sub.epsp_c01_mV.mean() if len(sub) else np.nan))
for r in csv.DictReader(open(os.path.join(ROOT, "ebner/ebner_targets.csv"))):
    pid, loc = r["protocol_id"], (r["synapse_location"] or "all").split()[0]
    m, s, n = float(r["mean_ratio"]), float(r["sem"]), r["n"]
    if pid.startswith("letzkus_1ap") or (pid.startswith("letzkus_") and loc == "pooled"):
        add(pid, pid, "all pairs", lambda d: d.ratio_full.notna(), m, s, n, loc)
    elif pid.startswith("letzkus_"):
        if loc == "distal":
            add(pid, pid, "rise>2.7ms", lambda d: d.letzkus_distal.astype(bool), m, s, n, loc)
        else:
            add(pid, pid, "rise<=2.7ms", lambda d: ~d.letzkus_distal.astype(bool), m, s, n, loc)
    elif pid == "sjostrom_50hz_dt+10ms" and loc == "distal":
        add(pid, pid, "S&H distal (path>200um, EPSP<1mV)", lambda d: d.sh_distal.astype(bool), m, s, n, loc)
# Sjostrom & Hausser 2006 Results: all unitary L2/3->L5 pairs 1.06 +- 0.09 (n=19)
pid = "sjostrom_50hz_dt+10ms"
for scope in ("all_protocols", "any"):
    sub = df[df.proto == pid]; sub = sub[sub.all_protocols.astype(bool)] if scope == "all_protocols" else sub
    a = pd.read_csv(os.path.join(ROOT, "ebner/cooker_l23l5_offline_records.csv")); a = a[a.proto == pid]
    a = a[a.all_protocols.astype(bool)] if scope == "all_protocols" else a
    out.append(dict(target=f"{pid}@all_L23L5_unitary", selection="all pairs (sjh06_02)", scope=scope, data=1.06,
                    data_sem=0.09, n_exp=19, offline=a.ratio.mean(), offline_n=len(a), full=sub.ratio_full.mean(),
                    full_sem=sub.ratio_full.sem(), n_full=len(sub), frac_ltp=(sub.ratio_full > 1.05).mean(),
                    frac_ltd=(sub.ratio_full < 0.95).mean(), offline_same_pairs=sub.ratio_offline.mean(),
                    baseline_epsp_mV=sub.epsp_c01_mV.mean()))
res = pd.DataFrame(out); res.to_csv(args.out + ".csv", index=False)
pd.set_option("display.width", 220)
print(res[["target", "selection", "scope", "data", "data_sem", "offline", "full", "full_sem", "n_full",
           "offline_same_pairs", "frac_ltp", "frac_ltd"]].round(3).to_string(index=False))
