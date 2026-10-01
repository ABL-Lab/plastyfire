"""Score a finished fit (fit_v2.py --save X -> X.json) once on a (usually larger) record set.

Fits run DE on a pair subset; this gives the chi2 and per-target table of the chosen parameters on all
pairs, plus the fit.py rules (theta_p > theta_d, active and potentiation-reachable fractions).

    python glusynapse_v2/eval_v2.py --fit glusynapse_v2/results/fit_v2_sub.json \
        --dirs glusynapse_v2/extracted/ebner_delta-prefire,glusynapse_v2/extracted/markram_delta-prefire-tr \
        --save glusynapse_v2/results/fit_v2_sub_allpairs
"""
import argparse, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fit_v2 as F                     # noqa: E402
from batch_v2 import BatchV2           # noqa: E402
from targets import load_targets       # noqa: E402


LOC_GROUPS = {"basal_only": (0.0, 0.0), "apical_lt_1/3": (1e-9, 1 / 3 - 1e-9), "apical_ge_1/3": (1 / 3 - 1e-9, 1.0)}


def by_loc(B, a, P, T, mode, save):
    """per-target tables on pair groups by apical synapse fraction (SONATA afferent_section_type 3 = apical)"""
    import pandas as pd
    st = np.load(os.path.join(HERE, "syn_section_type.npz")); typ = dict(zip(st["syn"].tolist(), st["section_type"]))
    frac = {r["pair"]: float(np.mean([typ[int(s)] == 3 for s in r["syn"]])) for r in B.recs}
    df = B.record_ratios(a, P, sorted({c for _, c in T}), mode)
    df["frac_apical"] = df.pair.map(frac)
    df.to_csv(save + "_records.csv", index=False)
    rows = []
    for g, (lo, hi) in LOC_GROUPS.items():
        sub = df[(df.frac_apical >= lo) & (df.frac_apical <= hi)]
        for (pid, cond), (m, s, n, src) in T.items():
            if pid.endswith("@distal"):
                continue        # as objective(): no tuft inputs in our L5-L5 pairs
            sel = sub[(sub.proto == pid.split("@")[0]) & (sub.condition == cond)].dropna(subset=["ratio"])
            if len(sel):
                pred = sel.ratio.mean()
                rows.append(dict(group=g, target=pid, condition=cond, target_mean=m, target_sem=s, pred=pred,
                                 pred_sem=sel.ratio.sem(), n_pairs=len(sel), z=(pred - m) / s, dpre=sel.dpre.mean()))
    out = pd.DataFrame(rows); out.to_csv(save + "_byloc.csv", index=False)
    print(out.pivot_table(index=["target", "condition"], columns="group", values="pred").round(2)
          .join(out.groupby(["target", "condition"]).target_mean.first()).to_string())
    for g, d in out.groupby("group"):
        print(f"{g}: chi2 {(d.z ** 2).sum():.2f} over {len(d)} targets, {d.n_pairs.max()} pairs")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True, help="json written by fit_v2.py")
    ap.add_argument("--dirs", required=True)
    ap.add_argument("--pairs", default=None, help="default: all pairs in --dirs")
    ap.add_argument("--save", required=True)
    ap.add_argument("--by-loc", action="store_true",
                    help="also score pair groups by apical synapse fraction (syn_section_type.npz): X_byloc.csv")
    args = ap.parse_args()
    fit = json.load(open(args.fit)); fa = fit["args"]
    conds = set(fa["conditions"].split(","))
    T = {k: v for k, v in load_targets(tuple(fa["groups"].split(","))).items() if k[1] in conds}
    fil = json.loads(fa["filters"])
    sig = () if fit["model"] == "post" else (("vdcc",) if fil.get("pre_drive") else ("shaft_cai",))
    if fil.get("no_drive") == 2:
        sig += ("cacr",)
    B = BatchV2(args.dirs.split(","), protocols=sorted({k[0].split("@")[0] for k in T}),
                pairs=set(args.pairs.split(",")) if args.pairs else None, fast=fa["rho"] == "fast", signals=sig)
    have = {r["proto"] for r in B.recs}
    T = {k: v for k, v in T.items() if k[0].split("@")[0] in have}
    P = {**json.loads(fa["filters"]), **fit["pre"]} if fit["model"] == "v2" else {}
    chi2, table = B.objective(fit["a"], P, T, mode="full" if fa["rho"] == "full" else "fast")
    td = np.concatenate([B.thetas(r, fit["a"])[0] for r in B.recs])
    tp = np.concatenate([B.thetas(r, fit["a"])[1] for r in B.recs])
    peak = np.concatenate([r["peak"] for r in B.recs])
    rule = dict(theta_p_gt_theta_d=bool(np.all(tp > td)), active=float((td < peak).mean()),
                pot_reachable=float((tp < peak).mean()), ok=bool(np.all(tp > td) and (td < peak).mean()
                >= F.MIN_ACTIVE and (tp < peak).mean() >= F.MIN_POT))
    table.to_csv(args.save + ".csv", index=False)
    out = dict(fit=args.fit, model=fit["model"], chi2=chi2, chi2_fit=fit["chi2"], n_targets=len(table),
               n_records=len(B.recs), n_pairs=len({r["pair"] for r in B.recs}), rules=rule)
    json.dump(out, open(args.save + ".json", "w"), indent=1, default=float)
    print(table.round(3).to_string(index=False))
    if args.by_loc:
        by_loc(B, fit["a"], P, T, "full" if fa["rho"] == "full" else "fast", args.save)
    print(f"chi2 {chi2:.2f} on {out['n_pairs']} pairs (fit subset: {fit['chi2']:.2f}); rules {rule}")


if __name__ == "__main__":
    main()
