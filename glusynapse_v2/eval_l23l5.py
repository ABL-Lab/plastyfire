"""Score a finished L5-L5 fit (fit_v2.py json) on the unitary L2/3->L5 TTPC records: targets.py paired_l23l5.

No refit: same parameters, new pathway (transfer test, L23_PATHWAYS.md). Per-target pair selection from
ebner/pair_geometry_L23PC_L5TTPC.csv (BatchV2.objective takes one loc map; these targets need two):
  letzkus_*@proximal / @distal  somatic uEPSP 10-90% rise < / > 2.7 ms (letzkus_distal; --split median uses
                                letzkus_distal_median)
  sjostrom_50hz_dt+10ms@distal  path_mean > 200 um and EPSP < 1 mV (sh_distal, S&H cooperativity)
  no suffix                     all pairs
Pairs: all_protocols pairs only (default) so every target sees the same cells.

    ANALYTICAL_BASIS_DIR=basis_results_edges_ebner_l23l5_delta python glusynapse_v2/eval_l23l5.py \
        --fit glusynapse_v2/results/reduced_gpu_subset_pl5_td2_s1.json --save glusynapse_v2/results/l23l5_td2_s1
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
os.environ.setdefault("ANALYTICAL_BASIS_DIR", os.path.join(ROOT, "basis_results_edges_ebner_l23l5_delta"))
sys.path.insert(0, HERE)
import batch_v2                         # noqa: E402
from batch_v2 import BatchV2           # noqa: E402
from targets import load_targets       # noqa: E402

GEOM = os.path.join(ROOT, "ebner", "pair_geometry_L23PC_L5TTPC.csv")
DIR = os.path.join(HERE, "extracted", "ebner_l23l5_delta-prefire")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True)
    ap.add_argument("--dirs", default=DIR)
    ap.add_argument("--save", required=True)
    ap.add_argument("--split", choices=["rise", "median"], default="rise")
    ap.add_argument("--all-pairs", action="store_true", help="also pairs missing a protocol")
    a_ = ap.parse_args()
    fit = json.load(open(a_.fit)); fa = fit["args"]
    batch_v2.BAP_GATE = fit.get("bap_gate")   # score with the gate the fit was made with
    T = load_targets(("paired_l23l5",))
    T = {k: v for k, v in T.items() if k[1] in set(fa["conditions"].split(","))}
    g = pd.read_csv(GEOM); g["pair"] = g.pregid.astype(str) + "-" + g.postgid.astype(str)
    if not a_.all_pairs:
        g = g[g.all_protocols]
    dist = dict(zip(g.pair, g["letzkus_distal" if a_.split == "rise" else "letzkus_distal_median"]))
    sh = dict(zip(g.pair, g.sh_distal))
    fil = json.loads(fa["filters"])
    sig = () if fit["model"] == "post" else (("vdcc",) if fil.get("pre_drive") else ("shaft_cai",))
    if fil.get("no_drive") == 2:
        sig += ("cacr",)
    B = BatchV2(a_.dirs.split(","), protocols=sorted({k[0].split("@")[0] for k in T}), pairs=set(g.pair),
                fast=fa.get("rho", "full") == "fast", signals=sig)
    P = {**fil, **fit["pre"]} if fit["model"] == "v2" else {}
    df = B.record_ratios(fit["a"], P, sorted({c for _, c in T}), "full" if fa.get("rho", "full") == "full" else "fast")
    df.to_csv(a_.save + "_records.csv", index=False)
    rows = []
    for (pid, cond), (m, s, n, src) in sorted(T.items()):
        proto, _, where = pid.partition("@")
        sel = df[(df.proto == proto) & (df.condition == cond)].dropna(subset=["ratio"])
        if where and proto.startswith("letzkus"):
            sel = sel[sel.pair.map(dist) == (where == "distal")]
        elif where == "distal":
            sel = sel[sel.pair.map(sh).fillna(False).astype(bool)]
        if sel.empty:
            continue
        pred = sel.ratio.mean()
        rows.append(dict(target=pid, condition=cond, target_mean=m, target_sem=s, n_exp=n, pred=pred,
                         pred_sem=sel.ratio.sem(), n_pairs=len(sel), z=(pred - m) / s, dpre=sel.dpre.mean(),
                         rho=sel.rho.mean()))
    out = pd.DataFrame(rows); out.to_csv(a_.save + ".csv", index=False)
    chi2 = float((out.z ** 2).sum())
    json.dump(dict(fit=a_.fit, chi2=chi2, n_targets=len(out), split=a_.split, n_pairs=int(df.pair.nunique()),
                   n_records=len(B.recs)), open(a_.save + ".json", "w"), indent=1)
    print(out.round(3).to_string(index=False))
    print(f"chi2 {chi2:.2f} over {len(out)} targets, {df.pair.nunique()} pairs ({a_.split} split)")


if __name__ == "__main__":
    main()
