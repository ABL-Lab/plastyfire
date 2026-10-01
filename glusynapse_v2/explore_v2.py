"""Score hand-made parameter sets on the GPU and print per-target predictions (model exploration).

    python glusynapse_v2/explore_v2.py --base glusynapse_v2/results/reduced_gpu_preview.json \
        --variants '[{"name": "no NO", "pre": {"A_NO": 0}}, {"name": "v2.4", "filters": {"no_drive": 2}, "pre": {...}}]'

Each variant starts from the base fit (a and pre), overrides "a" / "pre" entries, and may set fixed
"filters" (a new GPUModel is built per distinct filter set). GPU node only (jax).
"""
import argparse, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from batch_v2 import BatchV2
from targets import load_targets

PREVIEW = "180351-198084,181015-184976,182339-200396,182381-186941,182900-195902,183240-189246,183769-207759,186261-208490"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--variants", required=True, help="json list or @file")
    ap.add_argument("--dirs", default=",".join(os.path.join(HERE, "extracted", d) for d in ("ebner_preview", "markram_delta-cooker")))
    ap.add_argument("--pairs", default=PREVIEW)
    ap.add_argument("--show", default="", help="comma list of target substrings to print (default: |z| > 1.5)")
    args = ap.parse_args()
    import jax_v2 as JV
    base = json.load(open(args.base))
    V = json.load(open(args.variants[1:])) if args.variants.startswith("@") else json.loads(args.variants)
    T = load_targets(("markram", "nevian", "ebner"))
    B = BatchV2(args.dirs.split(","), protocols=sorted({k[0].split("@")[0] for k in T}), pairs=set(args.pairs.split(",")),
                fast=False, signals=("vdcc",))
    have = {r["proto"] for r in B.recs}; T = {k: v for k, v in T.items() if k[0].split("@")[0] in have}
    f0 = json.loads(base["args"]["filters"])
    groups = {}
    for v in V:
        f = {**f0, **v.get("filters", {})}
        groups.setdefault(json.dumps(f, sort_keys=True), []).append(v)
    res = []
    for fk, vs in groups.items():
        G = JV.GPUModel(B, T, json.loads(fk))
        A = [{**base["a"], **v.get("a", {})} for v in vs]
        for a in A:
            a.update(a20=a["a00"], a21=a["a01"], a30=a["a10"], a31=a["a11"])
        Ps = [{**base["pre"], **json.loads(fk), **v.get("pre", {})} for v in vs]
        rho, d = G.rho_dpre(A, Ps)
        chi, preds = G.chi2(rho, d, return_pred=True)
        for i, v in enumerate(vs):
            res.append((v["name"], chi[i], {k: p[i] for k, p in preds.items()}))
        del G
    show = [s for s in args.show.split(",") if s]
    keys = [k for k, _, _ in [(k, 0, 0) for k in T] if any(s in k[0] for s in show)] if show else None
    for name, c, p in res:
        print(f"\n=== {name}: chi2 {c:.2f}")
        for k, pr in p.items():
            m, s = T[k][0], T[k][1]; z = (pr - m) / s
            if (keys is None and abs(z) > 1.5) or (keys is not None and k in keys):
                print(f"  {k[0]:34s} {k[1]:10s} target {m:5.2f}  pred {pr:5.2f}  z {z:+5.2f}")


if __name__ == "__main__":
    main()
