"""Cooker (delta-cooker a-params, post-only rule) on the Ebner protocols: an early look at T12.

    python glusynapse_v2/preview_eval.py --dirs glusynapse_v2/extracted/ebner_preview
"""
import argparse, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from batch_v2 import BatchV2           # noqa: E402
from targets import load_targets       # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--dirs", required=True)
ap.add_argument("--a", default=os.path.join(os.path.dirname(HERE), "fit_results", "delta-cooker.json"))
args = ap.parse_args()
a = json.load(open(args.a))["params"]
B = BatchV2(args.dirs.split(","))
T = load_targets()
chi2, res = B.objective(a, {}, T, mode="full")
res = res.sort_values("target")
print(res[["target", "condition", "target_mean", "target_sem", "pred", "pred_sem", "n_pairs", "z", "rho"]]
      .round(3).to_string(index=False))
print(f"chi2 {chi2:.1f} over {len(res)} targets (cooker, post-only, delta-cooker a-params)")
out = os.path.join(HERE, "results", "preview_cooker_ebner.csv"); res.to_csv(out, index=False); print("->", out)
# shaft Ca statistics per protocol (for the v2 threshold boxes)
rows = {}
for r in B.recs:
    rows.setdefault(r["proto"], []).append(np.percentile(r["shaft_cai"].max(1), [50, 90]))
print("\nshaft Ca peak per synapse (uM), median / 90th pct over synapses, averaged over pairs:")
for p, v in sorted(rows.items()):
    v = np.array(v).mean(0) * 1e3
    print(f"  {p:28s} {v[0]:7.3f} {v[1]:7.3f}")
