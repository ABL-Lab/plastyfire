"""Why does the fit switch the NO pathway off? At the preview v2 fit, scan the NO parameters and print
mean dpre and ratio per protocol (which protocols NO potentiates, and by how much)."""
import json, os, sys, itertools
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from batch_v2 import BatchV2
PAIRS = "180351-198084,181015-184976,182339-200396,182381-186941,182900-195902,183240-189246,183769-207759,186261-208490".split(",")
fit = json.load(open(os.path.join(HERE, "results", "preview_fit_v2.json")))
a, P0 = fit["a"], fit["pre"]
B = BatchV2([os.path.join(HERE, "extracted", d) for d in ("ebner_preview", "markram_delta-cooker")], pairs=set(PAIRS))
grid = [dict(theta_NO=th, tau_NO=tn, tau_Z=tz) for th, tn, tz in itertools.product((5e-4, 1e-3, 2e-3, 4e-3), (5.0, 20.0), (30.0, 76.25))]
B.precompute([{**P0, **g} for g in grid], workers=16)
# NO-only K statistic: per protocol, mean total K (ms) per synapse -> how much NO overlap each protocol produces
rows = []
for g in grid:
    P = {**P0, **g}
    for (tT, K), r in zip(B.features(P), B.recs):
        rows.append(dict(**g, proto=r["proto"], K=float(K.sum(1).mean())))
k = pd.DataFrame(rows).groupby(["theta_NO", "tau_NO", "tau_Z", "proto"])["K"].mean().unstack("proto")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
print("mean total NO overlap K (ms) per synapse, by protocol:")
print(k.round(1).T.to_string())
