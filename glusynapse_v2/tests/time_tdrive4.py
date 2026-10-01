import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); V2 = os.path.join(HERE, ".."); sys.path.insert(0, V2)
from batch_v2 import BatchV2, _sig
from targets import load_targets
DIRS = [os.path.join(V2, "extracted", d) for d in ("sj03_delta-prefire-vca", "sj03r50_delta-prefire-vca")]
T = load_targets(("paired_l5",))
t0 = time.time()
B = BatchV2(DIRS, protocols=sorted({k[0].split("@")[0] for k in T}), pairs={"180351-198084", "181455-195199"}, signals=("vdcc",))
print(f"load {time.time()-t0:.1f}s, {len(B.recs)} records", flush=True)
for td in (3, 4):
    t0 = time.time()
    for r in B.recs:
        _sig(r, dict(t_drive=td), "T")
    print(f"_sig t_drive {td}: {time.time()-t0:.1f}s", flush=True)
for r in B.recs[:3]:
    import numpy as np
    n = np.isfinite(r["cev"]).sum(1)
    print(r["proto"], "arrivals", len(r["prespikes"]), "events/syn max", n.max(), "t len", len(r["t"]), flush=True)
