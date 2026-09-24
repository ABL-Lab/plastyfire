"""Rebuild cpre_cpost_cache/ion_channels_tau278.pkl from the extracted npz files.

extract.py copied c_pre/c_post out of the cache into every npz, so the cache is
recoverable exactly -- no re-simulation. c_pre/c_post are per-PAIR (protocol
independent), so the 7 protocols of a pair must agree; this asserts that rather
than trusting it.

Only valid while tau_effca is unchanged (278.3177658387), since that is the one
fitted parameter that shapes effcai and therefore c_pre/c_post.

    python rebuild_cache.py --out ../cpre_cpost_cache/ion_channels_tau278.pkl
"""
import argparse, glob, os, pickle
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

ap = argparse.ArgumentParser()
ap.add_argument("--extracted", default=os.path.join(HERE, "extracted"))
ap.add_argument("--out", default=os.path.join(HERE, "..", "cpre_cpost_cache",
                                             "ion_channels_tau278.pkl"))
a = ap.parse_args()

cache, seen = {}, {}
for f in sorted(glob.glob(os.path.join(a.extracted, "*.npz"))):
    pair = os.path.basename(f)[:-4].split("__")[0]
    pre, post = (int(x) for x in pair.split("-"))
    d = np.load(f)
    cp = {int(s): float(v) for s, v in zip(d["syn"], d["c_pre"])}
    cq = {int(s): float(v) for s, v in zip(d["syn"], d["c_post"])}
    key = (pre, post)
    if key in cache:                      # cross-check against another protocol
        for s in cp:
            assert abs(cache[key]["c_pre"][s]  - cp[s]) < 1e-12, f"c_pre disagrees {key} {s}"
            assert abs(cache[key]["c_post"][s] - cq[s]) < 1e-12, f"c_post disagrees {key} {s}"
        seen[key] += 1
    else:
        cache[key] = {"c_pre": cp, "c_post": cq}
        seen[key] = 1

n_syn = sum(len(v["c_pre"]) for v in cache.values())
print(f"pairs      : {len(cache)}")
print(f"synapses   : {n_syn}")
print(f"protocols/pair: min {min(seen.values())} max {max(seen.values())}")
allq = np.array([x for v in cache.values() for x in v["c_post"].values()])
allp = np.array([x for v in cache.values() for x in v["c_pre"].values()])
print(f"c_pre  med {np.median(allp):.4e}   c_post med {np.median(allq):.4e}")
print(f"c_post at floor(<1e-4): {(allq<1e-4).sum()}/{allq.size} "
      f"({100*(allq<1e-4).mean():.1f}%)")

out = os.path.abspath(a.out)
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, "wb") as f:
    pickle.dump(cache, f, -1)
print(f"\nwrote {out}  ({os.path.getsize(out)} bytes)")
