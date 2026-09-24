"""Where does the last ~0.01 of EPSP-ratio error come from?

Hypothesis: rho agrees to ~1e-3, but EPSP thresholds at rho>=0.5, so only the
handful of synapses sitting within ~1e-3 of 0.5 matter. Check that the flipped
synapses are exactly those, and that the affected pairs are the ones whose dt
does not reproduce exactly.
"""
import numpy as np
from constants import TIED_IC4
from batch import Batch

b = Batch(); td, tp = b.thetas(TIED_IC4)
rho = b.rho_full(td, tp); obs = b.rho_obs
flip = (rho >= .5) != (obs >= .5)

print(f"synapses           : {b.N}")
print(f"binary flips       : {flip.sum()}  ({100*flip.mean():.3f}%)")
print(f"mean |rho err|     : {np.abs(rho-obs).mean():.6f}")
print(f"max  |rho err|     : {np.abs(rho-obs).max():.6f}")

print(f"\n--- the flipped synapses ---")
print(f"{'obs rho':>10} {'pred rho':>10} {'|err|':>9} {'dist obs->0.5':>14}")
for i in np.where(flip)[0]:
    print(f"{obs[i]:>10.5f} {rho[i]:>10.5f} {abs(rho[i]-obs[i]):>9.5f} {abs(obs[i]-0.5):>14.5f}")

d = np.abs(obs - 0.5)
print(f"\n--- how close to the 0.5 cut? ---")
for lo, hi in [(0,.005),(.005,.02),(.02,.1),(.1,1)]:
    m = (d>=lo)&(d<hi)
    if m.sum(): print(f"  |obs-0.5| in [{lo:.3f},{hi:.3f}): n={m.sum():5d}  flipped {flip[m].sum():3d} ({100*flip[m].mean():5.1f}%)")

# which (pair,protocol) records contain a flip -> those are the dt that miss
off = 0; bad = {}
for m in b.meta:
    sl = slice(off, off+m["n"]); off += m["n"]
    if flip[sl].any(): bad.setdefault(m["dt"], []).append(m["pair"])
print(f"\n--- dt affected by >=1 flip ---")
for dt in sorted(bad): print(f"  dt={dt:>+4}: {len(bad[dt])} pair(s)  {bad[dt][:4]}")
print(f"  dt with NO flips reproduce the measured ratio exactly.")
