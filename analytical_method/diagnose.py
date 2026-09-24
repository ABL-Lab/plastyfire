"""Why is the offline curve biased -0.30? Three candidate causes, tested.

1. bias direction   - is predicted rho systematically BELOW observed?
2. who disagrees    - are the wrong-side synapses mostly pred-low (missed
                      potentiation) or pred-high?
3. decimation       - does the disagreement concentrate on synapses whose
                      effcai only BARELY crosses theta_p? Those are exactly the
                      ones a 2 ms sample can miss and CVODE root-finding cannot.
"""
import numpy as np
from constants import TIED_IC4
from batch import Batch

b = Batch()
td, tp = b.thetas(TIED_IC4)
rho = b.rho_full(td, tp)
obs = b.rho_obs

print(f"=== 1. bias direction (n={b.N}) ===")
print(f"  mean pred {rho.mean():.4f}   mean obs {obs.mean():.4f}   bias {rho.mean()-obs.mean():+.4f}")
print(f"  pred<obs on {100*(rho<obs).mean():.1f}% of synapses")
print(f"  by starting state:")
for r0 in (0.0, 1.0):
    m = b.rho0 == r0
    if m.sum(): print(f"    rho0={r0:.0f} n={m.sum():4d}  pred {rho[m].mean():.3f}  obs {obs[m].mean():.3f}  bias {rho[m].mean()-obs[m].mean():+.3f}")

pb, ob = rho >= 0.5, obs >= 0.5
print(f"\n=== 2. binary disagreement ===")
print(f"  agree                : {100*(pb==ob).mean():.1f}%")
print(f"  pred LOW  (obs 1, pred 0): {(ob & ~pb).sum():4d}  <- missed potentiation")
print(f"  pred HIGH (obs 0, pred 1): {(~ob & pb).sum():4d}")

print(f"\n=== 3. is disagreement concentrated near threshold? ===")
peak = b.effcai.max(axis=1)
marg = peak / tp                       # >1 means peak clears theta_p
bad  = pb != ob
for lo, hi, lab in [(0,0.9,'peak well BELOW tp'),(0.9,1.1,'peak NEAR tp'),(1.1,99,'peak well ABOVE tp')]:
    m = (marg>=lo)&(marg<hi)
    if m.sum(): print(f"  {lab:22s} n={m.sum():4d}   disagree {100*bad[m].mean():5.1f}%")
print("\n  -> if 'NEAR tp' dominates, the cause is sampling resolution, not the model.")

frac = (b.effcai > tp[:,None]).mean(axis=1)
print(f"\n  time above theta_p: median {np.median(frac)*100:.2f}% of the trace")
print(f"  synapses with <1% time above theta_p: {(frac<0.01).sum()} ({100*(frac<0.01).mean():.0f}%)")
