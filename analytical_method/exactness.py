"""How close can the offline model get to bluecellulab? Test each error source.

NEURON runs cvode=True: adaptive-step BDF with root-finding on the WATCH
conditions. We use fixed-step explicit Euler with crossings quantised to the
sample grid. Two independent error sources:

  (1) integrator truncation. With pot=1 one step moves rho by h*gp*(1-rho)
      ~5.7e-3 of the gap, so Euler error ~1.6e-5/step over ~600 active steps
      ~ 1e-2 -- the same order as the observed max error.
  (2) crossing placement. NEURON puts the theta crossing exactly; we snap it to
      a 2 ms sample, so time-above-threshold is quantised.

Variants:
  euler        explicit Euler                      (current)
  euler_sub    explicit Euler, N substeps/interval
  expint       exponential integrator: linear part solved exactly, cubic added
               explicitly. Exact when dep/pot are constant over the interval.
  expint_frac  expint + linear interpolation of the theta crossing inside each
               interval, so dep/pot become fractional weights in [0,1] --
               emulating cvode's root-finding.

    python exactness.py                  # 2 ms data
    ANALYTICAL_EXTRACTED=extracted_d0p25 python exactness.py
"""
import numpy as np
from constants import (RHO_STAR_GB, TAU_IND_GB, GAMMA_D_GB, GAMMA_P_GB, TIED_IC4)
from batch import Batch

MEASURED = {-50: 1.2848, -30: 1.1850, -10: 1.3329, 5: 1.6434,
            10: 1.6224, 30: 1.4710, 50: 1.3440}
T_ODE = 1e3 * TAU_IND_GB


def _cubic(rho):
    return -rho * (1 - rho) * (RHO_STAR_GB - rho) / T_ODE


def euler(b, td, tp, sub=1):
    h = (b.dt_ms / 1000.0) / TAU_IND_GB / sub
    rho = b.rho0.astype(np.float64).copy()
    E = b.effcaiT
    for k in range(b.T):
        e = E[k]; dep = e > td; pot = e > tp
        for _ in range(sub):
            rho += h * (-rho*(1-rho)*(RHO_STAR_GB-rho)
                        + pot*GAMMA_P_GB*(1-rho) - dep*(1-pot)*GAMMA_D_GB*rho)
            np.clip(rho, 0., 1., out=rho)
    return rho


def expint(b, td, tp, frac=False):
    """rho' = (a - c*rho)/T + cubic.  Linear part integrated exactly."""
    dt = b.dt_ms
    rho = b.rho0.astype(np.float64).copy()
    E = b.effcaiT
    prev = E[0]
    for k in range(b.T):
        e = E[k]
        if frac:
            # fraction of [k-1,k] spent above each threshold, by linear interp
            wp = _frac_above(prev, e, tp)
            wd = _frac_above(prev, e, td)
            prev = e
        else:
            wp = (e > tp).astype(np.float64)
            wd = (e > td).astype(np.float64)
        a = wp * GAMMA_P_GB
        c = wp * GAMMA_P_GB + wd * (1 - wp) * GAMMA_D_GB
        s = c * dt / T_ODE
        dec = np.exp(-s)
        eq = np.where(c > 0, a / np.maximum(c, 1e-300), 0.0)
        rho = np.where(c > 0, rho * dec + eq * (1 - dec), rho + a * dt / T_ODE)
        rho += dt * _cubic(rho)
        np.clip(rho, 0., 1., out=rho)
    return rho


def _frac_above(e0, e1, th):
    """Fraction of a linear segment e0->e1 lying above th."""
    a0, a1 = e0 > th, e1 > th
    both = (a0 & a1).astype(np.float64)
    nei  = (~a0 & ~a1)
    d = e1 - e0
    safe = np.where(np.abs(d) < 1e-30, 1e-30, d)
    x = (th - e0) / safe                       # crossing position in [0,1]
    x = np.clip(x, 0.0, 1.0)
    part = np.where(a1, 1.0 - x, x)            # rising: above after x
    return np.where(both, 1.0, np.where(nei, 0.0, part))


if __name__ == "__main__":
    b = Batch()
    td, tp = b.thetas(TIED_IC4)
    obs = b.rho_obs
    print(f"{b.N} synapses, dt={b.dt_ms:.3f} ms, T={b.T}\n")
    print(f"{'variant':>14} {'max|drho|':>10} {'binary%':>9} {'flips':>6} {'maxcurve':>9}")
    for name, fn in [("euler",       lambda: euler(b, td, tp, 1)),
                     ("euler_sub4",  lambda: euler(b, td, tp, 4)),
                     ("expint",      lambda: expint(b, td, tp, False)),
                     ("expint_frac", lambda: expint(b, td, tp, True))]:
        rho = fn()
        flips = int(((rho >= .5) != (obs >= .5)).sum())
        summ, _ = b.curve(rho)
        m = dict(zip(summ.dt, summ["mean"]))
        mc = max(abs(m[d] - MEASURED[d]) for d in MEASURED if d in m)
        print(f"{name:>14} {np.abs(rho-obs).max():>10.6f} "
              f"{100*((rho>=.5)==(obs>=.5)).mean():>8.2f}% {flips:>6} {mc:>9.5f}")
