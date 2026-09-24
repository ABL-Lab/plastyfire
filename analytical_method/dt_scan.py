"""Which rho-integration dt is the right accuracy/cost tradeoff?

No re-extraction needed: extracted_d0p25 (0.25 ms) strides down to any coarser
dt, so one load covers the whole scan. Accuracy is measured against the
bluecellulab rho of the tied-ic4 run; cost is the wall time of one evaluation,
which is what the optimiser inner loop actually pays.

    ANALYTICAL_EXTRACTED=extracted_d0p25 python dt_scan.py
"""
import time
import numpy as np
from constants import RHO_STAR_GB, TAU_IND_GB, GAMMA_D_GB, GAMMA_P_GB, TIED_IC4
from batch import Batch

MEASURED = {-50: 1.2848, -30: 1.1850, -10: 1.3329, 5: 1.6434,
            10: 1.6224, 30: 1.4710, 50: 1.3440}


def integrate(effT, dt_ms, td, tp, rho0):
    h = (dt_ms / 1000.0) / TAU_IND_GB
    rho = rho0.astype(np.float64).copy()
    for k in range(effT.shape[0]):
        e = effT[k]; dep = e > td; pot = e > tp
        rho += h * (-rho*(1-rho)*(RHO_STAR_GB-rho)
                    + pot*GAMMA_P_GB*(1-rho) - dep*(1-pot)*GAMMA_D_GB*rho)
        np.clip(rho, 0., 1., out=rho)
    return rho


if __name__ == "__main__":
    b = Batch()
    del b.effcaiT                      # free ~3 GB; rebuilt per stride below
    td, tp = b.thetas(TIED_IC4)
    obs = b.rho_obs
    base = b.dt_ms
    print(f"{b.N} synapses, base dt {base:.3f} ms, T={b.T}\n")
    print(f"{'dt (ms)':>8} {'T':>8} {'max|drho|':>10} {'mean|drho|':>11} "
          f"{'binary%':>8} {'flips':>6} {'maxcurve':>9} {'eval s':>8} {'RAM GB':>7}")

    rows = []
    for stride in [1, 2, 4, 8, 16, 32, 64]:
        dt = base * stride
        effT = np.ascontiguousarray(b.effcai[:, ::stride].T)
        ram = effT.nbytes / 2**30
        t0 = time.time()
        rho = integrate(effT, dt, td, tp, b.rho0)
        secs = time.time() - t0
        flips = int(((rho >= .5) != (obs >= .5)).sum())
        summ, _ = b.curve(rho)
        m = dict(zip(summ.dt, summ["mean"]))
        mc = max(abs(m[d] - MEASURED[d]) for d in MEASURED if d in m)
        print(f"{dt:>8.2f} {effT.shape[0]:>8} {np.abs(rho-obs).max():>10.6f} "
              f"{np.abs(rho-obs).mean():>11.6f} "
              f"{100*((rho>=.5)==(obs>=.5)).mean():>7.2f}% {flips:>6} "
              f"{mc:>9.5f} {secs:>8.1f} {ram:>7.2f}", flush=True)
        rows.append((dt, flips, mc, secs))
        del effT

    print("\n=== tradeoff ===")
    best = min(rows, key=lambda r: r[2])
    print(f"  most accurate : dt={best[0]:.2f} ms  maxcurve {best[2]:.5f}  {best[3]:.1f} s/eval")
    for dt, fl, mc, s in rows:
        # DE: 80 generations, 48 candidates, 60 cores -> ~80 sequential evals
        print(f"  dt={dt:>6.2f} ms : maxcurve {mc:.5f}  flips {fl:>2}  "
              f"DE(80 gen) ~ {80*s/60:.0f} min")
