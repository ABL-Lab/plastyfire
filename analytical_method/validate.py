"""Gate: does the offline model reproduce the bluecellulab run?

Run with the SAME a-params the BCL sims used (tied ic4) and compare against the
observed final rho and the measured STDP curve. If this fails, nothing built on
top of it means anything.

    python validate.py
"""
import argparse, os, sys
import numpy as np, pandas as pd
from constants import (TIED_IC4, DEFIT2, DEFIT2_GAMMA, GAMMA_D_GB, GAMMA_P_GB,
                       RHO_BINARY_THRESHOLD)
from batch import Batch

# (a-params, (gamma_d, gamma_p), curve measured by plot_stdp_ic_best.py)
REFS = {
    # tied run (hash cdf3a1e1db98), n100
    "tied_ic4": (TIED_IC4, (GAMMA_D_GB, GAMMA_P_GB),
                 {-50: 1.2848, -30: 1.1850, -10: 1.3329, 5: 1.6434,
                  10: 1.6224, 30: 1.4710, 50: 1.3440}),
    # DE fit #2 (hash b8c7ff3ecf0a), n120 Sabrina set, 806 runs
    "defit2_n120": (DEFIT2, DEFIT2_GAMMA,
                    {-50: 0.8522, -30: 0.8247, -10: 0.8698, 5: 1.1027,
                     10: 1.0496, 30: 0.9554, 50: 0.9044}),
}

def report(tag, rho, b):
    err = np.abs(rho - b.rho_obs)
    obs_b = (b.rho_obs >= RHO_BINARY_THRESHOLD)
    prd_b = (rho       >= RHO_BINARY_THRESHOLD)
    print(f"\n--- {tag} ---")
    print(f"  synapses            : {len(rho)}")
    print(f"  corr(pred, obs)     : {np.corrcoef(rho, b.rho_obs)[0,1]:.4f}")
    print(f"  mean |pred-obs|     : {err.mean():.4f}")
    print(f"  BINARY state agrees : {100*(obs_b==prd_b).mean():.1f}%   <-- what EPSP uses")
    return prd_b, obs_b

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", choices=sorted(REFS), default="tied_ic4")
    ref_a, (gd, gp), MEASURED = REFS[ap.parse_args().ref]
    stride = int(os.environ.get("ANALYTICAL_STRIDE", "1"))
    b = Batch(stride=stride)
    print(f"loaded {b.N} synapses over {len(b.meta)} (pair,protocol) records  "
          f"dt={b.dt_ms:.3f} ms")
    td, tp = b.thetas(ref_a)
    print(f"theta_d med {np.median(td):.4f}   theta_p med {np.median(tp):.4f}")

    rho_fu = b.rho_full(td, tp, gd, gp); report("full integration", rho_fu, b)
    rho_fa = b.rho_fast(td, tp, gd, gp); report("mean-field (fast)", rho_fa, b)
    print(f"\n  fast vs full: corr {np.corrcoef(rho_fa, rho_fu)[0,1]:.4f}  "
          f"mean|d| {np.abs(rho_fa-rho_fu).mean():.4f}  "
          f"binary agree {100*((rho_fa>=.5)==(rho_fu>=.5)).mean():.1f}%")

    print("\n=== STDP curve: offline vs bluecellulab ===")
    print(f"  {'dt':>5} {'measured':>10} {'full':>10} {'fast':>10}")
    sf, _ = b.curve(rho_fu); sa, _ = b.curve(rho_fa)
    mf = dict(zip(sf.dt, sf["mean"])); ma = dict(zip(sa.dt, sa["mean"]))
    for dt in sorted(MEASURED):
        print(f"  {dt:>5} {MEASURED[dt]:>10.4f} {mf.get(dt,float('nan')):>10.4f} "
              f"{ma.get(dt,float('nan')):>10.4f}")
    e = np.array([mf[dt]-MEASURED[dt] for dt in sorted(MEASURED) if dt in mf])
    print(f"\n  full-integration curve error: mean {e.mean():+.4f}  max|.| {np.abs(e).max():.4f}")
