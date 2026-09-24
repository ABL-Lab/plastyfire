"""STAGE 1: is my reconstructed effcai_GB right?

effcai_GB is not recorded, but it is recoverable. dep/pot are binary, so the
drive implied by the recorded rho_GB takes one of three values:

    drive = rho'*(1e3*tau_ind) + rho(1-rho)(rho*-rho)
          =  0                        dep=0 pot=0   (effcai < theta_d)
          = -gd*rho                   dep=1 pot=0   (theta_d < effcai < theta_p)
          =  gp(1-rho) - gd*rho       dep=1 pot=1   (effcai > theta_p)

Classifying each sample recovers NEURON's true threshold-crossing state, which
is a direct measurement of its effcai relative to theta_d/theta_p. Comparing
that against mine isolates stage 1 from the rho integration entirely.

rho moves ~1e-5 per ms, so differentiate on a 2 ms grid, not at 0.025 ms where
the step is near float32 resolution.

    python stage1_effcai.py [n_files]
"""
import glob, pickle, sys
import numpy as np
from scipy.signal import lfilter
from constants import (MIN_CA_CR, TAU_EFFCA_GB, RHO_STAR_GB, TAU_IND_GB,
                       GAMMA_D_GB, GAMMA_P_GB, TIED_IC4)

ROOT = "/lustre06/project/6077694/dhuruva/plastyfire"
SIMS = f"{ROOT}/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations"
GRID = 2.0

cache = pickle.load(open(f"{ROOT}/cpre_cpost_cache/ion_channels_tau278.pkl", "rb"))
cp, cq = {}, {}
for v in cache.values():
    for s, x in v["c_pre"].items():  cp[int(s)] = x
    for s, x in v["c_post"].items(): cq[int(s)] = x

nf = int(sys.argv[1]) if len(sys.argv) > 1 else 3
agg = []
for f in sorted(glob.glob(f"{SIMS}/*/10Hz_10ms/*/simulation_traces.pkl"))[:nf]:
    tr = pickle.load(open(f, "rb"))
    t = tr["t"]; dt = float(np.median(np.diff(t[:1000]))); step = int(round(GRID/dt))
    a = np.exp(-dt/TAU_EFFCA_GB); b = TAU_EFFCA_GB*(1-a)
    for sid, rr in tr["rho_GB"].items():
        sid = int(sid)
        if sid not in cp: continue
        rho = np.asarray(rr, dtype=np.float64)[::step]
        eff = lfilter([0., b], [1., -a],
                      np.asarray(tr["cai_CR"][sid], dtype=np.float64)-MIN_CA_CR)[::step]
        td = TIED_IC4["a00"]*cp[sid] + TIED_IC4["a01"]*cq[sid]
        tp = TIED_IC4["a10"]*cp[sid] + TIED_IC4["a11"]*cq[sid]

        drive = np.gradient(rho, GRID)*(1e3*TAU_IND_GB) + rho*(1-rho)*(RHO_STAR_GB-rho)
        # dep is gated off while pot is on: dep*(1-pot). So the third state is
        # potentiation ALONE, not both together.
        c = np.stack([np.zeros_like(rho), -GAMMA_D_GB*rho,
                      GAMMA_P_GB*(1-rho)])
        state = np.argmin(np.abs(c - drive), axis=0)      # 0,1,2
        # states are only identifiable while rho is moving; near rho=0 with pot=0
        # every candidate collapses to ~0, so restrict to samples with signal
        live = np.abs(c[2]-c[0]) > 1.0
        n_pot, n_dep = (state == 2), (state >= 1)
        m_pot, m_dep = (eff > tp), (eff > td)
        if live.sum() < 100: continue
        agg.append(dict(
            syn=sid, td=td, tp=tp,
            neuron_pot=100*n_pot[live].mean(), mine_pot=100*m_pot[live].mean(),
            neuron_dep=100*n_dep[live].mean(), mine_dep=100*m_dep[live].mean(),
            pot_agree=100*(n_pot[live] == m_pot[live]).mean(),
            # where NEURON says pot but I say not: how far below tp is my effcai?
            ratio_missed=float(np.median(eff[live & n_pot & ~m_pot]/tp))
                         if (live & n_pot & ~m_pot).sum() else np.nan,
            effmax_over_tp=float(eff.max()/tp)))

import pandas as pd
d = pd.DataFrame(agg)
print(f"=== STAGE 1: effcai threshold-crossing, {len(d)} synapses ===\n")
print(f"  % time pot ON   NEURON {d.neuron_pot.mean():6.2f}%   mine {d.mine_pot.mean():6.2f}%")
print(f"  % time dep ON   NEURON {d.neuron_dep.mean():6.2f}%   mine {d.mine_dep.mean():6.2f}%")
print(f"  pot state agreement (per sample): {d.pot_agree.mean():.1f}%")
print(f"\n  where NEURON pot=1 but mine=0, my effcai/theta_p median: "
      f"{np.nanmedian(d.ratio_missed):.3f}")
print(f"  (1.0 would mean my effcai only just misses; <<1 means it is far too low)")
print(f"\n  my effcai max / theta_p: median {d.effmax_over_tp.median():.3f}")
r = d.neuron_pot/d.mine_pot.replace(0, np.nan)
print(f"\n  NEURON pot-time / my pot-time: median {r.median():.2f}x")
