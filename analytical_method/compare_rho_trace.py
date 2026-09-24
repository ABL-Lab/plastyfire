"""Decisive test: my rho(t) vs NEURON's recorded rho_GB(t), same synapse.

The traces store rho_GB, so we can compare trajectories rather than endpoints.
Then invert NEURON's own trace to recover the (dep, pot) it actually applied:

    rho' * (1e3*tau_ind) = -rho(1-rho)(rho*-rho) + pot*gp*(1-rho) - dep*gd*rho

Given recorded rho(t) the left side is measurable, so the implied drive can be
compared against what my reconstructed effcai says dep/pot should have been.

  diverges immediately  -> my effcai is wrong (amplitude/units)
  tracks then drifts    -> integration or time-above-threshold is wrong
  NEURON drive > mine   -> NEURON had pot on more than my effcai implies

    python compare_rho_trace.py [n_files]
"""
import glob, os, pickle, sys
import numpy as np
from scipy.signal import lfilter
from constants import (MIN_CA_CR, TAU_EFFCA_GB, RHO_STAR_GB, TAU_IND_GB,
                       GAMMA_D_GB, GAMMA_P_GB, TIED_IC4)

ROOT="/lustre06/project/6077694/dhuruva/plastyfire"
SIMS=f"{ROOT}/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations"
cache=pickle.load(open(f"{ROOT}/cpre_cpost_cache/ion_channels_tau278.pkl","rb"))
cp={};cq={}
for v in cache.values():
    for s,x in v["c_pre"].items():  cp[int(s)]=x
    for s,x in v["c_post"].items(): cq[int(s)]=x

nf=int(sys.argv[1]) if len(sys.argv)>1 else 2
files=sorted(glob.glob(f"{SIMS}/*/10Hz_10ms/*/simulation_traces.pkl"))[:nf]

for f in files:
    pair=f.split("/simulations/")[1].split("/")[0]
    tr=pickle.load(open(f,"rb"))
    if "rho_GB" not in tr: print("no rho_GB in", f); continue
    t=tr["t"]; dt=float(np.median(np.diff(t[:1000])))
    a=np.exp(-dt/TAU_EFFCA_GB); b=TAU_EFFCA_GB*(1-a)
    print(f"\n=== {pair}  dt={dt:.4f} ms  n={len(t)} ===")
    shown=0
    for sid,rho_rec in tr["rho_GB"].items():
        sid=int(sid)
        if sid not in cp: continue
        rho_rec=np.asarray(rho_rec,dtype=np.float64)
        eff=lfilter([0.,b],[1.,-a], np.asarray(tr["cai_CR"][sid],dtype=np.float64)-MIN_CA_CR)
        td=TIED_IC4["a00"]*cp[sid]+TIED_IC4["a01"]*cq[sid]
        tp=TIED_IC4["a10"]*cp[sid]+TIED_IC4["a11"]*cq[sid]

        # my integration on the SAME grid
        h=dt/(1e3*TAU_IND_GB); rho=rho_rec[0]; mine=np.empty_like(rho_rec)
        dep=eff>td; pot=eff>tp
        for k in range(len(eff)):
            mine[k]=rho
            rho+=h*(-rho*(1-rho)*(RHO_STAR_GB-rho)+pot[k]*GAMMA_P_GB*(1-rho)-dep[k]*(1-pot[k])*GAMMA_D_GB*rho)
            rho=min(1.,max(0.,rho))

        # invert NEURON's trace for the drive it actually applied
        drive_rec=np.gradient(rho_rec,dt)*(1e3*TAU_IND_GB)+rho_rec*(1-rho_rec)*(RHO_STAR_GB-rho_rec)
        drive_mine=pot*GAMMA_P_GB*(1-mine)-dep*GAMMA_D_GB*mine
        print(f"  syn {sid}: rho0={rho_rec[0]:.2f} rec_final={rho_rec[-1]:.4f} mine={mine[-1]:.4f}")
        print(f"    theta_d={td:.4f} theta_p={tp:.4f}  effcai max={eff.max():.4f}")
        print(f"    time above tp: mine {100*pot.mean():.3f}%   |  integral drive: rec {drive_rec.sum()*dt:.1f}  mine {drive_mine.sum()*dt:.1f}")
        shown+=1
        if shown>=3: break
