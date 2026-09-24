"""Look at ONE synapse, one repetition, sample by sample.

Aggregate statistics have given contradictory answers (my effcai crosses
threshold MORE often, yet my rho ends LOWER). Both cannot hold, so stop
aggregating and read the actual traces.

Prints, over one 4 s repetition:
    cai_CR, my effcai, theta_d/theta_p, NEURON rho, my rho
plus where each crosses, so the divergence is visible rather than inferred.

    python inspect_one.py [file_index] [syn_index]
"""
import glob, pickle, sys
import numpy as np
from scipy.signal import lfilter
from constants import (MIN_CA_CR, TAU_EFFCA_GB, RHO_STAR_GB, TAU_IND_GB,
                       GAMMA_D_GB, GAMMA_P_GB, TIED_IC4)

ROOT="/lustre06/project/6077694/dhuruva/plastyfire"
SIMS=f"{ROOT}/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations"
fi = int(sys.argv[1]) if len(sys.argv)>1 else 0
si = int(sys.argv[2]) if len(sys.argv)>2 else 0

cache=pickle.load(open(f"{ROOT}/cpre_cpost_cache/ion_channels_tau278.pkl","rb"))
cp={};cq={}
for v in cache.values():
    for s,x in v["c_pre"].items():  cp[int(s)]=x
    for s,x in v["c_post"].items(): cq[int(s)]=x

f=sorted(glob.glob(f"{SIMS}/*/10Hz_10ms/*/simulation_traces.pkl"))[fi]
print("file:", f.split("/simulations/")[1])
tr=pickle.load(open(f,"rb"))
t=tr["t"]; dt=float(np.median(np.diff(t[:1000])))
sids=[int(s) for s in tr["rho_GB"] if int(s) in cp]
sid=sids[si]
print(f"synapse {sid}   dt={dt} ms   n={len(t)}")

rho_n=np.asarray(tr["rho_GB"][sid],dtype=np.float64)
cai  =np.asarray(tr["cai_CR"][sid],dtype=np.float64)
a=np.exp(-dt/TAU_EFFCA_GB); b=TAU_EFFCA_GB*(1-a)
eff=lfilter([0.,b],[1.,-a],cai-MIN_CA_CR)
td=TIED_IC4["a00"]*cp[sid]+TIED_IC4["a01"]*cq[sid]
tp=TIED_IC4["a10"]*cp[sid]+TIED_IC4["a11"]*cq[sid]

h=dt/(1e3*TAU_IND_GB); rho=rho_n[0]; mine=np.empty_like(rho_n)
for k in range(len(eff)):
    mine[k]=rho
    rho+=h*(-rho*(1-rho)*(RHO_STAR_GB-rho)+(eff[k]>tp)*GAMMA_P_GB*(1-rho)-(eff[k]>td)*(not (eff[k]>tp))*GAMMA_D_GB*rho)
    rho=min(1.,max(0.,rho))

print(f"theta_d={td:.5f}  theta_p={tp:.5f}   c_pre={cp[sid]:.5f} c_post={cq[sid]:.6f}")
print(f"cai_CR   : rest {cai[:100].mean():.3e}  max {cai.max():.3e}")
print(f"my effcai: max {eff.max():.5f}   ({eff.max()/tp:.2f} x theta_p)")
print(f"rho      : NEURON {rho_n[0]:.4f} -> {rho_n[-1]:.4f}    mine {mine[0]:.4f} -> {mine[-1]:.4f}")

print("\n=== trajectory at 10 checkpoints ===")
print(f"{'t (ms)':>9} {'NEURON rho':>11} {'my rho':>9} {'gap':>8} {'my effcai':>10} {'>td':>4} {'>tp':>4}")
for frac in np.linspace(0.05,1.0,10):
    k=int(frac*(len(t)-1))
    print(f"{t[k]:>9.0f} {rho_n[k]:>11.4f} {mine[k]:>9.4f} {mine[k]-rho_n[k]:>+8.4f} "
          f"{eff[k]:>10.5f} {str(eff[k]>td):>4} {str(eff[k]>tp):>4}")

print("\n=== first repetition, 950-1600 ms, every 25 ms ===")
print(f"{'t':>7} {'cai_CR':>11} {'my effcai':>10} {'>td':>4} {'>tp':>4} {'NEURON rho':>11} {'my rho':>9}")
for ms in range(950,1601,25):
    k=int(ms/dt)
    if k>=len(t): break
    print(f"{ms:>7} {cai[k]:>11.3e} {eff[k]:>10.5f} {str(eff[k]>td):>4} {str(eff[k]>tp):>4} "
          f"{rho_n[k]:>11.5f} {mine[k]:>9.5f}")
