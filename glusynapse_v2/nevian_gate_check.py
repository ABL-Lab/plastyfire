"""Analytic check for NEVIAN_EXTRACELLULAR.md: t_drive 2 eCB gate (paired_l5 td2 fits s1-s3) and the fast
post-AP trace S_M (tau 20 ms) at the pre arrival, on Nevian and Sjostrom spike timings. Seconds on the login node."""
import numpy as np
P={1:(0.5695,8.2072,0.2787,0.0058),2:(0.1348,7.136,3.6027,0.0055),3:(0.3114,8.803,2.9228,0.0098)}
dt=0.01; t=np.arange(-400,5,dt)
def gate(aps,th,tT,tg,tau=100.):
    S=np.zeros_like(t)
    for a in aps: S+= (t>=a)*np.exp(-(t-a)/tau)
    u=np.maximum(S-th,0); T=np.zeros_like(t)
    for k in range(1,len(t)): T[k]=T[k-1]+dt*(-T[k-1]/tT+u[k-1])
    i=np.searchsorted(t,0.0); return np.tanh(max(T[i]-tg,0))
def SM(aps,tau=20.): return sum(np.exp(a/tau) for a in aps if a<0)
b=lambda first,n,f: [first+k*1000./f for k in range(n)]
prot={ # post AP times relative to pre spike (EPSP onset) at 0; name: (aps, target)
 "N 3AP50 dt'-90":(b(-90,3,50),1.00),"N 3AP50 -10(dt'-50)":(b(-50,3,50),0.68),"N 3AP50 dt'-30":(b(-30,3,50),0.98),
 "N 3AP50 dt'-10":(b(-10,3,50),1.42),"N 3AP50 +10":(b(10,3,50),2.01),"N 1AP -10":([-10],0.80),"N 1AP +10":([10],1.04),
 "N 2AP50 -10":(b(-30,2,50),0.72),"N 3AP20 -10":(b(-110,3,20),0.72),"N 3AP100 -10":(b(-30,3,100),0.52),
 "Sj 1AP -10":([-10],None),"Sj 1AP -25":([-25],None),"Sj 1AP -50":([-50],None),"Sj 1AP -100":([-100],None),"Sj 1AP -120":([-120],None)}
print("%-22s %6s | %s | %s"%("protocol","target","gate s1 s2 s3 (t_drive2 fit)","S_M(tau20) at arrival"))
for k,(aps,tg) in prot.items():
    g=[gate(aps,*P[s][:3]) for s in (1,2,3)]
    print("%-22s %6s | %s | %.2f"%(k,tg,"  ".join("%.3f"%x for x in g),SM(aps)))
