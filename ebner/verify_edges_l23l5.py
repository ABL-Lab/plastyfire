"""T29 step 1 check: new edges vs old (L5-L5 identical) and L23->L5 rho0 stats before/after."""
import h5py, numpy as np, pandas as pd, os
R='/home/dhuruva/ws/plastyfire'
P='edges/S1nonbarrel_neurons__S1nonbarrel_neurons__chemical/0/'
F=['rho0_GB','Use_d_TM','Use_p_TM','gmax_d_AMPA','gmax_p_AMPA']
a=h5py.File(f'{R}/data/dhuruva_modified_edges.h5'); b=h5py.File(f'{R}/data/dhuruva_modified_edges_l23l5.h5')
n=a[P+'rho0_GB'].shape[0]
diff=np.zeros(n,bool)
for f in F:
    for s in range(0,n,50_000_000):
        diff[s:s+50_000_000]|=a[P+f][s:s+50_000_000]!=b[P+f][s:s+50_000_000]
idx=np.where(diff)[0]; print('changed synapses',len(idx))
E_=h5py.File('/project/rrg-emuller/datasets/SSCx_plasticity/O1_2023a_Ecker/S1nonbarrel_neurons__S1nonbarrel_neurons__chemical/edges.h5')
src=E_['edges/S1nonbarrel_neurons__S1nonbarrel_neurons__chemical/source_node_id']; tgt=E_['edges/S1nonbarrel_neurons__S1nonbarrel_neurons__chemical/target_node_id']
cond=E_[P+'conductance']
S='refitting_results/fitting/n120/seed20262009/'
def pairs(d): return {tuple(map(int,x.split('-'))) for x in os.listdir(f'{R}/{S}{d}/simulations') if '-' in x}
l5=pairs('Sabrina_L5TTPC_L5TTPC_STDP'); l23=pairs('Ebner2019_L23PC_L5TTPC')
ii=idx
s_,t_=src[ii[0]:ii[-1]+1][ii-ii[0]],tgt[ii[0]:ii[-1]+1][ii-ii[0]]
ch=set(zip(s_.tolist(),t_.tolist()))
print('changed pairs',len(ch),'in l23 set',len(ch&l23),'in l5 set',len(ch&l5))
# all synapses of l23 pairs
rows=[]
for pre,post in sorted(l23):
    pass
# use the per-synapse arrays of the changed set: find all synapses of l23 pairs via target index
ix=E_['edges/S1nonbarrel_neurons__S1nonbarrel_neurons__chemical/indices/target_to_source']
n2r=ix['node_id_to_ranges'][:]; r2e=ix['range_to_edge_id'][:]
sid=[]
for pre,post in sorted(l23):
    for s,e in r2e[n2r[post][0]:n2r[post][1]]:
        m=np.where((src[s:e]==pre)&(E_[P+'syn_type_id'][s:e]>=100))[0]+s; sid+=m.tolist()
sid=np.array(sorted(sid)); print('L23->L5 synapses',len(sid))
def get(h,f): return np.array([h[P+f][i] for i in sid]) if len(sid)<5000 else None
c=np.array([cond[i] for i in sid]); ro=get(a,'rho0_GB'); rn=get(b,'rho0_GB')
print('changed among l23 synapses',int((get(a,'gmax_d_AMPA')!=get(b,'gmax_d_AMPA')).sum()+0))
for nm,r in (('before',ro),('after',rn)):
    print(nm,'frac rho0=1 %.3f'%r.mean(),'mean cond rho0=1 %.3f rho0=0 %.3f nS'%(c[r==1].mean(),c[r==0].mean()))
sl5=[];
print('L5-L5 identical: changed set has no L5 pairs ->',len(ch&l5)==0)
