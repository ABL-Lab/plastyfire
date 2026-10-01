import pandas as pd, numpy as np
import sys; F=sys.argv[1]
S=pd.read_csv(f'syn_{F}.csv'); R=pd.read_csv(f'rec_{F}.csv')
S24=open('/project/rrg-emuller/dhuruva/plastyfire/glusynapse_v2/subset24_pairs.txt').read().strip().split(',')
S['d']=(S.rho_f>=.5).astype(int)-(S.rho0>=.5).astype(int); S['fl']=S.c_post<1e-3
S['c']=S.d*S.gain
k=['path','pair','proto']
A=S.groupby(k).apply(lambda x: pd.Series(dict(net=x.c.sum(), net_floor=x[x.fl].c.sum(), net_bap=x[~x.fl].c.sum(),
    up=x[x.d>0].c.sum(), dn=-x[x.d<0].c.sum(), nflip=x.d.sum())),include_groups=False).reset_index()
M=R.merge(A,on=k); M['cvf']=1+M.cv2.clip(upper=.25)
M['Rpost_chk']=(1+M.net)*M.cvf; print('check', (M.Rpost_chk-M.R_post).abs().max())
M['R_nofloor']=(1+M.net_bap)*M.cvf+(M.R-M.R_post)       # floor synapses frozen, pre part kept
M['R_nofloor_post']=(1+M.net_bap)*M.cvf
M['pre']=M.R-M.R_post
M['grp']=M.path; X=M[(M.path=='L5')&M.pair.isin(S24)].copy(); X['grp']='L5s24'; M=pd.concat([M,X])
cols=['R','R_post','pre','cvf','net','up','dn','net_floor','net_bap','nflip','R_nofloor','R_nofloor_post','R_allpot']
pd.set_option("display.width",300); pd.set_option("display.max_columns",20)
print(M.groupby(['proto','grp'])[cols].mean().round(3).T)
g=pd.read_csv('/project/rrg-emuller/dhuruva/plastyfire/ebner/pair_geometry_L23PC_L5TTPC.csv'); g['pair']=g.pregid.astype(str)+'-'+g.postgid.astype(str)
sh=dict(zip(g.pair,g.sh_distal)); ld=dict(zip(g.pair,g.letzkus_distal))
L=M[M.grp=='L23'].copy()
gs=S[S.d!=0].groupby('path').gain.mean(); k=gs['L5']/gs['L23']; print('gain per flip L23 %.3f L5 %.3f scale %.3f'%(gs['L23'],gs['L5'],k))
L['R_c']=1.013*(1+L.net*k)+L.pre; L['R_both']=1.013*(1+L.net_bap*k)+L.pre
for name,sel in [('all',L.proto!=''),('sh_distal',L.pair.map(sh).fillna(False).astype(bool)),('sh_prox',~L.pair.map(sh).fillna(False).astype(bool))]:
    print(name, L[sel].groupby('proto')[['R','R_post','pre','R_c','R_nofloor','R_both']].agg('mean').round(3).assign(n=L[sel].groupby('proto').size()).to_string())
print('rho0 mean', S.groupby('path').rho0.mean().round(3).to_dict(), 'floor frac', S.groupby('path').fl.mean().round(3).to_dict())
