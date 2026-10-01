import glob,pickle,collections
S="refitting_results/fitting/n120/seed20262009/Ebner2019_L23PC_L5TTPC/simulations"
c=collections.Counter(); bad=[]
for f in glob.glob(S+"/*/*/simulation_edges_delta-prefire-vseg-rs.pkl"):
    r=pickle.load(open(f,"rb")); k=(bool(r.get("guardrail_forced")), r["n_post_spikes"]==r["n_post_spikes_exp"]); c[k]+=1
    if k!=(False,True): bad.append(f.split("/")[-3]+"__"+f.split("/")[-2])
print(dict(c)); print(len(bad)); open("/tmp/claude-3117582/-lustre09-project-6070394-dhuruva/c979f8f7-556a-4034-950d-1a47fb951115/scratchpad/l23bad.txt","w").write("\n".join(bad))
