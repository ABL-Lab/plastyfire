import pickle, glob, numpy as np, pandas as pd
C = {"L5": "sabrina_n120_split1", "L23L5": "ebner_l23l5_split1_rs", "L23L23": "zilberter_l23l23_split1_rs"}
for p, c in C.items():
    cache = pickle.load(open(f"/scratch/dhuruva/split1/cache/{c}.pkl", "rb"))
    for f in glob.glob(f"/scratch/dhuruva/split1/cexp/parts/{p}/*.csv"):
        d = pd.read_csv(f); k = (int(d.pre_gid[0]), int(d.post_gid[0])); e = cache[k]
        rp = d.cpre / d.syn_id.map(e["c_pre"]); rq = d.cpost / d.syn_id.map(e["c_post"])
        print(p, k, "cpre/cache", rp.round(5).tolist(), "cpost/cache", rq.round(5).tolist())
        print("  apv<pre<mg0:", bool(((d.cpre_apv < d.cpre) & (d.cpre < d.cpre_mg0)).all()),
              "nmda_q_apv max", d.nmda_q_apv.abs().max(), "nmda_q med", d.nmda_q.median(), "vdcc_q_apv/vdcc_q", (d.vdcc_q_apv / d.vdcc_q).round(2).tolist())
        print(d[["loc", "dist_um", "cpre", "cpre_apv", "cpre_mg0", "cpost", "cpre_cai", "cpre_apv_cai", "cpre_mg0_cai", "cpost_cai"]].round(4).to_string())
