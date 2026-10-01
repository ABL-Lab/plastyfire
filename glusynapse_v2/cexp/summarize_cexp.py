import pickle, numpy as np, pandas as pd
O = "/scratch/dhuruva/split1/cexp"
C = {"L5L5": "sabrina_n120_split1", "L23L5": "ebner_l23l5_split1_rs", "L23L23": "zilberter_l23l23_split1_rs"}
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 50)
q = lambda x: f"{x.median():.4g} [{x.quantile(.25):.3g}-{x.quantile(.75):.3g}]"
for p, c in C.items():
    d = pd.read_csv(f"{O}/{p}.csv"); cache = pickle.load(open(f"/scratch/dhuruva/split1/cache/{c}.pkl", "rb"))
    d["nmda_share"] = 1 - d.cpre_apv / d.cpre
    d["r_apv"] = d.cpre_apv / d.cpre; d["r_mg0"] = d.cpre_mg0 / d.cpre; d["r_post"] = d.cpost / d.cpre
    d["r_apv_cai"] = d.cpre_apv_cai / d.cpre_cai; d["r_mg0_cai"] = d.cpre_mg0_cai / d.cpre_cai; d["r_post_cai"] = d.cpost_cai / d.cpre_cai
    k = list(zip(d.pre_gid, d.post_gid))
    d["cache_cpre"] = [cache[kk]["c_pre"].get(s, np.nan) if kk in cache else np.nan for kk, s in zip(k, d.syn_id)]
    d["cache_cpost"] = [cache[kk]["c_post"].get(s, np.nan) if kk in cache else np.nan for kk, s in zip(k, d.syn_id)]
    d["off_pre"] = d.cpre / d.cache_cpre; d["off_post"] = d.cpost / d.cache_cpost
    print(f"\n######## {p}: {len(d)} syn, {d.groupby(['pre_gid','post_gid']).ngroups} pairs; cpost NaN {d.cpost.isna().mean():.3f}; "
          f"pairs not in cache {int(d.cache_cpre.isna().sum())} syn; error col: {d.get('error', pd.Series(dtype=str)).fillna('').value_counts().to_dict()}")
    cols = ["cpre", "cpre_apv", "cpre_mg0", "cpost", "cpre_cai", "cpre_apv_cai", "cpre_mg0_cai", "cpost_cai", "vdcc_q", "nmda_q", "nmda_q_mg0",
            "nmda_share", "r_apv", "r_mg0", "r_post", "r_apv_cai", "r_mg0_cai", "r_post_cai", "dist_um"]
    for loc, g in [("all", d)] + list(d.groupby("loc")):
        print(f"-- {p} {loc} n={len(g)}  median [IQR]")
        for cn in cols:
            print(f"   {cn:14s} {q(g[cn].dropna())}")
    print(f"-- cache offset {p}: cpre/cache all: {q(d.off_pre.dropna())}  min {d.off_pre.min():.3f} max {d.off_pre.max():.3f}; cpost/cache: {q(d.off_post.dropna())} min {d.off_post.min():.3f} max {d.off_post.max():.3f}")
    g = d.groupby(["pre_gid", "post_gid"]).agg(n=("syn_id", "size"), off_pre=("off_pre", "median"), off_pre_min=("off_pre", "min"), off_pre_max=("off_pre", "max"),
                                                gmax=("gmax_NMDA_nS", "median"), cpre=("cpre", "median"))
    print("   pairs with |cpre/cache-1|>0.02:", int(((g.off_pre - 1).abs() > .02).sum()), "of", len(g))
    print(g.loc[(g.off_pre - 1).abs().sort_values(ascending=False).index[:8]].round(4).to_string())
    print("   corr(off_pre, dist):", round(d[["off_pre", "dist_um"]].corr().iloc[0, 1], 3), " corr(off_pre, cpre):", round(d[["off_pre", "cpre"]].corr().iloc[0, 1], 3),
          " off_pre by loc:", d.groupby("loc").off_pre.median().round(4).to_dict())
