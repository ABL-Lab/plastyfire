"""Old (cexp) vs new (cexp_ljp25) per pathway and loc: median of NMDA share, APV/pre, Mg0/pre, post/pre (free Ca above rest 0.07 uM, and effcai), vdcc_q_post."""
import pandas as pd, numpy as np
S = "/scratch/dhuruva/split1"; REST = 0.07  # uM = min_ca_CR
def stats(d):
    r = {}
    r["nmda_share"] = (1 - d.cpre_apv / d.cpre).median()
    r["apv/pre"] = (d.cpre_apv / d.cpre).median(); r["mg0/pre"] = (d.cpre_mg0 / d.cpre).median()
    r["post/pre_free"] = ((d.cpost_cai - REST) / (d.cpre_cai - REST)).median(); r["post/pre_eff"] = (d.cpost / d.cpre).median()
    r["vdcc_q_post"] = d.vdcc_q_post.median(); r["n"] = len(d)
    return pd.Series(r)
rows = []
for p in ("L5L5", "L23L5", "L23L23"):
    o = pd.read_csv(f"{S}/cexp/{p}.csv"); n = pd.read_csv(f"{S}/cexp_ljp25/{p}.csv")
    for loc in ("basal", "apical"):
        a = stats(o[o["loc"] == loc]); b = stats(n[n["loc"] == loc])
        for k in a.index: rows.append((p, loc, k, a[k], b[k]))
df = pd.DataFrame(rows, columns=["path", "loc", "metric", "old", "new"])
df["new/old"] = df.new / df.old
pd.set_option("display.width", 200)
print(df.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
