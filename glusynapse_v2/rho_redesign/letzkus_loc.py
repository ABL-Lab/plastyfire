"""Where do our L2/3 -> L5 TTPC synapses sit, compared with Letzkus 2006 (J Neurosci 26:10420)?
Per synapse of the 120 Ebner/Letzkus pairs: dendrite class (basal / oblique / trunk / tuft) and path distance.
Trunk = the apical path from its root to the main bifurcation (deepest common ancestor of the apical leaves in the
top 20% of apical height, same idea as morph_tool.apical_point); tuft = sections beyond it; oblique = other apical.
Also: our rise-time -> distance mapping (Letzkus Fig 1F: 2.7 ms ~ 400-450 um), how many of our pairs would pass
an apical-distance criterion, and all L2/3 PC contacts in the circuit onto the same posts (is the > 400 um pool there?).
Read-only on the circuit. usage: python letzkus_loc.py <out dir>"""
import os
import sys
import numpy as np
import pandas as pd
from bluepysnap import Circuit
from morphio import SectionType, IterType
from plastyfire.pipette import syn_path_distances

ROOT = "/project/rrg-emuller/dhuruva/plastyfire"
CIRCUIT = os.path.join(ROOT, "data", "dhuruva_delta_l23l5_circuit_config.json")   # edges: dhuruva_modified_edges_l23l5.h5
GEOM = os.path.join(ROOT, "ebner", "pair_geometry_L23PC_L5TTPC.csv")
NODE_POP, EDGE_POP = "S1nonbarrel_neurons", "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"
PROPS = ["@source_node", "afferent_section_id", "afferent_section_pos"]
CLASSES = ["soma", "basal", "oblique", "trunk", "tuft"]
BINS = [0, 100, 200, 300, 400, 500, 600, 800, 1000, 2000]
out_dir = sys.argv[1]
os.makedirs(out_dir, exist_ok=True)


def classify(morph):
    """morphio section id -> class; also the apical-point path distance (end of the trunk)"""
    cls = {s.id: ("basal" if s.type == SectionType.basal_dendrite else "axon") for s in morph.iter()}
    apical = [r for r in morph.root_sections if r.type == SectionType.apical_dendrite]
    if not apical:
        return cls, np.nan
    ap = apical[0]
    pts = np.vstack([s.points for s in ap.iter()])
    y0 = morph.soma.center[1]
    up = ap.points[:, 1].mean() > y0
    h = (pts[:, 1].max() - y0) if up else (y0 - pts[:, 1].min())
    leaves = [s for s in ap.iter() if not s.children and ((s.points[-1, 1] - y0) if up else (y0 - s.points[-1, 1])) >= 0.8 * h]
    common = None
    for lf in leaves:
        anc = [s.id for s in lf.iter(IterType.upstream)][::-1]          # root ... leaf
        common = anc if common is None else [a for a, b in zip(common, anc) if a == b]
    if len(leaves) == 1:
        common = common[:-1] or common
    trunk = set(common)
    end = morph.section(common[-1])
    for s in ap.iter():
        cls[s.id] = "trunk" if s.id in trunk else "oblique"
    for s in end.iter():
        if s.id != end.id:
            cls[s.id] = "tuft"
    d_end = syn_path_distances(morph, np.array([end.id + 1]), np.array([1.0]))[0]
    return cls, d_end


c = Circuit(CIRCUIT)
nodes, edges = c.nodes[NODE_POP], c.edges[EDGE_POP]
g = pd.read_csv(GEOM)
posts = sorted(g.postgid.unique())
mt_post = nodes.get(posts, ["mtype"])["mtype"]
rows, cells = [], []
for post in posts:
    morph = nodes.morph.get(int(post), transform=False, extension="asc")
    cls, d_end = classify(morph)
    df = edges.afferent_edges(int(post), PROPS)
    mt = nodes.get(df["@source_node"].unique(), ["mtype"])["mtype"]
    df["pre_mtype"] = df["@source_node"].map(mt).astype(str)
    df = df[df.pre_mtype.str.match(r"^L[23]_.*PC")].copy()            # every L2/3 PC afferent of this post
    sid = df.afferent_section_id.to_numpy().astype(int)
    df["path_dist"] = syn_path_distances(morph, sid, df.afferent_section_pos.to_numpy())
    df["cls"] = ["soma" if s == 0 else cls[s - 1] for s in sid]
    df["postgid"] = post
    rows.append(df[["@source_node", "postgid", "pre_mtype", "cls", "path_dist"]].rename(columns={"@source_node": "pregid"}))
    apd = [syn_path_distances(morph, np.array([s.id + 1]), np.array([1.0]))[0] for s in morph.iter()
           if s.type == SectionType.apical_dendrite and not s.children]
    cells.append(dict(postgid=post, mtype=mt_post[post], apical_point_um=d_end, apical_max_um=max(apd) if apd else np.nan))
allsyn = pd.concat(rows, ignore_index=True)
cells = pd.DataFrame(cells)
ours = allsyn.merge(g[["pregid", "postgid", "all_protocols", "rise_10_90_ms", "letzkus_distal", "epsp_mV"]],
                    on=["pregid", "postgid"])
ours.to_csv(os.path.join(out_dir, "ours_synapses.csv"), index=False)
cells.to_csv(os.path.join(out_dir, "post_cells.csv"), index=False)

L = []
p = L.append
p("posts: %d (%s); apical point (trunk end) path um: median %.0f [%.0f-%.0f]; apical max path um: median %.0f"
  % (len(cells), dict(cells.mtype.value_counts()), cells.apical_point_um.median(), cells.apical_point_um.min(),
     cells.apical_point_um.max(), cells.apical_max_um.median()))
p("our pairs: %d, synapses: %d (index pairs %d)" % (ours.groupby(["pregid", "postgid"]).ngroups, len(ours), len(g)))


def cls_table(d, by):
    t = pd.crosstab(d[by], d.cls).reindex(columns=CLASSES, fill_value=0)
    f = t.div(t.sum(1), axis=0).round(3)
    q = d.groupby(by).path_dist.describe(percentiles=[.25, .5, .75, .9])[["count", "25%", "50%", "75%", "90%", "max"]].round(0)
    return "counts\n%s\nfractions\n%s\npath_dist um\n%s" % (t.to_string(), f.to_string(), q.to_string())


p("\n== A. our synapses by Letzkus split (letzkus_distal = rise > 2.7 ms) ==")
p(cls_table(ours, "letzkus_distal"))
p("\n-- all_protocols pairs only --")
p(cls_table(ours[ours.all_protocols], "letzkus_distal"))
p("\n== A2. our synapses by distance bin x class ==")
ours["bin"] = pd.cut(ours.path_dist, BINS)
p(pd.crosstab(ours["bin"], ours.cls).reindex(columns=CLASSES, fill_value=0).to_string())

pair = ours.groupby(["pregid", "postgid"]).agg(
    n=("cls", "size"), path_mean=("path_dist", "mean"), path_max=("path_dist", "max"),
    n_basal=("cls", lambda s: (s == "basal").sum()), n_obl=("cls", lambda s: (s == "oblique").sum()),
    n_trunk=("cls", lambda s: (s == "trunk").sum()), n_tuft=("cls", lambda s: (s == "tuft").sum()),
    rise=("rise_10_90_ms", "first"), distal=("letzkus_distal", "first"), allp=("all_protocols", "first")).reset_index()
pair["n_tt"] = pair.n_trunk + pair.n_tuft
ap = ours[ours.cls.isin(["oblique", "trunk", "tuft"])]
for X in [300, 400, 450, 600]:
    s = ours.assign(f=(ours.cls.isin(["oblique", "trunk", "tuft"])) & (ours.path_dist > X)).groupby(["pregid", "postgid"]).f.mean()
    pair["frac_ap_gt%d" % X] = pair.set_index(["pregid", "postgid"]).index.map(s).to_numpy()
pair.to_csv(os.path.join(out_dir, "ours_pairs.csv"), index=False)

p("\n== B. rise time vs pair mean path distance (Letzkus Fig 1F: distance = a*rise + b, r 0.98, 2.7 ms ~ 400-450 um) ==")
a, b = np.polyfit(pair.rise, pair.path_mean, 1)
p("ours: distance = %.0f*rise %+.0f um; pearson r %.2f; distance at 2.7 ms = %.0f um; rise range %.2f-%.2f ms, median %.2f"
  % (a, b, np.corrcoef(pair.rise, pair.path_mean)[0, 1], a * 2.7 + b, pair.rise.min(), pair.rise.max(), pair.rise.median()))
for lo, hi in [(0, 2), (2, 2.7), (2.7, 3.5), (3.5, 5), (5, 99)]:
    q = pair[(pair.rise >= lo) & (pair.rise < hi)]
    p("  rise %.1f-%.1f ms: %3d pairs, path_mean median %4.0f um, frac syn trunk+tuft %.2f, oblique %.2f, basal %.2f"
      % (lo, hi, len(q), q.path_mean.median() if len(q) else np.nan, q.n_tt.sum() / max(q.n.sum(), 1),
         q.n_obl.sum() / max(q.n.sum(), 1), q.n_basal.sum() / max(q.n.sum(), 1)))

p("\n== C. our pairs that would pass an apical-distance criterion (all 120 / all_protocols) ==")
for X in [300, 400, 450, 600]:
    f = pair["frac_ap_gt%d" % X]
    for lab, m in [("all syn apical > X", f == 1), (">= half syn apical > X", f >= 0.5), ("any syn apical > X", f > 0)]:
        p("  X=%d %-24s %3d / %d   (all_protocols %3d; of them letzkus_distal %d)"
          % (X, lab, m.sum(), len(pair), (m & pair.allp).sum(), (m & pair.allp & pair.distal).sum()))
p("  pairs with >=1 trunk/tuft synapse: %d; >= half trunk/tuft: %d; any tuft: %d"
  % ((pair.n_tt > 0).sum(), (pair.n_tt >= pair.n / 2).sum(), (pair.n_tuft > 0).sum()))

p("\n== D. circuit: all L2/3 PC -> these %d posts (pool for re-selection) ==" % len(posts))
allsyn["bin"] = pd.cut(allsyn.path_dist, BINS)
p(pd.crosstab(allsyn["bin"], allsyn.cls).reindex(columns=CLASSES, fill_value=0).to_string())
p("synapses: %d; apical (obl/trunk/tuft) > 400 um: %d (%.1f%%); trunk+tuft > 400 um: %d; tuft: %d"
  % (len(allsyn), ((allsyn.cls != "basal") & (allsyn.path_dist > 400)).sum(),
     100 * ((allsyn.cls != "basal") & (allsyn.path_dist > 400)).mean(),
     (allsyn.cls.isin(["trunk", "tuft"]) & (allsyn.path_dist > 400)).sum(), (allsyn.cls == "tuft").sum()))
con = allsyn.groupby(["pregid", "postgid"]).agg(n=("cls", "size"), path_mean=("path_dist", "mean"))
fa = allsyn.assign(f=(allsyn.cls != "basal") & (allsyn.path_dist > 400)).groupby(["pregid", "postgid"]).f.mean()
con["f400"] = fa
p("connections: %d (%.1f syn/conn); all syn apical > 400 um: %d; >= half: %d; posts with >=1 such connection: %d"
  % (len(con), con.n.mean(), (con.f400 == 1).sum(), (con.f400 >= .5).sum(),
     con[con.f400 >= .5].reset_index().postgid.nunique()))
p("pre mtypes: %s" % dict(allsyn.drop_duplicates(["pregid", "postgid"]).pre_mtype.value_counts()))
pq = allsyn.groupby(["pregid", "postgid"]).pre_mtype.first()
for X in [400, 600]:
    m = allsyn.assign(f=(allsyn.cls != "basal") & (allsyn.path_dist > X)).groupby(["pregid", "postgid"]).f.mean() >= .5
    p("  >= half syn apical > %d um, by pre mtype: %s" % (X, dict(pq[m].value_counts())))
txt = "\n".join(L)
print(txt)
with open(os.path.join(out_dir, "summary.txt"), "w") as f:
    f.write(txt + "\n")
