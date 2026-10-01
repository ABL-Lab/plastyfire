"""Letzkus-calibrated geometric distal/proximal split of the L2/3 -> L5 pairs (LETZKUS_LOCATION.md).
distal = >= half of the pair's synapses on apical dendrites (oblique/trunk/tuft) at path distance > 430 um
(Letzkus 2006 Fig 1F: rise 2.7 ms ~ 400-450 um). Proximal = complement. Writes a copy of the pair-geometry csv with
letzkus_distal replaced (original column kept as letzkus_distal_rise_ms_split; not a *_distal name, so it is not read as a split).
Pairs with >= half basal synapses beyond 200 um are flagged (basal_far_flag) but NOT removed: dropping rows would
also drop them from the non-Letzkus L2/3 targets, which read the same csv.
usage: python letzkus_geo.py <ours_synapses.csv> <old geometry csv> <out csv>"""
import sys
import pandas as pd
syn = pd.read_csv(sys.argv[1]); old = pd.read_csv(sys.argv[2])
syn["ap_far"] = syn.cls.isin(["oblique", "trunk", "tuft"]) & (syn.path_dist > 430.)
syn["bas_far"] = (syn.cls == "basal") & (syn.path_dist > 200.)
g = syn.groupby(["pregid", "postgid"]).agg(frac_ap_far=("ap_far", "mean"), frac_bas_far=("bas_far", "mean")).reset_index()
out = old.merge(g, on=["pregid", "postgid"], how="left")
assert out.frac_ap_far.notna().all(), "pairs missing in ours_synapses.csv"
out["letzkus_distal_rise_ms_split"] = out["letzkus_distal"]
out["letzkus_distal"] = out.frac_ap_far >= 0.5
out["basal_far_flag"] = out.frac_bas_far >= 0.5
out.to_csv(sys.argv[3], index=False)
for name, m in [("all 120", out.all_protocols | ~out.all_protocols), ("all_protocols", out.all_protocols)]:
    o = out[m]
    print(f"{name}: n {len(o)}  NEW distal {o.letzkus_distal.sum()} proximal {(~o.letzkus_distal).sum()}  "
          f"(OLD rise>2.7 ms distal {o.letzkus_distal_rise_ms_split.sum()} proximal {(~o.letzkus_distal_rise_ms_split).sum()})  "
          f"basal_far_flag among proximal {(o.basal_far_flag & ~o.letzkus_distal).sum()}")
o = out[out.all_protocols]
print("crosstab all_protocols (rows OLD, cols NEW):"); print(pd.crosstab(o.letzkus_distal_rise_ms_split, o.letzkus_distal).to_string())
print("NEW distal path_mean median %.0f, n_basal mean %.2f; NEW proximal path_mean median %.0f"
      % (o[o.letzkus_distal].path_mean.median(), o[o.letzkus_distal].n_basal.mean(), o[~o.letzkus_distal].path_mean.median()))
