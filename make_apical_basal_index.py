"""
Split connection pairs by synapse location (apical / basal / mixed).

Writes a mapper CSV describing, per (pre_gid, post_gid), how many of that
connection's synapses land on apical vs basal dendrite, plus a `category`:

    apical  - every synapse on apical dendrite (afferent_section_type == 3)
    basal   - every synapse on basal
    mixed   - both present (the common case; captured with both counts)

Section type comes from edges.h5 afferent_section_type, the same source
simulator_edges._apply_theta_from_a_params uses to pick the a-param row, so the
split here matches what the plasticity model actually does.

Synapse membership per pair is read from the c_pre/c_post cache, whose keys are
(pre_gid, post_gid) and whose inner dicts are keyed by synapse global id -- the
same global ids that index afferent_section_type.

Also emits index_<label>_apical.csv / _basal.csv alongside the source simulation
index, filtered to pairs whose synapses are purely apical / purely basal, so the
evaluator's index-driven pair selection can be pointed at one population.

Usage:
    python make_apical_basal_index.py \
        --cache cpre_cpost_cache/ion_channels_tau278.pkl \
        --out apical_basal_map.csv
"""

import argparse
import os
import pickle

import h5py
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
EDGES_H5 = os.path.join(HERE, "data/dhuruva_modified_edges.h5")
EDGE_POP = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"
APICAL_SECTION_TYPE = 3          # matches _apply_theta_from_a_params
INDEX_DIR = os.path.join(HERE, "refitting_results/fitting/n100/seed19091997")
SRC_INDEX = os.path.join(INDEX_DIR, "index_L5TTPC_L5TTPC_STDP.csv")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache", default="cpre_cpost_cache/ion_channels_tau278.pkl")
    ap.add_argument("--out", default="apical_basal_map.csv")
    ap.add_argument("--write-index", action="store_true",
                    help="also write index_L5TTPC_L5TTPC_{apical,basal}.csv")
    args = ap.parse_args()

    cache = pickle.load(open(os.path.join(HERE, args.cache), "rb"))
    with h5py.File(EDGES_H5, "r") as f:
        sec = f["edges/%s/0/afferent_section_type" % EDGE_POP][:]

    rows = []
    for (pre, post), v in sorted(cache.items()):
        gids = sorted(v["c_pre"].keys())
        types = np.array([sec[g] for g in gids])
        n_ap = int((types == APICAL_SECTION_TYPE).sum())
        n_ba = int(len(types) - n_ap)
        cat = "apical" if n_ba == 0 else ("basal" if n_ap == 0 else "mixed")
        # mean c_post per location: the quantity that drives theta, useful when
        # asking whether one population can potentiate on its own.
        cq = v["c_post"]
        cp = v["c_pre"]
        ap_post = [cq[g] for g, t in zip(gids, types) if t == APICAL_SECTION_TYPE]
        ba_post = [cq[g] for g, t in zip(gids, types) if t != APICAL_SECTION_TYPE]
        ap_pre = [cp[g] for g, t in zip(gids, types) if t == APICAL_SECTION_TYPE]
        ba_pre = [cp[g] for g, t in zip(gids, types) if t != APICAL_SECTION_TYPE]
        rows.append(dict(
            pregid=pre, postgid=post, n_syn=len(gids),
            n_apical=n_ap, n_basal=n_ba,
            frac_apical=n_ap / len(gids) if gids else np.nan,
            category=cat,
            mean_cpre_apical=np.mean(ap_pre) if ap_pre else np.nan,
            mean_cpost_apical=np.mean(ap_post) if ap_post else np.nan,
            mean_cpre_basal=np.mean(ba_pre) if ba_pre else np.nan,
            mean_cpost_basal=np.mean(ba_post) if ba_post else np.nan,
            apical_syn_ids=";".join(str(g) for g, t in zip(gids, types)
                                    if t == APICAL_SECTION_TYPE),
            basal_syn_ids=";".join(str(g) for g, t in zip(gids, types)
                                   if t != APICAL_SECTION_TYPE),
        ))

    df = pd.DataFrame(rows)
    out = os.path.join(HERE, args.out)
    df.to_csv(out, index=False)

    print("cache      : %s  (%d pairs)" % (args.cache, len(df)))
    print("wrote      : %s\n" % out)
    print("Pairs by category:")
    print(df["category"].value_counts().to_string())
    print("\nSynapse totals: apical=%d  basal=%d  (%.1f%% apical)" % (
        df.n_apical.sum(), df.n_basal.sum(),
        100.0 * df.n_apical.sum() / (df.n_apical.sum() + df.n_basal.sum())))
    print("\nPer-pair apical fraction: median=%.2f  min=%.2f  max=%.2f" % (
        df.frac_apical.median(), df.frac_apical.min(), df.frac_apical.max()))
    print("\nMean c_pre / c_post by location (over pairs that have them):")
    print("  apical : c_pre=%.4e  c_post=%.4e" % (
        df.mean_cpre_apical.mean(), df.mean_cpost_apical.mean()))
    print("  basal  : c_pre=%.4e  c_post=%.4e" % (
        df.mean_cpre_basal.mean(), df.mean_cpost_basal.mean()))

    if args.write_index:
        src = pd.read_csv(SRC_INDEX)
        for cat in ("apical", "basal"):
            keep = {(int(r.pregid), int(r.postgid))
                    for r in df.itertuples() if r.category == cat}
            sub = src[[(int(a), int(b)) in keep
                       for a, b in zip(src.pregid, src.postgid)]]
            dst = os.path.join(INDEX_DIR, "index_L5TTPC_L5TTPC_%s.csv" % cat)
            sub.to_csv(dst, index=False)
            print("\nwrote %s  (%d rows, %d pairs)" % (dst, len(sub), len(keep)))


if __name__ == "__main__":
    main()
