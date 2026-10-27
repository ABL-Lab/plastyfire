"""Synapse location summary per pair of an index csv (path distance from soma, section type),
used to split L2/3 -> L5 pairs into proximal/distal as in Letzkus 2006.
With --epsp also the unitary somatic EPSP of each pair (prefire setup: GluSynapse, stochastic release, 30 spikes
at 1 Hz): mean single-sweep peak (epsp_mV), peak and 10-90% rise time of the averaged waveform. Splits:
Letzkus proximal/distal by rise time 2.7 ms (primary; path_mean is the cross-check), Sjostrom & Hausser 2006
distal = path_mean > 200 um and EPSP < 1 mV (distal LTD only below ~1 mV). all_protocols flags pairs that have
every protocol of the index (posts without a valid stimulus lack some).
usage: python ebner/pair_geometry.py <index csv> <out csv> [--epsp <sims dir> <workers>]"""
import sys
import numpy as np
import pandas as pd
from bluepysnap import Circuit
from plastyfire.pipette import afferent_table, SONATA_APICAL, SONATA_BASAL

CIRCUIT = "/project/rrg-emuller/dhuruva/plastyfire/data/dhuruva_delta_circuit_config.json"
NODE_POP, EDGE_POP = "S1nonbarrel_neurons", "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"

idx = pd.read_csv(sys.argv[1])
c = Circuit(CIRCUIT)
protos = set(idx["protocol_id"])
has_all = idx.groupby(["pregid", "postgid"])["protocol_id"].agg(lambda p: set(p) == protos)
rows = []
for (pre, post), _ in idx.groupby(["pregid", "postgid"]):
    df = afferent_table(c, NODE_POP, EDGE_POP, int(post))
    d = df.loc[df["@source_node"] == pre]
    ap = d.loc[d["afferent_section_type"] == SONATA_APICAL, "path_dist"]
    rows.append({"pregid": pre, "postgid": post, "all_protocols": bool(has_all[(pre, post)]), "n_syn": len(d),
                 "n_basal": int((d["afferent_section_type"] == SONATA_BASAL).sum()), "n_apical": len(ap),
                 "path_mean": d["path_dist"].mean(), "path_max": d["path_dist"].max(),
                 "n_apical_gt200": int((ap > 200).sum()), "n_apical_gt300": int((ap > 300).sum()), "n_apical_gt450": int((ap > 450).sum()), "n_apical_gt500": int((ap > 500).sum()),
                 "n_apical_gt600": int((ap > 600).sum()),
                 "post_apical_max": df.loc[df["afferent_section_type"] == SONATA_APICAL, "path_dist"].max()})
out = pd.DataFrame(rows)
if "--epsp" in sys.argv:
    import os, json
    from concurrent.futures import ThreadPoolExecutor
    from plastyfire.pipette import measure_epsp
    sims_dir, workers = os.path.abspath(sys.argv[sys.argv.index("--epsp") + 1]), int(sys.argv[sys.argv.index("--epsp") + 2])
    st = [1000. * (k + 1) for k in range(30)]

    def epsp(row):
        wd = os.path.join(sims_dir, "%i-%i" % (row.pregid, row.postgid))
        wd = os.path.join(wd, sorted(os.listdir(wd))[0])
        with open(os.path.join(wd, "prefire_simulation_config.json")) as f:
            cfg = json.load(f)
        cfg["node_sets_file"] = os.path.join(wd, "node_sets.json")
        cfg["run"]["tstop"] = st[-1] + 200.
        cfg["inputs"] = {}
        cfg["output"] = {"output_dir": os.environ.get("SLURM_TMPDIR", "/tmp")}  # workdir out/ may not exist yet
        tmp = os.path.join(os.environ.get("SLURM_TMPDIR", "/tmp"), "epsp_%i-%i.json" % (row.pregid, row.postgid))
        with open(tmp, "w") as f:
            json.dump(cfg, f)
        return measure_epsp(tmp, row.postgid, [row.pregid], NODE_POP, st,
                            "/project/rrg-emuller/dhuruva/DEES_cell_packages/")

    with ThreadPoolExecutor(workers) as pool:
        res = list(pool.map(epsp, out.itertuples()))
    out["epsp_mV"] = [r["epsp"] for r in res]
    out["epsp_avg_peak_mV"] = [r["peak_avg"] for r in res]
    out["rise_10_90_ms"] = [r["rise_10_90"] for r in res]
    out["letzkus_distal"] = out["rise_10_90_ms"] > 2.7
    out["letzkus_distal_median"] = out["rise_10_90_ms"] > out["rise_10_90_ms"].median()
    out["sh_distal"] = (out["path_mean"] > 200.) & (out["epsp_mV"] < 1.)
out.to_csv(sys.argv[2], index=False)
print(out.describe().T.to_string())
if "epsp_mV" in out:
    print("unitary EPSP: mean %.2f median %.2f max %.2f mV; pairs > 1 mV: %d / %d"
          % (out["epsp_mV"].mean(), out["epsp_mV"].median(), out["epsp_mV"].max(), (out["epsp_mV"] > 1.).sum(), len(out)))
print("all_protocols: %d / %d" % (out["all_protocols"].sum(), len(out)))
if "rise_10_90_ms" in out:
    ok = out["all_protocols"]
    print("rise 10-90: median %.2f ms; distal (>2.7 ms): %d / %d (all_protocols: %d)"
          % (out["rise_10_90_ms"].median(), out["letzkus_distal"].sum(), len(out), (out["letzkus_distal"] & ok).sum()))
    print("rise vs path_mean: spearman r = %.2f" % out["rise_10_90_ms"].corr(out["path_mean"], method="spearman"))
    print("S&H distal & < 1 mV: %d (all_protocols: %d)" % (out["sh_distal"].sum(), (out["sh_distal"] & ok).sum()))
for thr in [200, 300, 450, 500, 600]:
    print("pairs with >=1 apical syn beyond %d um: %d / %d" % (thr, (out["n_apical_gt%d" % thr] > 0).sum(), len(out)))
