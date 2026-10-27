"""Unitary somatic EPSP per pair on a chosen circuit, same readout as ebner/pair_geometry.py --epsp (prefire setup,
GluSynapse, stochastic release, 30 spikes at 1 Hz; mean single-sweep peak, averaged-waveform peak, 10-90% rise).
Only the "network" of each pair's prefire config is swapped, so base SSCx (original emodels_hoc + original edges)
and delta (og-delta emodels + dhuruva_modified_edges) are compared on identical pairs and stimulation.
usage: python ebner/epsp_base_sscx.py <index csv> <circuit config> <out csv> <workers>"""
import os, sys, json
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
from plastyfire.pipette import measure_epsp

NODE_POP = "S1nonbarrel_neurons"
idx_csv, circuit, out_csv, workers = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
pairs = pd.read_csv(idx_csv).groupby(["pregid", "postgid"]).path.first().reset_index()
st = [1000. * (k + 1) for k in range(30)]
tmpdir = os.environ.get("SLURM_TMPDIR", "/tmp")


def epsp(row):
    cfg_path = os.path.join(os.path.dirname(row.path), "prefire_simulation_config.json")
    with open(cfg_path) as f:
        cfg = json.load(f)
    cfg["network"] = circuit
    cfg["node_sets_file"] = os.path.join(os.path.dirname(cfg_path), "node_sets.json")
    cfg["run"]["tstop"] = st[-1] + 200.
    cfg["inputs"] = {}
    cfg["output"] = {"output_dir": tmpdir}
    tmp = os.path.join(tmpdir, "epsp_%i-%i.json" % (row.pregid, row.postgid))
    with open(tmp, "w") as f:
        json.dump(cfg, f)
    try:
        r = measure_epsp(tmp, row.postgid, [row.pregid], NODE_POP, st, "/project/rrg-emuller/dhuruva/DEES_cell_packages/")
    except RuntimeError as e:
        print("fail %i-%i: %s" % (row.pregid, row.postgid, str(e).splitlines()[0]), flush=True)
        return {"pregid": row.pregid, "postgid": row.postgid}
    return {"pregid": row.pregid, "postgid": row.postgid, "n_syn": r["n_syn"], "epsp_mV": r["epsp"],
            "epsp_avg_peak_mV": r["peak_avg"], "rise_10_90_ms": r["rise_10_90"], "spiked": r["spiked"]}


with ThreadPoolExecutor(workers) as pool:
    out = pd.DataFrame(list(pool.map(epsp, pairs.itertuples())))
out.to_csv(out_csv, index=False)
print("%s on %s: n=%d, EPSP mean %.3f median %.3f mV" % (os.path.basename(idx_csv), os.path.basename(circuit),
      out.epsp_mV.notna().sum(), out.epsp_mV.mean(), out.epsp_mV.median()))
