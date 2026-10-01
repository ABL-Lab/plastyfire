"""Per-protocol GPU cost of the v4 fit inputs: records, synapses, trace length (effcai header only, no data load)."""
import glob, os, zipfile, csv, sys
import numpy as np
V2 = "/lustre09/project/6070394/dhuruva/plastyfire/glusynapse_v2"; ROOT = os.path.dirname(V2); X = V2 + "/extracted"
L5B = ROOT + "/basis_results_edges_sabrina_n120_delta"; L23B = ROOT + "/basis_results_edges_ebner_l23l5_delta_rs"
pairs5 = set(open(V2 + "/subset24_pairs.txt").read().strip().split(","))
g = list(csv.DictReader(open(ROOT + "/ebner/pair_geometry_L23PC_L5TTPC.csv")))
pairs23 = {r["pregid"] + "-" + r["postgid"] for r in g if r["all_protocols"].strip().lower() in ("true", "1")}
def shape(f, key="effcai"):
    with zipfile.ZipFile(f) as z, z.open(key + ".npy") as fh:
        v = np.lib.format.read_magic(fh); return np.lib.format._read_array_header(fh, v)[0]
rows = []
for d, P, B in [(x, pairs5, L5B) for x in ("ebner_delta-prefire-vca", "markram_delta-prefire-vca", "sj03_delta-prefire-vca",
                "sj03r50_delta-prefire-vca", "sj07_delta-prefire-vca")] + [("ebner_l23l5_delta-prefire-vseg-rs", pairs23, L23B)]:
    agg = {}
    for f in sorted(glob.glob(f"{X}/{d}/*.npz")):
        pair, proto = os.path.basename(f)[:-4].split("__")
        if pair not in P or not os.path.isfile(f"{B}/basis_{pair.split('-')[0]}_{pair.split('-')[1]}.csv"):
            continue
        n, T = shape(f)
        a = agg.setdefault(proto, [0, 0, 0, 10**9, 0]); a[0] += 1; a[1] += n; a[2] += n * T; a[3] = min(a[3], T); a[4] = max(a[4], T)
    for p, a in sorted(agg.items()):
        rows.append([d, p] + a)
w = csv.writer(sys.stdout); w.writerow(["dir", "proto", "n_rec", "n_syn", "syn_steps", "T_min", "T_max"]); w.writerows(rows)
