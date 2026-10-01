"""Which term of the offline pre path makes the live (continuous GluSynapseV4, bin_GB 0) dpre differ from offline?
(equiv 22131696: dpre live > offline at L5 sj50 -10 and sj07 pair, both prefire / cooker, identical event counts.)

Per record, offline dpre_final (control, mglu_block = NO only, no_block = eCB only) with
  std   the fit's features (N driven by pos(bin mean of -ica_VDCC - theta_NOi), extract_v2.step_mean, 0.25 ms grid)
  rect  N driven by the bin mean of pos(-ica_VDCC - theta_NOi) from the full-resolution trace (what a continuous
        mechanism integrates), same grid and everything else
against the live prefire / cooker dpre of the 22131696 jsonl. rect ~ live -> the difference is the rectification of
the bin mean (NO path), not the eCB path.

    python glusynapse_v2/bcl_validation/diag_dpre_v4.py --fit <json> --live <jsonl>
"""
import argparse, json, os, pickle, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); V2 = os.path.dirname(HERE); ROOT = os.path.dirname(V2)
sys.path.insert(0, V2); sys.path.insert(0, os.path.join(V2, "rho_redesign")); sys.path.insert(0, HERE)
import batch_v2                       # noqa: E402
from batch_v2 import BatchV2, _sig    # noqa: E402
import model_v2 as MV                 # noqa: E402
from live_v4 import PATHS, v4_globals  # noqa: E402
from compare_v4 import BASIS          # noqa: E402

CASES = [("L5", "180351-198084", "sjostrom_50hz_dt-10ms", "delta-prefire-vseg"),
         ("L5", "180351-198084", "sjostrom07_step200ms_pair", "delta-prefire-vseg"),
         ("L23", "7471-187420", "sjostrom_50hz_dt+10ms", "delta-prefire-vseg-rs")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True)
    ap.add_argument("--live", required=True)
    a = ap.parse_args()
    fit = json.load(open(a.fit))
    _, P = v4_globals(fit, "control", "cooker")
    live = {}
    for line in open(a.live):
        d = json.loads(line)
        if d.get("ok") and d["phase"] == "prefire" and d["express"] == "cooker":
            live[(d["pair"], d["proto"], d["cond"])] = dict(zip(d["syn"], d["end"]["dpre_GB"]))
    th = P["theta_NOi"]
    for path, pair, proto, hsh in CASES:
        batch_v2.BASIS_DIR = BASIS[path]
        B = BatchV2(PATHS[path]["dirs"], protocols=[proto], pairs={pair}, fast=False, signals=("vdcc",), verbose=False)
        r = B.recs[0]
        tr = pickle.load(open(os.path.join(PATHS[path]["sims"], pair, proto, f"bluecellulab_results_{hsh}",
                                           "simulation_traces.pkl"), "rb"))
        t_raw0 = float(np.asarray(tr["t"])[0]); n_raw = len(next(iter(tr["ica_VDCC"].values())))
        raw_dt = float(r.get("raw_dt_ms", 0.025))
        keep = np.rint((np.asarray(r["t"]) - t_raw0) / raw_dt).astype(np.int64)
        ends = np.append(keep[1:], n_raw)
        rect = np.zeros_like(np.asarray(r["vdcc"], np.float64)); std_chk = np.zeros_like(rect)
        for i, s in enumerate(r["syn"]):
            q = -np.asarray(tr["ica_VDCC"][int(s)], np.float64)
            c = np.concatenate([[0.0], np.cumsum(np.maximum(q - th, 0.0))])
            rect[i] = (c[ends] - c[keep]) / (ends - keep) + th
            c2 = np.concatenate([[0.0], np.cumsum(q)])
            std_chk[i] = (c2[ends] - c2[keep]) / (ends - keep)
        del tr
        tT = MV.feature_T(_sig(r, P), r["t"], r["arr"], P)
        K_std = MV.feature_K(np.asarray(r["vdcc"], np.float64), r["t"], r["arr"], P)
        K_rect = MV.feature_K(rect, r["t"], r["arr"], P)
        print(f"\n{path} {pair} {proto}: n_syn {len(r['syn'])}, grid check max|vdcc - step_mean| "
              f"{np.max(np.abs(std_chk - r['vdcc'])):.3g} nA, sum K std {K_std.sum():.4g} rect {K_rect.sum():.4g}, "
              f"tT>0 at {int((tT > 0).sum())} arrivals, mean tT {tT[tT > 0].mean() if (tT > 0).any() else 0:.4f}")
        for cond, am, an in (("control", P["A_mglu"], P["A_NO"]), ("mglu_block", 0.0, P["A_NO"]),
                             ("no_block", P["A_mglu"], 0.0)):
            ds = MV.dpre_final(tT, K_std, am, an, P["dpre_min"], P["dpre_max"], P["dpre0"])
            dr = MV.dpre_final(tT, K_rect, am, an, P["dpre_min"], P["dpre_max"], P["dpre0"])
            lv = live.get((pair, proto, cond))
            msg = f"  {cond:10s} offline std {ds.mean():.4f}  rect {dr.mean():.4f}"
            if lv:
                dl = np.array([lv.get(int(s), np.nan) for s in r["syn"]])
                msg += (f"  live {np.nanmean(dl):.4f}  max|live-std| {np.nanmax(np.abs(dl - ds)):.4f}"
                        f"  max|live-rect| {np.nanmax(np.abs(dl - dr)):.4f}")
            print(msg, flush=True)


if __name__ == "__main__":
    main()
