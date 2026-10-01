"""Screen delta-bd candidates on real post cells (bluecellulab, the circuit's og-delta hoc, SK_E2 off as in the runs).

One cell per pair, instantiated once; each candidate rewrites the basal densities (bd.apply) and runs one timeline:

  blk  t (ms)   stimulus (soma TStim)                             measures
  1ap   300     1 x 2 ms @ A_k                                    bAP dV by distance, Ca_1
  f66   700     5 x 2 ms @ 66 Hz, A_k                             Kampa fig1g: distal Ca(f) / Ca(66 Hz)
  f100 1100     5 x 2 ms @ 100 Hz
  f150 1500     5 x 2 ms @ 150 Hz
  f200 1900     5 x 2 ms @ 200 Hz
  r3   2300     3 x 2 ms @ 200 Hz                                 Kampa fig1e R3 = Ca(3AP 200 Hz) / Ca(1AP)
  n100 2700     3 x 5 ms @ 100 Hz, nevian 100 Hz amp              somatic following (protocol pulses)
  l200 3100     3 x 2 ms @ 200 Hz, letzkus amp
  n50  3500     3 x 5 ms @ 50 Hz, nevian 50 Hz amp
  fi   3900..   3 x 600 ms steps at 0.3/0.5/0.7 x nevian 50 Hz amp   somatic f-I (compared with og-delta)
A_k = the pair's Letzkus 2 ms pulse amp (calibrated to fire). Pharmacology (basal Ka off = 4-AP, basal LVA off = Ni2+)
reruns blocks 1ap..r3 only. Ca = peak cai minus the value just before the block, per basal segment; groups by path
distance: prox < 75 um, mid 75-130, dist >= 130 (Kampa: distal >= 130 um).

    python glusynapse_v2/emodel_bd/screen_bd.py --pairs 186261-208490 --cands cands.json --out out.csv
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import bd

SIMS = ("/project/rrg-emuller/dhuruva/plastyfire/refitting_results/fitting/n120/seed20262009/"
        "Sabrina_L5TTPC_L5TTPC_STDP/simulations")
MECH = "/project/rrg-emuller/dhuruva/DEES_cell_packages/"
POP = "S1nonbarrel_neurons"
BLK = dict(ap1=300.0, f66=700.0, f100=1100.0, f150=1500.0, f200=1900.0, r3=2300.0, n100=2700.0, l200=3100.0,
           n50=3500.0)
FI0, FIW, FIG = 3900.0, 600.0, 300.0
FI_FR = (0.3, 0.5, 0.7)
T_PHARM = 2700.0
T_END = FI0 + 3 * (FIW + FIG)
GROUPS = dict(prox=(0, 75), mid=(75, 130), dist=(130, 1e9))


def pulse_amp(pair, proto):
    f = os.path.join(SIMS, pair, proto, "prefire_simulation_config.json")
    if not os.path.exists(f):
        return None
    pc = json.load(open(f))["inputs"]
    k = next(k for k, v in pc.items() if v.get("input_type") == "current_clamp")
    return pc[k]["amp_start"], pc[k]["width"]


def trains(amps):
    """[(t, width, amp)] for every somatic pulse, and the f-I steps."""
    A_k = amps["letzkus"][0]
    P = [(BLK["ap1"], 2.0, A_k)]
    for b, f, n in (("f66", 66, 5), ("f100", 100, 5), ("f150", 150, 5), ("f200", 200, 5), ("r3", 200, 3)):
        P += [(BLK[b] + i * 1000.0 / f, 2.0, A_k) for i in range(n)]
    for b, key, f in (("n100", "n100", 100), ("l200", "letzkus", 200), ("n50", "n50", 50)):
        a, w = amps[key]
        P += [(BLK[b] + i * 1000.0 / f, w, a) for i in range(3)]
    S = [(FI0 + i * (FIW + FIG), FIW, fr * amps["n50"][0]) for i, fr in enumerate(FI_FR)]
    return P, S


def build(pair):
    import bluecellulab
    from bluecellulab import CircuitSimulation
    h = bluecellulab.neuron.h
    cfgp = os.path.join(SIMS, pair, "nevian_3ap_50hz_dt+10ms", "prefire_simulation_config.json")
    cfg = json.load(open(cfgp))
    for k in ("reports", "inputs"):
        cfg.pop(k, None)
    cfg["run"]["tstop"] = T_END
    tmp = os.path.join(os.environ.get("SLURM_TMPDIR", "/tmp"), f"screen_bd_{pair}_{os.getpid()}.json")
    json.dump(cfg, open(tmp, "w"))
    sim = CircuitSimulation(tmp)
    post = int(pair.split("-")[1])
    sim.instantiate_gids([(POP, post)], add_synapses=False, add_minis=False, add_pulse_stimuli=False)
    cell = sim.cells[(POP, post)]
    amps = dict(letzkus=pulse_amp(pair, "letzkus_3ap_200hz_dt+10ms"), n100=pulse_amp(pair, "nevian_3ap_100hz_dt+10ms"),
                n50=pulse_amp(pair, "nevian_3ap_50hz_dt+10ms"))
    amps = {k: v or amps["n50"] for k, v in amps.items()}
    if amps["letzkus"] is amps["n50"]:                           # no Letzkus config: 2 ms at 2.5x the 5 ms amp
        amps["letzkus"] = (2.5 * amps["n50"][0], 2.0)
    P, S = trains(amps)
    keep = []
    for t, w, a in P:
        s = h.TStim(0.5, sec=cell.soma); s.train(t, w, a, 0.1, w); keep.append(s)
    for t, w, a in S:
        s = h.TStim(0.5, sec=cell.soma); s.train(t, w, a, 0.1, w); keep.append(s)
    cell.persistent.extend(keep)
    segs = bd.basal_segments(cell, h)
    # record every basal segment beyond 20 um (v, cai) and the soma, on a 0.1 ms grid
    rec_t = h.Vector(); rec_t.record(h._ref_t, 0.1)
    vs = h.Vector(); vs.record(cell.soma(0.5)._ref_v, 0.1)
    R = []
    for seg, d in segs:
        if d < 20:
            continue
        v = h.Vector(); v.record(seg._ref_v, 0.1); c = h.Vector(); c.record(seg._ref_cai, 0.1)
        R.append((d, v, c))
    return sim, cell, segs, dict(t=rec_t, vs=vs, R=R, keep=keep), amps


def spikes(t, v, t0, t1):
    m = (t >= t0) & (t < t1)
    tt, vv = t[m], v[m]
    return tt[1:][(vv[1:] > -10) & (vv[:-1] <= -10)]


def measure(rec, tstop, pharm=False):
    t = np.asarray(rec["t"]); vs = np.asarray(rec["vs"])
    d = np.array([x[0] for x in rec["R"]])
    V = np.array([np.asarray(x[1]) for x in rec["R"]]); C = np.array([np.asarray(x[2]) for x in rec["R"]])
    n = min(len(t), V.shape[1]); t, vs, V, C = t[:n], vs[:n], V[:, :n], C[:, :n]
    out = {}
    blocks = ["ap1", "f66", "f100", "f150", "f200", "r3"] + ([] if pharm else ["n100", "l200", "n50"])
    for b in blocks:
        t0 = BLK[b]; i0 = np.searchsorted(t, t0 - 1.0); m = (t >= t0) & (t < t0 + 250.0)
        ca = C[:, m].max(1) - C[:, i0]
        for g, (lo, hi) in GROUPS.items():
            sel = (d >= lo) & (d < hi)
            out[f"ca_{b}_{g}"] = float(np.median(ca[sel])) if sel.any() else np.nan
        out[f"nsp_{b}"] = int(len(spikes(t, vs, t0 - 1, t0 + 60)))
        if b == "ap1":
            dv = V[:, m].max(1) - V[:, i0]
            for g, (lo, hi) in dict(d50_100=(50, 100), d100_150=(100, 150), d200=(200, 1e9)).items():
                sel = (d >= lo) & (d < hi)
                out[f"bap_dv_{g}"] = float(np.median(dv[sel])) if sel.any() else np.nan
    if not pharm:
        for i, fr in enumerate(FI_FR):
            t0 = FI0 + i * (FIW + FIG)
            out[f"fi_{fr}"] = int(len(spikes(t, vs, t0, t0 + FIW)))
    return out


def derived(o, p=""):
    """Kampa-style ratios from one measure() dict."""
    r = {}
    for f in ("f100", "f150", "f200"):
        r[p + "R_" + f] = o[f"ca_{f}_dist"] / o["ca_f66_dist"]
    for g in GROUPS:
        r[p + "R3_" + g] = o[f"ca_r3_{g}"] / o[f"ca_ap1_{g}"]
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--cands", required=True, help="json list of {name, params}")
    ap.add_argument("--out", required=True)
    ap.add_argument("--no-pharm", action="store_true")
    args = ap.parse_args()
    import bluecellulab
    bluecellulab.neuron.load_mechanisms(MECH)
    cands = json.load(open(args.cands))
    rows = []
    for pair in args.pairs.split(","):
        sim, cell, segs, rec, amps = build(pair)
        for c in cands:
            t0 = time.time()
            bd.apply(segs, c["params"])
            sim.run(T_END, cvode=False, dt=0.025)
            o = measure(rec, T_END); row = dict(pair=pair, cand=c["name"], **c["params"], **o, **derived(o))
            if not args.no_pharm:
                for tag, kw in (("ka_off", dict(ka_off=True)), ("lva_off", dict(lva_off=True))):
                    bd.apply(segs, c["params"], **kw)
                    sim.run(T_PHARM, cvode=False, dt=0.025)
                    q = measure(rec, T_PHARM, pharm=True)
                    row.update({f"{tag}_{k}": v for k, v in q.items() if k.startswith(("ca_", "nsp_"))})
                    row.update(derived(q, tag + "_"))
            row["sec"] = time.time() - t0
            rows.append(row)
            print(pair, c["name"], {k: (round(row[k], 2) if isinstance(row[k], float) else row[k]) for k in
                  ("bap_dv_d50_100", "bap_dv_d200", "R_f100", "R_f200", "R3_prox", "R3_dist", "nsp_n100",
                   "nsp_l200", "nsp_n50", "fi_0.5", "sec")}, flush=True)
            pd.DataFrame(rows).to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
