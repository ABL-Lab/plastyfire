"""[ca_vs_data copy of cell_audit/bap_ca_map/bap_ca_map.py] Adds: --glob (GluSynapse GLOBALs, e.g. ljp25_g031), 3ap_100hz /
3ap_200hz protocols, and the 1/e decay time (from the peak) of spine cai_CR and shaft cai per site.

bAP-evoked Ca vs path distance along basal and apical dendrites (og-delta / antic-delta), spine (GluSynapse_v2 cai_CR) and shaft.

Real post cell (bluecellulab, delta circuit, SK_E2 off, as diag_burst.py/sims.py), somatic pulses only. A GluSynapse
point process with no NetCon is put every ~20 um path distance on selected basal terminal paths and on the apical
trunk + oblique paths. Two sims per protocol: with spines (v, shaft cai, cai_CR, ica_VDCC) and without (unperturbed shaft cai, v).

  python bap_ca_map.py --cell 181015-184976 --variant delta --out out_L5_181015-184976_delta.json
  python bap_ca_map.py --cell 143065 --l23 --regions basal ...
"""
import argparse, json, os, sys
import numpy as np

G = "/project/rrg-emuller/dhuruva/plastyfire/glusynapse_v2"
sys.path.insert(0, G)
import diag_burst as db

L5 = os.path.join(db.SIMS)
L23 = "/project/rrg-emuller/dhuruva/plastyfire/refitting_results/fitting/n120/seed20262009/Ebner2019_L23PC_L23PC_basal/simulations"
POP = "S1nonbarrel_neurons"
PROTOS = {"1ap": ([0.0], "nevian_1ap_dt+10ms"), "3ap_50hz": ([0.0, 20.0, 40.0], "nevian_3ap_50hz_dt+10ms"),
          "5ap_20hz": ([0.0, 50.0, 100.0, 150.0, 200.0], "nevian_3ap_20hz_dt+10ms"),
          "3ap_100hz": ([0.0, 10.0, 20.0], "nevian_3ap_100hz_dt+10ms"), "3ap_200hz": ([0.0, 5.0, 10.0], "letzkus_3ap_200hz_dt+10ms")}
GLOB = {}        # GluSynapse GLOBALs set in each run (--glob)


def tau_1e(t, x, w):
    """time from the peak (inside window w) to the first fall below peak/e; x is rest-subtracted."""
    i0 = np.where(w)[0][0]; ip = i0 + int(np.argmax(x[w])); pk = x[ip]
    if pk <= 0: return float("nan")
    j = np.where(x[ip:] < pk / np.e)[0]
    return float(t[ip + j[0]] - t[ip]) if len(j) else float("nan")
BIN = 20.0
VOL = 0.153      # um3, median volume_CR of the 68 real L5->L5 synapses (NMDA_VDCC_CALIBRATION.md)


def sec_path(h, sec, soma):
    out = []
    while sec is not None and sec != soma:
        out.append(sec); sr = h.SectionRef(sec=sec); sec = sr.parent if sr.has_parent() else None
    return out[::-1]


def terminals(h, secs):
    return [s for s in secs if h.SectionRef(sec=s).nchild() == 0]


def pick_sites(h, cell, secs, nmax, exclude_first=False):
    """Greedy paths (terminals) maximising new segments; then one segment per 20 um bin on each path."""
    soma = cell.soma
    dist = lambda seg: h.distance(soma(0.5), seg)
    tips = sorted(terminals(h, secs), key=lambda s: -dist(s(1.0)))
    paths, covered = [], set()
    for _ in range(nmax):
        best = max(tips, key=lambda s: sum(1 for p in sec_path(h, s, soma) if p not in covered) * 1.0 + 1e-6 * dist(s(1.0)))
        p = sec_path(h, best, soma)
        if not p or all(x in covered for x in p): break
        paths.append((best.name(), p)); covered.update(p)
    sites, seen = [], set()
    for name, p in paths:
        segs = [(seg, dist(seg)) for s in p for seg in s]
        for c in np.arange(BIN / 2, max(d for _, d in segs), BIN):
            seg, d = min(segs, key=lambda x: abs(x[1] - c))
            if abs(d - c) > BIN / 2 or (seg.sec.name(), seg.x) in seen: continue
            seen.add((seg.sec.name(), seg.x)); sites.append((seg, d, name))
    return sites


def run(cell_id, is23, variant, proto, regions, cfgp, amp, width, with_spines):
    import bluecellulab
    from bluecellulab import CircuitSimulation
    h = bluecellulab.neuron.h
    post = int(cell_id.split("-")[-1])
    sim = CircuitSimulation(cfgp)
    sim.instantiate_gids([(POP, post)], add_synapses=False, add_minis=False, add_pulse_stimuli=False)
    cell = sim.cells[(POP, post)]
    for sec in cell.basal:
        for k, val in db.VARIANTS[variant].items(): setattr(sec, k, val)
    h.cao_CR_GluSynapse = 2.0
    for k, val in GLOB.items():
        setattr(h, f"{k}_GluSynapse", val)
    sites = []
    if "basal" in regions: sites += [("basal",) + s for s in pick_sites(h, cell, list(cell.basal), 6)]
    if "apical" in regions: sites += [("apical",) + s for s in pick_sites(h, cell, list(cell.apical), 4)]
    ds = PROTOS[proto][0]
    for d in ds:
        s = h.TStim(0.5, sec=cell.soma); s.train(db.T0 + d, width, amp, 0.1, width); cell.persistent.append(s)
    rec = dict(t=h.Vector(), vs=h.Vector()); rec["t"].record(h._ref_t, 0.025); rec["vs"].record(cell.soma(0.5)._ref_v, 0.025)
    R = []
    for kind, seg, d, path in sites:
        r = {}
        refs = [("v", seg._ref_v), ("cai", seg._ref_cai)]
        if with_spines:
            syn = h.GluSynapse(seg); syn.volume_CR = VOL; cell.persistent.append(syn)
            refs += [("cacr", syn._ref_cai_CR), ("ica_vdcc", syn._ref_ica_VDCC)]
        for k, ref in refs:
            r[k] = h.Vector(); r[k].record(ref, 0.025)
        R.append((kind, seg.sec.name(), seg.x, d, path, r))
    sim.run(db.T0 + ds[-1] + 250.0, cvode=False, dt=0.025)   # 250 (orig 100): room for the shaft 1/e decay
    t = np.array(rec["t"]); vs = np.array(rec["vs"])
    up = np.where((vs[1:] > -10) & (vs[:-1] <= -10) & (t[1:] >= db.T0 - 5))[0]
    out = dict(soma_spikes=(t[1:][up] - db.T0).round(2).tolist(), sites=[])
    w = (t >= db.T0 - 1) & (t < db.T0 + ds[-1] + 60.0); pre = t < db.T0
    for kind, sn, x, d, path, r in R:
        a = {k: np.array(v) for k, v in r.items()}
        e = dict(kind=kind, sec=sn, x=x, dist=float(d), path=path, v_rest=float(a["v"][pre][-1]),
                 vpk=float(a["v"][w].max() - a["v"][pre][-1]), vmax=float(a["v"][w].max()),
                 cai_rest=float(a["cai"][pre][-1]), cai_pk=float(a["cai"][w].max() - a["cai"][pre][-1]))
        wt = t >= db.T0 - 1
        e["tau_cai"] = tau_1e(t, a["cai"] - a["cai"][pre][-1], wt[:len(a["cai"])])
        if with_spines:
            e.update(cacr_rest=float(a["cacr"][pre][-1]), cacr_pk=float(a["cacr"][w].max() - a["cacr"][pre][-1]),
                     vdcc_int=float(-np.trapz(a["ica_vdcc"][w], t[w])),
                     tau_cacr=tau_1e(t, a["cacr"] - a["cacr"][pre][-1], wt[:len(a["cacr"])]))
            # per-spike peaks of local V (bAP amplitude of each spike)
            e["vpk_spike"] = [float(a["v"][(t >= db.T0 + d0) & (t < db.T0 + d0 + 20)].max() - a["v"][pre][-1]) for d0 in ds]
        out["sites"].append(e)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", required=True, help="L5 pair '181015-184976' (post = 2nd) or L2/3 post gid with --l23")
    ap.add_argument("--l23", action="store_true")
    ap.add_argument("--variant", default="delta")
    ap.add_argument("--protos", default=",".join(PROTOS))
    ap.add_argument("--regions", default="basal,apical")
    ap.add_argument("--glob", default="{}", help="json {GluSynapse GLOBAL: value}, e.g. {\"ljp_VDCC\": 25, \"gca_bar_VDCC\": 0.02306}")
    ap.add_argument("--tag", default="", help="label stored in the output (e.g. ljp25_g031)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    GLOB.update(json.loads(a.glob))
    import bluecellulab
    bluecellulab.neuron.load_mechanisms(db.MECH)
    base = os.path.join(L23, f"pip0-{a.cell}") if a.l23 else os.path.join(L5, a.cell)
    cfgp0 = os.path.join(base, "nevian_3ap_50hz_dt+10ms", "prefire_simulation_config.json")
    pulses = {}
    for p in a.protos.split(","):
        pc = json.load(open(os.path.join(base, PROTOS[p][1], "prefire_simulation_config.json")))["inputs"]
        k = next(k for k, v in pc.items() if v.get("input_type") == "current_clamp")
        pulses[p] = (pc[k]["amp_start"], pc[k]["width"])
    cfg = json.load(open(cfgp0))
    for k in ("reports", "inputs"): cfg.pop(k, None)
    cfg["run"]["tstop"] = db.T0 + 400.0
    tmp = os.path.join(os.environ.get("SLURM_TMPDIR", os.path.dirname(os.path.abspath(a.out))), f".bapmap_{a.cell}_{os.getpid()}.json")
    json.dump(cfg, open(tmp, "w"))
    import multiprocessing as mp
    ctx = mp.get_context("fork")
    res = []
    for p in a.protos.split(","):
        r = {}
        for ws in (True, False):
            with ctx.Pool(1, maxtasksperchild=1) as pool:
                r[ws] = pool.apply(run, (a.cell, a.l23, a.variant, p, a.regions.split(","), tmp, *pulses[p], ws))
        for s, s0 in zip(r[True]["sites"], r[False]["sites"]):
            s["shaft_cai_pk_nospine"] = s0["cai_pk"]; s["vpk_nospine"] = s0["vpk"]; s["tau_cai_nospine"] = s0["tau_cai"]
        print(a.cell, a.variant, p, "pulse", pulses[p], "soma spikes", r[True]["soma_spikes"], "n sites", len(r[True]["sites"]), flush=True)
        res.append(dict(proto=p, pulse=pulses[p], soma_spikes=r[True]["soma_spikes"], sites=r[True]["sites"]))
    os.remove(tmp)
    json.dump(dict(cell=a.cell, l23=a.l23, variant=a.variant, glob=GLOB, tag=a.tag, results=res), open(a.out, "w"))


if __name__ == "__main__":
    main()
