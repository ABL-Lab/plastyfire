"""Letzkus 2006 audit: does the L5TTPC post cell give a distance-dependent burst Ca boost along the apical dendrite?

Same cell, circuit config and somatic pulses as diag_burst.py (post cell of an L5-L5 pair, no pre spikes), but the
recording sites are placed along the apical trunk/tuft path at chosen path distances. At each site a GluSynapse is
cloned from one of the pair's real synapses (all PARAMETERs copied), so `ica_VDCC` is the same spine VDCC signal the
fit uses. Also recorded: local v, cai, and the dendritic Ca currents of each channel present (Ca_HVA2, Ca_LVAst, ...).

    python glusynapse_v2/cell_audit/letzkus_ca.py --pair 181015-184976 --out glusynapse_v2/cell_audit/letzkus_ca_181015-184976.json
"""
import argparse, json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import diag_burst as db

TARGETS = [50, 150, 250, 350, 450, 550, 650, 750, 900]
PROTOS = ["1ap", "2ap_50hz", "3ap_200hz", "3ap_100hz", "3ap_50hz"]


def trunk_path(cell, h):
    """Segments on the path soma -> farthest apical tip, with path distance."""
    far = max(cell.apical, key=lambda s: h.distance(cell.soma(0.5), s(1.0)))
    secs = []
    sec = far
    while sec is not None and sec != cell.soma:
        secs.append(sec)
        sr = h.SectionRef(sec=sec)
        sec = sr.parent if sr.has_parent() else None
    segs = [(seg, h.distance(cell.soma(0.5), seg)) for s in secs[::-1] for seg in s]
    return segs


def clone_syn(h, tmpl, seg):
    new = h.GluSynapse(seg)
    ms = h.MechanismStandard("GluSynapse", 1)
    sref = h.ref("")
    for i in range(int(ms.count())):
        ms.name(sref, i)
        nm = sref[0]
        try:
            setattr(new, nm, getattr(tmpl, nm))
        except Exception:
            pass
    return new


def run(variant, proto, sim_cfg, pre, post, amp, width):
    import bluecellulab
    from bluecellulab import CircuitSimulation
    h = bluecellulab.neuron.h
    sim = CircuitSimulation(sim_cfg)
    pop = "S1nonbarrel_neurons"
    sim.instantiate_gids([(pop, post)], add_synapses=True, add_minis=False, add_pulse_stimuli=False,
                         intersect_pre_gids=[(pop, pre)])
    cell = sim.cells[(pop, post)]
    for sec in cell.basal:
        for k, val in db.VARIANTS.get(variant, {}).items():
            setattr(sec, k, val)
    tmpl = next(iter(cell.synapses.values())).hsynapse
    path = trunk_path(cell, h)
    stims = []
    for d in db.PROTOS[proto]:
        s = h.TStim(0.5, sec=cell.soma); s.train(db.T0 + d, width, amp, 0.1, width); stims.append(s)
    cell.persistent.extend(stims)
    rec = dict(t=h.Vector(), vs=h.Vector())
    rec["t"].record(h._ref_t, 0.025); rec["vs"].record(cell.soma(0.5)._ref_v, 0.025)
    sites, mechs = [], None
    for tgt in TARGETS:
        seg, dist = min(path, key=lambda p: abs(p[1] - tgt))
        if abs(dist - tgt) > 40:
            continue
        syn = clone_syn(h, tmpl, seg); cell.persistent.append(syn)
        r = {}
        for key, ref in (("v", seg._ref_v), ("cai", seg._ref_cai), ("ica_vdcc", syn._ref_ica_VDCC), ("ica", seg._ref_ica)):
            r[key] = h.Vector(); r[key].record(ref, 0.025)
        dm = list(seg.sec.psection()["density_mechs"].keys())
        mechs = mechs or {}
        mechs[round(dist)] = dm
        for m in dm:
            if m.startswith("Ca") and hasattr(seg, f"_ref_ica_{m}"):
                r[f"ica_{m}"] = h.Vector(); r[f"ica_{m}"].record(getattr(seg, f"_ref_ica_{m}"), 0.025)
        sites.append((dist, seg.sec.name(), r))
    for sid, sy in cell.synapses.items():                  # the pair's real synapses (L23->L5: basal .. distal tuft)
        hs = sy.hsynapse; x = hs.get_loc(); sec = h.cas(); h.pop_section(); seg = sec(x)
        r = {}
        for key, ref in (("v", seg._ref_v), ("cai", seg._ref_cai), ("ica_vdcc", hs._ref_ica_VDCC), ("ica", seg._ref_ica)):
            r[key] = h.Vector(); r[key].record(ref, 0.025)
        for m in ("Ca_HVA2", "Ca_LVAst"):
            if hasattr(seg, f"_ref_ica_{m}"):
                r[f"ica_{m}"] = h.Vector(); r[f"ica_{m}"].record(getattr(seg, f"_ref_ica_{m}"), 0.025)
        sites.append((h.distance(cell.soma(0.5), seg), "syn:" + sec.name(), r))
    sim.run(db.T0 + 200.0, cvode=False, dt=0.025)
    t = np.array(rec["t"]); vs = np.array(rec["vs"])
    w = t >= db.T0 - 5
    up = np.where((vs[w][1:] > -10) & (vs[w][:-1] <= -10))[0]
    spk = (t[w][1:][up] - db.T0).round(2).tolist()
    ds = db.PROTOS[proto]; ends = ds[1:] + [ds[-1] + 60.0]
    mw = (t >= db.T0 - 1) & (t < db.T0 + db.TWIN)
    out = dict(variant=variant, proto=proto, soma_spikes=spk, mechs=mechs, sites=[])
    for dist, name, r in sites:
        a = {k: np.array(v) for k, v in r.items()}
        vr = float(a["v"][t < db.T0][-1])
        per = []
        for d0, d1 in zip(ds, ends):
            m = (t >= db.T0 + d0) & (t < db.T0 + d1)
            per.append(dict(vpk=float(a["v"][m].max() - vr), ica_int=float(-np.trapz(a["ica_vdcc"][m], t[m]))))
        e = dict(dist=float(dist), sec=name, v_rest=vr, per_spike=per,
                 vdcc_int=float(-np.trapz(a["ica_vdcc"][mw], t[mw])),
                 ica_seg_int=float(-np.trapz(a["ica"][mw], t[mw])),
                 v_int=float(np.trapz(a["v"][mw] - vr, t[mw])),
                 cai_pk=float(a["cai"][mw].max()), cai_rest=float(a["cai"][t < db.T0][-1]))
        for k in a:
            if k.startswith("ica_C"):
                e[k + "_int"] = float(-np.trapz(a[k][mw], t[mw]))
        out["sites"].append(e)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pair", default="181015-184976")
    ap.add_argument("--variants", default="delta,antic_delta")
    ap.add_argument("--protos", default=",".join(PROTOS))
    ap.add_argument("--sims", default=db.SIMS)
    ap.add_argument("--base-proto", default="nevian_3ap_50hz_dt+10ms", help="workdir whose config/pulse is the fallback")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    import bluecellulab
    bluecellulab.neuron.load_mechanisms(db.MECH)
    base = os.path.join(args.sims, args.pair)
    cfgp = os.path.join(base, args.base_proto, "prefire_simulation_config.json")
    cfg = json.load(open(cfgp))
    pulses = {}
    for p in args.protos.split(","):
        f = os.path.join(base, db.PULSE_CFG[p], "prefire_simulation_config.json")
        pc = json.load(open(f if os.path.exists(f) else cfgp))["inputs"]
        k = next(k for k, v in pc.items() if v.get("input_type") == "current_clamp")
        pulses[p] = (pc[k]["amp_start"], pc[k]["width"])
    for k in ("reports", "inputs"):
        cfg.pop(k, None)
    cfg["run"]["tstop"] = db.T0 + 200.0
    tmp = os.path.join(os.environ.get("SLURM_TMPDIR", os.path.dirname(os.path.abspath(args.out))),
                       f".letzkus_{args.pair}_{os.getpid()}.json")
    json.dump(cfg, open(tmp, "w"))
    pre, post = (int(x) for x in args.pair.split("-"))
    import multiprocessing as mp
    ctx = mp.get_context("fork")
    res = []
    for v in args.variants.split(","):
        for p in args.protos.split(","):
            with ctx.Pool(1, maxtasksperchild=1) as pool:
                r = pool.apply(run, (v, p, tmp, pre, post, *pulses[p]))
            print(v, p, "soma", r["soma_spikes"], " ".join(
                f"{s['dist']:.0f}:{s['vdcc_int']:.3g}/{s['per_spike'][0]['vpk']:.0f}mV" for s in r["sites"] if not s["sec"].startswith("syn:")), flush=True)
            res.append(r)
    os.remove(tmp)
    json.dump(dict(pair=args.pair, pulses=pulses, results=res), open(args.out, "w"))


if __name__ == "__main__":
    main()
