"""Why does a bAP burst give no more synaptic Ca than one bAP in the delta emodel? Single-cell dissection.

The real post cell of a pair (bluecellulab, the pair's own prefire config, delta circuit, SK_E2 off as in the runs),
the real synapse locations, the same somatic pulses as the induction runs (config pulse amp/width). No pre spikes.

  soma:       spikes per burst (does the soma even fire 3 APs?)
  synapse:    per spike peak local v, peak and integral of GluSynapse ica_VDCC (= the fit's vdcc signal),
              dendritic cai (cad), kBK open prob / ik, NaTg ina, Ka_kampa gka at the synapse segment
  variants:   og-delta vs antic-delta basal Na/Ka (the only difference between the two hocs), singly and together

    python glusynapse_v2/diag_burst.py --pair 186261-208490 --out glusynapse_v2/results/diag_burst_186261.json
"""
import argparse, json, os, sys
import numpy as np

SIMS = ("/project/rrg-emuller/dhuruva/plastyfire/refitting_results/fitting/n120/seed20262009/"
        "Sabrina_L5TTPC_L5TTPC_STDP/simulations")
MECH = "/project/rrg-emuller/dhuruva/DEES_cell_packages/"
# og-delta (the circuit's emodel) vs antic-delta: the two hocs differ only in basal Na and Ka (Ca, BK identical)
ANTIC_BASAL = dict(gNaTgbar_NaTg=0.0143, gbar_Ka_kampa=0.025)          # antic_delta_values.json; og-delta 0.003, 0.002
VARIANTS = {
    "delta": {},
    "antic_na": dict(gNaTgbar_NaTg=ANTIC_BASAL["gNaTgbar_NaTg"]),
    "antic_ka": dict(gbar_Ka_kampa=ANTIC_BASAL["gbar_Ka_kampa"]),
    "antic_delta": dict(ANTIC_BASAL),
}
# each post protocol uses its own induction config's pulse amplitude/width
PULSE_CFG = {"1ap": "nevian_1ap_dt+10ms", "3ap_20hz": "nevian_3ap_20hz_dt+10ms", "3ap_50hz": "nevian_3ap_50hz_dt+10ms", "2ap_50hz": "nevian_2ap_50hz_dt+10ms",
             "3ap_100hz": "nevian_3ap_100hz_dt+10ms", "3ap_200hz": "letzkus_3ap_200hz_dt+10ms",
             "epsp": "nevian_3ap_50hz_dt+10ms"}
BD = {}                                                   # delta-bd candidates (--bd-cands): name -> params
SYNG = {}                            # spine variants (--syn-cands): name -> {"glob": {GluSynapse GLOBAL: value}, "bd": params}
PROTOS = {"1ap": [0.0], "3ap_20hz": [0.0, 50.0, 100.0], "3ap_50hz": [0.0, 20.0, 40.0], "2ap_50hz": [0.0, 20.0], "3ap_100hz": [0.0, 10.0, 20.0],
          "3ap_200hz": [0.0, 5.0, 10.0], "epsp": []}
# pairings: "<post proto>@<dt>" = one pre spike arriving dt ms before the first somatic pulse (dt>0: pre first)
for _p in ("1ap", "3ap_50hz", "3ap_100hz", "3ap_200hz"):
    for _dt in (10, -10):
        PROTOS[f"{_p}@{_dt:+d}"] = PROTOS[_p]
TWIN = 120.0
T0 = 1000.0
# optional per-proto extensions (set by wrappers, e.g. cell_audit/basal_nakv/grid.py); defaults keep old behaviour
PRE = {}        # proto -> pre spike times rel. T0 (a train); overrides the '@dt' single spike
TEND = {}       # proto -> sim length after T0 (ms); default 200, measurement window then runs to the end
PULSE_OVR = {}  # proto -> (amp nA, width ms)


def run(pair, variant, proto, sim_cfg, pre, post, amp, width):
    import bluecellulab
    from bluecellulab import CircuitSimulation
    h = bluecellulab.neuron.h
    sim = CircuitSimulation(sim_cfg)
    pop = "S1nonbarrel_neurons"
    base, _, dt = proto.partition("@")
    pre_t = None
    if proto == "epsp":
        pre_t = T0
    elif dt:
        pre_t = T0 - float(dt)
    pre_train = None if pre_t is None else [pre_t]
    if proto in PRE:
        pre_train = [T0 + x for x in PRE[proto]]; pre_t = pre_train[0]
    tend = TEND.get(proto, 200.0)
    sim.instantiate_gids([(pop, post)], add_synapses=True, add_minis=False, add_pulse_stimuli=False,
                         intersect_pre_gids=[(pop, pre)],
                         pre_spike_trains=None if pre_train is None else {(pop, pre): pre_train})
    cell = sim.cells[(pop, post)]
    if variant in SYNG:                                      # spine VDCC variants: GluSynapse GLOBALs (+ optional delta-bd)
        for k, val in SYNG[variant].get("glob", {}).items():
            setattr(h, f"{k}_GluSynapse", val)
        if SYNG[variant].get("bd"):
            BD[variant] = SYNG[variant]["bd"]
    if variant in BD:                                        # delta-bd: distance-dependent basal densities
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "emodel_bd"))
        import bd
        bd.apply(bd.basal_segments(cell, h), BD[variant])
    else:
        for sec in cell.basal:                               # the delta edits are basal-only
            for k, val in VARIANTS.get(variant, {}).items():
                setattr(sec, k, val)
    stims = []
    if pre_t is not None:                   # drop the axonal delay: the pre spike reaches every synapse at pre_t (+0.1)
        for c in cell.connections.values():
            if c.post_netcon is not None:
                c.post_netcon.delay = 0.1
    for d in PROTOS[proto]:
        s = h.TStim(0.5, sec=cell.soma); s.train(T0 + d, width, amp, 0.1, width); stims.append(s)
    cell.persistent.extend(stims)
    rec = dict(t=h.Vector(), vs=h.Vector())
    rec["t"].record(h._ref_t, 0.025); rec["vs"].record(cell.soma(0.5)._ref_v, 0.025)
    syns = []
    for sid, syn in cell.synapses.items():
        hs = syn.hsynapse; x = hs.get_loc(); sec = h.cas(); h.pop_section(); seg = sec(x)
        name = sec.name(); kind = "basal" if ".dend" in name else ("apical" if ".apic" in name else "other")
        r = {}
        for key, ref in (("v", seg._ref_v), ("ica_vdcc", hs._ref_ica_VDCC), ("cai", seg._ref_cai),
                         ("ica_nmda", hs._ref_ica_NMDA), ("cacr", hs._ref_cai_CR), ("effcai", hs._ref_effcai_GB)):
            r[key] = h.Vector(); r[key].record(ref, 0.025)
        for key, attr in (("bk_p", "_ref_p_kBK"), ("bk_ik", "_ref_ik_kBK"), ("ina", "_ref_ina_NaTg"),
                          ("gka", "_ref_gka_Ka_kampa")):
            if hasattr(seg, attr):
                r[key] = h.Vector(); r[key].record(getattr(seg, attr), 0.025)
        dist = h.distance(cell.soma(0.5), seg)
        syns.append((sid[1], kind, dist, r))
    sim.run(T0 + tend, cvode=False, dt=0.025)
    t = np.array(rec["t"]); vs = np.array(rec["vs"])
    w = t >= T0 - 5
    up = np.where((vs[w][1:] > -10) & (vs[w][:-1] <= -10))[0]
    spk = (t[w][1:][up] - T0).round(2).tolist()
    out = dict(pair=pair, variant=variant, proto=proto, soma_spikes=spk, syn=[])
    # windows around each pulse (spike): [d, d + gap) for per-spike peaks
    ds = PROTOS[proto] or [0.0]; ends = ds[1:] + [ds[-1] + 60.0]
    for sid, kind, dist, r in syns:
        a = {k: np.array(v) for k, v in r.items()}
        per = []
        for d0, d1 in zip(ds, ends):
            m = (t >= T0 + d0) & (t < T0 + d1)
            e = dict(vpk=float(a["v"][m].max()), ica_pk=float(-a["ica_vdcc"][m].min()),
                     ica_int=float(-np.trapz(a["ica_vdcc"][m], t[m])), cai_pk=float(a["cai"][m].max()))
            if "bk_p" in a: e["bk_p_pk"] = float(a["bk_p"][m].max()); e["bk_p_start"] = float(a["bk_p"][m][0])
            if "bk_ik" in a: e["bk_ik_pk"] = float(a["bk_ik"][m].max())
            if "ina" in a: e["ina_pk"] = float(-a["ina"][m].min())
            if "gka" in a: e["gka_pk"] = float(a["gka"][m].max())
            per.append(e)
        mt = (t >= T0) & (t < T0 + ends[-1])
        twin = TWIN if proto not in TEND else tend
        mw = (t >= T0 - 60) & (t < T0 + twin)              # whole event, both orders
        ext = dict(nmda_int=float(-np.trapz(a["ica_nmda"][mw], t[mw])), vdcc_int=float(-np.trapz(a["ica_vdcc"][mw], t[mw])),
                   cacr_pk=float(a["cacr"][mw].max()), effcai_end=float(a["effcai"][t < T0 + twin][-1]),
                   effcai_pk=float(a["effcai"][mw].max()), effcai_int=float(np.trapz(a["effcai"][mw], t[mw])),
                   vpk=float(a["v"][mw].max()))
        out["syn"].append(dict(sid=int(sid), kind=kind, dist=float(dist), v_rest=float(a["v"][t < T0][-1]),
                               ica_total=float(-np.trapz(a["ica_vdcc"][mt], t[mt])), per_spike=per, **ext))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pair", required=True)
    ap.add_argument("--variants", default=",".join(VARIANTS))
    ap.add_argument("--protos", default=",".join(PROTOS))
    ap.add_argument("--bd-cands", help="json list of delta-bd {name, params}; their names become variants")
    ap.add_argument("--syn-cands", help="json list of {name, glob: {GluSynapse GLOBAL: value}, bd: params|null}")
    ap.add_argument("--width", type=float, help="override the config pulse width (ms)")
    ap.add_argument("--amp", type=float, help="override the config pulse amplitude (nA)")
    ap.add_argument("--var-cands", help="json {name: {basal section attr: value}}: uniform basal variants added to VARIANTS")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    import bluecellulab
    bluecellulab.neuron.load_mechanisms(MECH)
    if args.bd_cands:
        BD.update({c["name"]: c["params"] for c in json.load(open(args.bd_cands))})
    if args.var_cands:
        VARIANTS.update(json.load(open(args.var_cands)))
    if args.syn_cands:
        SYNG.update({c["name"]: c for c in json.load(open(args.syn_cands))})
    cfgp = os.path.join(SIMS, args.pair, "nevian_3ap_50hz_dt+10ms", "prefire_simulation_config.json")
    cfg = json.load(open(cfgp))
    pulses = {}
    for base, d in PULSE_CFG.items():
        f = os.path.join(SIMS, args.pair, d, "prefire_simulation_config.json")
        pc = json.load(open(f if os.path.exists(f) else cfgp))["inputs"]
        k = next(k for k, v in pc.items() if v.get("input_type") == "current_clamp")
        pulses[base] = (args.amp or pc[k]["amp_start"], args.width or pc[k]["width"])
    # bluecellulab cannot build synapse reports; stimuli are added here by hand (as simulator_edges does)
    for k in ("reports", "inputs"):
        cfg.pop(k, None)
    cfg["run"]["tstop"] = T0 + max([200.0, *TEND.values()])
    tmp = os.path.join(os.environ.get("SLURM_TMPDIR", "/tmp"), f"diag_burst_{args.pair}_{os.getpid()}.json")
    json.dump(cfg, open(tmp, "w")); cfgp = tmp
    pre, post = (int(x) for x in args.pair.split("-"))
    res = []
    import multiprocessing as mp
    ctx = mp.get_context("fork")                 # a fresh NEURON state per run: cells are not freed otherwise
    for v in args.variants.split(","):
        for p in args.protos.split(","):
            with ctx.Pool(1, maxtasksperchild=1) as pool:
                r = pool.apply(run, (args.pair, v, p, cfgp, pre, post, *PULSE_OVR.get(p, pulses[p.split("@")[0]])))
            s = r["syn"]
            print(f"   nmda {np.median([x['nmda_int'] for x in s]):.3g} vdcc {np.median([x['vdcc_int'] for x in s]):.3g} "
                  f"cacr_pk {np.median([x['cacr_pk'] for x in s]):.3g} effcai_pk {np.median([x['effcai_pk'] for x in s]):.3g}", flush=True)
            ps = np.array([[e["vpk"] for e in x["per_spike"]] for x in s])
            ic = np.array([[e["ica_int"] for e in x["per_spike"]] for x in s])
            print(f"{args.pair} {v:9s} {p:9s} soma spikes {r['soma_spikes']}  syn vpk per spike (median) "
                  f"{np.round(np.median(ps, 0), 1).tolist()}  ica charge per spike rel. to 1st "
                  f"{np.round(np.median(ic / ic[:, :1], 0), 2).tolist()}  total {np.median([x['ica_total'] for x in s]):.3g}",
                  flush=True)
            res.append(r)
    json.dump(dict(pair=args.pair, pulses=pulses, results=res), open(args.out, "w"))


if __name__ == "__main__":
    main()
