"""Synapse-local T-drive test sims (LOCAL_T_DRIVE.md). Test-only; model code is unchanged.

diag_burst.py setup (real post cell, real synapses, pair's own pulse config, pre spike at every synapse at pre+0.1 ms),
but it keeps the full per-synapse traces of the allowed local signals, at 0.1 ms from T0-20 ms:
    v (local membrane V of the synapse segment), cai (shaft Ca of the segment), cacr (spine cai_CR), effcai_GB,
    ica_vdcc, ica_nmda (GluSynapse, nA)
Protocols (T0 = first post pulse; dt>0 = pre first). The eCB gate reads T at the pre arrival, so a post-only run
gives T for every post-before-pre single pairing and for the burst rows:
    ap1        one post AP, no pre (Sjostrom 2003 0.1 Hz, read at +10/+25/+50/+100/+120/+200 ms)
    burst      5 post APs at 20 Hz, no pre (read 120/200 ms after the last)
    epsp       one pre spike, no post (does the own EPSP trip the drive?)
    sj{10,20,50}@{+10,-10}  Sjostrom 2001: 5 pairings at f
    python glusynapse_v2/local_t/sims.py --pair P --variants og,s5 --out X.npz
"""
import argparse, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..")))
import diag_burst as D

SYNG = {"og": {}, "s5": {"ljp_VDCC": 5.0}}          # s5 = NMDA_VDCC_CALIBRATION.md recommended fix
PROTOS = {}                                          # name -> (post pulse times, pre times, tend, pulse cfg)
PROTOS["ap1"] = ([0.0], [], 300.0, "sjostrom_0.1hz_dt-25ms")
PROTOS["burst"] = ([0.0, 50.0, 100.0, 150.0, 200.0], [], 450.0, "sjostrom_burst5x20hz_r50_dt-120ms")
PROTOS["epsp"] = ([], [0.0], 200.0, "sjostrom_0.1hz_dt+10ms")
for f in ("10", "20", "50"):
    post = [k * 1000.0 / float(f) for k in range(5)]
    for dt in (10, -10):
        PROTOS[f"sj{f}@{dt:+d}"] = (post, [x - dt for x in post], post[-1] + 100.0, f"sjostrom_{f}hz_dt+10ms")
DS = 4                                               # 0.025 ms -> 0.1 ms


def run(pair, variant, proto, sim_cfg, pre, post, amp, width):
    import bluecellulab
    from bluecellulab import CircuitSimulation
    h = bluecellulab.neuron.h
    T0 = D.T0
    posts, pres, tend, _ = PROTOS[proto]
    sim = CircuitSimulation(sim_cfg)
    pop = "S1nonbarrel_neurons"
    train = [T0 + x for x in pres] or None
    sim.instantiate_gids([(pop, post)], add_synapses=True, add_minis=False, add_pulse_stimuli=False,
                         intersect_pre_gids=[(pop, pre)], pre_spike_trains=None if train is None else {(pop, pre): train})
    cell = sim.cells[(pop, post)]
    for k, val in SYNG[variant].items():
        setattr(h, f"{k}_GluSynapse", val)
    if train:
        for c in cell.connections.values():
            if c.post_netcon is not None:
                c.post_netcon.delay = 0.1
    stims = []
    for d in posts:
        s = h.TStim(0.5, sec=cell.soma); s.train(T0 + d, width, amp, 0.1, width); stims.append(s)
    cell.persistent.extend(stims)
    t = h.Vector(); t.record(h._ref_t, 0.025); vs = h.Vector(); vs.record(cell.soma(0.5)._ref_v, 0.025)
    syns = []
    for sid, syn in cell.synapses.items():
        hs = syn.hsynapse; x = hs.get_loc(); sec = h.cas(); h.pop_section(); seg = sec(x)
        name = sec.name(); kind = "basal" if ".dend" in name else ("apical" if ".apic" in name else "other")
        r = {}
        for key, ref in (("v", seg._ref_v), ("cai", seg._ref_cai), ("cacr", hs._ref_cai_CR), ("effcai", hs._ref_effcai_GB),
                         ("ica_vdcc", hs._ref_ica_VDCC), ("ica_nmda", hs._ref_ica_NMDA)):
            r[key] = h.Vector(); r[key].record(ref, 0.025)
        syns.append((int(sid[1]), kind, float(h.distance(cell.soma(0.5), seg)), r))
    sim.run(T0 + tend, cvode=False, dt=0.025)
    t = np.array(t); vs = np.array(vs)
    up = np.where((vs[1:] > -10) & (vs[:-1] <= -10))[0]
    spk = t[1:][up] - T0
    w = np.where(t >= T0 - 20)[0][0]
    out = dict(t=(t[w::DS] - T0).astype(np.float32), soma_spikes=spk[spk > -20], sid=np.array([s[0] for s in syns]),
               kind=np.array([s[1] for s in syns]), dist=np.array([s[2] for s in syns]))
    for key in ("v", "cai", "cacr", "effcai", "ica_vdcc", "ica_nmda"):
        out[key] = np.array([np.array(s[3][key])[w::DS] for s in syns], np.float32)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pair", required=True)
    ap.add_argument("--variants", default="og,s5")
    ap.add_argument("--protos", default=",".join(PROTOS))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    import bluecellulab
    bluecellulab.neuron.load_mechanisms(D.MECH)
    base = os.path.join(D.SIMS, a.pair, "nevian_3ap_50hz_dt+10ms", "prefire_simulation_config.json")
    cfg = json.load(open(base))
    for k in ("reports", "inputs"):
        cfg.pop(k, None)
    cfg["run"]["tstop"] = D.T0 + max(p[2] for p in PROTOS.values())
    tmp = os.path.join(os.environ.get("SLURM_TMPDIR", "/tmp"), f"local_t_{a.pair}_{os.getpid()}.json")
    json.dump(cfg, open(tmp, "w"))
    pre, post = (int(x) for x in a.pair.split("-"))
    res = {}
    import multiprocessing as mp
    ctx = mp.get_context("fork")
    for v in a.variants.split(","):
        for p in a.protos.split(","):
            pc = json.load(open(os.path.join(D.SIMS, a.pair, PROTOS[p][3], "prefire_simulation_config.json")))["inputs"]
            k = next(k for k, x in pc.items() if x.get("input_type") == "current_clamp")
            with ctx.Pool(1, maxtasksperchild=1) as pool:
                r = pool.apply(run, (a.pair, v, p, tmp, pre, post, pc[k]["amp_start"], pc[k]["width"]))
            print(f"{a.pair} {v} {p}: soma spikes {np.round(r['soma_spikes'], 1).tolist()} n_syn {len(r['sid'])} "
                  f"vpk median {np.median(r['v'].max(1)):.1f}", flush=True)
            for key, val in r.items():
                res[f"{v}|{p}|{key}"] = val
    np.savez_compressed(a.out, **res)


if __name__ == "__main__":
    main()
