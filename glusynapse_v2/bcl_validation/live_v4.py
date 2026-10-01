"""Live bluecellulab runs of a v4 fit (rho_redesign fit_v4 json) with the GluSynapseV4 mechanism.

One task = (path, pair, proto, cond, phase, express):
  path     L5 (Sabrina simulations, data/dhuruva_delta_circuit_config.json) or L23 (Ebner2019_L23PC_L5TTPC,
           data/dhuruva_delta_l23l5_circuit_config.json), the circuit configs of the fit's prefire traces.
  cond     control | mglu_block (A_mglu 0) | no_block (A_NO 0) | nmdar_block (bath APV: gmax_NMDA 0 on the pair's
           synapses from 500 ms before the induction to the snap, A_mglu = A_NO = 0, as batch_v2 nmdar_block).
  phase    prefire: the workdir's prefire_simulation_config.json (the fit's traces), final rho / dpre per synapse.
           full: C01 test pulses (4 min) + the same induction shifted by 4 min + C02 (4 min), built in memory from the
           prefire config exactly as simwriter writes simulation_config.json (C01_T = C02_T = the yaml T: 4000 ms
           Markram 10Hz_*, else 10000 ms). Snap 500 ms before the first C02 spike: rho -> 0/1, Use_GB =
           min(1, u0 (1 + dpre)), gmax_AMPA = gmax_d/p (the offline readout), then C02. EPSP ratio = mean C02 /
           mean C01 amplitude (ephysutils.get_epsp_vector, 100 ms window), as simulator_edges full_protocol.
  express  v4: expression follows the v4 rho_GB / dpre_GB during the induction (the live model).
           cooker: express_x_GB 1, expression follows the shadow v1 rule (delta-prefire cooker params,
           fit_results/cooker.json), so the induction reproduces the fit's prefire traces while the v4 rule runs
           passively (equivalence test). At the snap expression switches to the v4 state.
Per-synapse thresholds: rho_v4.thetas (cq = C_REF^(1-gamma) c_post^gamma) from the c_pre / c_post stored in the
fit's extracted record (= the prefire cache). rho0 is the circuit's (checked against the record's).

    python glusynapse_v2/bcl_validation/live_v4.py --fit <json> --tasks <tasks.json> --out <jsonl> --workers N
    python glusynapse_v2/bcl_validation/live_v4.py --fit <json> --from-fit [--pairs-l5 a,b] [--pairs-l23 c] [--no-loc-filter] ...
--from-fit: the task list is every (path, pair, proto, cond) record behind the fit's targets (tasks_from_fit), as
prefire / cooker runs (the fit's own traces; V4_BIN=0 gives the continuous live rule).
"""
import argparse, json, os, sys, tempfile, time, traceback
import numpy as np

ROOT = "/project/rrg-emuller/dhuruva/plastyfire"
V2 = os.path.join(ROOT, "glusynapse_v2")
DEES = "/project/rrg-emuller/dhuruva/DEES_cell_packages/"
V4LIB = os.path.join(V2, "mod_build_v4", "x86_64", "libnrnmech.so")
SIM_ROOT = os.path.join(ROOT, "refitting_results/fitting/n120/seed20262009")
PATHS = {
    "L5": dict(sims=os.path.join(SIM_ROOT, "Sabrina_L5TTPC_L5TTPC_STDP/simulations"),
               circuit=os.path.join(ROOT, "data/dhuruva_delta_circuit_config.json"),
               dirs=[os.path.join(V2, "extracted", f"{d}_delta-prefire-vca")
                     for d in ("markram", "sj03", "sj03r50", "sj07", "ebner")]),
    "L23": dict(sims=os.path.join(SIM_ROOT, "Ebner2019_L23PC_L5TTPC/simulations"),
                circuit=os.path.join(ROOT, "data/dhuruva_delta_l23l5_circuit_config.json"),
                dirs=[os.path.join(V2, "extracted", "ebner_l23l5_delta-prefire-vseg-rs")]),
}
NODE_POP = "S1nonbarrel_neurons"
EDGE_POP = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"
COOKER = os.path.join(ROOT, "fit_results/cooker.json")    # the delta-prefire* runs (run_de_fit2_pool --cooker)
OFFSET, C_MIN, EPSP_WINDOW, SNAP_BEFORE = 1000.0, 4.0, 100.0, 500.0
MECH = "GluSynapseV4"


GEOM = os.path.join(ROOT, "ebner", "pair_geometry_L23PC_L5TTPC.csv")
CONDS = ("control", "mglu_block", "no_block", "nmdar_block")      # the conditions v4_globals / run_task implement


def tasks_from_fit(fit, pairs5=None, pairs23=None, loc_filter=True, phase="prefire", express="cooker"):
    """Every record behind the fit's targets, its drop_targets (validation only) included. L5: load_targets(groups)
    on the fit's --pairs; L23 (--joint fits): paired_l23l5 on the all_protocols pairs of GEOM (as fit_v4 --joint).
    Conditions: the fit's that have a target. A pair enters a (proto, cond) if it has the extracted record; with
    loc_filter, a (proto, cond) whose L23 targets are all located (@proximal / @distal) only takes the pairs that the
    compare selections use (letzkus_* -> letzkus_distal, else sh_distal). pairs5 / pairs23 restrict (pilot)."""
    import pandas as pd
    sys.path.insert(0, V2)
    from targets import load_targets
    fa = fit["args"]; conds = set(fa["conditions"].split(","))
    g = pd.read_csv(GEOM); g["pair"] = g.pregid.astype(str) + "-" + g.postgid.astype(str); g = g[g.all_protocols]
    dist = dict(zip(g.pair, g.letzkus_distal.astype(bool))); sh = dict(zip(g.pair, g.sh_distal.astype(bool)))
    assert fa.get("pairs"), "fit without --pairs"
    spec = [("L5", tuple(fa["groups"].split(",")), list(fa["pairs"].split(",")), pairs5)]
    if fa.get("joint"):
        spec.append(("L23", ("paired_l23l5",), list(g.pair), pairs23))
    tasks = []
    for path, groups, pairs, only in spec:
        if only is not None:
            assert set(only) <= set(pairs), f"{path}: {set(only) - set(pairs)} not in the fit's pairs"
            pairs = [p for p in pairs if p in set(only)]
        need = {}
        for (pid, cond) in sorted(load_targets(groups)):
            if cond not in conds:
                continue
            assert cond in CONDS, f"condition {cond} not implemented live"
            proto, _, where = pid.partition("@")
            if path == "L5" and where:
                continue                                  # compare_v4 skips L5 located targets
            sel = pairs
            if loc_filter and path == "L23" and where:
                m = dist if proto.startswith("letzkus") else sh
                sel = [p for p in pairs if m.get(p, False) == (where == "distal")]
            need.setdefault((proto, cond), set()).update(sel)
        for (proto, cond), ps in sorted(need.items()):
            for p in sorted(ps):
                if find_record(path, p, proto):
                    tasks.append(dict(path=path, pair=p, proto=proto, cond=cond, phase=phase, express=express))
    return tasks


def find_record(path, pair, proto):
    for d in PATHS[path]["dirs"]:
        f = os.path.join(d, f"{pair}__{proto}.npz")
        if os.path.isfile(f):
            return f
    return None


def v4_globals(fit, cond, express):
    """{mod global (no suffix): value} for the fit (MV.DEFAULTS + filters + set + pre, as eval_v4)."""
    sys.path.insert(0, V2); sys.path.insert(0, os.path.join(V2, "rho_redesign"))
    import model_v2 as MV
    import rho_v4
    fa = fit["args"]
    P = {**MV.DEFAULTS, **json.loads(fa["filters"]), **json.loads(fa.get("set", "{}")), **fit["pre"]}
    gd, gp = rho_v4.rates(P); vmode, thV = rho_v4.vamp(P)
    assert int(P["t_drive"]) == 4 and int(P["pre_drive"]) == 1 and int(P["no_drive"]) == 0, "only t_drive 4 / VDCC"
    assert fit.get("bap_gate") is None and rho_v4.opts(P)[1] == 0.0 and rho_v4.opts(P)[2] == 0.0
    g = dict(tau_effca_GB=278.3177658387, gamma_d_GB=gd, gamma_p_GB=gp,
             pre_drive_GB=1, i_scale_GB=P["i_scale"], t_drive_GB=4, tau_E1_GB=P["tau_E1"],
             theta_Te_GB=P["theta_Te"], Te_scale_GB=P["Te_scale"], tau_T_GB=P["tau_T"], theta_Tg_GB=P["theta_Tg"],
             A_mglu_GB=P["A_mglu"], A_NO_GB=P["A_NO"], dpre_min_GB=P["dpre_min"], dpre_max_GB=P["dpre_max"],
             theta_NOi_GB=P["theta_NOi"], theta_Ti_GB=P.get("theta_Ti", P["theta_NOi"]), tau_NO_GB=P["tau_NO"],
             tau_Z_GB=P["tau_Z"], theta_Z_GB=P["theta_Z"], theta_N_GB=P["theta_N"], no_drive_GB=0, z_het_GB=0,
             K_ca_GB=P["K_ca"] * P["K_mult"], ca_hyst_GB=MV.CA_HYST, tau_glu_GB=P["tau_d_NMDA"],
             vamp_mode_GB=vmode, theta_VA_GB=thV, express_x_GB=1 if express == "cooker" else 0)
    ck = json.load(open(COOKER))["params"]
    g.update(gamma_dx_GB=ck["gamma_d_GB_GluSynapse"], gamma_px_GB=ck["gamma_p_GB_GluSynapse"])
    if cond in ("mglu_block", "nmdar_block"):
        g["A_mglu_GB"] = 0.0
    if cond in ("no_block", "nmdar_block"):
        g["A_NO_GB"] = 0.0
    return g, P


def build_config(wd, path, phase, proto, tmpdir):
    """-> (cfg dict, pre spike times, meta). Full phase: simwriter's construction (see module doc)."""
    from libsonata import SpikeReader
    cfg = json.load(open(os.path.join(wd, "prefire_simulation_config.json")))
    ns = json.load(open(os.path.join(wd, "node_sets.json")))
    pre_ids = [int(g) for g in ns["precell"]["node_id"]]
    assert len(pre_ids) == 1, "pairs only"
    pre_gid, post_gid = pre_ids[0], int(ns["postcell"]["node_id"][0])
    cfg["node_sets_file"] = os.path.join(wd, "node_sets.json")
    cfg["network"] = PATHS[path]["circuit"]
    cfg.setdefault("output", {})["output_dir"] = tmpdir
    cfg.pop("reports", None)
    glusyn = dict(cfg.get("conditions", {}).get("mechanisms", {}).pop("GluSynapse", {}))
    spk = SpikeReader(os.path.join(wd, "prefire_prespikes.h5"))[NODE_POP].get_dict()
    pre = np.asarray(spk["timestamps"], float)[np.asarray(spk["node_ids"]) == pre_gid]
    t_pf = float(cfg["run"]["tstop"])
    meta = dict(pre_gid=pre_gid, post_gid=post_gid, glusyn=glusyn, t_prefire=t_pf)
    if phase == "full":
        T = 4000.0 if proto.startswith("10Hz_") else 10000.0
        nb = int(C_MIN * 60000.0 / T); before = nb * T
        pairing = t_pf - OFFSET - 1000.0
        c01 = OFFSET + T * np.arange(nb)
        c02 = OFFSET + before + pairing + T * np.arange(nb)
        tstop = before + pairing + nb * T
        for k, inp in cfg["inputs"].items():
            if inp.get("module") == "pulse":
                inp["delay"] = float(inp["delay"]) + before
            elif inp.get("module") == "synapse_replay":
                inp["duration"] = tstop
        cfg["run"]["tstop"] = tstop
        pre = np.concatenate([c01, pre + before, c02])
        # bluecellulab always adds the config's synapse_replay and it REPLACES the pre_spike_trains connection
        # (cell.add_synapse_replay overwrites cell.connections[sid]), so the replay file must hold the full train
        # (22131696: the unshifted prefire file was replayed, no induction / C02 spikes -> NaN ratios).
        import h5py
        sf = os.path.join(tmpdir, "prespikes_full.h5")
        with h5py.File(sf, "w") as f:                    # simwriter.save_spikes format
            grp = f.require_group(f"spikes/{NODE_POP}")
            grp.create_dataset("timestamps", data=np.sort(pre)); grp["timestamps"].attrs["units"] = "ms"
            grp.create_dataset("node_ids", data=np.full(len(pre), pre_gid), dtype=int)
        n_rep = 0
        for inp in cfg["inputs"].values():
            if inp.get("module") == "synapse_replay":
                inp["spike_file"] = sf; n_rep += 1
        assert n_rep == 1, f"{n_rep} synapse_replay inputs"
        meta.update(before=before, pairing=pairing, c01=c01, c02=c02, tstop=tstop, T=T,
                    t_ind0=before - 500.0, t_snap=c02[0] - SNAP_BEFORE)
        old = os.path.join(wd, "simulation_config.json")
        if os.path.isfile(old):   # Markram / L2/3 workdirs: the written full config must agree
            o = json.load(open(old))
            dl = sorted(v["delay"] for v in o["inputs"].values() if v.get("module") == "pulse")
            dn = sorted(v["delay"] for v in cfg["inputs"].values() if v.get("module") == "pulse")
            meta["matches_written_full_config"] = bool(abs(float(o["run"]["tstop"]) - tstop) < 1e-6
                                                       and np.allclose(dl, dn))
    else:
        meta.update(tstop=t_pf)
    return cfg, pre, meta


def _patch_helper():
    """bluecellulab GluSynapse helper, creating GluSynapseV4 (mech_name stays 'GluSynapse', as the helper)."""
    import neuron
    from bluecellulab.synapse import synapse_types as st
    from bluecellulab.circuit import SynapseProperty
    from bluecellulab.rngsettings import RNGSettings

    def helper(self):
        self.mech_name = "GluSynapse"
        self.hsynapse = getattr(neuron.h, MECH)(self.hoc_args.location, sec=self.hoc_args.section)
        d = self.syn_description; hs = self.hsynapse
        hs.Use_d = d["Use_d_TM"] * d["u_scale_factor"]
        hs.Use_p = d["Use_p_TM"] * d["u_scale_factor"]
        hs.theta_d_GB = d["theta_d"]; hs.theta_p_GB = d["theta_p"]; hs.rho0_GB = d["rho0_GB"]
        hs.volume_CR = d["volume_CR"]
        hs.gmax_d_AMPA = d["gmax_d_AMPA"]; hs.gmax_p_AMPA = d["gmax_p_AMPA"]
        if hs.rho0_GB > getattr(neuron.h, f"rho_star_GB_{MECH}"):
            hs.gmax0_AMPA = hs.gmax_p_AMPA; hs.Use = hs.Use_p
        else:
            hs.gmax0_AMPA = hs.gmax_d_AMPA; hs.Use = hs.Use_d
        hs.gmax_NMDA = hs.gmax0_AMPA * d[SynapseProperty.CONDUCTANCE_RATIO]
        hs.tau_d_AMPA = d[SynapseProperty.DTC]
        hs.Dep = abs(d[SynapseProperty.D_SYN]); hs.Fac = abs(d[SynapseProperty.F_SYN])
        if d[SynapseProperty.NRRP] >= 0:
            hs.Nrrp = d[SynapseProperty.NRRP]
        self.randseed1 = self.post_gid
        self.randseed2 = 100000 + self.syn_id.sid
        self.randseed3 = RNGSettings.get_instance().synapse_seed + 200
        hs.setRNG(self.randseed1, self.randseed2, self.randseed3)
        hs.synapseID = self.syn_id.sid
    st.GluSynapse.use_glusynapse_helper = helper


def _spikes(t, v, thr=-30.0):
    i = np.flatnonzero((v[:-1] < thr) & (v[1:] >= thr))
    return t[i + 1]


def run_task(job):
    """One live run in a fresh process (spawn pool, maxtasksperchild 1). Returns a JSON-able dict."""
    task, fit_path = job
    t0 = time.time()
    tag = "{path}/{pair}/{proto}/{cond}/{phase}/{express}".format(**task)
    out = dict(task, ok=False)
    tmpdir = tempfile.mkdtemp(prefix="bclv4_")
    try:
        sys.path.insert(0, ROOT); sys.path.insert(0, V2); sys.path.insert(0, os.path.join(V2, "rho_redesign"))
        import bluecellulab
        from neuron import h
        bluecellulab.neuron.load_mechanisms(DEES)
        assert h.nrn_load_dll(V4LIB) == 1, f"cannot load {V4LIB}"
        _patch_helper()
        import rho_v4
        from plastyfire.simulator import _map_syn_idx
        from plastyfire.ephysutils import get_epsp_vector
        fit = json.load(open(fit_path)); ck = json.load(open(COOKER))["params"]
        path, pair, proto, cond, phase, express = (task[k] for k in ("path", "pair", "proto", "cond", "phase", "express"))
        wd = os.path.join(PATHS[path]["sims"], pair, proto)
        rec_f = find_record(path, pair, proto)
        rec = np.load(rec_f) if rec_f else None
        cfg, pre, meta = build_config(wd, path, phase, proto, tmpdir)
        cf = os.path.join(tmpdir, "sim.json"); json.dump(cfg, open(cf, "w"))
        g, P = v4_globals(fit, cond, express)
        # pre path on the offline grid (GluSynapseV4 bin mode): the record grid step (extract_v2 --decim 0.25), from
        # the start of the (shifted) prefire to its end (prefire) or the snap (full). V4_BIN=0: continuous v2 path.
        binw = float(os.environ.get("V4_BIN", "0.25"))
        if binw > 0 and rec is not None:
            assert abs(round(float(rec["dt_ms"]), 4) - binw) < 1e-9, f"record grid {float(rec['dt_ms'])} != bin {binw}"   # float32 dt_ms; batch_v2 rounds to 4 dp
        g.update(bin_GB=binw, bin_t0_GB=float(meta.get("before", 0.0)),
                 bin_t1_GB=float(meta["t_snap"] if phase == "full" else meta["tstop"]))
        gamma = rho_v4.opts(P)[0]

        sim = bluecellulab.CircuitSimulation(cf, base_seed=cfg["run"]["random_seed"])
        key = (NODE_POP, meta["post_gid"]); pk = (NODE_POP, meta["pre_gid"])
        sim.instantiate_gids([key], add_synapses=True, add_minis=False, add_pulse_stimuli=True,
                             intersect_pre_gids=[pk], pre_spike_trains={pk: pre})
        cell = sim.cells[key]
        for sec in cell.somatic + cell.axonal:          # simulator_edges fixhp
            sec.uninsert("SK_E2")
        gl = meta["glusyn"]                              # prefire config conditions, as simulator_edges
        for k, hk in (("cao_CR", "cao_CR"), ("init_depleted", "init_depleted"),
                      ("minis_single_vesicle", "minis_single_vesicle")):
            if k in gl:
                setattr(h, f"{hk}_{MECH}", float(gl[k]))
        for k, v in g.items():
            setattr(h, f"{k}_{MECH}", float(v))
        syns = list(cell.synapses.items())
        assert all(s.hsynapse.hname().startswith(MECH) for _, s in syns), "synapses are not GluSynapseV4"
        syn_idx = [sid[1] for sid, _ in syns]
        df = _map_syn_idx(cf, meta["post_gid"], syn_idx, EDGE_POP)
        gids = [int(df.loc[df["local_syn_idx"] == i].index[0]) for i in syn_idx]
        rmap = {}
        if rec is not None:
            rmap = {int(s): (float(a), float(b), float(r0)) for s, a, b, r0 in
                    zip(rec["syn"], rec["c_pre"], rec["c_post"], rec["rho0"])}
        a = fit["a"]; seed = int(cfg.get("run", {}).get("synapse_seed") or 0)
        missing, rho0_mis = [], 0
        for (sid, s), gid in zip(syns, gids):
            hs = s.hsynapse
            if gid in rmap:
                cp, cq, r0 = rmap[gid]
                rho0_mis += int(abs(hs.rho0_GB - r0) > 1e-9)
            else:
                cp = cq = 0.0; missing.append(gid)
            td, tp = rho_v4.thetas(np.array([cp]), np.array([cq]), a, gamma)
            hs.theta_d_GB = float(td[0]) if td[0] > 0 else -1.0
            hs.theta_p_GB = float(tp[0]) if tp[0] > 0 else -1.0
            tdx = ck["a00"] * cp + ck["a01"] * cq; tpx = ck["a10"] * cp + ck["a11"] * cq
            hs.theta_dx_GB = tdx if tdx > 0 else -1.0
            hs.theta_px_GB = tpx if tpx > 0 else -1.0
            hs.setRNG(meta["post_gid"], 100000 + int(sid[1]), 200 + seed)   # simulator_edges seeding
        nmda0 = [s.hsynapse.gmax_NMDA for _, s in syns]

        def state():
            return {k: [float(getattr(s.hsynapse, k)) for _, s in syns]
                    for k in ("rho_GB", "dpre_GB", "rhox_GB", "nev_GB", "Use_GB", "gmax_AMPA", "nbin_GB", "tTsum_GB")}

        def go(t_to, first):
            if first:
                sim.run(t_to, cvode=True)
            else:
                h.cvode_active(1); h.cvode.re_init(); h.continuerun(t_to)

        if phase == "prefire":
            if cond == "nmdar_block":
                for _, s in syns:
                    s.hsynapse.gmax_NMDA = 0.0
            go(meta["tstop"], True)
            out["end"] = state()
        else:
            first = True
            if cond == "nmdar_block":
                go(meta["t_ind0"], True); first = False
                for _, s in syns:
                    s.hsynapse.gmax_NMDA = 0.0
            go(meta["t_snap"], first)
            out["snap"] = state()
            for (_, s), n0 in zip(syns, nmda0):
                hs = s.hsynapse
                rb = 1.0 if hs.rho_GB >= 0.5 else 0.0
                u0 = hs.Use_p if rb else hs.Use_d
                hs.rho_GB = rb; hs.Use = u0; hs.Use_GB = min(1.0, u0 * (1.0 + hs.dpre_GB))
                hs.gmax_AMPA = hs.gmax_p_AMPA if rb else hs.gmax_d_AMPA; hs.gmax0_AMPA = hs.gmax_AMPA
                hs.gmax_NMDA = n0
            setattr(h, f"express_x_GB_{MECH}", 0.0)
            go(meta["tstop"], False)
            out["end"] = state()
            t = np.array(sim.get_time()); v = np.array(sim.get_voltage_trace(key))
            e01 = get_epsp_vector(t, v, meta["c01"], EPSP_WINDOW); e02 = get_epsp_vector(t, v, meta["c02"], EPSP_WINDOW)
            sp = _spikes(t, v)
            ind = (meta["t_ind0"], meta["c02"][0])
            out.update(epsp_c01=float(np.mean(e01)), epsp_c02=float(np.mean(e02)), n_c01=len(e01), n_c02=len(e02),
                       ratio=float(np.mean(e02) / np.mean(e01)),
                       ratio_last=float(np.mean(e02[-10:]) / np.mean(e01[-10:])),
                       n_post_test=int(np.sum((sp < ind[0]) | (sp > ind[1]))),
                       n_post_ind=int(np.sum((sp >= ind[0]) & (sp <= ind[1]))),
                       matches_written_full_config=meta.get("matches_written_full_config"))
            del t, v
        if phase == "prefire":
            t = np.array(sim.get_time()); v = np.array(sim.get_voltage_trace(key))
            out["n_post_ind"] = int(len(_spikes(t, v))); del t, v
        out.update(ok=True, syn=gids, rho0=[float(s.hsynapse.rho0_GB) for _, s in syns],
                   theta_d=[float(s.hsynapse.theta_d_GB) for _, s in syns],
                   theta_p=[float(s.hsynapse.theta_p_GB) for _, s in syns],
                   record=rec_f, missing_in_record=missing, rho0_mismatch=rho0_mis,
                   tstop=float(meta["tstop"]), globals=g)
    except Exception as e:                                           # report, keep the pool going
        out["error"] = f"{type(e).__name__}: {e}"; out["traceback"] = traceback.format_exc()[-3000:]
    finally:
        out["wall_s"] = time.time() - t0
        try:
            import shutil; shutil.rmtree(tmpdir, ignore_errors=True)
        except Exception:
            pass
    print(f"[live_v4] {tag} ok={out['ok']} ratio={out.get('ratio')} wall={out['wall_s']:.0f}s "
          f"{out.get('error', '')}", flush=True)
    return out


def _sub(task, fit, timeout):
    """Run one task in its own python process (a NEURON crash or hang costs one task, not the pool)."""
    import subprocess
    try:
        p = subprocess.run([sys.executable, os.path.abspath(__file__), "--fit", fit, "--one", json.dumps(task)],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=timeout)
        res = [l[len("RESULT "):] for l in p.stdout.splitlines() if l.startswith("RESULT ")]
        if res:
            return json.loads(res[-1])
        return dict(task, ok=False, error=f"exit {p.returncode}, no result", log_tail=p.stdout[-3000:])
    except subprocess.TimeoutExpired:
        return dict(task, ok=False, error=f"timeout {timeout} s")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True)
    ap.add_argument("--tasks", help="json list of {path, pair, proto, cond, phase, express}")
    ap.add_argument("--out", help="results jsonl (appended; done tasks are skipped)")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--timeout", type=float, default=7200.0, help="per task, s")
    ap.add_argument("--one", default=None, help="internal: run one task (json) and print 'RESULT <json>'")
    ap.add_argument("--from-fit", action="store_true", help="tasks = tasks_from_fit(--fit) (prefire / cooker)")
    ap.add_argument("--pairs-l5", default=None, help="--from-fit: only these L5 pairs (comma list)")
    ap.add_argument("--pairs-l23", default=None, help="--from-fit: only these L23 pairs (comma list)")
    ap.add_argument("--no-loc-filter", action="store_true", help="--from-fit: every (proto, cond) on every pair")
    a = ap.parse_args()
    if a.one:
        print("RESULT " + json.dumps(run_task((json.loads(a.one), os.path.abspath(a.fit)))), flush=True)
        return
    if a.from_fit:
        sp = lambda s: s.split(",") if s else None
        tasks = tasks_from_fit(json.load(open(a.fit)), sp(a.pairs_l5), sp(a.pairs_l23), not a.no_loc_filter)
        from collections import Counter
        print("[live_v4] from fit:", dict(Counter((t["path"], t["cond"]) for t in tasks)), flush=True)
    else:
        tasks = json.load(open(a.tasks))
    keyf = lambda d: tuple(d[k] for k in ("path", "pair", "proto", "cond", "phase", "express"))
    done = set()
    if os.path.isfile(a.out):
        for line in open(a.out):
            d = json.loads(line)
            if d.get("ok"):
                done.add(keyf(d))
    todo = [t for t in tasks if keyf(t) not in done]
    def _tstop(t):                                    # prefire sim length (s) of the workdir, for the ordering
        try:
            f = os.path.join(PATHS[t["path"]]["sims"], t["pair"], t["proto"], "prefire_simulation_config.json")
            return float(json.load(open(f))["run"]["tstop"])
        except Exception:
            return 0.0
    todo.sort(key=lambda t: (t["phase"] != "full", -_tstop(t)))     # long full runs first, then longest prefire
    print(f"[live_v4] {len(todo)} tasks ({len(tasks) - len(todo)} done), {a.workers} workers", flush=True)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    from concurrent.futures import ThreadPoolExecutor, as_completed
    with ThreadPoolExecutor(a.workers) as ex, open(a.out, "a") as fo:
        futs = [ex.submit(_sub, t, os.path.abspath(a.fit), a.timeout) for t in todo]
        for f in as_completed(futs):
            r = f.result()
            print(f"[live_v4] done {r.get('path')}/{r.get('pair')}/{r.get('proto')}/{r.get('cond')}/{r.get('phase')}/"
                  f"{r.get('express')} ok={r.get('ok')} ratio={r.get('ratio')} wall={r.get('wall_s', 0):.0f}s "
                  f"{r.get('error', '')}", flush=True)
            fo.write(json.dumps(r) + "\n"); fo.flush()


if __name__ == "__main__":
    main()
