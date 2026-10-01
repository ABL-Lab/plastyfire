"""Per-synapse experimental-style calcium measures on the delta-split1 emodel (no .mod edits).

Conditions (each in its own NEURON subprocess, same setup as simulator_edges._find_cpre_cpost: defit2 globals, fixhp,
per-synapse SONATA params, rho0=1 / Use_p / gmax_p, thresholds off):
  pre     1 pre spike at 1 s, 1 mM Mg                (Chindemi c_pre)
  apv     same, gmax_NMDA = 0 on every synapse       (APV: no NMDA)
  mg0     same, mg = 0 on every synapse              (Mg block removed, full NMDA)
  post    single bAP (soma TStim, amplitude/width from spike_threshold_finder), 1 mM Mg   (Chindemi c_post)
Recorded per synapse: effcai_GB peak (c_pre/c_post units), cai_CR peak (free spine Ca, uM), integral of ica_VDCC and
ica_NMDA from 990 ms to the end (inward charge, pC = nA ms, reported positive).
mg / gmax_NMDA are RANGE variables of GluSynapse (set per synapse after the SONATA params are applied).

  python measure_cexp.py --path L5 --start 0 --stop 3      # pairs [start, stop) of the sorted pair list -> parts
  python measure_cexp.py --path L5 --merge                 # parts -> /scratch/dhuruva/split1/cexp/<path>.csv
"""
import argparse, glob, json, logging, multiprocessing, os, sys, tempfile, time
import numpy as np, pandas as pd

ROOT = "/project/rrg-emuller/dhuruva/plastyfire"
sys.path.insert(0, ROOT)
S1 = "/scratch/dhuruva/split1"
SR = f"{S1}/refitting_results/fitting/n120/seed20262009"
OUT = os.environ.get("CEXP_OUT", f"{S1}/cexp")
PATHS = {
    "L5":     dict(name="L5L5",   sims=f"{SR}/Sabrina_L5TTPC_L5TTPC_STDP/simulations", circuit=f"{ROOT}/data/dhuruva_delta-split1_circuit_config.json",
                   cache=f"{S1}/cache/sabrina_n120_split1.pkl"),
    "L23L5":  dict(name="L23L5",  sims=f"{SR}/Ebner2019_L23PC_L5TTPC/simulations", circuit=f"{ROOT}/data/dhuruva_split1_l23l5_circuit_config.json",
                   cache=f"{S1}/cache/ebner_l23l5_split1_rs.pkl"),
    "L23L23": dict(name="L23L23", sims=f"{SR}/Zilberter2009_L23PC_L23PC/simulations", circuit=f"{ROOT}/data/dhuruva_split1_l23l23_circuit_config.json",
                   cache=f"{S1}/cache/zilberter_l23l23_split1_rs.pkl"),
}
NODE_POP = "S1nonbarrel_neurons"
EDGE_POP = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"
DEES = "/project/rrg-emuller/dhuruva/DEES_cell_packages/"
FIT = {"tau_effca_GB_GluSynapse": 278.3177658387, "gamma_d_GB_GluSynapse": 77.7558, "gamma_p_GB_GluSynapse": 299.9121}  # PARAM_PRESETS["defit2"] globals
if os.environ.get("GLUSYN_GLOBALS"):  # optional: JSON file (e.g. spine/delta_ljp25.json, key "globals" or flat dict) of extra *_GluSynapse HOC globals
    _j = json.load(open(os.environ["GLUSYN_GLOBALS"])); _g = {k: float(v) for k, v in _j.get("globals", _j).items()}
    assert _g and all(k.endswith("_GluSynapse") for k in _g), _g
    FIT.update(_g); print("GLUSYN_GLOBALS", os.environ["GLUSYN_GLOBALS"], _g, flush=True)
T_STIM, T_END, T_INT, DT_FAST = 1000.0, 1500.0, 990.0, 0.1
DIAG = os.environ.get("CEXP_DIAG") == "1"  # opt-in (A17): extra columns <cond>_vrest/_vpk (synapse segment v, mV), _cai0 (cai_CR at 989 ms, uM), _eff0 (effcai at 999 ms)
CONDS = ("pre", "apv", "mg0", "post")
# opt-in CEXP_SABATINI=1 (SABATINI_VALID.md): Chindemi 2022 Supp Fig 2 protocol. Conditions "syn" (SAB_N pre spikes at 0.2 Hz,
# stochastic release with the circuit's own Use/Nrrp/rho0, i.e. no full-RRP override) and "post" (one bAP). Per synapse and trial:
# peak free cai_CR in [t_k, t_k+SAB_WIN] minus cai_CR at t_k-1 ms (uM). Plus section name, diameter, branch order (primary = 1).
SAB = os.environ.get("CEXP_SABATINI") == "1"
SAB_N, SAB_ISI, SAB_WIN = int(os.environ.get("CEXP_SAB_NTRIALS", "10")), 5000.0, 200.0
SAB_THR = float(os.environ.get("CEXP_SAB_THR", "0.05"))  # uM; a trial with dCa above this counts as a release success


def pair_workdirs(sims):
    """{(pre, post): first workdir}"""
    out = {}
    for wd in sorted(glob.glob(f"{sims}/*-*/*")):
        if os.path.isfile(f"{wd}/prefire_simulation_config.json") and os.path.isfile(f"{wd}/prefire_prespikes.h5"):
            a, b = os.path.basename(os.path.dirname(wd)).split("-")
            out.setdefault((int(a), int(b)), wd)
    return out


def _child(conn, workdir, circuit, cond, stim):
    import bluecellulab
    from bluepysnap import Simulation
    from plastyfire.simulator import _map_syn_idx, _set_global_params, _set_local_params
    from plastyfire.simulator_edges import _read_syn_extra_params, _precell, spike_threshold_finder, C_POST_PULSE_WIDTHS_MS, \
        C_POST_MIN_AMP_NA, C_POST_MAX_AMP_NA, C_POST_AMP_LEVELS
    log = logging.getLogger("cexp")
    try:
        bluecellulab.neuron.load_mechanisms(DEES)
        h = bluecellulab.neuron.h
        sim_config = f"{workdir}/prefire_simulation_config.json"
        cfg = json.load(open(sim_config))
        ns = os.path.abspath(f"{workdir}/node_sets.json")
        if os.path.exists(ns):
            cfg["node_sets_file"] = ns
        cfg["network"] = os.path.abspath(circuit)
        cfg.pop("reports", None)
        cfg.get("conditions", {}).get("mechanisms", {}).pop("GluSynapse", None)
        tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".json", dir=tempfile.gettempdir(), delete=False)
        json.dump(cfg, tmp); tmp.close()
        snap = Simulation(tmp.name)
        pre_gid = _precell(snap.node_sets.content); post_gid = int(snap.node_sets.content["postcell"]["node_id"][0])
        pres = [pre_gid] if isinstance(pre_gid, int) else list(pre_gid)
        extra = _read_syn_extra_params(tmp.name, pre_gid, post_gid, EDGE_POP, log)
        res = dict(pre_gid=pre_gid, post_gid=post_gid)
        if cond == "post" and stim is None:                       # stimulus search, as _find_cpre_cpost
            for w in C_POST_PULSE_WIDTHS_MS:
                r = spike_threshold_finder(tmp.name, post_gid, 1, 0.1, w, 1000., C_POST_MIN_AMP_NA, C_POST_MAX_AMP_NA,
                                           C_POST_AMP_LEVELS, NODE_POP, True)
                if r is not None:
                    stim = (float(r["amp"]), float(r["width"])); break
            if stim is None:
                raise RuntimeError("no single-AP stimulus")
        kw = dict(pre_spike_trains={(NODE_POP, g): [T_STIM] for g in pres}) if cond != "post" else {}
        t_end = T_END
        if cond == "syn":
            tk = [T_STIM + k * SAB_ISI for k in range(SAB_N)]; t_end = tk[-1] + 500.0
            kw = dict(pre_spike_trains={(NODE_POP, g): list(tk) for g in pres})
        sim = bluecellulab.CircuitSimulation(tmp.name, base_seed=cfg["run"]["random_seed"])
        sim.instantiate_gids([(NODE_POP, post_gid)], add_synapses=True, add_minis=False, add_pulse_stimuli=False,
                             intersect_pre_gids=[(NODE_POP, g) for g in pres], **kw)
        cell = sim.cells[(NODE_POP, post_gid)]
        for sec in cell.somatic + cell.axonal:
            sec.uninsert("SK_E2")
        _set_global_params(FIT)
        if cond == "post":
            tstim = h.TStim(0.5, sec=cell.soma)
            tstim.train(T_STIM, stim[1], stim[0], 0.1, stim[1]); cell.persistent.append(tstim)
        syn_idx = [s[1] for s in cell.synapses]
        df = _map_syn_idx(sim_config, post_gid, syn_idx, EDGE_POP)
        rec = {}; geo = {}
        h.distance(0, 0.5, sec=cell.soma)
        for sid, syn in cell.synapses.items():
            hs = syn.hsynapse
            gid = int(df.loc[df["local_syn_idx"] == sid[1]].index[0])
            _set_local_params(syn, FIT, extra[gid])
            if cond not in ("post", "syn"):
                hs.rho0_GB = 1.0; hs.Use_p = 1.0; hs.Use = hs.Use_p; hs.gmax0_AMPA = hs.gmax_p_AMPA
            hs.theta_d_GB = -1.0; hs.theta_p_GB = -1.0
            if cond == "apv":
                hs.gmax_NMDA = 0.0
            elif cond == "mg0":
                hs.mg = 0.0
            x = hs.get_loc(); sec = h.cas(); h.pop_section(); seg = sec(x)
            dist = float(h.distance(1, seg))
            v = {}
            if cond == "syn":
                v["cai"] = h.Vector(); v["cai"].record(hs._ref_cai_CR, DT_FAST)
                geo[gid] = dict(Use=float(hs.Use), rho0=float(hs.rho0_GB))
            else:
                for k, ref, dt in (("eff", hs._ref_effcai_GB, 1.0), ("cai", hs._ref_cai_CR, DT_FAST),
                                   ("vdcc", hs._ref_ica_VDCC, DT_FAST), ("nmda", hs._ref_ica_NMDA, DT_FAST)):
                    v[k] = h.Vector(); v[k].record(ref, dt)
            if SAB:
                bo, s = 0, sec
                while True:
                    sr = h.SectionRef(sec=s)
                    if "soma" in s.name():
                        break
                    bo += 1
                    if not sr.has_parent():
                        break
                    s = sr.parent
                geo.setdefault(gid, {}).update(secname=sec.name(), diam_um=float(seg.diam), branch_order=bo)
            if DIAG:
                v["v"] = h.Vector(); v["v"].record(seg._ref_v, DT_FAST)
            rec[gid] = (v, extra[gid]["loc"], dist, float(hs.gmax_NMDA), float(hs.mg))
        sim.run(t_end, cvode=True)
        i1 = int(round(T_INT / 1.0)); i2 = int(round(T_INT / DT_FAST))
        out = {}
        for gid, (v, loc, dist, g, mg) in rec.items():
            if cond == "syn":
                cai = np.array(v["cai"]) * 1e3; d = []          # uM
                for t0 in tk:
                    ib = int(round((t0 - 1.0) / DT_FAST)); ie = int(round((t0 + SAB_WIN) / DT_FAST))
                    d.append(float(cai[ib:ie].max() - cai[ib]))
                out[gid] = dict(loc=loc, dist_um=dist, dca=d, **geo[gid]); continue
            eff = np.array(v["eff"]); cai = np.array(v["cai"]); vd = np.array(v["vdcc"]); nm = np.array(v["nmda"])
            out[gid] = dict(loc=loc, dist_um=dist, eff=float(eff[1000:].max() if cond != "post" else eff.max()),
                            cai=float(cai[i2:].max() * 1e3),                       # mM -> uM
                            vdcc_q=float(-np.trapz(vd[i2:], dx=DT_FAST)), nmda_q=float(-np.trapz(nm[i2:], dx=DT_FAST)),
                            gmax_NMDA=g, mg=mg)
            if DIAG:
                vv = np.array(v["v"])
                out[gid].update(vrest=float(vv[i2 - 10]), vpk=float(vv[i2:].max()), cai0=float(cai[i2 - 10] * 1e3),
                                eff0=float(eff[int(T_STIM) - 1]))
            if SAB:
                c0 = float(cai[int(round((T_STIM - 1.0) / DT_FAST))] * 1e3)
                out[gid].update(geo[gid], dca=float(cai[i2:].max() * 1e3 - c0), cai0_sab=c0)
        res.update(out=out, stim=stim)
        if cond == "post":
            t = np.array(sim.get_time()); vs = np.array(sim.get_voltage_trace((NODE_POP, post_gid)))
            res["post_vmax"] = float(vs.max())
        conn.send(res)
    except Exception as e:
        import traceback; traceback.print_exc()
        conn.send(None)
    finally:
        conn.close()


def run_cond(workdir, circuit, cond, stim=None):
    pc, cc = multiprocessing.Pipe(duplex=False)
    p = multiprocessing.Process(target=_child, args=(cc, workdir, circuit, cond, stim)); p.start(); cc.close()
    try:
        r = pc.recv()
    except EOFError:
        r = None
    p.join()
    if r is None:
        raise RuntimeError(f"{cond} failed for {workdir}")
    return r


def measure_pair(path, pair, wd):
    """A failing condition gives NaN for its columns and is named in the `error` column; the pair is still written."""
    P = PATHS[path]; R = {}; err = []
    for c in ("pre", "apv", "mg0", "post"):
        try:
            R[c] = run_cond(wd, P["circuit"], c)
        except Exception as e:
            err.append(c); print(f"  {c} failed for {pair}: {e}", flush=True)
    if "pre" not in R:
        raise RuntimeError(f"pre failed for {pair}")
    nan = dict(eff=np.nan, cai=np.nan, vdcc_q=np.nan, nmda_q=np.nan)
    g = lambda c, gid: R[c]["out"][gid] if c in R else nan
    rows = []
    for gid, a in R["pre"]["out"].items():
        r = dict(pathway=P["name"], pre_gid=pair[0], post_gid=pair[1], syn_id=gid, loc=a["loc"], dist_um=a["dist_um"])
        for k, c in (("cpre", "pre"), ("cpre_apv", "apv"), ("cpre_mg0", "mg0"), ("cpost", "post")):
            r[k] = g(c, gid)["eff"]
        for k, c in (("cpre", "pre"), ("cpre_apv", "apv"), ("cpre_mg0", "mg0"), ("cpost", "post")):
            r[k + "_cai"] = g(c, gid)["cai"]
        r["vdcc_q"] = a["vdcc_q"]; r["nmda_q"] = a["nmda_q"]
        for c in ("apv", "mg0", "post"):
            r["vdcc_q_" + c] = g(c, gid)["vdcc_q"]; r["nmda_q_" + c] = g(c, gid)["nmda_q"]
        if DIAG:
            for c in ("pre", "apv", "mg0", "post"):
                for k in ("vrest", "vpk", "cai0", "eff0"):
                    r[f"{c}_{k}"] = g(c, gid).get(k, np.nan)
        r["gmax_NMDA_nS"] = a["gmax_NMDA"]
        r["stim_nA"], r["stim_ms"] = R["post"]["stim"] if "post" in R else (np.nan, np.nan)
        r["error"] = ",".join(err)
        rows.append(r)
    return pd.DataFrame(rows)


def measure_pair_sab(path, pair, wd):
    """CEXP_SABATINI=1: per synapse dCa (uM, baseline subtracted) for SAB_N stochastic synaptic events and one bAP."""
    P = PATHS[path]; R = run_cond(wd, P["circuit"], "syn"); err = []
    try:
        B = run_cond(wd, P["circuit"], "post")
    except Exception as e:
        B = None; err.append("post"); print(f"  post failed for {pair}: {e}", flush=True)
    rows = []
    for gid, a in R["out"].items():
        d = np.array(a["dca"]); s = d > SAB_THR
        r = dict(pathway=P["name"], pre_gid=pair[0], post_gid=pair[1], syn_id=gid, loc=a["loc"], dist_um=a["dist_um"],
                 secname=a["secname"], diam_um=a["diam_um"], branch_order=a["branch_order"], Use=a["Use"], rho0=a["rho0"],
                 syn_dca_mean=float(d.mean()), syn_dca_succ=float(d[s].mean()) if s.any() else np.nan,
                 syn_psucc=float(s.mean()), syn_dca_trials=";".join(f"{x:.4f}" for x in d),
                 bap_dca=B["out"][gid]["dca"] if B else np.nan, bap_cai0=B["out"][gid]["cai0_sab"] if B else np.nan,
                 stim_nA=B["stim"][0] if B else np.nan, stim_ms=B["stim"][1] if B else np.nan, error=",".join(err))
        rows.append(r)
    return pd.DataFrame(rows)


def sab_summary(d):
    """Chindemi 2022 Supp Fig 2 pools: L5 basal, diam < 2 um, 2 <= branch order <= 4; synaptic path < 150 um, bAP < 60 um."""
    ref = dict(syn=((0.67, 0.44), (0.7, 0.4)), bap=((1.4, 0.6), (1.7, 0.6)))   # (Chindemi model, Sabatini 2002 data), uM mean, SD
    bas = d[d["loc"] == "basal"]
    acc = bas[(bas.diam_um < 2.0) & bas.branch_order.between(2, 4)]
    print(f"synapses {len(d)}, basal {len(bas)}, Sabatini-accepted basal {len(acc)}; trials/synapse {SAB_N}, success thr {SAB_THR} uM")
    for lab, col, key in (("synaptic all-trial mean", "syn_dca_mean", "syn"), ("synaptic successes only", "syn_dca_succ", "syn"),
                          ("bAP", "bap_dca", "bap")):
        for pool, p in (("accepted", acc), ("all basal", bas)):
            for lim in ((60.0, 150.0) if key == "syn" else (60.0,)):
                x = p.loc[p.dist_um < lim, col].dropna()
                (cm, cs), (sm, ss) = ref[key]
                ok = "PASS" if len(x) and abs(x.mean() - sm) <= ss else "FAIL"
                print(f"{lab:25s} {pool:9s} <{lim:.0f} um: {x.mean():.2f} +- {x.std():.2f} uM (n={len(x)}) | Chindemi {cm} +- {cs}"
                      f" | Sabatini {sm} +- {ss} | {ok} (within 1 SD of Sabatini)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--path", required=True, choices=list(PATHS))
    ap.add_argument("--start", type=int, default=0); ap.add_argument("--stop", type=int, default=10 ** 6)
    ap.add_argument("--merge", action="store_true")
    a = ap.parse_args()
    logging.basicConfig(level=logging.WARNING)
    parts = f"{OUT}/parts/{a.path}"; os.makedirs(parts, exist_ok=True)
    if a.merge:
        d = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{parts}/*.csv"))], ignore_index=True)
        d.to_csv(f"{OUT}/{PATHS[a.path]['name']}.csv", index=False)
        print(a.path, "pairs", d.groupby(["pre_gid", "post_gid"]).ngroups, "synapses", len(d))
        if SAB:
            sab_summary(d)
        return
    pw = pair_workdirs(PATHS[a.path]["sims"]); keys = sorted(pw)
    print(a.path, "pairs total", len(keys), "this task", keys[a.start:a.stop], flush=True)
    failed = []
    for k in keys[a.start:a.stop]:
        t0 = time.time()
        f = f"{parts}/{k[0]}-{k[1]}.csv"
        if os.path.exists(f):
            continue
        try:
            d = (measure_pair_sab if SAB else measure_pair)(a.path, k, pw[k]); d.to_csv(f, index=False)
        except Exception as e:
            print(f"FAILED PAIR {k}: {e}", flush=True); failed.append(k); continue
        print(f"pair {k}: {len(d)} syn, {time.time() - t0:.0f} s, error='{d.error.iloc[0]}'", flush=True)
    print("failed pairs:", failed, flush=True)


if __name__ == "__main__":
    main()
