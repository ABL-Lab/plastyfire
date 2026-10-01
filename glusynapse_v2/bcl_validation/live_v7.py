"""Live bluecellulab runs of a v7 fit (rho_redesign gpu_v7_rho.py through fit_v6: objective gpu_v6_rho, v5_mode 2,
set veto_T / ecb_ref) with the GluSynapseV7 mechanism (mod/GluSynapseV7.mod, built by compile_v7.sh). live_v5.py is
reused unchanged (imported and patched here, not edited); the differences:
  * mechanism GluSynapseV7 from mod_build_v7 (live_v5.MECH / V5LIB);
  * globals = live_v5.v5_globals (its v5 asserts hold for v7) + veto_T_GB = the fit's veto_T (0 = V5 rule);
  * per-synapse uE_GB (theta_eCB,i = theta_eCB uE_i), set before the run:
      ecb_ref 0 : uE 1 (absolute theta_eCB, V7v);
      ecb_ref 2 : uE = 1e5 x vdcc_q_post from the A4 cexp tables (<cexp_dir>/L5L5.csv, L23L5.csv; keyed by pre_gid,
                  post_gid, syn_id = the edge id) and NaN -> fit["v6"]["seed_medians"]["uE"], the pooled median over all
                  pathways that fit_v6 used (same rule, same number). A synapse missing from the table is an error.
      ecb_ref 1 is not implemented (not uniform across pathways, LTD_DIAG V7).
  * the result adds v7 = {nveto, nvpend, uE} per synapse (necb_GB is in the saved state as in v5); nvpend must be 0.
  * joint fits with three pathways (W2): the L5 dirs add l5extra_<E>-prefire-vca (sims in the Sabrina folder), and a
    path L23L23 (fit --extra l23l23:GROUPS:...) with the Zilberter + L23L23extra (Egger, Banerjee) records, sims
    BCL_L23L23_SIMS (':'-list, the task's workdir is taken from the first that has it), circuit BCL_L23L23_CIRCUIT
    (default data/dhuruva_split1_l23l23_circuit_config.json: rho0-split edges, as l23l5). Tasks: every (proto, cond)
    of the extra groups + VAL_GROUPS (validation only, BCL_NO_VAL=1 drops them) on every pair with the record;
    BCL_PAIRS_L23L23 (comma list, "none") restricts (shards). Located (@) targets are skipped as for L5.
ECB_REF (env, optional) must equal the fit's ecb_ref (guard against a wrong json in the run scripts).

    python glusynapse_v2/bcl_validation/live_v7.py <live_v5 arguments>
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import live_v5 as L                                    # noqa: E402

V7LIB = os.path.join(L.V2, "mod_build_v7", "x86_64", "libnrnmech.so")
CEXP_DIR = "/scratch/dhuruva/split1/cexp"
CEXP_FILE = {"L5": "L5L5.csv", "L23": "L23L5.csv", "L23L23": "L23L23.csv"}
VAL_GROUPS = {"L23L23": ("paired_l23l23_val",)}       # Banerjee 2014: scored, never fitted
_v5_globals, _run_task, _tasks_from_fit = L.v5_globals, L.run_task, L.tasks_from_fit
_SR = "/scratch/dhuruva/split1/refitting_results/fitting/n120/seed20262009"
_X = os.path.join(L.V2, "extracted")
_d = os.path.join(_X, f"l5extra_{L._E}-prefire-vca")
if os.path.isdir(_d) and _d not in L.PATHS["L5"]["dirs"]:
    L.PATHS["L5"]["dirs"].append(_d)
_s23 = (os.environ.get("BCL_L23L23_SIMS") or
        f"{_SR}/Zilberter2009_L23PC_L23PC/simulations:{_SR}/L23L23extra/simulations").split(":")
L.PATHS["L23L23"] = dict(
    sims=_s23[0], sims_list=_s23,
    circuit=os.environ.get("BCL_L23L23_CIRCUIT") or os.path.join(L.ROOT, "data/dhuruva_split1_l23l23_circuit_config.json"),
    dirs=[os.path.join(_X, f"{d}_{L._E}-prefire-vseg-rs") for d in ("zilberter_l23l23", "l23l23extra")])


def extra_groups(fit, name="l23l23"):
    """GROUPS of the fit's --extra NAME:GROUPS:DIRS:BASIS entry (() if none) and its basis dir."""
    for e in fit["args"].get("extra") or []:
        f = e.split(":")
        if f[0] == name:
            return tuple(f[1].split(",")), f[3]
    return (), ""


def tasks_from_fit(fit, pairs5=None, pairs23=None, loc_filter=True, phase="prefire", express="cooker"):
    """live_v5.tasks_from_fit (L5, L23) + the L23L23 tasks of the fit's l23l23 extra model and VAL_GROUPS."""
    tasks = _tasks_from_fit(fit, pairs5, pairs23, loc_filter, phase, express)
    groups, _ = extra_groups(fit)
    if not groups:
        return tasks
    from targets import load_targets
    conds = set(fit["args"]["conditions"].split(","))
    P = L.PATHS["L23L23"]
    pairs = sorted({f.split("__")[0] for d in P["dirs"] if os.path.isdir(d) for f in os.listdir(d) if f.endswith(".npz")})
    only = os.environ.get("BCL_PAIRS_L23L23")
    if only:
        pairs = [p for p in pairs if p in set(only.split(","))]
    gsets = [groups] + ([] if os.environ.get("BCL_NO_VAL") else [VAL_GROUPS["L23L23"]])
    need = {}
    for gs in gsets:
        for (pid, cond) in sorted(load_targets(gs)):
            proto, _, where = pid.partition("@")
            if cond not in conds or where:
                continue
            assert cond in L.CONDS, f"condition {cond} not implemented live"
            need.setdefault((proto, cond), set()).update(pairs)
    for (proto, cond), ps in sorted(need.items()):
        for p in sorted(ps):
            if L.find_record("L23L23", p, proto):
                tasks.append(dict(path="L23L23", pair=p, proto=proto, cond=cond, phase=phase, express=express))
    return tasks


def fit_P(fit):
    """The fit's parameter dict as live_v5.v5_globals builds it (MV.DEFAULTS + v5 defaults + filters + set + pre)."""
    sys.path.insert(0, L.V2)
    import model_v2 as MV
    fa = fit["args"]
    return {**MV.DEFAULTS, "v5_mode": 0, "A_eCB": 0.0, "theta_eCB": 0.0, **json.loads(fa["filters"]),
            **json.loads(fa.get("set", "{}")), **fit["pre"]}


def ecb_ref(fit):
    P = fit_P(fit); e = int(P.get("ecb_ref", 0))
    assert e in (0, 2), f"ecb_ref {e}: live v7 implements 0 (absolute) and 2 (cexp vdcc_q_post) only"
    v6 = fit.get("v6") or {}
    assert not v6.get("scale_V") and not v6.get("theta_pre"), "scale_V / theta_pre not implemented live"
    if e == 2:
        assert v6.get("scale_E") == {"vdcc_q_post": 1e5}, f"ecb_ref 2 needs scale_E vdcc_q_post x 1e5: {v6.get('scale_E')}"
    else:
        assert not v6.get("scale_E"), "ecb_ref 0 with a scale_E fit"
    env = os.environ.get("ECB_REF")
    assert env is None or int(env) == e, f"ECB_REF {env} but the fit has ecb_ref {e}"
    return e


def ue_map(fit, path):
    """{(pre_gid, post_gid, syn_id): uE} for ecb_ref 2 (fit_v6 lincomb + NaN rule), None for ecb_ref 0."""
    if ecb_ref(fit) != 2:
        return None
    import numpy as np, pandas as pd
    v6 = fit["v6"]
    df = pd.read_csv(os.path.join(v6.get("cexp_dir") or os.environ.get("CEXP_DIR") or CEXP_DIR, CEXP_FILE[path]))
    u = 1e5 * df.vdcc_q_post.to_numpy(np.float64)
    u = np.where(np.isfinite(u), u, float(v6["seed_medians"]["uE"]))
    assert np.all(u >= 0), "negative uE"
    return {(int(a), int(b), int(c)): float(x) for a, b, c, x in zip(df.pre_gid, df.post_gid, df.syn_id, u)}


def v7_globals(fit, cond, express):
    g, P = _v5_globals(dict(fit, objective="gpu_v5_rho"), cond, express)   # fit_v6 writes objective gpu_v6_rho
    assert int(P["v5_mode"]) == 2, "v7 fits are v5_mode 2 (v5c + options)"
    g["veto_T_GB"] = float(P.get("veto_T", 0.0))
    return g, P


def run_task(job, hook=None):
    task, fit_path = job
    fit = json.load(open(fit_path))
    P = L.PATHS[task["path"]]
    if "sims_list" in P:                          # L23L23: the workdir is in one of the sims folders (one task per process)
        P["sims"] = next((s for s in P["sims_list"] if os.path.isdir(os.path.join(s, task["pair"], task["proto"]))),
                         P["sims_list"][0])
    um = ue_map(fit, task["path"])

    def hook7(task, sim, cell, syns, gids, meta, h):
        assert getattr(h, f"bin_GB_{L.MECH}") == 0.0, "the veto reads qca_GB: bin_GB must be 0"
        miss = []
        for (_, s), gid in zip(syns, gids):
            if um is None:
                s.hsynapse.uE_GB = 1.0
            else:
                k = (*map(int, task["pair"].split("-")), int(gid))     # cexp key = fit_v6 (pair, syn)
                if k not in um:
                    miss.append(gid); continue
                s.hsynapse.uE_GB = um[k]
        assert not miss, f"{len(miss)} synapses not in the cexp table (first {miss[:3]})"
        fin0 = hook(task, sim, cell, syns, gids, meta, h) if hook is not None else None

        def fin(out):
            out["v7"] = {k: [float(getattr(s.hsynapse, f"{k}_GB")) for _, s in syns] for k in ("nveto", "nvpend", "uE")}
            out["v7"]["veto_T"] = float(getattr(h, f"veto_T_GB_{L.MECH}"))
            if fin0 is not None:
                fin0(out)
        return fin

    return _run_task(job, hook7)


L.MECH = "GluSynapseV7"
L.V5LIB = V7LIB
L.v5_globals = v7_globals
L.run_task = run_task
L.tasks_from_fit = tasks_from_fit
L.__file__ = os.path.abspath(__file__)        # live_v5._sub re-launches this file per task (--one)

if __name__ == "__main__":
    if "--one" not in sys.argv:                  # pool process: check the fit / ECB_REF once before the tasks
        _f = sys.argv[sys.argv.index("--fit") + 1]
        print(f"[live_v7] {_f}: ecb_ref {ecb_ref(json.load(open(_f)))}, mechanism {L.MECH} ({V7LIB})", flush=True)
    L.main()
