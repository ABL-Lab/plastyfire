"""Live bluecellulab runs of a v8 fit (rho_redesign gpu_v8_rho.py through fit_v6) with GluSynapseV8 (mod/GluSynapseV8.mod,
built by compile_v8.sh). live_v7.py (and through it live_v5.py) is reused unchanged (imported and patched here); the
differences:
  * mechanism GluSynapseV8 from mod_build_v8;
  * globals = live_v7.v7_globals + the gate: gate_src_GB = the fit's gate_src (0 = V7 rule; 2 = G1 shaft gate, the only
    gpu_v8 source implemented live, with gate_win > 0), gate_win_GB = gate_win (ms), theta_sh_GB = fit["v8"]["theta_G"]
    (uM -> mM);
  * gate_src 2: t_rest_GB = min(pre spikes, post spikes of the task's record) - 1 ms, so each synapse's cai_rest_GB is
    its own segment cai just before the stimulus (gpu_v8_rho._shaft_rest: last grid sample before that time);
  * the result adds v8 = {gate_src, gate_win, theta_sh_uM, t_rest, cai_rest (mM), nsh (licence openings),
    tsh (licence open time, ms, the open interval at tstop included), gsh_end} per synapse.
The fit's t_exact (eCB timing) needs no switch: live always reads W at the exact arrival (= offline t_exact 1).

    python glusynapse_v2/bcl_validation/live_v8.py <live_v5 arguments>
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import live_v7 as L7                                   # noqa: E402  (patches live_v5 for V7)
L = L7.L

V8LIB = os.path.join(L.V2, "mod_build_v8", "x86_64", "libnrnmech.so")


def gate(fit):
    """(gate_src, gate_win ms, theta_G uM) of a fit (gate_src 0 for v7 fits without a v8 block)."""
    P = L7.fit_P(fit)
    src = int(P.get("gate_src", 0))
    if src == 0:
        return 0, 0.0, 0.0
    assert src == 2, f"gate_src {src}: live v8 implements 0 (Vg, V7) and 2 (shaft G1) only"
    win = float(P.get("gate_win", 100.0))
    assert win > 0, "gate_win 0 (G2 low-pass) is not implemented live"
    v8 = fit.get("v8") or {}
    assert int(v8.get("gate_src", -1)) == 2, "gate_src 2 fit without its v8 block"
    return src, win, float(v8["theta_G"])


def v8_globals(fit, cond, express):
    g, P = L7.v7_globals(fit, cond, express)
    src, win, thG = gate(fit)
    g.update(gate_src_GB=float(src), gate_win_GB=win if src else 100.0, theta_sh_GB=thG * 1e-3 if src else 1.5e-4,
             t_rest_GB=-1.0)
    return g, P


def run_task(job, hook=None):
    task, fit_path = job

    def hook8(task, sim, cell, syns, gids, meta, h):
        src = getattr(h, f"gate_src_GB_{L.MECH}")
        if src > 1.5:
            import numpy as np
            rf = L.find_record(task["path"], task["pair"], task["proto"])
            assert rf, "gate_src 2 needs the task's record (rest time)"
            d = np.load(rf)
            ev = np.concatenate([np.asarray(d[k], float).ravel() for k in ("prespikes", "postspikes") if k in d.files])
            t0 = float(np.asarray(d["t_ms"]).ravel()[0]) if "t_ms" in d.files else 0.0
            setattr(h, f"t_rest_GB_{L.MECH}", max(float(ev.min()) - 1.0, t0) if ev.size else t0)
        fin0 = hook(task, sim, cell, syns, gids, meta, h) if hook is not None else None

        def fin(out):
            hs = [s.hsynapse for _, s in syns]
            out["v8"] = dict(gate_src=float(src), gate_win=float(getattr(h, f"gate_win_GB_{L.MECH}")),
                             theta_sh_uM=1e3 * float(getattr(h, f"theta_sh_GB_{L.MECH}")),
                             t_rest=float(getattr(h, f"t_rest_GB_{L.MECH}")),
                             cai_rest=[float(x.cai_rest_GB) for x in hs], nsh=[float(x.nsh_GB) for x in hs],
                             tsh=[float(x.tsh_GB + (h.t - x.ton_GB if x.gsh_GB > 0.5 else 0.0)) for x in hs],
                             gsh_end=[float(x.gsh_GB) for x in hs])
            if fin0 is not None:
                fin0(out)
        return fin

    return L7.run_task(job, hook8)


L.MECH = "GluSynapseV8"
L.V5LIB = V8LIB
L.v5_globals = v8_globals
L.run_task = run_task
L.__file__ = os.path.abspath(__file__)        # live_v5._sub re-launches this file per task (--one)

if __name__ == "__main__":
    if "--one" not in sys.argv:
        _f = sys.argv[sys.argv.index("--fit") + 1]
        _fit = json.load(open(_f))
        print(f"[live_v8] {_f}: ecb_ref {L7.ecb_ref(_fit)}, gate (src, win, theta_G uM) {gate(_fit)}, mechanism {L.MECH} "
              f"({V8LIB})", flush=True)
    L.main()
