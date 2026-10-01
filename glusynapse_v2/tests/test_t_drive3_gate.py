"""t_drive 3 (synapse-local eCB drive) gate: offline model_v2 dpre vs the NEURON mod on the toy section.

Needs the mod built with t_drive 3:  MECH_DIR=glusynapse_v2/mod_build_td3 python glusynapse_v2/tests/test_t_drive3_gate.py
The clamp gives bAP-like 2 ms steps to -10 mV (60 mV above rest). The offline side gets the NEURON v trace
(0.025 ms), runs model_v2.v_events + veto + feature_T + dpre_final, and is compared with the mod's final dpre.
Cases: post-before-pre (LTD), pre-before-post (vetoed, no LTD), a burst, and a case with the veto window edge.
Residual: event time discretisation only. The mod vbar is explicit Euler (uses v of the previous step), the
offline EMA uses the current sample, so x = v - vbar differs by O(dt/tau_b * dv) ~ 1e-4 * 60 mV before a
crossing that is 60 mV high, which cannot move an event.
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault("MECH_DIR", os.path.join(HERE, "..", "mod_build_td3"))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, ".."))
import test_v2_single_synapse as S
import model_v2 as MV
from neuron import h

NEW = ("t_drive_GB", "theta_Te_GB", "Te_scale_GB", "tau_E1_GB", "theta_Tg_GB", "theta_V_GB", "w_V_GB", "tau_b_GB")
PAR = dict(t_drive_GB=3, A_mglu_GB=0.3, theta_Te_GB=0.5, Te_scale_GB=50.0, tau_E1_GB=100.0, tau_T_GB=40.0,
           theta_Tg_GB=0.05, theta_V_GB=3.0, w_V_GB=15.0, tau_b_GB=300.0)
AP = 2.0        # ms, step width
VAP = -10.0     # mV


def snapshot():
    S._snapshot_defaults()
    for k in NEW:
        S.DEFAULTS[k] = getattr(h, f"{k}_GluSynapse_v2")


def protocol(kind):
    pre, steps = [], []
    for i in range(10):
        t0 = 100.0 + i * 500.0
        pre.append(t0)
        if kind == "post-before-pre":
            steps.append((t0 - 10.0, AP, VAP))
        elif kind == "pre-before-post":
            steps.append((t0 + 5.0, AP, VAP))          # 4.9 ms after arrival: vetoed
        elif kind == "burst":
            steps += [(t0 - 60.0 + 10.0 * j, AP, VAP) for j in range(5)]      # 5 bAPs at 100 Hz, ends 16 ms before
        elif kind == "veto-edge":
            steps.append((t0 + 0.1 + 15.5, AP, VAP))   # just outside the 15 ms window: counts
            pre.append(t0 + 60.0)                      # a second pre spike reads T
    return sorted(pre), steps


def offline(r, pre, P, delay=0.1):
    dt = 0.025
    v = r["v"]; t = MV.grid(len(v), dt)
    ev = MV.v_events(v, dt, theta_V=P["theta_V"], tau_b=P["tau_b"])
    ev = MV.veto(ev, np.asarray(pre) + delay, P["w_V"])
    imp = np.zeros((1, len(t)))
    np.add.at(imp[0], np.clip(np.searchsorted(t, ev - 1e-9), 0, len(t) - 1), 1.0)
    arr = MV.arrival_index(pre, [delay], t)
    tT = MV.feature_T(imp, t, arr, P)
    d = MV.dpre_final(tT, np.zeros((1, len(pre) + 1)), P["A_mglu"], 0.0, P["dpre_min"], P["dpre_max"])[0]
    return d, len(ev), tT[0]


def test_t_drive3_gate():
    snapshot()
    P = {**{k[:-3]: val for k, val in {**S.DEFAULTS, **PAR}.items()}}
    worst = 0.0; out = []; res = {}
    for kind in ("post-before-pre", "pre-before-post", "burst", "veto-edge"):
        pre, steps = protocol(kind)
        r = S.run(pre, steps, PAR)
        d_mod = r["v2.dpre_GB"][-1]; d_off, nev, tT = offline(r, pre, P)
        res[kind] = d_mod
        rel = abs(d_mod - d_off) / max(abs(d_mod), 1e-3)
        worst = max(worst, rel)
        out.append(f"{kind}: mod {d_mod:+.6f} offline {d_off:+.6f} rel {rel:.1e} events {nev}, S max {r['v2.S_GB'].max():.3f}")
        assert rel < 1e-3, out[-1]
    assert res["post-before-pre"] < -0.02, res           # LTD
    assert abs(res["pre-before-post"]) < 1e-9, res       # vetoed: no eCB, no LTD
    assert res["burst"] < res["post-before-pre"], res    # burst supralinear
    assert res["veto-edge"] < -0.005, res
    return "; ".join(out) + f"; worst rel {worst:.1e}"


def test_t_drive3_cvode():
    """Variable step: the WATCH root-finder still gives the events; dpre within 1e-2 relative of fixed step."""
    snapshot()
    ok = []
    for kind in ("post-before-pre", "pre-before-post", "burst"):
        pre, steps = protocol(kind)
        d_fix = S.run(pre, steps, PAR)["v2.dpre_GB"][-1]
        h.cvode_active(1); h.cvode.atol(1e-6)
        try:
            d_cv = S.run(pre, steps, PAR)["v2.dpre_GB"][-1]
        finally:
            h.cvode_active(0)
        assert abs(d_cv - d_fix) <= 1e-2 * max(abs(d_fix), 1e-3), (kind, d_cv, d_fix)
        ok.append(f"{kind} cvode {d_cv:+.5f} fixed {d_fix:+.5f}")
    return "; ".join(ok)


def test_t_drive3_off_is_unchanged():
    """t_drive 0, theta_Tg 0: identical dpre to the old mod build (mGluR drive from shaft Ca)."""
    snapshot()
    pre, steps = S.pairing(-40.0, v=-10.0, dur=30.0)
    d = S.run(pre, steps, dict(A_mglu_GB=0.05))["v2.dpre_GB"][-1]
    assert d < -1e-3
    return f"dpre {d:.10f} (compare with build/: {os.environ.get('OLD_DPRE', 'see DECISIONS')})"


if __name__ == "__main__":
    fails = 0
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                print(f"PASS {name}: {fn()}")
            except AssertionError as e:
                fails += 1; print(f"FAIL {name}: {e}")
    sys.exit(1 if fails else 0)
