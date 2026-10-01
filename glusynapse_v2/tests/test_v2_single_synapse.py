"""Single-synapse tests for GluSynapse_v2 (tsk T14).

A single dendrite-like section with Ca_HVA + CaDynamics_DC0 carries one GluSynapse (v1) and
one GluSynapse_v2 at the same location, driven by the same presynaptic spikes and the same
Random123 streams. A voltage clamp supplies postsynaptic depolarisation, and with it shaft Ca.

  test_v1_equivalence  v2 with A_mglu = A_NO = 0 reproduces v1 exactly (rho, Use, gmax, effcai)
  test_pre_ltd         post-before-pre depolarisation, mGluR on  -> dpre < 0
  test_pre_ltp         pre spike 10 ms before a strong depolarisation, NO on -> dpre > 0
  test_no_pre          depolarisation without presynaptic spikes -> dpre stays 0
  test_pre_only        presynaptic spikes at rest -> dpre stays 0

    python glusynapse_v2/tests/test_v2_single_synapse.py
"""
import os, sys
import numpy as np
from neuron import h

HERE = os.path.dirname(os.path.abspath(__file__))
# MECH_DIR overrides the compiled mechanism directory (e.g. glusynapse_v2/mod_build_td3), default build/
h.nrn_load_dll(os.path.join(os.environ.get("MECH_DIR") or os.path.join(HERE, "..", "build"), "x86_64", "libnrnmech.so"))
h.load_file("stdrun.hoc")


def make_cell():
    sec = h.Section(name="dend")
    sec.L, sec.diam, sec.nseg = 20.0, 1.0, 1
    sec.insert("pas"); sec.g_pas = 3e-5; sec.e_pas = -70.0
    sec.insert("Ca_HVA2"); sec.gCa_HVAbar_Ca_HVA2 = 1e-3
    sec.insert("CaDynamics_DC0")
    sec.cm = 1.0
    return sec


def make_syn(cls, sec, seed):
    s = getattr(h, cls)(sec(0.5))
    s.gmax0_AMPA = 1.0; s.gmax_d_AMPA = 1.0; s.gmax_p_AMPA = 2.0
    s.Use = 0.5; s.Use_d = 0.3; s.Use_p = 0.7; s.Nrrp = 2
    s.rho0_GB = 0.0
    s.theta_d_GB = 0.006; s.theta_p_GB = 0.012
    s.setRNG(seed, 7, 11)
    return s


def run(pre_times, clamp_steps, v2_params=None, tstop=None):
    """clamp_steps: list of (t_on, dur, v). Returns recorded traces for v1 and v2."""
    sec = make_cell()
    syn1 = make_syn("GluSynapse", sec, 1)
    syn2 = make_syn("GluSynapse_v2", sec, 1)
    for k, v in (v2_params or {}).items():
        setattr(h, f"{k}_GluSynapse_v2", v)

    vs = h.VecStim(); tv = h.Vector(sorted(pre_times)); vs.play(tv)
    ncs = []
    for s in (syn1, syn2):
        nc = h.NetCon(vs, s); nc.weight[0] = 1.0; nc.delay = 0.1; ncs.append(nc)

    clamp = h.SEClamp(sec(0.5)); clamp.rs = 1.0
    # piecewise clamp through a played vector: rest -70, steps as given
    tstop = tstop or (max([p for p in pre_times] + [a + d for a, d, _ in clamp_steps] + [0]) + 500.0)
    tvec = [0.0]; vvec = [-70.0]
    for a, d, v in sorted(clamp_steps):
        tvec += [a, a, a + d, a + d]; vvec += [-70.0, v, v, -70.0]
    tvec.append(tstop); vvec.append(-70.0)
    clamp.dur1 = 1e9
    vplay = h.Vector(vvec); tplay = h.Vector(tvec)
    vplay.play(clamp._ref_amp1, tplay, 1)

    rec = {}
    for name, s in (("v1", syn1), ("v2", syn2)):
        for var in ("rho_GB", "Use_GB", "gmax_AMPA", "effcai_GB", "cai_CR"):
            r = h.Vector(); r.record(getattr(s, f"_ref_{var}")); rec[f"{name}.{var}"] = r
    for var in ("dpre_GB", "T_GB", "N_GB", "Z_GB", "Use_tgt_GB", "ica_VDCC"):
        r = h.Vector(); r.record(getattr(syn2, f"_ref_{var}")); rec[f"v2.{var}"] = r
    cai = h.Vector(); cai.record(sec(0.5)._ref_cai); rec["shaft_cai"] = cai
    vv = h.Vector(); vv.record(sec(0.5)._ref_v); rec["v"] = vv
    for var in ("S_GB", "vbar_GB"):          # t_drive 3 states, only in mod_build_td3 and later
        if hasattr(syn2, f"_ref_{var}"):
            r = h.Vector(); r.record(getattr(syn2, f"_ref_{var}")); rec[f"v2.{var}"] = r
    t = h.Vector(); t.record(h._ref_t)

    h.dt = 0.025; h.steps_per_ms = 40; h.celsius = 34.0
    h.finitialize(-70.0)
    h.continuerun(tstop)
    out = {k: np.array(v) for k, v in rec.items()}
    out["t"] = np.array(t)
    # restore defaults so tests don't leak into each other
    for k in (v2_params or {}):
        setattr(h, f"{k}_GluSynapse_v2", DEFAULTS[k])
    return out


DEFAULTS = {}
def _snapshot_defaults():
    for k in ("A_mglu_GB", "A_NO_GB", "theta_T_GB", "tau_T_GB", "theta_NO_GB",
              "tau_NO_GB", "tau_Z_GB", "theta_Z_GB", "z_het_GB", "dpre_min_GB", "dpre_max_GB",
              "ca_scale_GB", "ca_sh_rest_GB", "pre_drive_GB", "i_scale_GB", "theta_Ti_GB", "theta_NOi_GB",
              "no_drive_GB", "theta_NOe_GB", "e_scale_GB", "theta_NOc_GB", "c_scale_GB", "theta_N_GB"):
        DEFAULTS[k] = getattr(h, f"{k}_GluSynapse_v2")


# ---------------------------------------------------------------------------
PAIRING = dict(n=20, period=200.0)   # 20 pairings at 5 Hz


def pairing(dt_post_minus_pre, v=0.0, dur=30.0):
    """pre at t0; clamp step onset at t0 + dt (dt < 0: post before pre)."""
    pre, steps = [], []
    for i in range(PAIRING["n"]):
        t0 = 100.0 + i * PAIRING["period"]
        pre.append(t0)
        steps.append((t0 + dt_post_minus_pre, dur, v))
    return pre, steps


def test_v1_equivalence():
    pre, steps = pairing(+10.0)
    r = run(pre, steps)
    for var in ("rho_GB", "Use_GB", "gmax_AMPA", "effcai_GB", "cai_CR"):
        a, b = r[f"v1.{var}"], r[f"v2.{var}"]
        assert np.array_equal(a, b), f"{var}: max |diff| {np.max(np.abs(a - b))}"
    assert np.all(r["v2.dpre_GB"] == 0.0)
    assert np.ptp(r["v1.rho_GB"]) > 0 or np.ptp(r["v1.effcai_GB"]) > 0, "protocol did nothing"
    return f"identical over {len(r['t'])} samples; effcai max {r['v1.effcai_GB'].max():.4f}, " \
           f"shaft cai max {r['shaft_cai'].max():.2e} mM"


def test_pre_ltd():
    # post (clamp to -10 mV, 30 ms) starts 40 ms BEFORE the pre spike and ends 10 ms before it
    pre, steps = pairing(-40.0, v=-10.0, dur=30.0)
    r = run(pre, steps, dict(A_mglu_GB=0.05))
    d = r["v2.dpre_GB"][-1]
    assert d < -1e-3, f"dpre {d}"
    return f"dpre {d:+.4f}, T max {r['v2.T_GB'].max():.3f}"


def test_pre_ltp():
    # pre 10 ms before a strong depolarisation (0 mV, 30 ms) -> shaft Ca spike -> NO
    pre, steps = pairing(+10.0, v=0.0, dur=30.0)
    r = run(pre, steps, dict(A_NO_GB=500.0, tau_Z_GB=20.0))
    d = r["v2.dpre_GB"][-1]
    assert d > 1e-3, f"dpre {d}"
    return f"dpre {d:+.4f}, N max {r['v2.N_GB'].max():.3f}, shaft cai max {r['shaft_cai'].max():.2e}"


def test_no_pre():
    _, steps = pairing(+10.0, v=0.0, dur=30.0)
    r = run([], steps, dict(A_NO_GB=500.0, A_mglu_GB=0.05))
    assert np.all(r["v2.dpre_GB"] == 0.0), r["v2.dpre_GB"].min()
    return "dpre 0 (no presynaptic activity)"


def test_pre_only():
    pre, _ = pairing(0.0)
    r = run(pre, [], dict(A_NO_GB=500.0, A_mglu_GB=0.05))
    d = r["v2.dpre_GB"][-1]
    assert abs(d) < 1e-4, f"dpre {d}"
    return f"dpre {d:+.2e}, shaft cai max {r['shaft_cai'].max():.2e} (EPSPs alone)"


if __name__ == "__main__":
    _snapshot_defaults()
    fails = 0
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                print(f"PASS {name}: {fn()}")
            except AssertionError as e:
                fails += 1; print(f"FAIL {name}: {e}")
    sys.exit(1 if fails else 0)
