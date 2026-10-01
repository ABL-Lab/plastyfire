"""Smoke test of the compiled GluSynapseV8 (run in compile_v8.sh), loaded next to the DEES library and GluSynapseV7, on a
single passive compartment (v5_mode 2, theta_eCB 10, veto_T 25, pot / dep thresholds low so rho moves):
  A) gate_src 0: the same run with GluSynapseV7 and GluSynapseV8 gives bit-identical traces (t, v, rho, effcai) and
     final states (rho, rhox, dpre, necb, nveto, W, Vg, qca, Use, gmax_AMPA, cai_CR);
  B) gate_src 2, segment cai set by events (no Ca mechanism writes it): rest frozen at t_rest_GB, licence opens at the
     crossing, stays open gate_win after the last fall (a second crossing inside the window extends it, a 0.1 uM bump
     does not open it), nsh / tsh counts;
  C) gate_src 2 replaces the Vg licence: flat shaft Ca -> pot is never licensed (rho falls below rho0 although Vg >
     theta_VA); shaft Ca held above theta -> rho rises above rho0.
nvpend_GB is 0 at the end of every run."""
import numpy as np
from neuron import h
import bluecellulab
bluecellulab.neuron.load_mechanisms("/project/rrg-emuller/dhuruva/DEES_cell_packages/")
V2 = "/project/rrg-emuller/dhuruva/plastyfire/glusynapse_v2"
assert h.nrn_load_dll(f"{V2}/mod_build_v7/x86_64/libnrnmech.so") == 1
assert h.nrn_load_dll(f"{V2}/mod_build_v8/x86_64/libnrnmech.so") == 1
h.load_file("stdrun.hoc")
s = h.Section(name="s"); s.insert("pas"); s.L = s.diam = 10
base = dict(t_drive_GB=4, vamp_mode_GB=1, theta_VA_GB=1e-6, express_x_GB=1, i_scale_GB=1e-5, A_mglu_GB=0.0,
            A_NO_GB=0.0, pre_drive_GB=1, dpre_min_GB=-0.29, A_eCB_GB=1.0, theta_eCB_GB=10.0, tau_glu_GB=70.0,
            v5_mode_GB=2, bin_GB=0.0, veto_T_GB=25.0)
ns = h.NetStim(); ns.noise = 0; ns.start = 50; ns.number = 5; ns.interval = 20
ic = h.IClamp(0.5, sec=s); ic.delay = 40; ic.dur = 3; ic.amp = 0.3
ic2 = h.IClamp(0.5, sec=s); ic2.delay = 60; ic2.dur = 3; ic2.amp = 0.3
KEYS = ("rho_GB", "rhox_GB", "dpre_GB", "necb_GB", "nveto_GB", "nvpend_GB", "W_GB", "Vg_GB", "qca_GB", "Use_GB",
        "gmax_AMPA", "cai_CR", "effcai_GB")
TH = 1e-7                 # theta_d = theta_p (us/liter): pot and dep on whenever effcai > TH
CA0 = 5e-5                # cai0_ca_ion (mM)


def run(mech, tstop=400, glob=None, rho0=0.5, cai_ev=(), trec=False, probes=()):
    for k, v in {**base, **(glob or {})}.items():
        setattr(h, f"{k}_{mech}", v)
    syn = getattr(h, mech)(0.5, sec=s)
    assert len(list(s(0.5).point_processes())) == 3, "a previous synapse is still on the section"
    syn.theta_d_GB = syn.theta_p_GB = TH; syn.theta_dx_GB = syn.theta_px_GB = 1e9; syn.rho0_GB = rho0
    syn.setRNG(1, 2, 3)
    nc = h.NetCon(ns, syn); nc.weight[0] = 1; nc.delay = 1
    rec = {}
    if trec:
        for k, ref in (("t", h._ref_t), ("v", s(0.5)._ref_v), ("rho", syn._ref_rho_GB), ("eff", syn._ref_effcai_GB)):
            rec[k] = h.Vector(); rec[k].record(ref)
    pr = []

    def setca(x):
        s(0.5).cai = x; h.cvode.re_init()

    h.cvode_active(1); h.finitialize(-65)
    for tt, x in cai_ev:
        h.cvode.event(tt, lambda x=x: setca(x))
    for tt in probes:
        h.cvode.event(tt, lambda: pr.append((h.t, syn.gsh_GB)))
    h.continuerun(tstop)
    r = {k: getattr(syn, k) for k in KEYS}
    if mech == "GluSynapseV8":
        r.update({k: getattr(syn, k) for k in ("gsh_GB", "nsh_GB", "tsh_GB", "ton_GB", "cai_rest_GB")})
    assert r["nvpend_GB"] == 0, f"pending veto decisions at the end: {r}"
    tr = {k: np.array(v) for k, v in rec.items()}
    del nc, syn
    s(0.5).cai = CA0
    return r, tr, pr


# A) gate_src 0 = V7, bit for bit
a7, t7, _ = run("GluSynapseV7", trec=True)
a8, t8, _ = run("GluSynapseV8", glob=dict(gate_src_GB=0, t_rest_GB=-1), trec=True)
print("A) V7:", {k: a7[k] for k in KEYS})
print("A) V8 gate_src 0:", {k: a8[k] for k in KEYS})
diff = [k for k in KEYS if a7[k] != a8[k]] + [k for k in t7 if t7[k].shape != t8[k].shape or np.any(t7[k] != t8[k])]
assert not diff, f"gate_src 0 differs from V7 in {diff}"
assert abs(a7["rho_GB"] - 0.5) > 1e-6 and a7["necb_GB"] + a7["nveto_GB"] >= 1, "test does not exercise rho / eCB"
print(f"A) gate_src 0 == V7 bit for bit: {len(KEYS)} states, {len(t7['t'])} trace samples x 4; rho 0.5 -> "
      f"{a7['rho_GB']:.6f}, necb {a7['necb_GB']:.0f}, nveto {a7['nveto_GB']:.0f}")

# B) licence timing (theta 0.15 uM, win 100): rest 6e-5 mM set at t 20 < t_rest 30
R = 6e-5
ev = [(20, R), (100, R + 2e-4), (110, R), (250, R + 1e-4), (260, R), (300, R + 3e-4), (301, R), (330, R + 3e-4),
      (331, R)]
pt = (99, 101, 109, 205, 209, 211, 255, 299, 302, 332, 430, 432, 480)
want = (0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0)     # open 100-210 (fall 110 + 100), 0.1 uM bump 250 shut, 300-431
h.cvode.maxstep(0.1)      # B / C only (A runs both mechanisms with the defaults): crossing detection lag <= 0.1 ms
b, _, pr = run("GluSynapseV8", tstop=500, glob=dict(gate_src_GB=2, gate_win_GB=100, theta_sh_GB=1.5e-4, t_rest_GB=30),
               cai_ev=ev, probes=pt)
print("B) probes (t, gsh):", pr, "| nsh", b["nsh_GB"], "tsh", round(b["tsh_GB"], 3), "rest", b["cai_rest_GB"])
assert abs(b["cai_rest_GB"] - R) < 1e-15, f"rest {b['cai_rest_GB']} != {R}"
assert tuple(int(g) for _, g in pr) == want, f"licence timeline {pr} != {want}"
assert b["nsh_GB"] == 2 and abs(b["tsh_GB"] - (110 + 131)) < 0.5 and b["gsh_GB"] == 0, "nsh / tsh"
print("B) licence timeline OK")

# C) the shaft licence replaces Vg: flat cai -> no LTP, cai held above theta -> LTP
g2 = dict(gate_src_GB=2, gate_win_GB=100, theta_sh_GB=1.5e-4, t_rest_GB=30)
c0, _, _ = run("GluSynapseV8", glob=g2)
c1, _, _ = run("GluSynapseV8", glob=g2, cai_ev=[(35, CA0 + 1e-3)])
print(f"C) rho0 0.5 -> V7 {a7['rho_GB']:.6f}, gate_src 2 flat shaft {c0['rho_GB']:.6f} (nsh {c0['nsh_GB']:.0f}), "
      f"shaft held +1 uM {c1['rho_GB']:.6f} (nsh {c1['nsh_GB']:.0f}, gsh {c1['gsh_GB']:.0f})")
assert c0["nsh_GB"] == 0 and c0["rho_GB"] < 0.5 - 1e-6, "flat shaft Ca licensed LTP (Vg licence not replaced?)"
assert c1["nsh_GB"] == 1 and c1["gsh_GB"] == 1 and c1["rho_GB"] > 0.5 + 1e-6, "shaft licence did not license LTP"
print("smoke v8 OK")
