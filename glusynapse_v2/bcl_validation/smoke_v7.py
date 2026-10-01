"""Smoke test of the compiled GluSynapseV7 (run in compile_v7.sh): loads next to the DEES library; on a single passive
compartment (v5_mode 2, theta_eCB 10 pool units):
  1) veto_T 0: an own arrival 8 ms after a post depolarisation steps at once (= V5);
  2) veto_T 25, a second depolarisation 9 ms after that arrival: the step is vetoed (nveto 1, necb 0);
  3) veto_T 25, arrivals 18-58 ms after the depolarisation, nothing after them: the deferred steps are applied;
  4) uE_GB 1e3 (theta_eCB,i = 1e4 > W): no trigger at all;
and nvpend_GB is 0 at the end of every run."""
from neuron import h
import bluecellulab
bluecellulab.neuron.load_mechanisms("/project/rrg-emuller/dhuruva/DEES_cell_packages/")
assert h.nrn_load_dll("/project/rrg-emuller/dhuruva/plastyfire/glusynapse_v2/mod_build_v7/x86_64/libnrnmech.so") == 1
h.load_file("stdrun.hoc")
M = "GluSynapseV7"
s = h.Section(name="s"); s.insert("pas"); s.L = s.diam = 10
syn = h.GluSynapseV7(0.5, sec=s)
base = dict(t_drive_GB=4, vamp_mode_GB=1, theta_VA_GB=1e-6, express_x_GB=1, i_scale_GB=1e-5, A_mglu_GB=0.0,
            A_NO_GB=0.0, pre_drive_GB=1, dpre_min_GB=-0.29, A_eCB_GB=1.0, theta_eCB_GB=10.0, tau_glu_GB=70.0,
            v5_mode_GB=2, bin_GB=0.0)
syn.theta_d_GB = syn.theta_p_GB = 1e9; syn.theta_dx_GB = syn.theta_px_GB = 1e9
syn.setRNG(1, 2, 3)
ns = h.NetStim(); ns.noise = 0
nc = h.NetCon(ns, syn); nc.weight[0] = 1; nc.delay = 1
ic = h.IClamp(0.5, sec=s); ic.delay = 40; ic.dur = 3; ic.amp = 0.3
ic2 = h.IClamp(0.5, sec=s); ic2.delay = 60; ic2.dur = 3
KEYS = ("necb_GB", "nveto_GB", "nvpend_GB", "dpre_GB", "W_GB", "Vg_GB", "qca_GB", "bglu_GB")


def run(veto, start, number, interval, amp2, uE, tstop):
    for k, v in {**base, "veto_T_GB": veto}.items():
        setattr(h, f"{k}_{M}", v)
    syn.uE_GB = uE
    ns.start = start; ns.number = number; ns.interval = interval; ic2.amp = amp2
    h.cvode_active(1); h.finitialize(-65); h.continuerun(tstop)
    r = {k: getattr(syn, k) for k in KEYS}
    assert r["nvpend_GB"] == 0, f"pending veto decisions at the end: {r}"
    return r


a = run(0.0, 50, 1, 10, 0.3, 1.0, 150)
print("1) veto_T 0, post -> pre -> post:", {k: round(v, 6) for k, v in a.items()})
assert a["necb_GB"] == 1 and a["nveto_GB"] == 0 and abs(a["dpre_GB"] + 0.29) < 1e-12, "veto_T 0 is not V5"
b = run(25.0, 50, 1, 10, 0.3, 1.0, 150)
print("2) veto_T 25, post -> pre -> post:", {k: round(v, 6) for k, v in b.items()})
assert b["nveto_GB"] == 1 and b["necb_GB"] == 0 and b["dpre_GB"] == 0.0, "bAP 9 ms after the arrival did not veto"
c = run(25.0, 60, 3, 20, 0.0, 1.0, 150)
print("3) veto_T 25, post -> pre x3:", {k: round(v, 6) for k, v in c.items()})
assert c["necb_GB"] >= 1 and c["nveto_GB"] == 0 and abs(c["dpre_GB"] + 0.29) < 1e-12, \
    "no deferred step without a bAP in the window (resting VDCC charge above theta_eCB?)"
d = run(25.0, 50, 1, 10, 0.3, 1e3, 150)
print("4) uE 1e3:", {k: round(v, 6) for k, v in d.items()})
assert d["necb_GB"] == 0 and d["nveto_GB"] == 0, "uE_GB does not scale theta_eCB"
print("smoke v7 OK")
