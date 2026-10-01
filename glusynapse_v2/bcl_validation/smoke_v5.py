"""Smoke test of the compiled GluSynapseV5 (run in compile_v5.sh): loads next to the DEES library; on a single passive
compartment the v5_mode 2 eCB step fires at an own pre arrival after a post depolarisation (W > theta_eCB ->
dpre = dpre_min), is off at v5_mode 0, and W stays below Vg when the own arrivals precede the depolarisation
((1 - min(b, 1)) weight)."""
from neuron import h
import bluecellulab
bluecellulab.neuron.load_mechanisms("/project/rrg-emuller/dhuruva/DEES_cell_packages/")
assert h.nrn_load_dll("/project/rrg-emuller/dhuruva/plastyfire/glusynapse_v2/mod_build_v5/x86_64/libnrnmech.so") == 1
h.load_file("stdrun.hoc")
M = "GluSynapseV5"
s = h.Section(name="s"); s.insert("pas"); s.L = s.diam = 10
syn = h.GluSynapseV5(0.5, sec=s)
base = dict(t_drive_GB=4, vamp_mode_GB=1, theta_VA_GB=1e-6, express_x_GB=1, i_scale_GB=1e-5, A_mglu_GB=0.0,
            A_NO_GB=0.0, pre_drive_GB=1, dpre_min_GB=-0.29, A_eCB_GB=1.0, theta_eCB_GB=1e-6, tau_glu_GB=70.0)
syn.theta_d_GB = syn.theta_p_GB = 1e9; syn.theta_dx_GB = syn.theta_px_GB = 1e9
syn.setRNG(1, 2, 3)
ns = h.NetStim(); ns.noise = 0
nc = h.NetCon(ns, syn); nc.weight[0] = 1; nc.delay = 1
ic = h.IClamp(0.5, sec=s); ic.delay = 40; ic.dur = 3; ic.amp = 0.3
KEYS = ("necb_GB", "dpre_GB", "W_GB", "Vg_GB", "bglu_GB", "nev_GB", "rho_GB", "rhox_GB")


def run(mode, start, number, interval, tstop):
    for k, v in {**base, "v5_mode_GB": mode}.items():
        setattr(h, f"{k}_{M}", v)
    ns.start = start; ns.number = number; ns.interval = interval
    h.cvode_active(1); h.finitialize(-65); h.continuerun(tstop)
    return {k: getattr(syn, k) for k in KEYS}


# 1) post (40-43 ms) before pre (arrivals 61, 81, 101): eCB at the first arrival
a = run(2, 60, 3, 20, 150)
print("smoke v5_mode 2 post->pre:", {k: round(v, 8) for k, v in a.items()})
assert a["necb_GB"] >= 1 and abs(a["dpre_GB"] + 0.29) < 1e-12, "v5_mode 2: no eCB step"
# 2) the same with v5_mode 0 (V4 behaviour): no eCB step
b = run(0, 60, 3, 20, 150)
print("smoke v5_mode 0:", {k: round(v, 8) for k, v in b.items()})
assert b["necb_GB"] == 0 and b["dpre_GB"] == 0.0, "v5_mode 0 changed dpre"
# 3) own arrivals 31-43 ms (b >= 1 during the clamp) then the clamp: W << Vg
c = run(2, 30, 5, 3, 60)
print("smoke v5_mode 2 pre->post:", {k: round(v, 8) for k, v in c.items()})
assert c["Vg_GB"] > 0 and c["W_GB"] < 0.5 * c["Vg_GB"], "W is not (1 - b)-weighted"
print("smoke v5 OK")
