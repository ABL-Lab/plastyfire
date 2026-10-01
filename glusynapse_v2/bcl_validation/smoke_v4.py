"""Smoke test of the compiled GluSynapseV4 (run in compile_v4.sh): loads next to the DEES library, t_drive 4 events,
C1 gate and shadow rule are live on a single passive compartment."""
from neuron import h
import bluecellulab
bluecellulab.neuron.load_mechanisms("/project/rrg-emuller/dhuruva/DEES_cell_packages/")
assert h.nrn_load_dll("/project/rrg-emuller/dhuruva/plastyfire/glusynapse_v2/mod_build_v4/x86_64/libnrnmech.so") == 1
h.load_file("stdrun.hoc")
s = h.Section(name="s"); s.insert("pas"); s.L = s.diam = 10
syn = h.GluSynapseV4(0.5, sec=s)
for k, v in dict(t_drive_GB=4, vamp_mode_GB=1, theta_VA_GB=0.5, express_x_GB=1, i_scale_GB=1e-5, A_mglu_GB=0.4,
                 theta_Tg_GB=0.0, theta_Te_GB=0.0, pre_drive_GB=1).items():
    setattr(h, f"{k}_GluSynapseV4", v)
syn.theta_d_GB = syn.theta_p_GB = 1e9; syn.theta_dx_GB = syn.theta_px_GB = 1e9
syn.setRNG(1, 2, 3)
ns = h.NetStim(); ns.start = 10; ns.number = 5; ns.interval = 20; ns.noise = 0
nc = h.NetCon(ns, syn); nc.weight[0] = 1; nc.delay = 1
ic = h.IClamp(0.5, sec=s); ic.delay = 40; ic.dur = 3; ic.amp = 0.3
h.cvode_active(1); h.finitialize(-65); h.continuerun(200)
vals = {k: getattr(syn, k) for k in ("nev_GB", "S_GB", "T_GB", "bglu_GB", "Vg_GB", "gv_GB", "rho_GB", "rhox_GB", "dpre_GB")}
print("smoke GluSynapseV4:", {k: round(v, 6) for k, v in vals.items()})
assert vals["nev_GB"] >= 1, "no t_drive 4 event"
print("smoke OK")
# bin mode (offline-grid pre path): one mGluR step per arrival at the next 0.25 ms grid sample, dpre only from bin_step
for k, v in dict(bin_GB=0.25, bin_t0_GB=0.0, bin_t1_GB=150.0, A_NO_GB=5000.0, theta_NOi_GB=1e-9).items():
    setattr(h, f"{k}_GluSynapseV4", v)
h.finitialize(-65); h.continuerun(200)
vb = {k: getattr(syn, k) for k in ("nbin_GB", "tTsum_GB", "Nd_GB", "Td_GB", "dpre_GB", "nev_GB")}
print("smoke bin mode:", {k: round(v, 6) for k, v in vb.items()})
assert int(round(vb["nbin_GB"])) == 5, "bin mode: expected 5 mGluR steps"
print("smoke bin OK")
