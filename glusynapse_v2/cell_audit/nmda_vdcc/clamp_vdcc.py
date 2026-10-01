"""Mod-only test: clamp a 1-compartment spine-bearing segment to a bAP waveform and read GluSynapse spine Ca (VDCC only).
Waveform: rest -72 mV, peak Vp, alpha-like rise 0.5 ms, exponential decay (tau_d) -> half-width ~ 0.7*tau_d + 0.5.
Output: peak cai_CR - rest (uM) vs Vp, tau_d, spine volume, and GLOBAL variants."""
import json, sys, itertools
import numpy as np
from neuron import h
h.nrn_load_dll("/project/rrg-emuller/dhuruva/DEES_cell_packages/x86_64/.libs/libnrnmech.so")
h.load_file("stdrun.hoc"); h.celsius = 34
sec = h.Section(name="d"); sec.L = sec.diam = 5; sec.insert("pas")
syn = h.GluSynapse(sec(0.5))
clamp = h.SEClamp(sec(0.5)); clamp.rs = 1e-3; clamp.dur1 = 1e9
dt = 0.005
def wave(vp, td, rest=-72.0, t0=20.0, T=120.0):
    t = np.arange(0, T, dt); v = np.full_like(t, rest); m = t >= t0; s = t[m] - t0
    shape = (1 - np.exp(-s / 0.25)) * np.exp(-s / td); shape /= shape.max(); v[m] = rest + (vp - rest) * shape
    return t, v
def run(vp, td, vol, glob):
    for k, x in glob.items(): setattr(h, f"{k}_GluSynapse", x)
    syn.volume_CR = vol
    t, v = wave(vp, td); tv, vv = h.Vector(t), h.Vector(v); vv.play(clamp._ref_amp1, tv, 1)
    ca = h.Vector(); ca.record(syn._ref_cai_CR); vr = h.Vector(); vr.record(sec(0.5)._ref_v)
    h.dt = dt; h.finitialize(-72); h.continuerun(t[-1])
    ca = np.array(ca); vr = np.array(vr)
    hw = dt * np.sum(vr > (vr.min() + vr.max()) / 2)
    return float((ca.max() - 70e-6) * 1e3), float(vr.max()), float(hw)
base = dict(gca_bar_VDCC=0.0744, ljp_VDCC=0.0)
out = []
for name, glob in (("mod", base), ("g1.65", dict(base, gca_bar_VDCC=0.0744 * 1.65)), ("g2", dict(base, gca_bar_VDCC=0.1488)),
                   ("s10g3", dict(base, gca_bar_VDCC=0.2232, ljp_VDCC=10.0))):
    for vp, td, vol in itertools.product((-50, -40, -30, -20, -10, 0, 10, 20, 30), (1.0, 2.0, 3.0), (0.087, 0.153)):
        ca, vpk, hw = run(vp, td, vol, glob)
        out.append(dict(var=name, vp=vp, td=td, vol=vol, ca=ca, vpk=vpk, hw=hw))
for k, x in base.items(): setattr(h, f"{k}_GluSynapse", x)
json.dump(out, open("/project/rrg-emuller/dhuruva/plastyfire/glusynapse_v2/cell_audit/nmda_vdcc/clamp_vdcc.json", "w"))
for r in out:
    if r["var"] == "mod" and r["td"] == 2.0: print(r)
