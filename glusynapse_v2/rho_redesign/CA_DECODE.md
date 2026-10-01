# Can a spine decode the Ca source from one pool? (CA_DECODE)

Question (user, 2026-10-01): "the spine doesnt know which calcium source its gets, is there a way to decode this".
The current C1 potentiation gate reads c_VDCC, a separate pool fed only by the synapse's own spine VDCC current
(dc_VDCC/dt = -c_VDCC/tauE1 + (-ica_VDCC)/i_scale, tauE1 100 ms; pot = H(c* - theta_p) H(c_VDCC - theta_VDCC)).
Allowed replacements: own spine Ca (cai_CR / effcai), local shaft Ca, local Ca currents, own pre spikes, own state.
Diagnostic: `diag_ca_decode.py` (fit `results/v4_A0g_s3.json`, same groups as `diag_l23_ltp.py`).

## 1. Literature

All references checked against PubMed (DOI given). "Shown" = experiment; "modelled" = model only.

**(a) Nanodomain sensors tethered at the channel.** A spine is one diffusional compartment for *bulk* Ca, but not
for Ca within ~10-100 nm of an open channel. Sensors anchored there read the channel's own influx.
- Shown, spines: CaV2.3 (R-type) channels sit on CA1 spines, and their Ca selectively activates SK channels in the
  same spine; NMDAR Ca does not substitute. The authors call these functional Ca microdomains inside the spine
  (Bloodgood & Sabatini 2007, Neuron, doi:10.1016/j.neuron.2006.12.017). The SK-NMDAR loop is in Ngo-Anh et al. 2005
  (Nat Neurosci, doi:10.1038/nn1449).
- Shown, non-spine: CaMKII tethers to the CaV1.2 C-terminus and then acts as a local integrator of that channel's Ca
  (Hudmon et al. 2005, J Cell Biol, doi:10.1083/jcb.200505155). CaMKII near L-type channels encodes channel *opening
  frequency* rather than integrated Ca flux (Wheeler et al. 2008, J Cell Biol, doi:10.1083/jcb.200805048). For the
  same bulk [Ca]i rise, CaV1 nanodomain Ca is about 10x more effective than CaV2 Ca in signalling to CREB (Wheeler et
  al. 2012, Cell, doi:10.1016/j.cell.2012.03.041). These are SCG neurons and cultured hippocampal neurons, not
  neocortical spines.
- Shown, the NMDA counterpart: CaMKII binds GluN2B (NR2B) and is locked active there (Bayer et al. 2001, Nature,
  doi:10.1038/35081080). This is a sensor tethered to the *NMDA* source.
- Modelled: in an MCell spine model, CaM activation in the PSD depends on local Ca gradients near the channels and on
  pre/post order, beyond what the volume-averaged Ca shows (Keller et al. 2008, PLoS One,
  doi:10.1371/journal.pone.0002045).
- STDP link (shown): in L2/3 the LTD and LTP arms use different Ca sources, VDCC + mGluR for LTD and NMDAR for LTP
  (Nevian & Sakmann 2006, J Neurosci, doi:10.1523/JNEUROSCI.1749-06.2006; Bender et al. 2006, J Neurosci,
  doi:10.1523/JNEUROSCI.0176-06.2006). In both papers volume-averaged spine Ca did not predict the sign.

*Verdict on c_VDCC.* Read as "occupancy of a sensor tethered at the spine's VDCCs (nanodomain), integrated with the
sensor's ~100 ms memory", c_VDCC is defensible. Its input, -ica_VDCC, is what sets nanodomain Ca at the channel
mouth: nanodomain Ca follows the unitary current with sub-ms lag and does not mix with the bulk pool. The 100 ms
tauE1 has to be the sensor/kinase memory, not a Ca decay. Caveats: (i) direct evidence for a VDCC-tethered *plasticity*
sensor in spines is only the SK result. The CaMKII-CaV1 coupling comes from other preparations and is extrapolated.
(ii) In real cells the potentiation sensor is usually NMDA-tethered (GluN2B-CaMKII), so with this reading the gate is
a "VDCC-nanodomain coincidence" condition, not the LTP kinase itself. (iii) It is still a channel-labelled signal. If
the requirement is a strictly single-pool readout, the kinetic route (b) is the alternative, and section 2 tests it.

**(b) Kinetic decoding from one pool.**
- Shown: CaM binds Ca faster than the other spine buffers, so CaM directly intercepts incoming Ca (Faas et al. 2011,
  Nat Neurosci, doi:10.1038/nn.2746). The two lobes differ: the N-lobe is fast with low affinity and the C-lobe is
  slower with high affinity. C-lobe-only Ca2-CaM already binds and partially activates CaMKII (Shifman et al. 2006,
  PNAS, doi:10.1073/pnas.0606433103). [The exact lobe rate constants were not re-extracted here; use Faas 2011
  Table 1 if a parameterised sensor is built.]
- Modelled: frequency dependence of CaMKII activation arises from Ca-CaM-CaMKII binding kinetics in fluctuating Ca
  (Pepke et al. 2010, PLoS Comput Biol, doi:10.1371/journal.pcbi.1000675). An allosteric CaM model activates PP2B vs
  CaMKII differentially by Ca level (Stefan et al. 2008, PNAS, doi:10.1073/pnas.0804672105).
- Shown: in spines, CaMKIIa decodes input frequency/number supralinearly while calcineurin counts inputs (Fujii et al.
  2013, Cell Rep, doi:10.1016/j.celrep.2013.03.033).
- STDP models that read the Ca *time course*: the Rubin et al. detector model (P, A, V, D detectors with different
  kinetics; a veto on LTD by the right Ca time course; J Neurophysiol 2005, doi:10.1152/jn.00803.2004). Those authors
  also show that detectors based only on postsynaptic Ca are fragile to small bAP/synaptic changes. The
  biochemical CaMKII/PP1 bistable synapse is Graupner & Brunel 2007 (PLoS Comput Biol,
  doi:10.1371/journal.pcbi.0030221). Graupner & Brunel 2012 (PNAS, doi:10.1073/pnas.1109359109) and Shouval et al.
  2002 (PNAS, doi:10.1073/pnas.152343099) are amplitude/threshold rules (the current GluSynapse family) and do not
  decode the source.

*Verdict on kinetics.* Biologically this is a real mechanism: fast low-affinity sensors respond to brief, steep
influx (VDCC/bAP), and slow high-affinity sensors integrate slow NMDA Ca. Whether it *works in this model* depends on
whether the spine Ca trace (cacr, volume-averaged, 0.25 ms grid) still carries a kinetic signature that separates the
groups. Section 2 tests this directly.

## 2. Diagnostic: single-pool features vs vd_int

Features come from the synapse's own spine Ca (cacr = step mean of cai_CR - min, recovered exactly from effcai) on the
extracted grid (0.25 ms in the spike windows, 5 ms coarse elsewhere; all pairing windows lie in the fine part).
Pairing windows are [first pre - 30, last pre + 150] ms. Each feature is the mean over pairings.
- `dca_*`: dCa/dt peak; integral of pos(dCa/dt - k), k = 0 and k = ca_ref / {1, 3, 10} ms.
- `s_n{n}_kd{lo,mid,hi}_to{1,3,10}_*`: sensor dS/dt = kon Ca^n (1 - S) - koff S, Kd = {0.5, 1, 2} x ca_ref
  (half-saturation), tau_off = 1/koff. Readouts: `pk` peak S, `int` integral S, `tab` time S > 0.5, `R` peak of a
  100 ms low-pass of S (as tauE1), `Ratp` the same low-pass at samples where effcai > theta_p (gate coincidence).
- `hp{10,30,100}_*`: Ca - lowpass(Ca, tau): peak and positive integral.
- `sharp`: pairing peak / pairing integral.
- References: `vd_int`, `ca_pk`, `ca_int`, `sh_pk`, `m_p` (from the diag_l23_ltp code path), `cv_pk`/`cv_atp`
  (100 ms low-pass of own VDCC current, i.e. c_VDCC, its peak and its value while effcai > theta_p).

AUCs are P(L2/3 should-depress > L5 should-potentiate). Targets: auc_up <= ~0.25 (vd_int 0.20), auc_3ap (3AP+10
distal > proximal) <= ~0.3 (vd_int 0.21). A value >= 0.75 also separates, with the gate sign inverted.

<!-- TABLE -->
ca_ref = 0.003537 mM; n up-flippers dep 235, pot 112; 3AP distal 218, prox 142. Top 25 single-pool features by |auc_up - 0.5|, then references. Full table: `results/ca_decode_A0g_s3_auc.csv`.

| feat | family | auc_up | auc_all | auc_3ap | rho_vd_all | rho_vd_up | med_dep_up | med_pot_up |
|---|---|---|---|---|---|---|---|---|
| s_n1_kdlo_to10_int | s_n1 | 0.188 | 0.171 | 0.256 | 0.64 | 0.66 | 55.5 | 77.9 |
| s_n1_kdlo_to3_int | s_n1 | 0.191 | 0.172 | 0.261 | 0.63 | 0.648 | 54.9 | 75.9 |
| s_n1_kdlo_to1_int | s_n1 | 0.191 | 0.173 | 0.262 | 0.628 | 0.644 | 54.7 | 75.4 |
| sharp | sharp | 0.795 | 0.697 | 0.203 | 0.289 | 0.0329 | 0.0178 | 0.0154 |
| s_n1_kdmid_to10_int | s_n1 | 0.241 | 0.197 | 0.248 | 0.655 | 0.669 | 35.4 | 47.9 |
| s_n1_kdmid_to3_int | s_n1 | 0.242 | 0.197 | 0.251 | 0.646 | 0.657 | 34.9 | 46.9 |
| s_n1_kdmid_to1_int | s_n1 | 0.242 | 0.197 | 0.252 | 0.643 | 0.654 | 34.8 | 46.6 |
| s_n1_kdhi_to3_int | s_n1 | 0.283 | 0.214 | 0.247 | 0.659 | 0.666 | 20.9 | 26.8 |
| s_n1_kdhi_to1_int | s_n1 | 0.284 | 0.214 | 0.248 | 0.657 | 0.663 | 20.8 | 26.7 |
| s_n1_kdhi_to10_int | s_n1 | 0.285 | 0.214 | 0.245 | 0.666 | 0.676 | 21.4 | 27.2 |
| hp10_pint | hp | 0.328 | 0.227 | 0.223 | 0.817 | 0.812 | 0.0284 | 0.0452 |
| s_n2_kdlo_to10_int | s_n2 | 0.334 | 0.294 | 0.24 | 0.669 | 0.656 | 44.7 | 56.1 |
| dca_pos_k0 | dca | 0.345 | 0.218 | 0.218 | 0.865 | 0.843 | 0.005 | 0.00892 |
| s_n2_kdlo_to3_int | s_n2 | 0.346 | 0.296 | 0.242 | 0.653 | 0.622 | 42.6 | 52 |
| s_n2_kdlo_to1_int | s_n2 | 0.35 | 0.297 | 0.243 | 0.648 | 0.609 | 42.3 | 50.9 |
| hp30_pint | hp | 0.367 | 0.262 | 0.226 | 0.762 | 0.74 | 0.0549 | 0.0714 |
| s_n2_kdhi_to3_Ratp | s_n2 | 0.626 | 0.466 | 0.282 | 0.58 | 0.457 | 0.038 | 0.0262 |
| s_n2_kdhi_to10_Ratp | s_n2 | 0.626 | 0.467 | 0.282 | 0.58 | 0.456 | 0.0383 | 0.0263 |
| s_n2_kdhi_to1_Ratp | s_n2 | 0.626 | 0.466 | 0.282 | 0.579 | 0.456 | 0.0376 | 0.026 |
| s_n1_kdmid_to3_tab | s_n1 | 0.618 | 0.458 | 0.287 | 0.621 | 0.405 | 6.35 | 0.5 |
| s_n2_kdhi_to10_pk | s_n2 | 0.618 | 0.399 | 0.228 | 0.731 | 0.445 | 0.145 | 0.103 |
| s_n2_kdmid_to3_tab | s_n2 | 0.617 | 0.459 | 0.287 | 0.626 | 0.414 | 6.52 | 0.6 |
| s_n2_kdmid_to10_tab | s_n2 | 0.614 | 0.465 | 0.308 | 0.571 | 0.367 | 5.39 | 0 |
| s_n1_kdmid_to10_tab | s_n1 | 0.613 | 0.468 | 0.314 | 0.555 | 0.345 | 3.83 | 0 |
| s_n2_kdmid_to1_Ratp | s_n2 | 0.611 | 0.466 | 0.281 | 0.576 | 0.459 | 0.111 | 0.0814 |
| vd_int | ref | 0.201 | 0.164 | 0.205 | 1 | 1 | 0.000232 | 0.00244 |
| cv_pk | ref | 0.251 | 0.202 | 0.204 | 0.992 | 0.987 | 1.59e-06 | 1.07e-05 |
| cv_atp | ref | 0.279 | 0.402 | 0.273 | 0.65 | 0.97 | 8.74e-07 | 7.37e-06 |
| ca_int2 | ref | 0.353 | 0.239 | 0.242 | 0.68 | 0.66 | 0.188 | 0.228 |
| ca_int | ref | 0.353 | 0.239 | 0.242 | 0.68 | 0.66 | 0.188 | 0.228 |
| sh_pk | ref | 0.42 | 0.249 | 0.228 | 0.943 | 0.803 | 0.000586 | 0.000839 |
| ca_pk2 | ref | 0.559 | 0.336 | 0.22 | 0.814 | 0.595 | 0.00362 | 0.00328 |
| ca_pk | ref | 0.559 | 0.336 | 0.22 | 0.814 | 0.595 | 0.00362 | 0.00328 |
| m_p | ref | 0.519 | 0.383 | 0.211 | 0.625 | 0.577 | 1.22 | 1.23 |

Useful (|auc_up - 0.5| >= 0.25 and 3AP on the same side): s_n1_kdlo_to10_int, s_n1_kdlo_to3_int, s_n1_kdlo_to1_int, vd_int, s_n1_kdmid_to10_int, s_n1_kdmid_to3_int, s_n1_kdmid_to1_int

## Window check (pilot 22134502): the integral result is a confound
The pairing-window length alone gives AUC 0.182 (235 dep / 112 pot up-flippers), the same as the S_kd0.5 integral (0.191). Normalised by window length (mean S, AUC 0.591) or over a fixed 200 ms window (0.487), the sensor no longer separates the up-flippers. The integral was counting how long the protocol lasts, not which source the Ca came from. A single sensor on total spine Ca cannot decode the source here. The source-specific VDCC pool (nanodomain coupling, as for CaV2.3 → SK in Bloodgood & Sabatini 2007) stays the gate.
