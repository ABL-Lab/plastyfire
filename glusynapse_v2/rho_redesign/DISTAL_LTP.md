# Distal t-LTP for Letzkus 3AP 200 Hz dt -10 (data 1.42 +- 0.09)

Diagnosis job 22115348 (diag_distal.py, A0_joint_s4 parameters, no refit). Outputs: results/distal_A0_joint_s4_{syn,rec,sep}.csv and figs/fig5_distal.png.

## Conclusion

The rule is not what stops distal LTP. In these simulations, no synapse-local signal is larger at distal 3AP -10 than at
3AP +10 while also reversing that preference at proximal synapses. The cause is the local voltage waveform, i.e. the
morphology and emodel. The current Ca-threshold rule does give -10 LTP in the one cell whose waveform has the right
shape (see below). So no change to the rule is recommended for this target.

## What Letzkus 2006 reports (J Neurosci 26:10420, Figs 4-8)

- **dt sign.** dt -10 means pre comes 10 ms after the burst onset (post before pre). The model protocols match this:
  pre minus post onset is +10.0 ms in the dt-10 files and -10.0 ms in the dt+10 files.
- **Pharmacology.**
  - Distal -10 LTP is blocked by APV (0.97) and by 100 uM NiCl2 (0.99), so it needs NMDARs and T/R-type dendritic Ca
    electrogenesis.
  - Distal +10 LTD is blocked by APV (1.03) but not by NiCl2 (0.81).
- **Mechanism (their model, Fig 8).**
  - At distal sites the dendritic depolarisation peaks after the 3rd AP, as a slow Ca-spike plateau lasting about
    20-30 ms.
  - At -10 the EPSP lands on that peak, giving maximal NMDA activation; at +10 the EPSP arrives before it.
  - A plain threshold on integrated NMDA conductance then reproduces both the proximal and the distal timing rules.
    That is the same class of rule as ours.
- **Synapse location.** Their distal synapses are on the apical dendrite: dendritic recordings at 430-660 um, Ca spike
  generated in the tuft. Ours are basal/oblique pairs classified by EPSP rise time (> 2.7 ms).

## Diagnosis (distal: 139 synapse-protocol rows at -10 and +10, 531 at 1AP +10)

Each per-pairing window is aligned to the synapse's own pre arrival. AUC(-10 > +10) is computed over synapses. A usable
signal needs AUC > 0.5 at distal synapses and < 0.5 at proximal ones (proximal data: +10 gives LTP 1.30, -10 gives LTD 0.89).

| signal | AUC distal | AUC proximal | reading |
|---|---|---|---|
| ca_bef (spine Ca in the 30 ms before own arrival) | 0.90 | 1.00 | favours -10 at both locations: it just means the AP came first |
| vd_int / sh_int (VDCC current or shaft Ca over the window) | 0.55 / 0.55 | 0.52 / 0.53 | paired difference about +5%, also at both locations |
| ca_late (spine Ca from arr+10 to arr+60, the NMDA window) | 0.38 | 0.20 | favours +10 at both locations |
| m_p (pairing effcai peak / theta_p) | 0.43 | 0.36 | favours +10 |
| n_cev_aft, b_cev (VDCC events after own arrival; NMDA-bound b at those events) | 0.14, 0.28 | 0.02, 0.06 | at -10 every VDCC event comes before arrival (b = 7e-7) |

- **Distal synapses with bAP Ca.**
  - Median ca_late: 0.024 at -10 vs 0.038 at +10 (mM ms).
  - P_up: 0.03 at -10 vs 0.31 at +10.
  - So in our cells a distal synapse behaves like a weaker proximal one, with the same timing preference rather than
    the reversed one.
- **What fig5 panels a-c show.** The burst's VDCC current and Ca are finished about 5 ms after the -10 arrival. There
  is no post-burst plateau for the EPSP to coincide with.
- **The exception: pair 20877-180143** (post cell 180143, all 13 synapses at the c_post floor).
  - At -10 the VDCC current after own arrival is 0.012, against 3e-6 at 1AP +10.
  - That gives ca_late 0.9-1.5 (40-60x the median) and m_p 9-14.
  - The record ratio is 2.08, versus 0.84 at 1AP +10.
  - This is the Letzkus mechanism, a regenerative local depolarisation triggered by burst plus EPSP, appearing in one
    cell. The A0 rule potentiates it without any change. This matches the earlier Letzkus Ca audit, which found a Ca
    spike in 1 of 8 cells.

## Candidates (synapse-local, Ca-based, uniform across pathways, <= 1 parameter, reduce to A0 when off)

1. **Post-before-pre Ca credit.** Potentiation reads effcai + k * (spine Ca trace at own pre arrival), with k = 0 = A0.
   - *For:* ca_bef has AUC 0.90 at distal synapses.
   - *Against:* AUC 1.00 at proximal synapses, so it would potentiate proximal -10 (data 0.89) and the L5 Sjostrom -10
     protocols (LTD).
   - **Rejected.**
2. **Slow VDCC / shaft-Ca potentiation drive** (an analogue of the T/R-type NiCl2 dependence). Potentiation reads
   effcai + k * filtered -ica_VDCC (or shaft Ca), with k = 0 = A0.
   - *Against:* the effect is +5% (AUC 0.55) and the same at proximal synapses (0.52). It cannot produce 1.42 without
     potentiating everything.
   - **Rejected.**
3. **No rule change; fix the preparation (recommended).** Letzkus distal needs synapses where the burst produces a
   dendritic plateau that lasts past the EPSP.
   - (a) Until then, treat `letzkus_3ap_200hz_dt-10ms@distal` (and its APV/NiCl2 rows) as unreachable with this
     placement, and leave it out of the joint chi2 or report it as validation only.
   - (b) Place the Letzkus distal synapses on apical obliques/tuft at > 400 um, in an emodel with BAC firing and Ca
     spikes (check og-delta vs antic-delta for tuft Ca-spike generation), then re-extract the npz files and rescore A0
     unchanged.
   - Evidence: pair 20877-180143 above, and Letzkus Fig 8, whose own model needs nothing beyond an NMDA-Ca threshold.

The +10 over-potentiation (distal 3AP +10 0.99 vs 0.79; 1AP +10 0.88 distal and 0.82 proximal vs 0.72) is a separate
problem and is not addressed here.
