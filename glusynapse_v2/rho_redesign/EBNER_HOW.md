# How Ebner et al. 2019 reproduced the STDP data, and what we are missing (2026-10-02)

Question: how did Ebner, Clopath, Jedlicka & Cuntz 2019 reproduce so many STDP results, including location dependence?
What do they have, and what do the source papers of our targets have, that GluSynapse_v2 lacks?
Literature and code reading only, no compute. "Stated" = in the paper or code; "inferred" = ours.

Exact title: **"Unifying Long-Term Plasticity Rules for Excitatory Synapses by Modeling Dendrites of Cortical Pyramidal
Neurons"**, Cell Reports 29:4295-4307 (2019), PMID 31875541, doi 10.1016/j.celrep.2019.11.068. Code: ModelDB 251493
(github.com/ModelDBRepository/251493: syn_4p.mod, Fig2B.hoc, Fig3.hoc, Fig4B.hoc, models/L5PCbiophys4.hoc).

## Short answer
- **Ebner did not fit everything with one parameter set.** They used **three parameter sets**, one per dataset:
  set 1 for Sjostrom 2001, set 2 for Nevian 2006, set 3 for Letzkus 2006. The 23 shape parameters are shared. The 4 pathway
  amplitudes were re-tuned per dataset, by hand, over up to a 13x range. The paper says so: "we had to readjust plasticity
  amplitudes ... to reproduce different experiments" (Discussion).
- **The rule is voltage-based, not Ca-based.** It extends Clopath 2010. Every drive is the local membrane voltage at the
  synapse compartment plus presynaptic event traces. There is no Ca variable, and the NMDA conductance is not used by the
  rule. That breaks our no-voltage rule as written, but each term has a Ca analogue.
- **Location dependence comes almost entirely from the neuron model.** They used the Hay 2011 L5b cell (biophys4, adult
  P36 morphology). It was fitted to a **45 +- 10 mV bAP at 620 um and 36 +- 9 mV at 800 um**, and has a 100x Ca_LVA hot zone at
  700-900 um. Their 200 Hz bursts evoke a dendritic Ca spike delayed by about 20 ms at 669 um. The rule has no
  distance-dependent term. Our emodel gives 2-13 mV distal bAPs (DISTAL_BAP_LIT.md), so the same rule could not work on it.
- **What we miss, in order:** (1) distal bAP and burst Ca electrogenesis (emodel); (2) an LTD arm that integrates long,
  intermediate post signals (the duration hypothesis), plus a high-frequency LTP arm that wins over eCB-LTD instead of a
  veto; (3) no cap on LTP size; (4) protocol details that differ between papers (age, Ca_o, pulse shape, readout window).

## 1. Ebner's model (stated, from the paper's STAR Methods and the ModelDB code)
**Rule (syn_4p.mod).** It has four pathways and weight w = w_pre * w_post. The input u is the local membrane voltage `u = v`.
- **pre-LTD (eCB-like).** u is rectified above -60 mV and low-pass filtered (tau 10 ms) to give T. A delta pulse at each
  pre spike reads T: LTD_pre = A * D * T. This is a pure post-before-pre detector with no veto. Later bAPs do not
  cancel it; at high frequency it is simply outweighed by the LTP arms.
- **pre-LTP (NO-like).** u is rectified above -30 mV (L-type-like) and filtered twice (tau 7.5 ms, then 30 ms), with
  saturations. The product N_alpha*N_beta minus a threshold of 0.2 is multiplied by a presynaptic trace Z (rise 1 ms,
  decay 15 ms). This arm is frequency dependent through summation. Blocking it "substantially reduced the amount of LTP"
  at 40/50 Hz (Results, Fig 2).
- **post-LTD (NMDA/phosphatase-like).** C = G * [u + 68]+, where G is a glutamate trace (rise 2 ms, decay 50 ms). P is a
  quadratic bump of C between theta_C- = 15 and theta_C+ = 35, normalised to 1. The amplitude is kept small, so brief
  APs give little LTD while long intermediate plateaus (dendritic spikes) accumulate it. The authors call this the
  "duration hypothesis" (Discussion).
- **post-LTP (CaMKII-like).** C above theta_C+ drives K_alpha. It is filtered (tau 15 ms), then filtered again (tau 20 ms), with
  a negative-feedback limiter. LTP = K_alpha*K_beta*K_gamma, so it is frequency dependent.
- **Bounds.** w_pre is in [0,1] (init 0.5) and w_post in [0,5] (init 2). One synapse can therefore grow up to 5x. There is no
  bistability and no population cap.

**Parameters.** 23 shape parameters (time constants, thresholds, slopes) are shared, plus 4 amplitudes per set
(Fig2B.hoc / Fig3.hoc / Fig4B.hoc):

| amplitude | set 1 (Sjostrom 2001) | set 2 (Nevian 2006) | set 3 (Letzkus 2006) | range |
|---|---|---|---|---|
| A_LTD_pre | 3.0e-3 | 2.8e-3 | 1.5e-3 | 2x |
| A_LTP_pre | 3.3e-3 | 1.3e-3 | 2.5e-4 | 13x |
| A_LTD_post | 3.6e-4 | 3.6e-4 | 7.5e-4 | 2x |
| A_LTP_post | 0.20 | 0.57 | 0.078 | 7x |

That makes about 35 hand-tuned values. All three sets also share the same synaptic weight (3.5 nS, AMPA:NMDA 0.5:0.5).
Fitting was a "manual search in a two-stage process": first shape parameters, to match all datasets qualitatively, then
amplitudes per protocol. No error measure or chi2 is reported. A sensitivity analysis says the model is "relatively
sensitive" to thresholds (STAR Methods).

**Neuron.** Hay et al. 2011 L5b, biophys4 (axonal AP initiation), NEURON 7.4, dt 0.025 ms, no inhibition.
- **Apical dendrites:** uniform NaTs2_t 21.5 mS/cm2, low SKv3_1 (1.8 mS/cm2), no A-type K, Ca_HVA uniform, Ca_LVA 100x at
  700-900 um, cm 2.
- **Basal dendrites:** passive plus Ih.
- **bAP targets:** Hay fitted the bAP at 620 um and 800 um (Hay 2011 Table 1, adult P36 Wistar, 33-35 C). Hay notes one
  model variant attenuates more strongly than the data (Fig S6 legend). Which variant that is, is not resolved here.

**Synapse placement.** One synapse, at hand-picked compartments:
- Sjostrom and Letzkus: apic[1] at 90 um ("L5->L5") and apic[50] at 669 um ("L2/3->L5");
- Nevian: dend[80] at 55 um. This is an L5 basal dendrite standing in for L2/3 basal.

**Protocols.**
- **Post spikes:** somatic pulses of 5 ms at 2.7 nA (Sjostrom), 5 ms at 2.1 nA (Nevian), and **2 ms at 5.5 nA (Letzkus, chosen
  "to reliably produce distal dendritic Ca spikes")**.
- **Weight update:** the change is computed from **one sweep and multiplied linearly** by the number of pairings (10/15, 60,
  100), then clipped. Pairing repetition rate therefore has no effect in Ebner.

**Reproduced (stated).**
- Sjostrom 2001 Fig 1D/7B, +-10 ms at 0.1-50 Hz, proximal: LTD below about 30 Hz and LTP above (Fig 2B).
- Distal (669 um): only one data point, S&H 2006 50 Hz +10. Ebner took it from an exponential fit of inset data, not a
  measured mean (ebner_targets.md). Elsewhere distal is a prediction: no LTP, slight LTD above 20 Hz.
- Nevian 2006: all 14 conditions, said to be "faithfully captured" (Fig 3), but on an L5 basal dendrite.
- Letzkus 2006 Fig 5: all four 3AP 200 Hz conditions (proximal +10 LTP, -10 LTD; distal +10 LTD, -10 LTP) (Fig 4B).
- Qualitative: Weber 2016 cluster plasticity (Fig 5).

**Not attempted (inferred from the paper's scope):**
- Markram 1997;
- Letzkus 1AP +10 distance-dependent LTD (0.72), the BAC-pairing conversion, APV/Ni arms;
- Sjostrom 2001 cooperativity and the timing curve;
- Sjostrom 2003/2004/2007 pharmacology (only in-silico pathway knock-outs);
- Zilberter 2009, Sjostrom & Hausser 2006 measured unitary distal values;
- any dt-0 or random-firing data.

No dataset is fitted together with another under one parameter set.

## 2. What actually produces the location dependence, and does it break our rules?

| ingredient (Ebner) | role | our rules |
|---|---|---|
| Hay cell: distal bAP 45 +- 10 mV at 620 um, adult morphology | Distal 1AP/5AP pairings give intermediate post signal -> distal LTD and no LTP (Fig 2C) | Allowed as a new emodel (new name, ion_fitter). Our split2 passes 0/30 on bap_620 (DISTAL_BAP_LIT.md). |
| Burst-evoked distal Ca spike, about 20 ms late, long plateau | Distal 3AP -10 coincides with the plateau -> LTP. Distal +10 sees a decayed glutamate trace -> intermediate -> LTD (Results, Fig 4) | The Ca hot zone (Ca_LVA x100 at 700-900 um) **conflicts with the user's no-hotspot decision** (2026-10-01). A uniform-Ca variant that still spikes at 200 Hz would be allowed. |
| 5.5 nA, 2 ms somatic pulses | Guarantees the Ca spike | Must follow the paper (Letzkus: 2 ms, 3-5 nA). Over-driving would be unanchored. |
| Post-LTD as a band with small amplitude, integrated over time | Long intermediate plateaus give LTD, brief APs do not | **Allowed.** It is the Chindemi depression band (theta_d < c < theta_p with gamma_d). Ebner's per-unit-time LTD is weak relative to LTP; ours sits at the gamma_d bound (IMPROVE_PLAN). |
| G decay 50 ms (glutamate/NMDA-like) times local depolarisation | Timing relative to the late dendritic spike | **Allowed.** Our NMDA-GHK spine Ca already is glutamate x Mg-unblock. That is closer to biology than G x linear voltage. |
| pre-LTD read from 10 ms low-pass of V > -60 at the pre spike, no veto | eCB LTD at post-pre timings; overwhelmed at high frequency | **Allowed if driven by local VDCC Ca**, not V. This replaces our 25 ms bAP veto (which kills 20 Hz LTD) with a "LTP wins" balance, as in Sjostrom 2001's own model 3. |
| pre-LTP from twice-filtered V > -30 x presyn trace | High-frequency LTP, including the 40/50 Hz -10 switch | **Allowed if driven by local high-threshold VDCC Ca** (our NO arm, S2C/3N). Rates must be anchored (Sjostrom 2007 NO block). |
| Voltage as the rule input | all four arms | **Violates** the synapse-local Ca rule; translate every term to a Ca signal. |
| Per-dataset amplitudes | Fits Sjostrom, Nevian and Letzkus separately | **Violates** uniform-across-pathways. Ebner's success is not evidence that one uniform set exists. |
| One synapse at 90/669 um; linear sweep extrapolation; 5x weight range | No population averaging, no saturation within a protocol, no rho0 cap | Not adoptable (we use real contacts and rho0 splits). But it explains why Ebner never meets our LTP cap (section 4, item 4). |

Inferred: with set 3 (A_LTP_post 0.078, A_LTP_pre 2.5e-4), Sjostrom 50 Hz LTP would be far smaller than with set 1.
With set 1, Letzkus distal +10 would probably flip. The Letzkus and Sjostrom tension we see is therefore also present in
Ebner's model; it is hidden by the per-dataset amplitudes. This has not been tested by running their code.

## 3. Experimental conditions in the source papers (stated unless marked)

| | Markram 1997 | Sjostrom 2001 | Letzkus 2006 | Zilberter 2009 | (S&H 2006, ref.) |
|---|---|---|---|---|---|
| pathway, area | L5->L5, S1 | L5->L5, V1 | L2/3->L5, S1 | L2/3->L2/3, V1 | L2/3/L5->L5, V1 |
| age, strain | P14-16 Wistar | P12-21 Long-Evans | **3-6 weeks** Wistar | P14-21 SD | P14-21 |
| temperature | 32-34 C (companion J Physiol) | 32-34 C | 34-35 C | 32-34 C | 32-35 C |
| Ca / Mg (mM) | not in report | **2.5 / 1** | 2 / 1 | 2 / 1 | 2 / 1 |
| GABA-A block | not stated | **no** ("not blocked") | not stated (none listed) | none for pairs; gabazine had no effect on extracellular LTP | never blocked |
| post spike pulse | 5 ms, 1 nA | 5 ms, 0.8-1.5 nA | **2 ms, 3-5 nA** | not extracted | 5 ms, 1.0-1.8 nA |
| induction | 5 APs at 10 Hz, 10 (Fig 2) or 10-15 (Fig 3) bursts every 4 s | 0.1 Hz: 50 pairings; >=10 Hz: 15 x 5 APs every 10 s | **100-200 pairings at 1 Hz** | **40 pairings every 5 s** | 15 x 5 APs at 50 Hz every 10 s |
| readout | **maximum deviation** in a 10-50 min window | mean from 10 min to end (>= 40 min) | mean at 20-30 min | mean from **5 min** to end | mean from 10 min to end |
| washout control | pairing about 13-15 min after break-in | all cells broken in "in quick succession to prevent unequal LTP washout" | 10 min baseline | +15 min diffusion with D890 | quick succession |
| post Vm / Rin | -60 mV | -64.7 mV (no LJP corr.); Rin 142-191 MOhm | -66.3 mV (no LJP corr.); **Rin 20.9 MOhm** | not extracted | - |

The differences most likely to matter for one uniform rule (inferred):
1. **Age and Rin.** Letzkus' L5 cells (3-6 weeks, Rin 21 MOhm, full adult-like apical tree) are not P14 cells (Sjostrom Rin
   about 190 MOhm at P14). bAP propagation and burst Ca electrogenesis mature over this window. Hay's cell is P36. Our
   emodels are juvenile SSCx: they match Markram, Sjostrom and Zilberter, not Letzkus. Part of the "distal bAP 2-13 vs
   45-65 mV" gap is an age mismatch.
2. **Ca_o 2.5 vs 2.0 mM.** Sjostrom 2001 has about 25 % more Ca_o. Through the Ca_o dependence of release and NMDA Ca, that
   raises LTP drive in exactly the targets where our LTP is capped (50 Hz -10, 1.70). Check that the yaml uses 2.5 for
   Sjostrom and 2.0 for Letzkus and Zilberter.
3. **Readout.** Markram's maximum deviation inflates magnitude relative to the window means used elsewhere. Zilberter's
   5-min start includes early potentiation. Our model reads a steady state.
4. **Pairing rate and count.** Letzkus pairs 100-200 times at 1 Hz, Zilberter 40 times at 0.2 Hz, Sjostrom 0.1 Hz bursts.
   Our rule has slow integrators (eCB T, rho), so the repetition rate matters. In Ebner it does not (linear extrapolation).
   ebner_targets.md and PROTOCOLS.md say Letzkus now uses 100 at 1 Hz. Use the paper values everywhere, never Ebner's.
5. **Spike pulse shape.** 5 ms pulses (Markram, Sjostrom, S&H) add a somatic depolarising plateau. Sjostrom 2001 Fig 5D/E
   shows LTP grows sigmoidally with residual depolarisation. 2 ms high-amplitude pulses (Letzkus) do not add one.
6. **GABA-A:** probably not a cause. No paper blocked it, and unitary stimulation recruits little inhibition, so our
   no-inhibition sims are consistent.

Items PROTOCOLS.md does not show, so our implementation probably ignores them (verify in configs/*.yaml):
- per-paper Ca_o (2.5 vs 2.0);
- per-paper post pulse width and amplitude (2 ms vs 5 ms);
- the Letzkus age/emodel mismatch;
- the Markram maximum-deviation readout;
- Zilberter's 5-min start of the readout window.

PROTOCOLS.md records the Letzkus pairing rate and count (100 at 1 Hz) and that inhibition is absent; nothing else.

## 4. Changes we could adopt, ranked by expected payoff

1. **New L5 emodel with a data-like distal bAP** (45 +- 10 mV at 620 um, 36 +- 9 at 800 um; Letzkus Fig 6D 0.45-0.65 of
   the 100 um value) and a 200 Hz burst Ca plateau at 450-700 um. Hay's recipe: uniform apical Na about 21 mS/cm2, low
   apical K, no A-type, cm 2, an older morphology. This fixes the input to Letzkus distal +-10, S&H distal and the distal
   eCB trigger at once. **Needs user decision** (ion_fitter work; new emodel name). The Ca hot zone part **violates** the
   no-hotspot decision unless the user lifts it for this purpose.
2. **eCB arm without a veto.** Read a 10-ms low-pass of the synapse's own VDCC Ca at each pre spike, and let the NO and rho
   LTP arms outweigh it at 20-50 Hz. This is Ebner pre-LTD translated to Ca, and the "LTP wins" model of Sjostrom 2001.
   It would directly address 20 Hz -10 LTD and the 20/40 Hz dt0 sign errors. **Allowed.** The time constant must be anchored,
   e.g. by the Sjostrom 2003 timing data (-25 LTD, -120/-200 none at 0.1 Hz). It competes with the planned v11 veto window
   (IMPROVE_PLAN step 2); compare both in the same round.
3. **High-frequency presynaptic LTP arm driven by summed high-threshold VDCC Ca** (NO arm). In Ebner this is what makes
   40/50 Hz -10 switch to LTP. **Allowed.** It is already in S2C/3N. Anchor its rates on Sjostrom 2007 L-NAME/cPTIO and keep it
   uniform.
4. **Diagnose the LTP ceiling before changing the rule.** Ebner's single synapse can grow 5x. Our population mean is bounded
   by the rho0 split and the potentiated/depressed basis ratio. Compute the maximum achievable ratio for the 50 Hz -10
   records with all rho -> 1 and dpre at its maximum. If it is below 1.70, the cap is structural (rho0/basis), not a gamma
   issue. **Allowed** (diagnostic). Any change to the rho0 file **needs a user decision**.
5. **Per-paper conditions in the protocols:** Ca_o 2.5 for Sjostrom 2001, 2.0 for the others; post pulses 5 ms at paper
   amplitudes vs Letzkus 2 ms at 3-5 nA; paper pairing counts and rates. **Allowed** (anchored in the papers).
6. **Post-LTD band with weak gamma_d plus seeding a00 > 1** (DISTAL_BAP_LIT / IMPROVE_PLAN step 5), so that long
   intermediate plateaus give LTD and brief ones do not. **Allowed** (Chindemi parameters only). It pays off only after item 1.
7. **Letzkus on an age-matched emodel** (3-6 weeks), same rule. This follows the "match the preparation" practice.
   **Needs user decision**: it changes the emodel for one pathway, not the rule.
8. **Markram readout:** score Markram as a window mean if the time course can be digitised (Fig 3C), or note that its
   magnitude is a maximum deviation. **Needs user decision** (Markram +-10 is a must-hit target).
9. **Not to adopt:**
   - per-dataset amplitudes (**violates** uniform rule);
   - local membrane voltage as the rule input (**violates** synapse-local Ca rule);
   - hand-placed single synapses and linear sweep extrapolation (conflicts with our real-contact, rho0-based design);
   - over-driven 5.5 nA pulses (unanchored).
10. **Zilberter CB1-independent, mGluR-dependent LTD:** Ebner offers nothing here (their post-LTD is NMDA-like and
    Zilberter's LTD is NMDAR-independent). Keep the MVD arm (IMPROVE_PLAN step 4). Zilberter states that the unitary L2/3
    contacts are proximal (36.5 +- 5.4 um), so location is not the explanation.

## Open points
- No figure was digitised: Ebner's model values per condition (Figs 2B, 3, 4B) are known only qualitatively. Running
  ModelDB 251493 in a small Slurm job (1 CPU) would give their numbers, and would test the inferred claim that set 1 fails
  Letzkus and set 3 fails Sjostrom.
- Ca_o, pulse shapes and readout windows in our yaml were not checked (outside the brief); verify before acting on item 5.
- Markram 1997 bath Ca/Mg is not in the Science report. The companion paper (Markram et al. 1997 J Physiol) gives 32-34 C
  and -60 mV; its ACSF was not extracted.
- Zilberter's post-spike pulse parameters were not extracted.

## References (PubMed DOIs)
- Ebner, Clopath, Jedlicka & Cuntz 2019, Cell Rep 29:4295, doi 10.1016/j.celrep.2019.11.068 (PMC6941234; ModelDB 251493).
- Hay, Hill, Schurmann, Markram & Segev 2011, PLoS Comput Biol 7:e1002107, doi 10.1371/journal.pcbi.1002107 (Table 1 bAP targets).
- Markram, Lubke, Frotscher & Sakmann 1997, Science 275:213, doi 10.1126/science.275.5297.213 (notes 6, 9; Figs 2, 3).
- Markram, Lubke, Frotscher, Roth & Sakmann 1997, J Physiol 500:409, doi 10.1113/jphysiol.1997.sp022031 (abstract: P14-16, 32-34 C).
- Sjostrom, Turrigiano & Nelson 2001, Neuron 32:1149, doi 10.1016/s0896-6273(01)00542-6 (Experimental Procedures pp. 1161-1162).
- Letzkus, Kampa & Stuart 2006, J Neurosci 26:10420, doi 10.1523/JNEUROSCI.2650-06.2006 (Methods).
- Zilberter et al. 2009, Cereb Cortex 19:2308, doi 10.1093/cercor/bhn247 (Materials and Methods).
- Sjostrom & Hausser 2006, Neuron 51:227, doi 10.1016/j.neuron.2006.06.017 (Experimental Procedures).
- Clopath, Busing, Vasilaki & Gerstner 2010, Nat Neurosci 13:344, doi 10.1038/nn.2479 (base rule; cited, not re-read).
