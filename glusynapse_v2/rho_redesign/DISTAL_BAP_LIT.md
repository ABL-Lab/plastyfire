# Distal bAP Ca at L2/3->L5 synapses: real or an emodel artefact? (agent D, 2026-10-02)

Question: at distal L2/3->L5 synapses (Letzkus distal, >~130 um up to the tuft), the model's bAP Ca drive is ~6x below
L5->L5 and the shaft-Ca licence never opens. Is that physiological? Is distal LTD a low-Ca LTD regime?
Literature only; no compute run. "Stated" = in the paper; "inferred" = ours. PubMed DOIs at the end.

## Short answer
- **The under-drive is a real mismatch, and it is much larger than 6x.** The 6x is a median over all L2/3->L5 synapses
  (uE 17 vs 109; S1C_L23_FAILURE section 1), and a large share of that is window-current pedestal (SPINE_CA_DIAG).
  - At the Letzkus distal sites the delta-family model gives 2-13 mV bAPs. The data give 45-65 mV.
  - Model spine bAP Ca there is 3e-4 to 3e-3 of proximal (cell_audit/BAP_CA_MAP, LETZKUS_CA). The data put it at roughly 0.1-0.3.
- **Distal LTD in Letzkus is a moderate-NMDA-Ca LTD driven by a weak but nonzero bAP.** The authors say so
  themselves (p. 10427). It is NMDAR-dependent, Ni-insensitive, absent with EPSPs alone, and deepens with distance.
  - This is the Chindemi depression band. It needs the bAP to add some NMDA/VDCC Ca to c* at distal sites, which the model does not do.
- **Sjostrom & Hausser distal LTD is a different mechanism (CB1/eCB).** A weak bAP lets it happen, but it is not rho-LTD.
- **Verdict: an emodel fix first, plus a rule-side condition that costs no new parameter.** No rule change can make
  distal 1AP +10 differ from EPSP-alone while the local bAP is ~2 mV.

## Literature vs model

| quantity | literature (stated unless marked) | model (source file) |
|---|---|---|
| bAP amplitude, L5 apical trunk | ~100 mV at 100 um, ~80 at 200, ~70 at 300, ~60-65 at 400, ~50-55 at 600, ~45 at 650 um. Fig 6D read by eye, 17 sites/9 cells, unchanged by Ni (Letzkus 2006) | og-delta trunk: 73 / 50 / 23 / 8.9 / 3.6 mV at 0-100 / 100-250 / 250-450 / 450-700 / 700-1000 um (BAP_CA_MAP). Real L2/3->L5 synapses: 74 / 13 / **2 mV** at <200 / 200-450 / >450 um (LETZKUS_CA) |
| bAP at 620 um (ion_fitter target) | 45 +- 10 mV | split2: **0/30 cells pass** bap_620 (delta_caspike WIDE_REPORT_v5). Best feasible ~8 mV at 0.497 L_ap; stated there as a model-structure limit (ION_CHANNELS_DECISIONS 2026-10-01; BAP_REPORT) |
| distal/proximal bAP voltage ratio (400-650 vs 100 um) | **0.45-0.65** (Letzkus Fig 6D) | 0.03-0.15 |
| single-AP Ca along the L5 trunk | Largest proximally: 128 +- 25 nM (fura-2), 5 APs at 10-12 Hz 710 +- 214 nM. Declines gradually up to the bifurcation; tuft profile depends on frequency (Schiller 1995, abstract). Na APs give little Ca near the bifurcation in vivo (Helmchen 1999). Single bAPs give "very small depolarization" and "little Ca2+ influx" in adult tuft (Larkum & Zhu 2002, discussion). Subcritical trains "failed to make large changes" distally (Larkum 1999, Fig 5B) | spine 1AP Ca 0.74 / 0.071 / 0.0020 / 0.00044 uM at 0-100 / 100-250 / 250-450 / 450-700 um (BAP_CA_MAP, og-delta) |
| distal/proximal single-bAP Ca ratio | No direct L5 number found. L2/3 analogue: 10 % of max at 65 % of the soma-pia distance, where the bAP is ~50 mV (Waters 2003, Fig 2f). 5-AP 50 Hz bAP Ca still visible at 557 and 816 um (Sjostrom & Hausser 2006, Fig 7C). **Inferred: ~0.1-0.3 shaft at 400-650 um, lower in the tuft** | spine **3e-4 to 3e-3** at >250 um. VDCC charge >250 um is 300-1000x below proximal, at leak level (LETZKUS_CA) |
| basal reference (L5->L5 sites) | spine 1.7 +- 0.6 uM, CA1 <150 um (Sabatini 2002). Basal bAP ~37 % of soma at 140 um (Nevian 2007, as quoted by Acker & Antic 2009). Basal Ca flat to ~250 um: >130/<130 um 0.95 +- 0.20 (Kampa & Stuart 2006) | spine 1.50 uM at <60 um after T32. bAP 77 +- 17 mV at <60 um (T31). >130/<130 um ratio 0.32 (DECISIONS: a 3x under-drive) |
| 3AP 200 Hz burst at distal apical sites | Ca electrogenesis: peak depolarisation follows AP 3 at 660 um (Fig 6A). Control/Ni voltage integral ~1.0 at 100-150 um, ~1.5 at 450-660 um (Fig 6C). Above the critical frequency (98 +- 6 Hz), distal Ca "rose to levels comparable to the soma", with the peak 400-700 um from the soma (Larkum 1999). Basal >150 um: 4.6 +- 0.5-fold supralinear (Kampa, Letzkus & Stuart 2006) | 3AP200/1AP is 2.5 / 1.6 / 1.2 at 450 / 650 / 900 um and **falls** with distance. 1 of 8 cells spikes (LETZKUS_CA). Shaft burst median **1.2 nM** vs licence 0.15 uM (S1C_L23_FAILURE) |
| bAP + EPSP supralinearity | L2/3 basal spines: 3AP 50 Hz at +10 is 1.8 +- 0.1-fold (Nevian & Sakmann 2006, Fig 5D). Peak spine Ca linear in AP number (Fig 5E). Distal L5 pairing with small unitary EPSPs (<1 mV): -2 +- 8 % (no supralinearity). With large EPSPs: 100 +- 30 % (Sjostrom & Hausser 2006, text and Fig 7E) | distal c* = EPSP-alone c* (nopost 0.98), so the supralinearity is ~0. L2/3 basal 3AP50/1AP = 1.0 (BAP_CA_MAP) |
| distal EPSP boosts the bAP | 3-4 fold at distal sites with well-timed EPSPs (Stuart & Hausser 2001); needs dendritic depolarisation or large EPSPs (S&H 2006 Fig 7) | not tested |

## Is distal LTD a low-Ca LTD? Paper by paper

**Letzkus 2006 (S1 L2/3->L5, 3-6 wk, 1 Hz): yes, this is the authors' own reading.**
- 1AP +10 gives LTD at every distance, and it is deeper distally: 0.72 +- 0.03, r = -0.69 against rise time (Fig 2D).
- The LTD is associative: EPSP alone 0.98, offset 500 ms 0.99 (Fig 2C).
- Subthreshold dendritic depolarisation (12 mV at 429 um) leaves it unchanged. BAC firing turns proximal LTD into LTP
  but leaves distal LTD (r = -0.95, Fig 3C).
- Distal 3AP +10 LTD is APV-sensitive (1.03) and Ni-insensitive (0.81) (Fig 7A).
- Their explanation (p. 10427): "moderate NMDA receptor activation" because the EPSP arrives before the dendritic
  depolarisation. Single bAPs, "associated with minimal calcium electrogenesis, also caused distance-dependent LTD".
- Their model, a threshold on integrated g_NMDA (strong = LTP, moderate = LTD, weak = none; Fig 8, suppl. Fig 5), is the Chindemi rule class.
- Inferred: "moderate" NMDA activation needs Mg unblock by a 45-65 mV distal bAP. A 2-13 mV bAP gives the EPSP-alone state.

**Sjostrom & Hausser 2006 (V1 L2/3->L5 and L5->L5, P14-21, 50 Hz): a weak bAP plus eCB, not rho-LTD.**
- Distal pairing LTD: 80.2 +- 4.5 %, n = 28 (Fig 3).
- EPSPs alone also give LTD (71 +- 6 %), and AM251 abolishes it (99 +- 3 %) (Fig 5B). Pairing LTD was not drug-tested.
- They propose that failing bAPs "may still depolarize the distal dendrites sufficiently to help evoke LTD" (p. 231).
- Small distal EPSPs give no supralinear Ca.
- So this LTD rides on low post Ca, but its effector is CB1. That matches the existing eCB arm, which at L2/3->L5 fires in 42 of 2688 rows because its VDCC trigger sees no distal bAP.

## Numbers the model needs (proximal = L5->L5 basal/proximal sites, spine bAP Ca ~1.5 uM)
1. **Single-bAP voltage ratio, distal (400-650 um apical) / proximal (100 um): 0.45-0.65**, i.e. 45-65 mV (Letzkus Fig 6D; bap_620 45 +- 10).
2. **Single-bAP Ca ratio, distal / proximal: ~0.1-0.3**, central ~0.2.
   - The shaft range is from the data (Schiller, Waters, Larkum; L5 is not digitised).
   - Inferred for the rule's spine VDCC (ljp25: V1/2 ~ -31 mV, k 9.5, m^2): a bAP peaking at -15 to -25 mV gives m_inf^2 ~0.5-0.7 of a full
     bAP, over a broader waveform. So, with V matched, spine VDCC Ca would come out ~0.2-0.5 by construction.
   - Today the model is 1e-3 to 3e-4: **two to three orders of magnitude too low, not 6x.**
3. **Distal spine bAP Ca ~0.15-0.5 uM per single bAP.** That is below the 0.15 uM shaft licence for single bAPs (the licence stays shut
   at 1AP +10, as LTD needs), but above it for 3AP 200 Hz bursts in most cells (Larkum 1999; Letzkus Fig 6C).
4. **Distal 3AP200/1AP must rise with distance.** Data: Ni share of the voltage integral 0 -> ~33 % at 450-660 um. Model: it falls.
5. **For the band:** at c_post -> 0, the band is (a10 - a00) c_pre. Letzkus needs a00 > 1 (EPSP alone 0.98 stays below theta_d) and
   1AP +10 distal c* = c_pre (1 + s_d) > theta_d, i.e. **a00 - 1 < s_d**. s_d is the distal bAP-NMDA increment, which is ~0 now.
   From data-like voltage it is inferred at ~0.2-0.8 (smaller than the 1.8 proximal L2/3 burst value; small unitary EPSPs, S&H Fig 7E).

## Verdict
- **Emodel fix (ion_fitter): required.** No new target is needed; the existing ones are what fail:
  - (i) apical trunk bAP 60-65 mV at 400 um and 45 +- 10 mV at 620 um, i.e. 0.45-0.65 of the 100 um value (Letzkus Fig 6D; bap_620 already listed);
  - (ii) Letzkus Fig 6C: control/Ni burst integral ~1.5 at 450-660 um;
  - (iii) single-AP shaft Ca ratio 0.1-0.3 at 400-650 um vs 100 um, as validation.

  ion_fitter has recorded (i) as a structural limit: cable bistability, apical cm 2, thin P14 trunks (BAP_REPORT). Its
  "apical cm 1" and spine-cm routes are the levers named there. This is a user/ion_fitter decision.
- **Rule side (no new parameter, Chindemi params only):** seed refits with a00 slightly > 1 and a10 - a00 > 0 (S1C_L23_FAILURE 5.2b: a01 ~ 2,
  a00 ~ 1.2), so the band is open at low c_post. This cannot produce distal LTD while s_d ~ 0. It only pays off after the emodel fix.
- **Do not** add a distal Ca floor, a distance-dependent gca or a rise-time-based drive. Each would be unanchored, and some would be non-local.
- **Until then:** treat the distal Letzkus +10 rows (1AP +10 distal-weighted, 3AP +10 distal) like 3AP -10 distal, i.e. as
  emodel-limited. Report them as validation, or accept their chi2. Count S&H distal as an eCB-trigger question, not a band question.
- **User decision:** (a) ask ion_fitter for an apical-bAP variant (new emodel name), or accept the limit and demote the distal
  Letzkus +10 rows; (b) whether the eCB trigger may read synaptic (NMDA) Ca at distal sites. S&H EPSP-alone LTD is CB1-dependent without any bAP.

## Open points
- The model numbers for the distal apical sites come from og-delta/delta (BAP_CA_MAP, LETZKUS_CA), not from split2. split2 only has the
  bap_620 fail. A cheap R-V/R-VS run at the real L2/3->L5 synapse sites on split2 (1 CPU, size from 22064266: 1.0 GB) would confirm.
- No L5 apical single-bAP Ca-vs-distance curve was digitised. Schiller 1995 (scanned PMC1156647) and Larkum 1999 Fig 5 are the sources to digitise.
- Our "distal" set is closer than Letzkus' (median 390 vs ~420+ um, 13 % basal; LETZKUS_LOCATION). Use the geometric split when scoring.

## References (PubMed DOIs)
Letzkus, Kampa & Stuart 2006 J Neurosci 10.1523/JNEUROSCI.2650-06.2006 (local PDF, Figs 2, 3, 6, 7, 8) ·
Sjostrom & Hausser 2006 Neuron 10.1016/j.neuron.2006.06.017 (local PDF, Figs 3, 5, 7, 8) ·
Nevian & Sakmann 2006 J Neurosci 10.1523/JNEUROSCI.1749-06.2006 (local PDF, Figs 3, 5) ·
Schiller, Helmchen & Sakmann 1995 J Physiol 10.1113/jphysiol.1995.sp020902 · Larkum, Kaiser & Sakmann 1999 PNAS
10.1073/pnas.96.25.14600 · Larkum & Zhu 2002 J Neurosci 10.1523/JNEUROSCI.22-16-06991.2002 · Larkum, Zhu & Sakmann 2001
J Physiol 10.1111/j.1469-7793.2001.0447a.x · Helmchen et al. 1999 Nat Neurosci 10.1038/14788 · Stuart, Schiller & Sakmann
1997 J Physiol 10.1111/j.1469-7793.1997.617ba.x · Stuart & Hausser 2001 Nat Neurosci 10.1038/82910 · Waters et al. 2003
J Neurosci 10.1523/JNEUROSCI.23-24-08558.2003 · Svoboda et al. 1999 Nat Neurosci 10.1038/4569 · Sabatini, Oertner & Svoboda
2002 Neuron 10.1016/S0896-6273(02)00573-1 · Nevian et al. 2007 Nat Neurosci 10.1038/nn1826 · Kampa & Stuart 2006 J Neurosci
10.1523/JNEUROSCI.3062-05.2006 · Kampa, Letzkus & Stuart 2006 J Physiol 10.1113/jphysiol.2006.111062
