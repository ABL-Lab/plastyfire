# Emodel constraints from the plasticity model (handover to the ion-channel session)

Written 2026-10-01 for the og-delta L5 TTPC refit. The refit is in split apical/basal mode with five current
targets: BAC ratio, Kampa 2006 basal Ca, Antic 2003 TTX/4-AP, adult somatic, and critical frequency. These
constraints come next, after those five.

Sources:
- the ion session's files, read-only: delta_caspike/TARGETS.md, DESIGN.md and ION_CHANNELS_DECISIONS.md;
- the plasticity docs: DECISIONS.md, rho_redesign/{DLTD_DIAG, CA_DECODE, ROUND2, DISTAL_LTP}.md;
- the per-target fit results: results/v4_C1Ajn_s5{,_l23}.csv;
- primary papers. FT = full text read; ABS = abstract only; "by eye" = read off the published figure, not
  digitised.

Statistics are given as the paper states them (SD or SEM). Use the SD, or a population range, as the tolerance.

---

## 0. How the plasticity rule reads the emodel (read this first)

GluSynapse v2 (`glusynapse_v2/mod/GluSynapse_v2.mod` and `GluSynapseV4.mod`) is a point process on a dendritic
segment. The **"spine" has no electrical compartment of its own: its voltage is the segment voltage `v`**. The
rule reads the following:

| rule input | computed from | what the emodel controls |
|---|---|---|
| spine Ca `cai_CR` (drives effcai, so rho) | own NMDA Ca + own spine VDCC Ca, in a 0.087 um3 head (gamma 0.04, tau 12 ms) | **local v(t) only**: Mg unblock and VDCC gating |
| spine VDCC current `ica_VDCC` (the C1 potentiation gate c_VDCC; the eCB/T event drive, t_drive 4) | R-type-like HVA: g = 0.0744 nS/um2 x head area (~0.95 um2, ~19 channels of 3.72 pS), m^2 h; m: V1/2 -5.9 mV, k 9.5, tau 1 ms; h: V1/2 -39 mV, k -9.2, tau 27 ms | **local v(t) only** |
| shaft Ca (emodel `cai`) | emodel Ca_HVA2 / Ca_LVAst | not read by the current rule (C1Ajn); diagnostic only |

Consequences for the emodel:
1. **What the emodel delivers to the rule is the local voltage waveform at each synapse site**:
   - bAP peak and width;
   - burst plateaus;
   - depolarisation during steps;
   - local EPSP size, which sets Mg unblock.

   The emodel's own Ca channels matter only through their effect on v. Shaft-Ca data are still useful because in
   real cells shaft and spine bAP Ca track each other: spine/shaft 1.1 in L5 basal spines (Koester & Sakmann
   1998) and about 2.7 in L5 obliques (Cornelisse 2007).
2. The spine VDCC is high-threshold: m^2 at V1/2 -5.9 mV. Spine VDCC Ca therefore depends steeply on the local
   bAP peak above about -30 mV and on the bAP width. An attenuated bAP (peak below about -30 mV) gives almost no
   spine VDCC Ca, and a 250 ms step to -52 mV gives about 1/250 of a bAP's current (DLTD_DIAG).
3. ROUND2 found that L2/3->L5 synapses get 10-100x less own VDCC charge than L5->L5 synapses:
   - all synapses: vd_int medians 3.4e-6 vs 2.6e-4;
   - up-flipping synapses only: 2.3e-4 vs 2.4e-3.

   Measured bAP-Ca attenuation in basal and oblique dendrites out to about 150-230 um is only about 1-3x
   (C1, C2 below). So this deficit is realistic only for sites far out on the apical tree. **Whether it is real
   is the main question these constraints should settle.**

### Readouts used below (all cheap, no synapses in the simulation)

- **R-VS ("virtual spine")**:
  1. Record `v(t)` at the synapse-site segments.
  2. Integrate the GluSynapse VDCC + spine-Ca ODEs offline (GluSynapse_v2.mod lines 381-446, NMDA term set to
     0) to get the peak `cai_CR - min_ca_CR` and ∫(-ica_VDCC)dt.
  3. Report **ratios** (site/proximal reference, condition/control). The absolute spine Ca is set by synapse
     parameters, not by the emodel.

  R-VS is the readout the rule actually sees. Alternative: insert inert GluSynapse instances (no NetCon) and set
  calcium_current_flag = 0, so that VDCC current does not feed back on v.
- **R-SH (shaft dye)**: the peak dye-bound Ca at the segment.
  - OGB-1: the existing `cadogb` (Kd 0.2 uM, coupled), as in dendfit2/stages/kampa.py.
  - Fluo-5F (S&H 2006; Kd about 0.8 uM in cells, Sabatini 2002): low affinity, so peak Δcai with Kd 0.8 uM
    binding is enough.
  - Report ratios, as the papers do.
- **R-V (voltage)**: the local v peak, the width at half amplitude, and the integral above rest.
- Sites:
  - Generic: basal and apical-oblique segments in 25 um path-distance bins, plus trunk points at 100-700 um.
  - Plasticity-specific: the actual synapse segments of the 123 plasticity posts. The plasticity side can export
    (gid, section, x, path distance, pathway) from the edges files on request.

---

## 1. Where the plasticity fit misses (C1Ajn_s5, best joint fit)

| target (data) | model | z | emodel feature it depends on |
|---|---|---|---|
| Letzkus 3AP 200 Hz -10 distal LTP, L2/3->L5 (1.42 +- 0.09, n7) | about 0.9-1.0; **dropped** from the fit (DROPT) | - | tuft/apical Ca plateau after the burst; og-delta has one in 1 of 8 cells (DISTAL_LTP) |
| Letzkus 3AP +10 proximal LTP (1.30 +- 0.10) | 1.16 | -1.4 | bAP and burst amplitude at proximal-to-mid apical sites |
| Sjostrom 20 Hz -10 LTD (0.65 +- 0.09) | 0.92 | +3.0 | bAP amplitude and VDCC event size during 20 Hz trains (eCB drive) |
| Sjostrom 20 Hz -10 mglu_block (1.02 +- 0.07) | 1.20 | +2.6 | the same: post LTP that should be absent |
| Sjostrom 10 Hz -10 LTD (0.57 +- 0.11) | 0.78 | +1.9 | as 20 Hz |
| Sjostrom 50 Hz -10 LTP (1.70 +- 0.19) | 1.22 | -2.5 | frequency-dependent dendritic boosting at 50 Hz |
| Sjostrom 0.1 Hz +10 (0.97 +- 0.04) | 1.06 | +2.3 | single-pairing NMDA/bAP supralinearity |
| Sjostrom 2007 step200ms pair, no_block (1.36 +- 0.07) | 1.24 | -1.7 | bAP boosting by a 200 ms subthreshold step |
| Sjostrom 2004 dLTD, 250 ms subthreshold step (0.69 +- 0.04), not in this fit | 0.89 (C1Aj4) | - | step depolarisation at the synapse; spine VDCC type (synapse side) |

L2/3->L5 now fits within |z| < 1.5 once Letzkus -10 distal is dropped. The L5 misses are mostly in the rule
(eCB and post-side gate; DLTD_DIAG). They depend on the emodel only through the bAP waveform in trains (C5) and
under depolarisation (C3).

**Location mismatch to flag (plasticity side, not an emodel job).** Letzkus 2006 Fig 1F: the uEPSP rise time is
linear in the anatomical contact distance (r = 0.98). Contacts are on proximal obliques, the main apical trunk
and secondary tuft branches. Their proximal/distal split, at 2.7 ms rise time, lies at about 400-450 um (by eye),
so "distal" means trunk or tuft at 400-900 um. Our "distal" L2/3->L5 synapses are basal/oblique, classified only
by rise time. Constraint C4 helps Letzkus -10 only if those synapses are re-placed at more than 400 um on the
apical tree (DISTAL_LTP option 3b).

---

## 2. Constraints

### C1. Single-bAP (and 5-AP 50 Hz) Ca vs distance in basal and apical-oblique dendrites, at synapse sites

**Data.**

| source | prep | value | stats, n |
|---|---|---|---|
| Kampa & Stuart 2006 J Neurosci 26:7424 (FT summary) | rat L5 TTPC, 3-4 wk, 35 C, OGB-1 200 uM, 2 ms somatic pulses | basal dF/F 1 AP: <130 um 18.3 +- 2.5 % (n12); >130 um 17.3 +- 2.7 % (n18): ratio about 0.95 | SEM |
| Krieger, de Kock & Frick 2017 Front Cell Neurosci 11:194 (FT) | rat L5B TT basal, P26-36, 32-36 C, OGB-1 200 uM | 3 APs: proximal (<10 um) 0.71 +- 0.13 (n7); distal 57-168 um 0.80 +- 0.41 at 40 Hz. Significant small decline over 0-125 um (p = 0.001); time x dF/F proximal 131 +- 50 vs distal 98 +- 61, ratio about 0.75 | SD |
| Nevian & Sakmann 2004 J Neurosci 24:1689 (FT) | rat L4 spiny stellate, **P13-15**, 34-36 C, OGB-1 200 uM, spines 30-210 um | spine bAP Ca falls with distance (r = -0.44); **>150 um is 41 % lower than <50 um** (ratio about 0.59). Single bAP 0.44 +- 0.20 dF/F (n36) | SD |
| Koester & Sakmann 1998 PNAS 95:9596 (FT) | rat L5 basal spines 20-80 um, **P13-14**, 32-34 C, Calcium Green-1 100 uM (Kd 190 nM) | spine AP dF/F 1.01 +- 0.31 (n13). **Shaft = 89 +- 36 % of spine** (n11). AP Ca comparable at all spines up to 80 um | SD |
| Cornelisse et al. 2007 PLoS One 2:e1073 (FT) | mouse L5 (visual), P6-15, 33-35 C, OGB-1 33-100 uM, oblique/secondary dendrites about 100 um | zero-dye Δ[Ca] per AP: spine 1.05 uM (95% CI 0.59-7.98), dendrite 0.38 uM (0.23-1.3); kappa_E spine 19, dendrite 62. Spine decay 91 +- 13 ms vs dendrite 201 +- 20 ms (100 uM OGB-1, n22) | CI / SEM |
| Sabatini, Oertner & Svoboda 2002 Neuron 33:439 (FT) | rat CA1, P14-20, 34 C, sites <150 um (basal and oblique) | zero-buffer spine Δ[Ca]AP 1.7 +- 0.6 uM (single-spine loading, n5); population 1.1 uM (0.6-8.2); small dendrite 1.5 +- 0.5 / 0.7 uM; kappa_E about 20; tau 12 ms | SD / CI |
| Grewe, Bonnan & Frick 2010 Front Cell Neurosci 4:13 (in TARGETS 3) | L5A (not TTPC), P14 and P25 | basal single-AP Ca L1/e 92 um (0.55 of basal length); no Ca in the distal 30 % of the tree; P14 attenuates more than P25 | SEM |

**Protocol.**
- Somatic 2 ms pulse, 1 AP (threshold x 1.1).
- Also 5 APs at 50 Hz (2 ms pulses), which matches S&H 2006 and the Sjostrom 50 Hz protocols.
- celsius 34, v_init at rest.

**Readout.**
- Profile: R-VS peak spine Ca (and ∫-ica_VDCC) vs path distance on basal and on oblique trees, normalised to
  the 20-50 um bin.
- Also R-SH dF/F (OGB-1 200 uM) as the cross-check against Kampa/Krieger.
- Score the distal/proximal ratio for 130-230 um vs <50 um:
  - **target 0.75-0.95 (adult: Kampa, Krieger), hard floor about 0.5** (P13-15: Nevian & Sakmann 0.59);
  - R-VS and R-SH should differ by no more than about 1.5x at the same site.

**Constrains:**
- every L2/3->L5 and L5->L5 target through the size of the VDCC events (the C1 gate, eCB events);
- the ROUND2 finding that L2/3 synapse sites get about 75x less VDCC charge.

If the og-delta R-VS ratio at the actual L2/3 synapse sites (<250 um) is far below 0.5, the emodel is outside the
data there.

**Priority 1. Difficulty: low** (the bap stage already exists; this adds the R-VS readout and the oblique sites).

**No synapses needed.**

Notes:
- TARGETS objective 6 already uses the emodel's shaft Δcai ratio. C1 adds the spine-VDCC readout, which is the
  one the rule sees and is much more sensitive to bAP peak, and the oblique tree.
- Basal and apical are fitted separately (split mode), so C1 can be met by basal parameters alone.

### C2. Hyperpolarisation sensitivity of bAP Ca vs distance (bAP safety margin) — Sjostrom & Hausser 2006 Fig 8B

**Data** (FT, PDF Fig 8A-B; values by eye):
- Prep: rat **P14-21**, thick-tufted L5, 32-35 C. Fluo-5F 180 uM + Alexa 594 10-25 uM, dG/R peak in a 10-25 ms
  window. Mean +- SEM. Collected from 5 cells.
- Protocol: APs alone vs APs paired with a somatic hyperpolarisation of -0.34 nA for 200 ms (-23 mV peak
  deflection). APs were evoked by 5 ms, 1.0-1.8 nA somatic steps. The AP number in the imaging runs is not stated
  in the legend; the plasticity runs used 5 APs at 50 Hz (Fig 8C), so simulate 5 at 50 Hz.
- Result: bAP Ca is reduced more the farther the site is from the soma, in basal (n25), apical (n12) and
  apical-oblique (n7) dendrites. The axon (n6) is unaffected (about 100-115 %).
- Fit-curve values (by eye): about 100 % at 0 um, about 80 % at 100 um, about 62 % at 200 um, about 50 % at
  300 um, about 40 % at 400 um, about 32 % at 600 um.
- Points:
  - basal 20-80 um: 85-100 %;
  - basal 100-130 um: 70-90 %;
  - basal 150-170 um: 30-80 %;
  - obliques 100-200 um: 70-90 %;
  - apical 380-600 um: 15-75 % (scattered).
- Example: basal 21 um is unchanged; basal 92 um is clearly reduced.
- Functional link: the same hyperpolarisation abolishes L5->L5 LTP (Fig 8D: 50 Hz +10 control 140 +- 6 %, n34;
  hyperpolarised about 90-100 %, n6).

**Readout.**
- R-SH (Fluo-5F, linear) and R-VS: peak Ca with hyperpolarisation / without, per site, vs distance.
- Score per 50 um bin against the fit curve, tolerance about +-20 % (the point scatter).
- Digitise Fig 8B before using it as an objective; until then, use the curve values above as validation.

**Constrains:**
- How close the basal and oblique bAP is to propagation failure. That margin decides whether a somewhat
  attenuated bAP still opens the spine VDCC: the L2/3 vs L5 vd_int split and the C1 gate at L2/3 synapses.
- The Sjostrom/Markram L5 LTP magnitudes.
- An emodel with too much basal Ka, or too little Na, shows all-or-none failure, i.e. too steep a profile.
- A near-passive tree shows almost no sensitivity, i.e. too flat a profile.

**Priority 1** (age-matched P14-21; directly tests the quantity the rule thresholds). **Difficulty: low-medium.**

**No synapses needed.**

### C3. Depolarisation-induced bAP boosting vs distance — Sjostrom & Hausser 2006 Fig 7A, B, E

**Data** (FT, PDF; values by eye; mean +- SEM; P14-21, 32-35 C, Fluo-5F 180 uM).

- **Protocol:** 0.4 nA subthreshold dendritic current for 200 ms, injected through the dendritic pipette at the
  recording site (dendritic recordings up to 450 um). APs at 50 Hz are delivered during the depolarisation.
- **Fig 7A** (one cell, 301 um):

  | measure | soma | dendrite |
  |---|---|---|
  | AP height | about 100-105 % | about 140 % |
  | AP width | about 100-105 % | about 140 % |
  | AP area | about 100-105 % | about 190 % |

- **Fig 7B** (n = 13 cells; triangles = current injection, circles = synaptic activation). Dendritic AP area
  boost vs distance, exponential fit:

  | distance (um) | 100 | 200 | 300 | 400 | 450-500 |
  |---|---|---|---|---|---|
  | boost (% of control) | about 100-110 | about 125-140 | about 160-200 | about 200-215 | about 300 |

  Points at 250-270 um span 105-220 %.
- **Fig 7E**: supralinear Ca (pairing APs with the 0.4 nA step, 18 positions in 4 cells; or with large EPSPs, 37
  positions in 10 cells), vs distance:
  - about 100-130 % at <300 um;
  - about 150-300 % at 400-500 um;
  - about 300-600 % at 600-800 um;
  - up to about 1000 % at 900 um.

  Text values for the 50 Hz paired supralinearity:
  - APs + dendritic depolarisation: 268 +- 68 % (n19);
  - APs + large EPSPs: 100 +- 30 % (n37);
  - APs + small distal EPSPs (<1 mV): -2 +- 8 % (n9).
- **Thresholds:** the supralinear-Ca threshold is 1.1 +- 0.32 mV somatic EPSP (n5); the LTP threshold is about
  1.0 mV.
- **Related data:**
  - Letzkus 2006 Fig 3: 100 ms injection at 429 +- 14 um. Subthreshold for BAC it gave 12.1 +- 0.9 mV (n6);
    suprathreshold, 24.1 +- 2.5 mV. Already in TARGETS 1g.
  - Williams & Stuart 2000 J Neurosci 20:8238 (ABS; PMC6773172 available): dendritic depolarisation amplifies
    single distal (>600 um) apical APs 4-7 fold.

**Readout.**
- R-V AP area and R-VS spine Ca at apical trunk/oblique sites 100-450 um: 5 APs at 50 Hz with vs without a
  0.4 nA, 200 ms injection at the site.
- Score the area-boost curve (Fig 7B) and the Ca-boost curve (Fig 7E).
- Digitise both before using them as objectives.

**Constrains:**
- S&H 50 Hz +10 distal (0.86 unboosted; S&H also report LTP 163 +- 7 % with dendritic depolarisation at
  247 +- 22 um, n5);
- Sjostrom 2007 step200ms pair (a 200 ms step; no_block z -1.7);
- the Letzkus BAC experiments.

Boosting through Na recruitment and Ka inactivation is the mechanism S&H propose for the distal LTD-to-LTP
switch. It is consistent with the approved distance-dependent apical Ka gradient.

**Priority 2. Difficulty: medium** (needs a dendritic IClamp and per-site runs).

**No synapses needed.**

### C4. Apical burst Ca electrogenesis vs distance and single-bAP amplitude — Letzkus, Kampa & Stuart 2006 Fig 6

**Prep** (FT, local PDF `plastyfire/papers/Letzkus et al. - 2006 ...pdf`): Wistar rats 3-6 weeks, 34-35 C,
simultaneous soma + apical-dendrite whole-cell recordings. Mean +- SEM.

**Protocol:**
- Bursts: 3 APs at 200 Hz from 2 ms, 3-5 nA somatic pulses.
- Pharmacology: NiCl2 100 uM, bath.
- An EPSP was evoked near the dendritic electrode (<30 um away).

**Fig 6A** (traces): at 100 um the first bAP is the largest. At 660 um the first bAP is the smallest and is
followed by dendritic Ca spikes on APs 2-3. Peak dendritic depolarisation comes after the third AP.

**Fig 6C** (15 sites, 9 neurons): local voltage-integral ratio, control / NiCl2, vs distance (by eye):

| distance (um) | control/Ni area | Ca-channel share of the burst integral (1 - 1/ratio) |
|---|---|---|
| 100-150 | about 1.00 | 0 |
| 200-250 | about 1.15-1.20 | about 15 % |
| 300-400 | about 1.30-1.45 | about 25-30 % |
| 450-660 | about 1.35-1.70 (mean about 1.5) | about 33 % |

**Fig 6D** (17 sites, 9 neurons): single-bAP amplitude vs distance, identical with and without Ni (by eye):

| distance (um) | bAP amplitude (mV) |
|---|---|
| 100 | about 100 |
| 200 | about 80-85 |
| 300 | about 70-75 |
| 400 | about 60-65 |
| 500 | about 55-60 |
| 600 | about 50-55 |
| 650 | about 45-48 |

This is consistent with TARGETS 2 (45 +- 10 mV at 620 um).

**Readout.**
- R-V: the integral of (v - v_rest) over the burst window (0-60 ms) at trunk and oblique sites 100-700 um,
  control vs "Ni".
- Ni proxy: set apical gCa_LVAst = gCa_HVA2 = 0 as the upper bound of the block, and LVA alone = 0 as the lower
  bound. 100 uM Ni blocks T-type and R-type but not most HVA, so the target should lie between the two.
- Also the single-bAP peak vs distance.
- Express distance as a fraction of the main-bifurcation distance, since the P14 tree is about 15 % shorter
  (TARGETS 0).

**Constrains:**
- Letzkus 3AP -10 distal LTP (1.42): the burst plateau the EPSP has to coincide with. Distal LTP is APV- and
  Ni-sensitive: 0.97 +- 0.06 and 0.99 +- 0.10; Fig 7B.
- Letzkus 3AP +10 proximal (z -1.4).
- Letzkus 3AP +10 distal LTD, which is Ni-insensitive (0.81 +- 0.01). So the plateau must not remove the +10
  LTD: the EPSP at +10 arrives before the plateau.

Today's model: og-delta gives a tuft plateau in 1 of 8 cells. TARGETS 5 already lists this as optional
objective 8 ("fit after digitising Fig 6C").

**Priority 1 for the Letzkus pathway** (it is the one dropped target), but it only pays off once the distal
L2/3->L5 synapses sit on the apical tree (section 1). **Difficulty: medium-high**: a graded profile has to come
from uniform apical Ca plus a Ka gradient.

**No synapses needed.**

Do not copy Letzkus's own NEURON model: it used gCaT 0.001 and HVA 1.25 pS/um2 beyond 600 um only, a distal
zone that the no-hotspot rule excludes. The graded Fig 6C profile is the data; how the model reaches it is free.

### C5. bAP train frequency dependence at synapse sites (10, 20, 50, 200 Hz)

**Data:**
- Kampa & Stuart 2006 (FT summary), basal: 3-AP Ca is supralinear only distally (>about 130-150 um, a sigmoid
  in distance, n31).
  - Critical frequency about 100 Hz; 5.0 +- 1.2-fold (n4) for supracritical vs 66 Hz trains.
  - VSD 3rd/1st AP ratio 2.11 +- 0.22 distally (Cd-sensitive).
  - Already a current target.
- Kampa, Letzkus & Stuart 2006 J Physiol 574:283 (FT summary; 3-4 wk, 35 C, OGB-1 200 uM, basal >150 um):
  - 200 Hz burst Ca is 4.6 +- 0.5-fold the single AP (n5); Ni 100 uM leaves 3.5 +- 0.6;
  - a 50 Hz burst gives no LTP (97 +- 5 %).
- Krieger 2017 (FT): distal basal 3-AP Ca at 40 Hz is 0.80 +- 0.41; with a 200 Hz component it is 0.95-1.00,
  i.e. about a 20 % boost (SD).
- Williams & Stuart 2000 J Neurosci 20:8238 (ABS, DOI 10.1523/JNEUROSCI.20-22-08238.2000; FT in PMC6773172,
  numbers not yet extracted):
  - physiological spike trains back-propagate 3-4x more effectively than mean-rate trains into the distal
    apical dendrite (>600 um);
  - maximal for 80-300 Hz components;
  - reduced by dendritic hyperpolarisation and by Na- or Ca-channel block.
- Larkum 2001 Fig 8C (last/first bAP in 20 Hz trains): figure-only (TARGETS 6a).
- Cornelisse 2007: at 5 APs and 50 Hz, the spine signal tracks each AP while the dendritic signal builds up
  (n9).

**Readout:**
- R-VS per-AP spine VDCC peak (the event size for the t_drive 4 eCB drive) and R-V per-AP bAP peak, for 5-AP
  trains at 10, 20, 50 and 100 Hz, at basal and oblique sites.
- Score: AP5/AP1 ratios, plus the Kampa 3-AP 200/66 Hz ratio already in the fit.

**Constrains:** the L5 frequency series:
- 10/20 Hz -10 LTD too weak (z +1.9, +3.0);
- 20 Hz -10 mglu_block too high (z +2.6);
- 50 Hz -10 LTP missing (z -2.5).

These are mostly rule misses (DLTD_DIAG, ROUND2). The emodel should only make sure that per-AP bAP and VDCC
events neither collapse nor grow unrealistically within 10-50 Hz trains.

**Priority 2. Difficulty: low** (same runs as C1). The only basal 10-50 Hz train data are Krieger (40 Hz) and
Kampa (66 Hz); extract Williams & Stuart 2000 for the apical side.

**No synapses needed.**

### C6. Spine Ca supralinearity for single EPSP + single bAP pairings (needs one synapse)

**Data:**
- Koester & Sakmann 1998 PNAS 95:9596 (FT).
  - Prep: rat L5 basal spines 20-80 um, **P13-14**, 32-34 C, 1 mM Mg / 2 mM Ca, mean +- SD. Extracellular theta
    pipette within 10 um of the dendrite; somatic EPSP 0.5-12 mV; AP from a 10 ms, 700 pA step.
  - Single responses:

    | stimulus | Calcium Green-1 100 uM | OG-BAPTA-2 500 uM |
    |---|---|---|
    | EPSP | 1.15 +- 0.41 dF/F (n13) | 0.69 +- 0.32 (n11) |
    | AP | 1.01 +- 0.31 (n13) | 0.62 +- 0.25 (n14) |

  - Paired at 50 ms:

    | pairing | Calcium Green-1 | OG-BAPTA-2 |
    |---|---|---|
    | EPSP->AP | 1.72 +- 0.39 | 1.21 +- 0.5 |
    | AP->EPSP | 1.44 +- 0.36 | 0.89 +- 0.40 |
    | **ratio (EPSP->AP)/(AP->EPSP)** | **121 +- 13 % (n14)** | **143 +- 32 % (n8)** |

  - OG-2 decomposition (n8): the AP component preceded by an EPSP is 157 +- 62 % of control; the EPSP component
    preceded by an AP is 72 +- 19 % of control.
  - Pharmacology: APV 100 uM leaves EPSP Ca at 2.9 +- 6.6 % (n8); NBQX 2 uM leaves 76 +- 15 % (n4).
- Nevian & Sakmann 2004 J Neurosci 24:1689 (FT).
  - Prep: L4 spiny stellate spines, **P13-15**, 34-36 C, OGB-1 200 uM, mean +- SD.
  - Nonlinearity factor (paired / linear sum of time-shifted single fits):

    | dt (ms) | +10 | +20 | +50 | -50 |
    |---|---|---|---|---|
    | factor | **1.7 +- 0.2 (n13)** | 1.5 +- 0.3 (n13) | 1.3 +- 0.1 (n5) | 0.86 +- 0.08 (n7) |

  - The bAP component boosted by a preceding EPSP is 2.3 +- 0.5 at +10, decaying with tau 32 ms.
  - The factor is independent of spine distance (r = 0.1).
  - APV abolishes the supralinearity; 0 Mg + APV gives 1.0 +- 0.2.
- Nevian & Sakmann 2006 J Neurosci 26:11001 (FT; L2/3 basal spines, P13-15, OGB-6F 500 uM, mean +- SEM;
  extracellular stimulation, so validation only):
  - nonlinearity 1.8 +- 0.1 at +10 (n39);
  - AP-before-EPSP pairings sum linearly.
- Kampa & Stuart 2006: EPSP + 200 Hz burst gives 1.9 +- 0.2-fold distal (>150 um) vs 1.4 +- 0.1 proximal
  (APV-sensitive).

**Readout.**
- One GluSynapse at a basal site (20-80 um, and 150 um).
- Runs: EPSP alone, bAP alone, then pairings at +10, +20, +50 and -50 ms.
- Spine Ca through the dye model: Calcium Green-1 Kd 0.19 uM at 100 uM, or OGB-1 at 200 uM, applied to
  `cai_CR`.
- Score: the factor paired / (EPSP alone + bAP alone, time-shifted).

**Constrains:**
- the Sjostrom 0.1 Hz +10 over-potentiation (z +2.3);
- Letzkus 1AP +10;
- the timing asymmetry of every single-pairing target.

The emodel part is how much the EPSP boosts the bAP at the site (the 157 % and 2.3x components: dendritic Na/Ka)
and how much the bAP unblocks NMDA. The NMDA and synapse parameters carry the rest. **Hand over only as a joint
synapse + emodel check.**

**Priority 2. Difficulty: medium.** Needs one synapse, but no network.

### C7. Subthreshold 250 ms step (Sjostrom 2004 dLTD)

**Data** (FT, local PDF):
- Prep: Long-Evans rat **P13-19**, visual cortex, 32-34 C. Unitary L5->L5 pairs whose synapses are proximal,
  presumably basal (EPSP 20-80 % rise 1.3 +- 0.14 ms, n38).
- Protocol:
  - a 250 ms somatic step at threshold minus 20 pA: **141 +- 12.4 pA**;
  - Vrest -67.8 +- 0.46 mV (junction potential not corrected); **peak Vm -52.1 +- 7.7 mV** (paper states SE; 7.7
    looks like an SD);
  - Fig 4 example traces plateau at -46 and -50 mV;
  - 50-60 pairings at 0.1-0.2 Hz.
- Result: LTD 70 +- 4 % (n14 bursts), 65 +- 14 % (n3 singles).
- Timing (Fig 4):

  | pre timing relative to the step | plasticity | n |
  |---|---|---|
  | 25-125 ms before | about 104 % | 7 |
  | coincident | about 72 % | 9 |
  | 50-400 ms after | about 68 % | 8 |

- Pharmacology: blocked by AM251 and by ifenprodil.
- **No Ca was measured.** The paper only suggests LVA Ca (citing Markram & Sakmann 1994, PNAS 91:5207, ABS:
  subthreshold EPSPs raise apical dendritic Ca through LVA channels; no numbers retrievable).
- Related:
  - Seong, Behnia & Carter 2014 J Neurophysiol 111:1960 (FT; mouse PFC L5, P21-28, 32-34 C, Fluo-4FF): spine
    uncaging Ca rises from 3.5 +- 0.4 % (-70 mV) to 6.0 +- 0.8 % (-50 mV), mostly through NMDA. A VGCC cocktail
    including mibefradil removes about 22 %, the same at -70 and at -50 mV.
  - Almog & Korngreen 2009 PLoS One 4:e4841 (FT; P12-16 L5 soma, room temperature): **no T-type found**; the
    low-voltage-inactivating component is R-type.

**Readout:** R-V local v at the basal synapse sites during the step: steady depolarisation and time to plateau.
This is nearly passive in basal dendrites, so local v is about the somatic v of -52 mV. Optionally R-SH shaft Ca
(LVA share).

**Constrains:** Sjostrom 2004 dLTD (0.69). The emodel can only fix the depolarisation, which should already be
right if Rin is right (about 110 MOhm implied by 16 mV / 141 pA at P13-19).

The missing eCB drive is a rule/synapse issue:
- the spine VDCC is HVA, so a -52 mV step gives about 1/250 of a bAP's current;
- the event logic counts the whole step as one event (DLTD_DIAG).

**Priority 3 for the emodel** (validation: step response at the soma and at basal sites). **Difficulty: low.**

**No synapses needed.** Do not raise emodel LVA to "fix" dLTD: the rule does not read shaft Ca, and the P12-16
data show no somatic T-type.

### C8. Basal dendritic Na/NMDA spikes and plateaus (validation only)

| source | prep | value |
|---|---|---|
| Schiller, Major, Koester & Schiller 2000 Nature 404:285 (ABS) | L5 basal, clustered input | somatic spike component 5.9 +- 1.5 mV, half-width 64.4 +- 19.8 ms; NMDA >= 80 % of charge; Ca confined to the activated segment |
| Major, Polsky, Denk, Schiller & Tank 2008 J Neurophysiol 99:2584 (ABS) | L5 basal | somatic plateau about 3 mV for distal vs about 23 mV for proximal inputs; high-Ca zone 10-40 um around the input |
| Nevian et al. 2007 Nat Neurosci 10:206 (ABS; closed) | L5 basal, direct patch | EPSP attenuation >30-fold to the soma; distal uEPSP >8 mV locally; Na and NMDA spikes but no Ca spikes over 75 % of basal length; bAPs "significantly attenuated" (per-distance values not obtained) |
| Milojkovic et al. 2004 J Physiol 558:193; 2005 J Neurosci 25:3940 (FT) | mPFC L5, P21-42, 29-34 C, SD | glutamate-evoked basal plateau: soma 17.05 +- 3.69 mV, half-width 361 +- 92 ms (n24); amplitude saturates; duration grows with glutamate |

**Readout:** needs clustered synapses (NMDA): a local EPSP size at 100-200 um > 8 mV, and no Ca spike in basal.

**Constrains:** nothing in the current fit, since no plasticity protocol drives clustered input. It keeps the
cluster-competition aim (project goal) realistic later.

**Priority 3. Difficulty: medium.** Needs synapses. The basal Ca-spike absence is the only part to watch now:
basal split parameters must not create basal Ca spikes, consistent with Nevian 2007.

### C9. Spine VDCC type and number (synapse side, NOT for the ion session)

- Bloodgood & Sabatini 2007 Neuron 53:249 (ABS): spine VDCCs are CaV2.3 (R-type), absent from the shaft. They
  damp uncaging spine Ca and the uEPSP through an SK loop; SNX-482 boosts both about 2x (Giessel & Sabatini 2011,
  PMC3111456).
- Sabatini & Svoboda 2000 Nature 408:589 (ABS): 1-20 VDCCs per spine, Po about 0.5 per AP.
- Yasuda, Sabatini & Svoboda 2003 Nat Neurosci 6:948 (ABS): L-type channels add no measurable spine influx.

Our mod has about 19 channels of R-type-like kinetics: consistent. The absolute zero-buffer spine Δ[Ca] per bAP
is 1.1-1.7 uM (Sabatini 2002; Cornelisse 2007, L5: 1.05 uM). **The plasticity side** should check the model's
`cai_CR` peak per proximal bAP against that value. The emodel sets only the voltage.

---

## 3. Ranked shortlist to hand over first

| rank | constraint | target (tolerance) | age match | cost |
|---|---|---|---|---|
| 1 | **C2** S&H 2006 Fig 8B: hyperpolarisation (-0.34 nA, 200 ms) sensitivity of bAP Ca vs distance, basal + oblique + trunk | fit curve 100 / 80 / 62 / 50 / 40 % at 0 / 100 / 200 / 300 / 400 um (+-20 %); digitise first | P14-21, 32-35 C | 2 runs/cell, no synapses |
| 2 | **C1** bAP Ca distal/proximal at synapse sites with the spine-VDCC readout (R-VS) + shaft dF/F, basal and oblique, 1 AP and 5 AP at 50 Hz | 130-230 / <50 um ratio 0.75-0.95, floor 0.5 | adult + P13-15 | reuses the bap stage |
| 3 | **C4** Letzkus 2006 Fig 6C/D: Ni-sensitive share of the 3-AP 200 Hz voltage integral vs distance + single-bAP amplitude vs distance on the apical trunk | control/Ni about 1.0 (<=150 um), 1.15-1.2 (200-250), 1.3-1.45 (300-400), about 1.5 (450-660); bAP about 100 / 80 / 70 / 60 / 50 mV at 100 / 200 / 300 / 400 / 600 um | 3-6 wk | 2-4 runs/cell |
| 4 | **C3** S&H 2006 Fig 7B/E: AP-area and Ca boost by a 0.4 nA, 200 ms dendritic step, 100-450 um | area about 125-140 % at 200 um, about 200 % at 400 um; Ca supralinearity 268 +- 68 % (n19) | P14-21 | per-site dendritic IClamp |
| 5 | **C5** per-AP bAP / spine-VDCC event size in 10-50 Hz trains at synapse sites | no collapse within 5 APs at 10-50 Hz; Kampa 200/66 Hz already fitted; Krieger 40 Hz about 0.8 vs 200 Hz component about 1.0 | adult | same runs as C1 |

C6 (pairing supralinearity, Koester & Sakmann 1998 121-143 %; Nevian & Sakmann 2004 1.7 at +10) is a joint
synapse + emodel check. Run it after a refit, not as an ion-session objective.

## 4. Check against the ion session's user decisions

- **No apical Ca hotspot (user, 2026-10-01).** None of C1-C5 needs one. C4 is a graded profile and must be met
  with uniform apical Ca. Its rise with distance can come from the falling bAP amplitude, the rising local input
  impedance, and the approved Ka gradient. Letzkus's own zone model must not be copied (see C4). If uniform Ca
  cannot reach C4, record it as a known limit, not as a reason for a zone.
- **Apical ≠ basal allowed.**
  - C1, C2 (basal arm) and C5 are basal-parameter objectives.
  - C2 (apical/oblique arm), C3 and C4 are apical.
  - No constraint requires shared values.
- **Distance-dependent apical Ka gradient approved (Hoffman 1997).** It fits C2 and C3: their distance-graded
  hyperpolarisation sensitivity and depolarisation boosting are the classic signatures of a distal Ka increase.
  C4 needs the gradient not to suppress the 200 Hz burst plateau at 450-660 um.
- **Existing targets.**
  - C1 extends TARGETS objective 6, which uses the same Kampa data; add the R-VS readout and the oblique sites,
    with no duplicate scoring.
  - C4 is TARGETS optional objective 8.
  - C5's 200 Hz part is the current Kampa target.
  - Nothing here conflicts with BAC ratio, Antic, adult somatic or critical frequency.
- **Age.** C2, C3 and C7 (and C6) are P13-21, which matches the P14-16 emodel. C1 mixes adult and P13-15 data;
  C4 is adult (3-6 wk). Weight age-matched data higher where they disagree. P14 basal attenuation is stronger
  (Grewe 2010), so the C1 floor of 0.5 is set for P14.

## 5. Gaps (need full text or digitising)

- Digitise S&H 2006 Fig 7B, 7E and 8B, and Letzkus 2006 Fig 1F and 6C-D. The PDFs are in plastyfire/papers/.
  The values above are by eye.
- Nevian et al. 2007 (basal bAP amplitude per distance, train dependence) is closed; the PDF has to come from a
  browser.
- Williams & Stuart 2000 (PMC6773172) train numbers are not yet extracted.
- No primary data were found for dendritic Ca during a long subthreshold somatic step with a T-type blocker in
  L5 (C7). Markram & Sakmann 1994 is a scan; no numbers.
- No quantitative bAP Ca vs distance was found for L5 TTPC apical obliques, apart from the S&H 2006 Fig 8B
  oblique points (n7).

## References (DOI)

- Sjostrom & Hausser 2006 Neuron 51:227 — 10.1016/j.neuron.2006.06.017 (PMC7616902)
- Letzkus, Kampa & Stuart 2006 J Neurosci 26:10420 — 10.1523/JNEUROSCI.2650-06.2006
- Sjostrom, Turrigiano & Nelson 2004 J Neurophysiol 92:3338 — 10.1152/jn.00376.2004
- Kampa & Stuart 2006 J Neurosci 26:7424 — 10.1523/JNEUROSCI.3062-05.2006
- Kampa, Letzkus & Stuart 2006 J Physiol 574:283 — 10.1113/jphysiol.2006.111062
- Krieger, de Kock & Frick 2017 Front Cell Neurosci 11:194 — 10.3389/fncel.2017.00194
- Koester & Sakmann 1998 PNAS 95:9596 — 10.1073/pnas.95.16.9596
- Nevian & Sakmann 2004 J Neurosci 24:1689 — 10.1523/JNEUROSCI.3332-03.2004
- Nevian & Sakmann 2006 J Neurosci 26:11001 — 10.1523/JNEUROSCI.1749-06.2006
- Sabatini, Oertner & Svoboda 2002 Neuron 33:439 — 10.1016/S0896-6273(02)00573-1
- Cornelisse et al. 2007 PLoS One 2:e1073 — 10.1371/journal.pone.0001073
- Grewe, Bonnan & Frick 2010 Front Cell Neurosci 4:13 — 10.3389/fncel.2010.00013
- Williams & Stuart 2000 J Neurosci 20:8238 — 10.1523/JNEUROSCI.20-22-08238.2000
- Williams & Stuart 1999 J Physiol 521:467 — 10.1111/j.1469-7793.1999.00467.x
- Nevian, Larkum, Polsky & Schiller 2007 Nat Neurosci 10:206 — 10.1038/nn1826
- Schiller, Major, Koester & Schiller 2000 Nature 404:285 — 10.1038/35005094
- Major et al. 2008 J Neurophysiol 99:2584 — 10.1152/jn.00011.2008
- Milojkovic et al. 2004 J Physiol 558:193 — 10.1113/jphysiol.2004.061416; 2005 J Neurosci 25:3940 — 10.1523/JNEUROSCI.5314-04.2005
- Seong, Behnia & Carter 2014 J Neurophysiol 111:1960 — PMC4044337
- Almog & Korngreen 2009 PLoS One 4:e4841 — 10.1371/journal.pone.0004841
- Markram & Sakmann 1994 PNAS 91:5207 — 10.1073/pnas.91.11.5207
- Bloodgood & Sabatini 2007 Neuron 53:249 — 10.1016/j.neuron.2006.12.017
- Sabatini & Svoboda 2000 Nature 408:589 — 10.1038/35046076
- Yasuda, Sabatini & Svoboda 2003 Nat Neurosci 6:948 — 10.1038/nn1112
- Larkum, Zhu & Sakmann 2001 J Physiol 533:447 — 10.1111/j.1469-7793.2001.0447a.x
