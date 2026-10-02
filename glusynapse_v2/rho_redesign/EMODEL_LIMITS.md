# Emodel limits that stop GluSynapse_v2 from fitting the plasticity data (2026-10-02)

**Question.** Which deficits of the delta-split2 emodels (L5TPC, L2TPC, L3TPC) block the plasticity fit, how much of the misfit does each explain, which misses would remain with a perfect emodel, and what emodel changes (new names) would fix them?

**Method.** Research only. I read the three hocs and `delta-split2_values.json` (grep/sed), the RS docs, and the glusynapse_v2 and cexp docs named in the brief, and checked the literature on PubMed. Nothing was run or submitted. "Stated" means the paper or project file says it. "Inferred" is my reading.

## Short answer
- **Three emodel deficits matter.** Ranked:
  1. **Almost no distal apical bAP in the L5 TTPC** (D1). This directly blocks Letzkus 1AP +10 and 3AP +10 distal, and it forces the Markram vs L2/3-LTD trade-off.
  2. **bAP Ca too low beyond ~100 um in L5 basal and oblique dendrites** (D2). It closes the shaft licence at about 30 % of L5 potentiation candidates, which partly explains the LTP cap at about 1.3-1.4.
  3. **The L2/3 cells carry unfitted L5 dendritic parameters** (D3). Train Ca is unphysiological (shaft up to ~16 uM), which matters for the Zilberter train rows.
- **Burst fidelity** (D4) and **kBK** (D6) are second order. **Ih** (D7) and the **spine NMDA/VDCC balance** (D5) are not emodel blockers.
- **Is there enough local Ca for LTP of 1.7?** In magnitude, yes. With the licence removed (S1A), the same Ca records give 1.59-1.82 at 40-50 Hz. The cap comes from the rule: the licence, plus the knife-edge between 20 Hz -10 and 50 Hz -10. Real dendrites add no supralinear Ca between 20 and 50 Hz, because the critical frequency is ~100 Hz. The emodel share of the cap is the set of synapses that never see bAP Ca (D2).
- **Rule limits that a perfect emodel would not fix:**
  - the eCB timing rows (20/40 Hz dt 0, 20 Hz +/-25);
  - Zilberter +5-first vs +4-last;
  - Zilberter CB1-independent train LTD under AM251;
  - Zilberter 1AP +/-10 against the L5 0.1 Hz nulls;
  - 50 Hz -10 vs 20 Hz -10;
  - Sj01 vs Markram 10 Hz -10.

## 1. Deficits

### D1. L5 distal apical bAP and Ca (direct)
**Model.**
- Apical NaTg is uniform 0.00998 S/cm2 with vshiftm 3.59 mV. Ka_kampa is uniform 0.00488, kBK 0.005, cm 2 (`cADpyr_L5TPC.hoc` 279-291). The decaying Na gradient is commented out (l. 348), and `ka_gradient` is null in the json.
- bAP amplitude at real L2/3->L5 synapses: 74 / 13 / 2 mV at <200 / 200-450 / >450 um. The og-delta trunk gives 8.9 mV at 450-700 um.
- split2 fails bap_620 in 30/30 cells. The best feasible value is ~8 mV.
- Distal/proximal single-bAP spine Ca: 3e-4 to 3e-3.
- 3AP200/1AP falls with distance (2.5 / 1.6 / 1.2 at 450 / 650 / 900 um), and only 1 of 8 cells makes a tuft plateau.
- Sources: DISTAL_BAP_LIT table; HARD_TARGETS l. 134-147.

**Data.**
- Letzkus 2006 Fig 6D (J Neurosci, doi:10.1523/JNEUROSCI.2650-06.2006): bAP ~100 / 80 / 70 / 60-65 / 50-55 / 45 mV at 100-650 um.
- Fig 6C: control/Ni burst integral ~1.5 at 450-660 um, rising with distance.
- Active backpropagation along the trunk: Stuart & Sakmann 1994 (Nature, doi:10.1038/367069a0).
- Burst Ca above a ~100 Hz critical frequency: Larkum, Kaiser & Sakmann 1999 (PNAS, doi:10.1073/pnas.96.25.14600).
- BAC firing: Larkum, Zhu & Sakmann 1999 (Nature, doi:10.1038/18686).
- Inferred distal/proximal Ca ratio: ~0.1-0.3. So the model is 2-3 orders of magnitude low, not 6x.

**Targets blocked.**
- Direct: Letzkus 1AP +10 (0.72 +- 0.03) and 3AP +10 distal (0.79). No threshold can separate "1AP +10 distal" from EPSP-alone (nopost 0.98) while s_d ~ 0 (DISTAL_BAP_LIT, numbers section, item 5).
- The models that hit 1AP over-depress 3AP distal (0.48-0.53), because their LTD comes from the basal or near "distal" synapses.
- The dropped 3AP -10 distal LTP (1.42) needs the tuft plateau.
- S&H distal eCB: its trigger fires in 42/2688 L2/3->L5 rows.

**Indirect effect.** With s_d ~ 0, the only way to get distal LTD is a global open band (a10 - a00 > 0). That band removes L5 LTP (r1E: Markram +10 0.96). A data-like distal bAP would let 1 < a00 < 1 + s_d give LTD only where c_post is small (S1C_L23_FAILURE 5.2). This is inferred. It would relax part of the Markram vs L2/3->L5 axis (trade-off index -0.90 with Letzkus 3AP +10 proximal).

**Mismatch with the project's Ka plan.** The approved "Hoffman 1997" Ka gradient is a CA1 result. In L5 the trunk data show:
- I_A rises only weakly with distance, +2.3 pA per 100 um, while I_K falls (Bekkers 2000, J Physiol, doi:10.1111/j.1469-7793.2000.t01-2-00611.x);
- total K density decreases distally (Korngreen & Sakmann 2000, J Physiol, doi:10.1111/j.1469-7793.2000.00621.x).

A steep increasing Ka gradient would make D1 worse. The apical Ka of og-delta (= today's L2/3 apical: 0.0228, equal to its Na 0.0229) also failed distally (8.9 mV at 450-700 um).

### D2. L5 basal (>~100 um) and oblique bAP Ca too low (partial: LTP cap; L5 single pairings)
**Model.**
- Basal Ca_HVA2 = Ca_LVAst = 0.0023897659621395684 and kBK 0.0427 (l. 320-322). These are og-delta values, never refitted: json `og_delta_basal` (l. 39) equals split2 `basal` (l. 21) except NaTg 0.003 -> 0.0143 and Ka 0.002 -> 0.025.
- After the ljp25 spine VDCC (E1), spine bAP Ca passes Sabatini at <60 um (1.22-1.31 uM against 1.7 +- 0.6; SABATINI_VALID) and Koester at <80 um (SPINE_CA_DIAG 3).
- Beyond 100 um the bAP peaks at -24 to -54 mV, and spine bAP/EPSP falls to 0.02-0.5.
- The basal >130/<130 um ratio is 0.32. Kampa & Stuart 2006 measured 0.95 +- 0.20 (J Neurosci, doi:10.1523/JNEUROSCI.3062-05.2006).
- Single-bAP shaft dCa at L5->L5 synapses: median 0.062 uM, proximal 0.53 uM (S1C_L23_FAILURE table). Physiological proximal dendritic transients are 150-300 nM per AP (Helmchen 1996, Biophys J, doi:10.1016/S0006-3495(96)79653-4).
- L5->L5 apical synapses (48 of 191) have cpost/cpre 0.002 (CEXP). L5 oblique spine data give ~1 uM per AP at ~100 um (Cornelisse 2007, doi:10.1371/journal.pone.0001073).

**Targets blocked.**
- Under the S1C shaft licence (0.15 uM), 509 of 1725 L5->L5 potentiation-candidate rows are blocked (S1C_L23_FAILURE table), so higher frequency recruits no new synapses (S1C_L5_FAILURE 3c).
- Inferred upper bound: licensing those rows adds up to ~0.42x the present potentiation gain. That would lift 40/50 Hz +10 from 1.37-1.40 towards ~1.52-1.57 (data 1.53 / 1.57) and the S07 pair towards its licence-off value of 1.57 (data 1.62). Not all blocked rows are distal, so this is an upper bound.
- No effect on 50 Hz -10: it stays at ~1.2 with the licence off (HARD_TARGETS l. 98).
- The eCB trigger and veto read own VDCC events, so the C2/C5 profile also sets which L5 synapses can carry eCB LTD.

**Remedy.** Untie basal HVA from LVA and fit them, with a low-threshold (Ni-sensitive) share as in Kampa, Letzkus & Stuart 2006 (burst 4.6x, 3.5x in Ni; J Physiol, doi:10.1113/jphysiol.2006.111062). Leave the bAP voltage alone: the model's 48 mV at 100-150 um already matches Nevian 2007's ~37 % at 140 um (doi:10.1038/nn1826; quoted via Acker & Antic 2009; BAP_CA_DISTANCE_LIT).

### D3. L2/3 dendrites are unfitted L5 og-delta values (partial: Zilberter train rows)
**Model.**
- L2TPC and L3TPC apical NaTg 0.0229, Ka 0.0228 and kBK 0.00356 (l. 286-290), and basal NaTg 0.003, Ka 0.002, HVA = LVA = 0.00239, kBK 0.0427 (l. 318-322). These are byte-for-byte `og_delta_apical` / `og_delta_basal` of the L5 json.
- The L2/3 hocs are identical in og-delta, split2 and split2-bk25 (cmp). L2 and L3 differ in their dendrites only in SKv3_1.
- The basal tree is nearly passive (Na 0.003). Rest at L2/3 synapses is -82 mV (e_pas -80.7 / -86.3, l. 276).

**Consequences seen.**
- Single-bAP shaft dCa 0.115 uM (proximal 0.23), but **trains reach ~16 uM** in thin L2/3 dendrites (S1C_L23_FAILURE table; C_MVD_CALIB cue e). Helmchen 1996 shows a plateau linear in AP frequency, built from 150-300 nM per AP with decay <100 ms, so a 50 Hz train should plateau near 1 uM (inferred).
- Spine bAP Ca is 1-2 uM to 40 um but 0.2-0.6 uM beyond 60 um (SPINE_CA_DIAG 3). Koester & Sakmann 2000 report only a "slight" decrease in L2/3 basal spines (doi:10.1111/j.1469-7793.2000.00625.x).
- The model's L2/3 VDCC pools are as large as or larger than L5's (C_MVD_CALIB 2). This is the opposite of what the Zilberter LTD rows need under a Ca-monotone rule (HARD_TARGETS l. 165).

**Data to fit (L2/3, not L5):**
- Waters 2003 (doi:10.1523/JNEUROSCI.23-24-08558.2003): bAP ~75 mV at 150 um, attenuation length 317 um; Ca at 10 % of maximum at 65 % of the soma-pia distance.
- Larkum 2007 (J Neurosci, doi:10.1523/JNEUROSCI.1717-07.2007): little sag (low Ih); critical frequency 130 Hz; dendritic spikes of 2-4 ms half-width.
- Zilberter 2009 (Cereb Cortex, doi:10.1093/cercor/bhn247): L2/3-L2/3 contacts are proximal.

**Targets.**
- Partial for train4 +4 (0.76), 5ap 10 Hz (0.76) and post-only (1.03). Their misses come from the licence being open in every train and from c_post summing over theta_p, and the excess train Ca feeds both.
- None for train -10 AM251, 1AP +/-10 and +5 first: the cue search found no local separation (C_MVD_CALIB 2), and correct L2/3 dendrites at proximal contacts would look like L5 basal.

### D4. Burst/doublet behaviour (partial: Letzkus 3AP rows, by sample)
- split2 L5TTPC fails 200 Hz following in 201 of 300 burst tasks (S1C_L23_FAILURE 3.5, DECISIONS T39). Early workdirs showed 296/300 spikes (L23_PATHWAYS 2).
- Doublets at threshold remove 18/120 L2/3->L5 posts and 2/24 L5 posts (DOUBLET_BIAS 2).
- Letzkus 2006 delivered 3 APs at 200 Hz from 2 ms, 3-5 nA pulses in every cell. Larkum 1999 PNAS drove trains up to 200 Hz.
- Doublets are allowed (user), but missed burst spikes are not physiological. They bias the 3AP rows towards the more excitable cells. This is the soma/AIS (Na recovery, K_Pst, SK_E2), not the dendrites.

### D5. Spine NMDA/VDCC balance (none as an emodel issue)
- The EPSP VDCC share is 0.3-0.6 % proximally (SPINE_CA_DIAG 2). Cortical spine EPSP Ca is "almost entirely" NMDAR (Kovalchuk 2000, doi:10.1523/JNEUROSCI.20-05-01791.2000; Nevian & Sakmann 2004, doi:10.1523/JNEUROSCI.3332-03.2004).
- bAP/EPSP passes Koester at <80 um. The CA1 10-30 % APV-resistant share would need a spine-neck voltage, which is a synapse-side choice (E2).
- The emodel enters only through rest (-78 / -82 mV; 1/B = 109 / 145) and the local EPSP (+9-12 mV). That is relevant only if the papers held cells near -60 to -70 mV (Markram 1997: -60 +- 2). This is a protocol check, not a refit.

### D6. kBK (partial, mostly settled for L5)
- split2 basal kBK 0.0427 is 8.5x apical (0.005). bk25 sets basal to 0.0025 for L5 only, and the user has locked it in.
- BK density is homogeneous from soma to 850 um on the apical dendrite (Benhassine & Berger 2005, doi:10.1111/j.1460-9568.2005.03934.x; no basal data), so the 8.5x basal excess had no anchor.
- L2/3 basal still has 0.0427. With ~16 uM shaft Ca, BK will narrow later bAPs in trains (inferred). The direction of the effect on Zilberter rows is unknown; check it in K1.

### D7. Ih (none to minor)
- L5 Ih uses the Kole-type exponential (length constant 323 um, ~37x at 1000 um; l. 347-350), consistent with the exponential channel gradient of Kole et al. 2006 (doi:10.1523/JNEUROSCI.3664-05.2006). The same law is applied to basal dendrites by path distance; I checked no basal data.
- L2TPC Ih is ~0 (2.8e-8 S/cm2 base) and L3TPC is 7.2e-6. That fits the "little sag" of Larkum 2007.
- Ih could cap the tuft plateau (D1) through the low distal input resistance. Letzkus' rise-time split is unaffected by ZD7288 (supp Fig 3). No current target depends on Ih directly.

## 2. Ranking and the rule limits
1. **D1, distal apical bAP.**
   - Direct on 2 L2/3->L5 rows (Letzkus 1AP SEM 0.03, so z up to +8 when missed).
   - It enables the no-new-parameter band route that would relax the Markram vs L2/3->L5 axis.
   - Highest value, highest cost.
2. **D2, L5 basal/oblique Ca beyond ~100 um.**
   - Partial on the LTP-cap rows (40/50 Hz +10, 40Hz_5ms, S07 pair and its AM251 lane).
   - Cheaper: basal parameters only, and split mode already exists.
3. **D3, L2/3 dendritic fit.** Partial on 3-4 Zilberter rows. It also removes an unanchored pathway difference: today any L2/3 vs L5 Ca contrast is an artefact of L5 parameters on L2/3 morphologies.
4. **D4, 200 Hz fidelity.** Sample bias on the Letzkus 3AP rows.
5. **D6 for L2/3, D5, D7.** Checks only.

**Misses that a perfect emodel would not fix** (HARD_TARGETS, S1C_L5_FAILURE 4):
- 20/40 Hz dt 0 and 20 Hz +/-25. Under any Ca-amplitude rule, f(0) would have to be <0.15 f(+15). They need a timing veto.
- Zilberter +5 first against +4 last, and post-only. These need coincidence-reading potentiation.
- Zilberter train -10 AM251. This is mGluR, CB1-independent; no local cue exists.
- Zilberter 1AP +/-10 against the L5 0.1 Hz and AM251 nulls. The contacts are proximal, so the same local state.
- 50 Hz -10 (1.70) against 20 Hz -10 (0.65) and its AM251 arm. Both protocols integrate the train in c*. Basal and apical critical frequencies (~100 Hz; 130 Hz in L2/3) mean real dendrites supply no 20-vs-50 Hz step, so the emodel must not be tuned to make one.
- Sj01 10 Hz -10 against Markram -10 (data conflict, z^2 ~3.7). The 0.1 Hz NMDAR-block lane.

## 3. Emodel fixes (ion_fitter; all new names; uniform apical Ca, no hotspot; doublets allowed)
**F1. `delta-split3-ap` (L5TPC, from split2-bk25).**
- Free parameters: apical Na and K only. NaTg density and shift, Ka with at most a weak gradient (Bekkers: +2.3 pA/100 um), and the "apical cm 1 on the trunk" lever named by ion_fitter (BAP_REPORT).
- Targets:
  - Letzkus 2006 Fig 6D: bAP 60-65 mV at 400 um and 45 +- 10 at 620 um (ratio 0.45-0.65 to 100 um);
  - Fig 6C: control/Ni ~1.15-1.2 at 200-250 um, ~1.3-1.45 at 300-400 um, ~1.5 at 450-660 um;
  - keep BAC 8/8 and the long plateau (split2 note);
  - validation: single-AP shaft Ca ratio 0.1-0.3 at 400-650 um.
- If uniform channels cannot reach Fig 6D, log it as a known limit, not a reason for a zone. Any new channel model needs a new .mod and SUFFIX (e.g. `NaTg_ap`).
- Cost: an ion_fitter stage refit, plus the full plasticity chain. That is about 308 CPU-h for L5->L5 + L2/3->L5 (EMODEL_PIPELINE table: total ~379 less ~71 for L2/3->L2/3), logged first.

**F2. `delta-split3-bca` (L5TPC; can be merged with F1 into one package to save a pipeline pass).**
- Basal HVA and LVA untied and free (no new channel unless an R-type is needed: then `CaR_dend`, new SUFFIX).
- Targets:
  - Kampa & Stuart 2006: 1-AP >130/<130 um = 0.95 +- 0.20, critical frequency ~100 Hz, VSD 3rd/1st 2.11 +- 0.22;
  - Kampa, Letzkus & Stuart 2006: 200 Hz burst 4.6 +- 0.5x, Ni 3.5 +- 0.6x;
  - Nevian 2007: no basal Ca spikes;
  - the plasticity-side R-VS readout at the real synapse sites, with floor 0.5 at 130-230 um (EMODEL_CONSTRAINTS C1);
  - S&H 2006 Fig 8B hyperpolarisation curve as validation (C2).
- Cost: basal stage only, the same pipeline pass as F1.

**F3. `l23delta-d1` (new dir `SSCx-AAD-l23delta-d1-emodels`, L2TPC + L3TPC).**
- Own apical/basal split fit on L2/3 data: Waters 2003 profile, Larkum 2007 (critical frequency 130 Hz, low Ih), Koester & Sakmann 2000 basal spine ~ shaft with a slight distance decline.
- Shaft train Ca plateau: linear in frequency, about 1 uM at 50 Hz, from Helmchen 1996. This is inferred from 150-300 nM per AP and tau < 100 ms; treat it as a guard, not a fit value.
- Basal kBK set like bk25 only if K1 shows it matters.
- Cost: ion_fitter fit of 2 e-types, plus the L2/3->L2/3 chain, ~71 CPU-h.

**F4. 200 Hz fidelity objective in F1/F2.** 3 APs for 3 pulses (2 ms, 3-5 nA) at 200 Hz in every cell, from Letzkus 2006 Methods. This is a soma/AIS objective; doublets at threshold stay allowed. It adds no extra pipeline pass.

**Not proposed:**
- a distal Ca floor, distance-dependent gca, or a rise-time drive (unanchored);
- a hotspot zone;
- raising LVA for Sjostrom 2004 dLTD (EMODEL_CONSTRAINTS C7);
- emodel tuning for the 20-vs-50 Hz step.

## 4. Simulation checks before any refit (outlines only; nothing submitted)
All of these are 1 CPU on CPU nodes, with runtime overrides of existing hocs as diagnostics (not models). Pilot one pair first, then size the arrays from the pilot's seff.

- **K1. Site audit, split2-bk25 L5 and the current L2/3 hocs, at the real synapse sites of the extracted pairs.**
  - Pairs: 24 L5->L5, 120 L2/3->L5, 120 L2/3->L2/3.
  - Stimuli: 1 AP; 5 AP at 10/20/50 Hz; 3 AP at 200 Hz; 10 AP at 50 Hz (Zilberter train); 1 AP and 5 AP at 50 Hz with -0.34 nA for 200 ms (S&H Fig 8B).
  - Record: v peak and width, shaft cai (cad), and R-VS spine Ca and VDCC charge (offline ODE, ljp25 globals). Report per-AP ratios and distance bins <60 / 60-130 / 130-230 / 230-450 / >450 um.
  - Variants: Ih = 0, L2/3 basal kBK 0.0025, split2 (kBK 0.0427).
  - Confirms D2, D3, D6, D7 and gives the D1 numbers on split2 (DISTAL_BAP_LIT open point 1).
  - Sizing: from the cexp arrays (40 s per pair for 4 conditions, MaxRSS 2.09 GB) -> ~70 s per pair for 7 conditions. CHUNK 12 is ~14 min, so request **2.7G, 0:30** per task. Pilot: 1 pair, 2.7G, 0:15.
- **K2. Letzkus Fig 6 trunk profile on 8-10 L5 posts.**
  - 1 AP and 3 AP at 200 Hz, recorded on the trunk at 100-700 um. Control, LVA = 0, and LVA = HVA = 0 as the Ni bounds.
  - Sizing: from 22064266 (1.0 GB) -> **1.3G, 0:15** pilot. Elapsed not recorded in my files, so take time from the pilot.
- **K3. 200 Hz fidelity scan on all 144 post cells.**
  - 3 x 2 ms pulses at 3/4/5 nA; count spikes; flag doublets separately.
  - Sizing: 2.7G (cexp basis); pilot 0:15 on 12 cells, then one task.
- **K4. Licence-blocked rows by distance (offline, existing records and the s1C_s rescore).**
  - Count blocked-pot rows per distance bin and per protocol.
  - If the kernel can take a per-synapse licence mask: rescore s1C_s with the licence forced open only where data say bAP Ca is flat (basal <250 um, obliques <200 um). This is an emulation of D2's upper bound, not a model.
  - Sizing: CPU count like 22296646 (1:28, 4.2 GB) -> **5.3G, 0:15**. GPU rescore like the L5-all rescores (7:31, 37.6 GB) -> 3g.40gb, **47G, 0:15**.

## 5. Open points
- The cad parameters (shell depth, decay) behind the shaft licence were not checked. The 16 uM L2/3 train value may be partly cad, not channel density.
- No digitised L5 apical single-AP Ca-vs-distance curve exists (Schiller 1995, Larkum 1999 Fig 5).
- The holding potentials of Zilberter and Letzkus have not been extracted (D5).
- **User decisions:**
  - (a) commission F1+F2 (one L5 package) knowing ion_fitter has called the distal bAP a structural limit;
  - (b) commission F3, a separate L2/3 dendritic fit;
  - (c) until then, treat the Letzkus distal +10 rows as emodel-limited validation.
