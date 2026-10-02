# Anchors for the gpu_v11 terms: MVD (mGluR-gated VDCC depression), eCB veto window, graded A_eCB (2026-10-02)

**Question.** Which measurements fix or bound (1) the MVD threshold(s), reference pool and own-glutamate window,
(2) the veto window (t0, Tv] and its peak threshold, (3) a graded eCB step A_eCB?
**Sources.** Local PDFs, read as text plus figure pages: Zilberter 2009, Sjöström 2003, Nevian & Sakmann 2006.
Full texts from PMC: Bender 2006, Marcaggi 2009, Birtoli & Ulrich 2004. Abstracts from PubMed: Sabatini & Svoboda 2000,
Czarnecki 2007. Sjöström 2001 numbers come from LTD_LIT.md and S1C_L5_FAILURE.md (the publisher page returned 403).
**S** = the paper states it, **I** = my inference, "fig" = read from a figure (about +-0.01 ΔR, +-0.02 for ratios).

## Short answer
1. **The 0.37 is not an MVD threshold.**
   - It is the 10-AP-train Ca under D890 divided by the control train Ca. With D890 the train-LTP protocol gives LTD.
   - What the data show is that LTP loses its VDCC licence. They do not show an upper LTD edge. Under APV the full-Ca
     train protocol gives LTD, and train-LTD (Z15) is LTD at the full train pool.
   - **theta_MVD_hi: none.** If the kernel needs one, set it ≥ the 10-AP pool (≈ 4.8 P1).
   - **theta_MVD_lo: (pre-only pool, 1.0] x P1** (I). A single bAP gives LTD; pre alone does not.
   - The 0.37 belongs to the LTP licence: θ_VDCC lies in (3.0, 4.4] x P1 (I). This agrees with Kampa's (3.0, 4.6] (ANCHORS_V9).
2. **No uniform theta_MVD works in P1_i units.**
   - Zilberter's L2/3 LTD occurs at pools ≤ 1-1.8 P1.
   - The L5 AM251 rows show no ρ-LTD at pools up to ≈ 2.3 P1.
   - Only an absolute reference (θ_VDCC charge units) can separate the two pathways, and only if the model's L2/3
     proximal pool is larger in absolute terms. This is checked by calibration C-MVD (below).
   - Do not use the own control-train pool as the reference: it is a per-protocol quantity, not a rule input.
3. **Own-glutamate window: W = 50 ms**, from mGluR1 deactivation τ 51.3 ± 1.9 ms (Marcaggi 2009). Range 12-200 ms.
   On the τ_glu = 70 ms trace this means b/b_peak > exp(-W/70) = 0.49, so **b > 0.5 is anchored**.
4. **Veto: t0 = 2 ms (1-2.5), Tv = 17 ms (15-19 at the synapse), threshold 0.5 x own single-bAP VDCC peak.**
   - The chi2 is piecewise constant in both t0 and Tv, because the protocol lags are discrete.
   - t0 = 3 may miss the Zilberter +4 ms bAP.
   - Outside the Sjöström targets, the support is Bender 2006 (+5 LTP dominates; the isolated eCB-LTD reaches +25) and
     Sjöström 2003 (ifenprodil leaves 50 Hz +10 unchanged).
5. **A_eCB: no dose curve exists.** Saturation gives a weak lower bound, A_eCB ≥ ~0.05-0.08 (I).

## 1. MVD: Zilberter 2009 (rat V1 L2/3→L2/3 pairs, P14-21, 32-34 °C, 40 pairings at 0.2 Hz)
| measurement | value | where (S) |
|---|---|---|
| D890 200 µM in post pipette: 10-AP 50 Hz train Ca (oblique dendrite, fura-2 100 µM) | 0.37 ± 0.04 of control (n 4) | Fig 5A, p.2313 |
| train-LTP protocol (pre 4 ms before the 10th AP) with D890 | LTD 0.57 ± 0.07 (n 5); time control 1.34 ± 0.11 (n 5) | Fig 5B, p.2314 |
| Ca transient peaks, Fig 5C (fig, scale bar 0.1 ΔR) | 1 bAP ≈ 0.09; 4 ≈ 0.21; 8 ≈ 0.24; 10 ≈ 0.24. 1-bAP decay τ ≈ 250 ms (dye) | Fig 5C (fig) |
| AP number at +4 before the last AP | 4 APs 0.76 ± 0.07 (4); 8 APs 1.15 ± 0.08 (8); 10 APs 1.49 ± 0.12 (11) | Fig 5D, p.2314 |
| post BAPTA 0.01 / 0.05 / 0.25 mM, train-LTP | ≈ 0.97 / ≈ 0.68 / ≈ 0.82 (3-11 per point) | Fig 5E (fig) |
| APV 50 µM | train-LTP → LTD 0.70 ± 0.07 (4); train-LTD 0.73 ± 0.08 (7) | Fig 6D, p.2315 |
| CPCCOEt 25 µM + EGLU 50 µM | train-LTD 1.02 ± 0.05 (4); train-LTP 1.23 ± 0.12 (4, not shown) | Fig 6D, p.2316-7 |
| AM251 | train-LTD 0.73 ± 0.07 (7). **L5→L5 5x5 10 Hz -10: 1.07 ± 0.08 (3)** | Fig 6C, p.2315 |
| single pairs, controls | +10 0.64 ± 0.07 (6); -10 0.56 ± 0.06 (4); pre alone 0.98 ± 0.08 (5); train alone 1.03 ± 0.04 (4) | p.2311-2 |
| authors | "A simple peak Ca2+ concentration threshold model ... does not explain the induction of LTD with the preceding 10 AP train, when [Ca2+]post is high" | Discussion p.2317 |

**Mapping to "VDCC pool / own single-bAP pool P1" (I).**
- The model pool (τ_E1 100 ms) after k bAPs at 50 Hz, assuming no attenuation, is P_k/P1 = (1-e^(-0.2k))/(1-e^(-0.2)).
  That gives **P4 3.04, P8 4.40, P10 4.77**.
- At 10 Hz, 5 bAPs give 1.57; at 20 Hz, 5 bAPs give 2.33.
- The D890 train is then ≈ 0.37 x 4.77 = **1.77 P1**.
- Fura-2 (Kd ≈ 0.2 µM) saturates. Fig 5C gives P10/P1 ≈ 2.6 in ΔR, which puts the D890 train at ≈ **0.96 P1**. The true
  value lies between these, ≈ 1.0-1.8 P1.
- The low-affinity, linear dye in Nevian (OGB-6F) shows peaks linear in AP number (r² > 0.98, Fig 5E). This supports the
  linear model conversion as the upper estimate.

**Edges.**
- **Lower edge.** A single bAP inside the window gives LTD (Z1, Z3), so theta_lo ≤ 1 P1. Pre alone gives nothing (Z4),
  so theta_lo is above the own EPSP-only VDCC pool.
- **The lower edge is not anchored to a point.** Use 0.5 P1 as the structural midpoint, and test 0.2 and 0.9 for sensitivity.
- **There is no upper edge.**
  - Under APV the train-LTP protocol (pool ≈ P10) gives LTD.
  - Train-LTD (Z15, pre 5-12 ms after the last AP, pool ≈ 4.3 P1) gives LTD.
  - Z11 (LTP) and Z15 (LTD) have the same pool and differ only in a bAP after the pre spike. No Ca band can split them;
    NMDAR coincidence LTP does.
  - The D890 switch is LTP failing: the 4 → 8 AP flip puts the LTP licence at **θ_VDCC ∈ (3.04, 4.40] P1**, with
    D890 at 1.77 below it.
- BAPTA 0.25 mM blocks LTD although the time integral of free Ca does not depend on the added buffer (I). So the real
  sensor reads fast Ca. A charge pool cannot show this, and Z25 stays out of the fit.

**Uniformity (I).**
- With a P1_i reference, the L5 AM251 rows require θ_lo > the pools they reach: 0.1 Hz -10 ≈ 1, Markram 10 Hz -10 ≈ 1.6,
  Sj 20 Hz -10 ≈ 2.3 P1. These are 1.01 ± 0.04, 1.07 ± 0.08 and 1.02 ± 0.07 (Sjöström 2003; Zilberter p.2315).
- Zilberter needs θ_lo ≤ 1-1.8. The two ranges do not overlap.
- **C-MVD calibration** (1 CPU, short job, from extracted c_VDCC only, no plasticity fit):
  - Take the absolute pools of Z 1AP ±10, Z train4 and 0.37 x Z P10 against L5 0.1 Hz -10, Markram 10 Hz -10 and Sj 20 Hz -10.
  - Put theta_MVD,abs in the gap.
  - If there is no gap, MVD cannot be uniform with this emodel, and the user must decide.

**Is there an mGluR/VDCC postsynaptic LTD at L5?**
- **At L5→L5 pairs: no.**
  - Every LTD is presynaptic by CV, with φ ≥ 0 in 68/68 experiments (Sjöström 2007 Fig 6).
  - AM251 abolishes it (Sjöström 2003; Zilberter).
  - Sjöström 2004 finds "no evidence for a postsynaptic component".
  - Sjöström 2003 tested mGluR block only on ACEA-LTD (LY341495 ≈ 0.58, n 6, fig), never on tLTD.
- **Nevian & Sakmann 2006 is L2/3.** Its LTD is mGluR-dependent (MCPG 1.06 ± 0.16) but CB1-dependent, so it is
  presynaptic eCB, not MVD.
- **At L2/3→L5 inputs: yes, but with a burst trigger** (extracellular, so validation only). Birtoli & Ulrich 2004
  (3-4 wk, 35 °C, 0.2 Hz):
  - EPSP + 3-4 APs at 150-300 Hz at +10 gives 59 ± 9 % (n 15). Single AP at +10 gives LTP, 123 ± 9 % (n 11).
  - MCPG 0.25 mM 101 ± 6 % (n 6); Ni2+ 100 µM 93 ± 10 % (n 7); APV 48 ± 9 % (n 7). It needs ≥ ~70 pairings (10-50 ineffective).
  - Window -100 to +200 ms.
  - Czarnecki 2007 shows mGluR1 → PLC → PKC → ER Ca → AMPAR endocytosis, i.e. postsynaptic expression.
  - So in L5 the same arm exists with θ above one bAP (≈ 3-4 P1 burst), which fits the L5 AM251 constraint, not Zilberter's single pairings (I).

**Own-glutamate window.**
- Marcaggi 2009 (mGluR1β FRET, PC12 cells, temperature not stated in the text read):
  - activation τ reaches a minimum of 10.1 ± 0.7 ms (n 7, Fig 2G);
  - **deactivation τ 51.3 ± 1.9 ms (n 41)**, independent of concentration (Fig 2H);
  - sensitization τ 363 ± 54 ms (Fig 3F);
  - EC50 for transients decaying in 10 / 1 / 0.1 ms: 60 µM / 262 µM / 1.46 mM (Fig 4H);
  - summation for trains above 20 Hz (Fig 4J).
- Bounds on W:
  - lower ≥ ~12 ms, so that the +10 single bAP falls inside (Z1);
  - upper: Bender's mGluR-eCB LTD window, -100 to +25 ms (~125 ms), and Birtoli's +200 ms.
- Mapping: b_c = exp(-W/τ_glu) x b_peak. W 50 → 0.49; W 12 → 0.84; W 125 → 0.17; W 200 → 0.06.
- Caveat (I): perisynaptic mGluRs respond poorly to one brief transient. Zilberter tested mGluR dependence only with trains.

## 2. eCB veto window (own spine VDCC peak in (t_a+t0, t_a+Tv])
**Data (S).**
- **Sjöström 2001 (lags pre→next post / previous post→pre).**
  - 20 Hz: +10/-40 → LTP 1.31; +25/-25 → LTD 0.65/0.70; 0/±50 → 0.76.
  - 40 Hz: +15/-10 → 1.51; 0/±25 → 0.93.
  - 50 Hz: 0/±20 → 0.92 (validation).
- **Sjöström 2003.** Ifenprodil leaves 50 Hz +10 tLTP unchanged (≈ 1.55, n 7, vs ≈ 1.6, n 16; fig), so the veto is full at +10.
  AM251 raises 200 ms step-pairing LTP to 2.13 ± 0.22 (vs 1.62 ± 0.07, Sjöström 2007), so a sustained depolarisation does not veto.
- **Bender 2006** (rat S1 L4→L2/3, P16-23, 22-24 °C, 100 pairings at 0.2-0.25 Hz):
  - +5 gives 1.27 ± 0.05 (n 20); -25 gives 0.72 ± 0.03 (n 40);
  - isolated t-LTP (AM251) appears only at +5 to +30;
  - **isolated eCB t-LTD (internal MK-801) spans -100 to +25**, so between +5 and +25 both arms run and LTP dominates.
- **Nevian & Sakmann 2006** (L2/3 basal spines, 3 APs at 50 Hz, Fig 5D, fig):
  - nonlinearity +10: 1.8 ± 0.1 (n 39); +50: ≈ 1.2 (significant); -10: ≈ 1.35; -30: ≈ 1.15 (not significant);
  - peak ΔG/R: +10 ≈ 0.18, +50 ≈ 0.095;
  - the excess falls with τ ≈ 29 ms (I).
  - So the Ca coincidence signal has no edge near 17 ms. A hard Tv is a model construct set by the L5 rows.

**Values (I).**
- The scored lags are 0, +4 (Z11), +5 (Z12), +10, +15 (40 Hz -10), +20 (50 Hz 0) and +25 (40 Hz 0, 20 Hz +25).
  So the chi2 is flat between neighbouring lags.
- At the synapse, a lag L arrives at about L - d_syn + d_bAP after t_a, with d_syn ≈ 1-2 ms and d_bAP ≲ 0.5 ms.
- **t0** must exclude the tail of the dt-0 bAP and include the +4/+5 bAP (Bender +5, Zilberter Z11): **(≈1, ≈2.5) ms; use 2**.
  t0 = 3 is borderline for +4.
- **Tv** must include +15 and exclude +20: **[≈14, ≈19) ms at the synapse; 17 is fine.** Check both against the
  model's t_a convention (pre spike vs EPSP onset).
- The upper limit of 20 rests only on the 50 Hz dt-0 validation row. The fitted rows give only < 25 (40 Hz 0, 20 Hz +25)
  and ≥ 15 (40 Hz -10). Bender caps it at ≤ 30 at room temperature.

**Peak threshold.**
- θ_veto = **0.5 x own single-bAP VDCC current peak** (structural; range: above the step/EPSP-only peak ratio, below the
  weakest bAP peak in a 50 Hz train, ≈ 0.6).
- Sjöström 2004 step: steady -I_VDCC ≈ 3.4 K ≈ 0.7 % of the median bAP peak (MECH_NECESSITY §6). The EPSP-only peak must
  be read from the pre-only extraction.
- Weak support: proximal CA1 spines hold 1-20 VSCCs, which open with p ≈ 0.5 per dendritic AP (Sabatini & Svoboda 2000).
  Half the mean bAP influx is therefore a typical single-trial value.

## 3. Graded A_eCB
- **No tLTD-vs-pairing-number curve exists** in Sjöström 2003 (checked: Figs 1, 8, 9 and Methods) or in Bender 2006.
  Bender reports only a correlation of tLTD with ΔPPR over 100-300 pairings (R² 0.30).
- I: tLTD is frequency-independent at 0.1-20 Hz and sits at the agonist ceiling: 0.73 ± 0.03 (n 60) vs ACEA 0.71 ± 0.03
  (n 14); 50 pairings at 0.1 Hz give 0.69 ± 0.07.
- So (1-A)^50 ≤ ~0.1 and **A ≥ 0.045**.
- Burst -120/-200 ≈ 0.79 (n 7, fig) after 15 pairings would need A ≥ ~0.08 if fully gated. This is weak: a partial gate
  confounds it.
- Birtoli's ≥ 70-pairing need is for the postsynaptic arm. It is a validation check on MVD, which runs at γ_d with no new rate.

## Table
| parameter | value | range | source | how mapped |
|---|---|---|---|---|
| theta_MVD_lo (P1_i ref) | 0.5 | (pre-only pool, 1.0] | Zilberter Z1, Z3, Z4 | single bAP in window → LTD; conflicts with L5 AM251 (> 2.3) |
| theta_MVD_lo (absolute) | C-MVD gap | L5 -10 pools < θ ≤ Z pools | Zilberter; Sjöström 2003 AM251 | calibration on extracted c_VDCC |
| theta_MVD_hi | none (∞) | ≥ P10 ≈ 4.8 P1 | Zilberter APV, Z15 | LTD at the full train pool |
| θ_VDCC (LTP licence), D890 use | (3.04, 4.40] P1 | D890 1.77 below | Zilberter Fig 5B,D; Kampa 2006 | 4 vs 8 bAPs, linear 100 ms pool |
| W_glu / b cutoff | 50 ms / 0.49 b_peak | 12-200 ms / 0.84-0.06 | Marcaggi 2009 Fig 2H | exp(-W/70) |
| t0 | 2 ms | 1-2.5 | Sj01 dt 0; Bender +5; Z11 +4 | between the dt-0 and +4 bAP peaks |
| Tv | 17 ms | 14-19 (protocol 15-20) | Sj01 40 Hz -10, 50 Hz 0 | between the +15 and +20 bAPs |
| θ_veto (peak) | 0.5 x P1peak_i | (EPSP/step ratio, ≈ 0.6) | Sj03 ifenprodil; Sj04; Sabatini & Svoboda 2000 | structural |
| A_eCB | not anchored | ≥ 0.05-0.08 | Sjöström 2003 Figs 4, 8B, 9C | saturation over 50 pairings |

## Open points and user decisions
1. **Decision:** MVD with a P1_i reference cannot hold under the uniform rule (section 1). Run C-MVD with the absolute
   reference. If the gap is empty, either drop Zilberter 1AP ±10 / 5x5 10 Hz to validation, or accept about 20 chi2.
2. **Decision:** the kernel's theta_MVD_hi should stay off. The D890 data constrain θ_VDCC, not MVD.
3. Check the t_a timing convention before fixing t0 (2 vs 3 ms).
4. Not read: Sjöström 2001 full text (403), Czarnecki 2007 full text, and the Bender Fig 7 per-delay values (not in the PMC text).

## References (DOIs; via PubMed where noted)
Zilberter 2009 Cereb Cortex 19:2308, 10.1093/cercor/bhn247 · Sjöström 2001 Neuron 32:1149, 10.1016/s0896-6273(01)00542-6 (PubMed) ·
Sjöström 2003 Neuron 39:641, 10.1016/S0896-6273(03)00476-8 · Sjöström 2004 J Neurophysiol, 10.1152/jn.00376.2004 ·
Sjöström 2007 Neuropharmacology 52:176, 10.1016/j.neuropharm.2006.07.021 · Nevian & Sakmann 2006 J Neurosci 26:11001,
10.1523/JNEUROSCI.1749-06.2006 · Bender 2006 J Neurosci 26:4166, 10.1523/JNEUROSCI.0176-06.2006 (PubMed/PMC) ·
Marcaggi 2009 PNAS 106:11388, 10.1073/pnas.0901290106 (PubMed/PMC) · Birtoli & Ulrich 2004 J Neurosci 24:4935,
10.1523/JNEUROSCI.0795-04.2004 (PubMed/PMC) · Czarnecki 2007 J Physiol 578:471, 10.1113/jphysiol.2006.123588 (PubMed) ·
Sabatini & Svoboda 2000 Nature 408:589, 10.1038/35046076 (PubMed) · Kampa 2006 J Physiol, 10.1113/jphysiol.2006.111062.
