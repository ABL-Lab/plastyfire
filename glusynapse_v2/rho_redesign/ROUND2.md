# Rho redesign round 2: L2/3->L5 pre-before-post (+10) over-potentiation

Diagnosis job 22120179 (diag_l23_ltp.py, A0g_s3 parameters, no refit). Outputs: results/l23ltp_A0g_s3_{decomp,med,sep,syn,rec}.csv
and figs/fig6_l23_ltp.png. The "full" column reproduces v4_A0g_s3.csv and l23_v4_A0g_s3.csv exactly.

## Conclusion

- **The miss is on the post side.** The pre pathway is silent at every L2/3 +10 target (mean dpre 0.000 to 0.013), so
  the model ratio equals the rho-only ratio.
- **Unweighting the eCB drive is not the fix.** With 1 - b set to 1:
  - The Letzkus targets do not change, because a single VDCC event per pairing never reaches theta_Tg = 3.15.
  - The S&H 50 Hz targets drop by about 0.32.
  - The L5 +10 LTP targets drop by the same amount.
- **What separates the two groups is local VDCC / shaft Ca, not effcai.**
  - L2/3 synapses that potentiate wrongly reach theta_p with almost no VDCC Ca, so the crossing is NMDA-dominated.
  - L5 synapses that potentiate (correctly) get 10-100x more own -ica_VDCC and shaft Ca.
  - Between the two groups of up-flipping synapses, effcai peak / theta_p is identical (AUC 0.52).
- **Letzkus 1AP +10 (0.72) cannot be reached by any potentiation gate.** Removing all potentiation still leaves 0.88,
  because most of these synapses never reach theta_d (median pairing peak / theta_d = 0.37).

## 1. Ratio decomposition (record means)

Column key:
- pre only: rho is held at rho0.
- pot off: rho_f is capped at rho0.
- eCB unw: the t_drive 4 events are unweighted (1 - b -> 1).
- P_up: fraction of rho0 = 0 synapses that end at rho_f >= 0.5.
- P_down: fraction of rho0 = 1 synapses that end at rho_f < 0.5.

| target | data | full | post only (dpre 0) | pre only | pot off | eCB unw | dpre / dpre unw | P_up / P_down |
|---|---|---|---|---|---|---|---|---|
| L2/3 1AP +10 all | 0.72 | 1.041 | 1.041 | 1.043 | 0.883 | 1.041 | 0 / 0 | 0.15 / 0.21 |
| L2/3 3AP +10 distal | 0.79 | 1.134 | 1.134 | 1.032 | 0.900 | 1.134 | 0 / 0 | 0.19 / 0.18 |
| L2/3 3AP +10 proximal (must stay LTP) | 1.30 | 1.439 | 1.439 | 1.043 | 0.874 | 1.439 | 0 / 0 | 0.58 / 0.25 |
| L2/3 S&H 50 Hz +10 distal | 0.86 | 1.331 | 1.319 | 1.049 | 0.944 | **1.017** | 0.013 / -0.21 | 0.35 / 0.14 |
| L2/3 S&H 50 Hz +10 all | 1.06 | 1.431 | 1.418 | 1.051 | 0.966 | **1.082** | 0.011 / -0.22 | 0.47 / 0.12 |
| L5 Markram 10 Hz +10 | 1.20 | 1.168 | 1.168 | 1.013 | 0.886 | **0.865** | 0 / -0.26 | 0.22 / 0.13 |
| L5 Sjostrom 0.1 Hz +10 | 0.97 | 1.093 | 1.093 | 1.013 | 1.011 | 1.093 | 0 / 0 | 0.07 / 0.01 |
| L5 Sjostrom 10 Hz +10 | 1.16 | 1.107 | 1.107 | 1.013 | 0.810 | **0.819** | 0 / -0.26 | 0.26 / 0.20 |
| L5 Sjostrom 20 Hz +10 | 1.31 | 1.312 | 1.312 | 1.013 | 0.916 | **0.969** | 0 / -0.26 | 0.37 / 0.13 |
| L5 Sjostrom 50 Hz +10 | 1.57 | 1.611 | 1.398 | 1.224 | 1.190 | **1.250** | 0.185 / -0.13 | 0.37 / 0.06 |

**eCB-unweighted test, plainly.** Removing the 1 - b weight does not separate the pathways.
- It lowers L2/3 S&H 50 Hz by 0.31-0.35 and every L5 +10 burst target by 0.30-0.36, almost the same amount.
- It does not touch Letzkus 1AP or 3AP at all: the unweighted gate tTu_max is 0 there.
- So the 1 - b weighting is not why L2/3 +10 LTD is missing.

## 2. Local signals (fig 6c-f, l23ltp_A0g_s3_sep.csv)

AUC(L2/3 should-depress > L5 should-potentiate). Groups:
- L2/3: 1AP +10, 3AP +10 distal, S&H 50 Hz distal and all; 1690 synapse rows.
- L5: Markram 10 Hz +10 and Sjostrom 10/20/50 Hz +10; 764 rows.

| signal | all: L2/3 / L5 median, AUC | up-flipping synapses only: L2/3 / L5 median, AUC |
|---|---|---|
| m_p (pairing effcai peak / theta_p) | 0.53 / 0.66, 0.38 | 1.22 / 1.23, **0.52** |
| t_p (s per pairing above theta_p) | 0 / 0.003, 0.46 | 0.143 / 0.168, 0.49 |
| vd_int (own -ica_VDCC integral per pairing) | 3.4e-6 / 2.6e-4, **0.16** | 2.3e-4 / 2.4e-3, **0.20** |
| sh_pk (shaft Ca peak) | 1.5e-6 / 2.0e-4, 0.25 | 5.9e-4 / 8.4e-4, 0.42 |
| n_cev (VDCC events per pairing) | 1 / 5, 0.10 | 4 / 5, 0.19 |
| b_cev (NMDA-bound b at those events) | 0.94 / 0.99, 0.12 | 0.97 / 0.98, 0.13 |
| c_post | 1.2e-4 / 3.2e-4, 0.47 | 0.008 / 0.012, 0.48 |
| section type (L2/3 3, L5 2) | AUC 0.78 | AUC 0.65 |

- **Inside one protocol (3AP +10).** Distal synapses (LTD 0.79) differ from proximal ones (LTP 1.30) by:
  - vd_int: AUC 0.21, medians 2e-6 vs 7e-5.
  - sh_pk: AUC 0.23.
  - m_p: AUC 0.21.
  - c_post: AUC 0.24.
- **n_cev and n_pre** separate mostly by protocol (one spike vs five per pairing), not by synapse.
- **Letzkus 1AP vs L5 0.1 Hz +10** (data 0.72 vs 0.97; both one spike per pairing):
  - Nearly the same m_d (0.37 vs 0.31) and n_cev (1 vs 1).
  - Different vd_int (2e-6 vs 8e-5) and sh_pk (1e-6 vs 1.4e-4).
  - Pairing rate differs (1 Hz x 100 vs 0.1 Hz x 50). That is a protocol difference and is not used.
- **Reading.** A potentiating crossing at an L2/3 synapse is NMDA Ca on an un-boosted bAP: m_p is the same as at L5, but
  the VDCC current is 10x smaller. This is the regime that Sjostrom & Hausser 2006 found gives LTD.

## 3. Literature (PubMed)

- **Sjostrom & Hausser 2006, Neuron 51:227** ([DOI](https://doi.org/10.1016/j.neuron.2006.06.017)).
  - Five AP-EPSP pairings at 50 Hz, +10 ms, give LTP at proximal L5-L5 synapses (1.40) but not at L2/3-L5 (1.06).
    Distal unitary inputs depress.
  - Distal LTP needs bAP boosting: a dendritic depolarisation, or an EPSP above about 1 mV. This gives supralinear
    dendritic Ca with the same threshold as LTP.
  - Unboosted bAP plus small EPSP gives no supralinear Ca and no LTP.
  - Distal 50 Hz LTD without post spikes is blocked by AM251, so it is eCB/CB1-dependent. The paired distal LTD was not
    tested pharmacologically.
  - Hyperpolarisation abolishes proximal LTP without giving LTD.
  - Sign is set by bAP efficacy at the synapse.
- **Letzkus et al. 2006, J Neurosci 26:10420** ([DOI](https://doi.org/10.1523/JNEUROSCI.2650-06.2006)). Abstract via
  PubMed; pharmacology as in DISTAL_LTP.md.
  - Low-frequency 1AP +10 gives a distance-dependent shift to LTD. Proximal LTD turns into LTP only with bursts or BAC
    depolarisation.
  - Distal 3AP +10 LTD is blocked by APV (1.03) but not by NiCl2 (0.81), so it needs NMDARs and not T/R-type VDCCs.
  - Distal -10 LTP needs both NMDARs and NiCl2-sensitive channels.
- **Location dependence (both papers).** Both attribute it to the local bAP / dendritic voltage at the synapse, not to
  the input pathway. Neither reports a presynaptic mechanism for paired +10 LTD.

## 4. Candidates

Each candidate is synapse-local, uniform across pathways, has at most 1 new parameter, and reduces to A0 when off.
They are ranked by L2/3 excess removed and L5 +10 risk.

### C1 (recommended to test first): potentiation needs own VDCC Ca. Parameter theta_V; 0 gives A0.

- **Rule.** V' = -V / tau_E1 + (-ica_VDCC) / i_scale.
  - tau_E1 = 100 ms and i_scale are the existing t_drive 4 constants.
  - V is the synapse's own spine VDCC current, filtered.
  - pot = [effcai > theta_p] * [V > theta_V]. dep is unchanged, i.e. dep * (1 - pot).
  - So a crossing above theta_p that lacks VDCC Ca drives depression instead of potentiation.
- **How C1 differs from candidate D (rho_v4d, L5 chi2 312.8).**
  - *D's gate is binary and fixed.* It opens on any K_ca crossing in the last 69 ms. Each crossing counts 1, whatever
    the current's size. K_ca is 0.2% of the median bAP peak, so it does not tell an attenuated bAP from a full one.
  - *D has no tunable level and nothing to scan.* It blocks potentiation only where there is no event at all, which
    hits L5 crossings that have no VDCC event close by, e.g. NMDA-only Ca build-up during the bursts.
  - *Under D, a gated-out crossing becomes depression.* In C1 that is the C1/C2 choice.
  - *C1 gates on the VDCC charge.* It integrates the amplitude of -ica_VDCC over 100 ms and thresholds it with theta_V.
  - *C1 is continuous in theta_V.* theta_V = 0 is A0 exactly, so the scan can find a level between the L2/3 and L5
    distributions, if one exists.
  - *Both groups have events; the difference is size.* Up-flippers have about 4 vs 5 events per pairing in the two
    groups, but the charge differs 10x. A count gate cannot separate them; an amplitude gate can.
- **Scan.** Job 22120612 (scan_vgate_amp.py; results/vgate_amp_A0g_s3*.csv; figs/fig7_vgate_amp.png).
  - It covers 12 theta_V values from 0 to 700 for C1 and C2, and scores the 29 L5 and 9 L2/3 targets with no refit.
  - _vmax.csv answers whether one theta_V can sit above the VDCC level of the L2/3 up-flippers and below the L5 one.
    It gives V at the A0 theta_p crossings, as quantiles per group.
- **Evidence.**
  - vd_int separates L2/3 from L5 up-flippers (AUC 0.20; medians 2.3e-4 vs 2.4e-3) where m_p does not (0.52).
  - It also separates distal from proximal 3AP +10 inside one protocol (AUC 0.21).
  - It matches the bAP-boosting / supralinear-Ca LTP condition of S&H 2006 and the NiCl2 sensitivity of distal LTP in
    Letzkus.
- **Size.**
  - Upper bound is the "pot off" column: L2/3 excess removed is 49% (1AP), 68% (3AP distal) and 82% (S&H distal); S&H
    all reaches 0.97.
  - The conversion to depression adds some LTD beyond that.
  - Expected under a threshold between the two medians: most of the S&H and 3AP-distal excess.
- **Risks.**
  - Any L5 target whose LTP synapses have low VDCC current: Markram 10 Hz +5 and +10 (already 1.17 vs 1.20), Sjostrom
    10 Hz +10 (already 1.107 vs 1.16), sj07 step-pair and its mglu_block row.
  - L2/3 3AP +10 proximal (1.30) if theta_V sits above its vd_int (7e-5).
  - Conversion to depression may over-depress L5 0.1 Hz -10, which has VDCC Ca from the AP before pre, so it is probably
    safe.
- **Next step.** Use a CPU numba loop in a new file (not the GPU kernel). Scan theta_V on the A0g_s3 parameters without
  a refit, then do one seeded refit with theta_V free.

### C2: the C1 gate without the depression conversion. Same parameter; 0 gives A0.

- **Rule.** Where effcai > theta_p and V < theta_V, rho gets neither drive: dep * (1 - [effcai > theta_p]).
- **Size and risk.**
  - Removes at most the "pot off" amount, so no extra LTD.
  - Safer for the L5 -10 and step targets.
  - Leaves Letzkus 1AP at about 0.9 and 3AP distal at about 0.95.
- **Next step.** Score it in the same scan as C1. It is the fallback if the C1 conversion costs L5 LTD targets.

### C3 (rejected): a partial glutamate weight on the eCB drive, w = 1 - k b with k = 1 as now.

- **Evidence against.** The k = 0 extreme (section 1) cuts L5 +10 LTP by as much as L2/3 S&H, and it never reaches
  Letzkus.
- **Size.** At most about 0.3 at S&H, 0 at Letzkus.
- **Risk.** Every L5 burst +10 target, including Markram 10 Hz +10 (0.865 vs 1.20).
- **Use.** Only as a re-balancing term after C1, if S&H stays high.

## 5. Residual (likely beyond a local rule)

Letzkus 1AP +10 (0.72) is bounded at about 0.88 by any change on the potentiation side, because most of its synapses
never reach theta_d (median m_d 0.37). Lowering theta_d alone would also depress L5 0.1 Hz +10, whose m_d (0.31) is about the same and
whose data is 0.97. C1's conversion only reaches synapses that cross theta_p.

The separating local quantity between these two targets is again VDCC / shaft Ca (about 40x). The remaining gap is
therefore in the same place as distal -10: the waveform from the emodel and morphology.

Keep it as a scored target and accept a residual.
