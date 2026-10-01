# LTD diagnosis (v5c rule, 2026-10-01)

## Findings
1. **L5->L5 LTD is all eCB in the model (rho-LTD ~0), as LTD_LIT requires.** The shortfall is trigger reliability: at Sj01 0.1 Hz -10
   only 38 % of synapses step (basal 45 %, apical 16 %). d_min -0.29 is right: with every synapse stepping, the Sj01/Sj03 LTD rows
   drop from z^2 36.1 to 5.9. Sj01 10 / 20 Hz -10 are only 1.3 / 0.8 SEM under the 0.72 ceiling, so they need no rho-LTD.
2. **The trigger pool is graded.** At own arrivals 10-17 % of synapses are < 0.01 theta_eCB and 45-58 % are at 0.01-1, so W misses
   theta_eCB by a factor of up to 100 at most synapses. Its size follows bAP-VDCC amplitude, not a basal/apical rule: L5-L5 apical
   contacts fail more, L2/3-L2/3 apical contacts trigger more.
3. **One absolute threshold cannot give the timing.** -10 steps 38 % but -120 still steps 22 % (data 1.05). A burst at -120 steps only 31 %
   (data 0.79). Lowering theta_eCB recruits synapses (LTD rows 39 -> 25) but puts eCB into every +10 LTP train (LTP rows 34 -> 150),
   because the (1 - b) weighting stops being enough once the threshold is low.
4. **A bAP veto (own spine VDCC influx within 25 ms after the own pre spike vetoes that eCB step) is the only veto that works.**
   It brings chi2 at the fitted params from 87.8 to 85.6 (39 targets) and 298.7 to 287.4 (54). L2/3->L5 improves 14.7 -> 8.8 (Letzkus
   3AP -10 proximal). It also removes most of the LTP damage of a low threshold (thE x0.1: LTP rows 77.5 -> 51.8). Vetoes read from c*
   (influx or a theta_p event) never fire within 25 ms: c* is too slow. Tv 15 / 25 / 35 give the same result.
5. **A trigger relative to the synapse's own bAP is the natural fix for points 2-3.** With W > k x (own single-bAP VDCC pool jump) and
   k ~ 0.5, a bAP 10-25 ms before the pre spike triggers at every synapse, and one 120 ms before does not (exp(-1.2) = 0.30 < k).
   It is now in the v7 kernel (ecb_ref) and being refit; it has no offline test yet.
6. **Data conflict, irreducible about 3.7:** Markram 10 Hz -10 (0.79 +- 0.03) and Sj01 10 Hz -10 (0.57 +- 0.11) are the same protocol
   in the model.
7. **Other errors that are not LTD-trigger errors:**
   - Sj01 20 Hz -10 has rho-LTP at 8 % of synapses (mGluR-block lane 1.113 vs 1.02).
   - Sj07 pre-only has rho-LTD from pre spikes alone (0.90 vs 1.00), because theta_d is low relative to c_pre (a00).
8. **L2/3->L2/3 (Zilberter) LTD needs rho-LTD, but its depression band is closed.** theta_d / theta_p = 0.96 there vs 0.65 at L5
   (a01 ~ a11 and c_post dominates). At 1AP -10 the c* peak is 0.78 theta_d, the same as L5 -10 (0.75), which must not depress.
   A uniform c* theta_d move cannot separate the two. It needs a depression band read on spine VDCC Ca (LTD_LIT 8), which is a user decision.
## 1. LTD and no-change targets in three fits
Preds from the fit csvs. S1 = v5_S1_V5c_s6_39 (delta-split1 emodel, chi2 87.8); og = v5_V5c_s7 (og emodel, 84.7);
C1 = v4_C1Ajn_s5 (graded eCB + NO, og emodel, 52.0). z^2 in brackets. Shortfall = S1 pred - data (> 0: too little LTD).

| target (L5-L5 unless marked) | data ± SEM | S1 | og | C1 | S1 shortfall |
|---|---|---|---|---|---|
| Markram 10Hz -10 | 0.79 ± 0.03 | 0.815 (0.8) | 0.817 (0.9) | 0.783 (0.1) | +0.02 |
| Sj01 0.1Hz -10 | 0.69 ± 0.07 | 0.882 (7.5) | 0.848 (5.1) | 0.799 (2.4) | +0.19 |
| Sj01 10Hz -10 | 0.57 ± 0.11 | 0.855 (6.7) | 0.776 (3.5) | 0.779 (3.6) | +0.29 |
| Sj01 20Hz -10 | 0.65 ± 0.09 | 0.940 (10.4) | 0.861 (5.5) | 0.918 (8.9) | +0.29 |
| Sj03 0.1Hz -25 | 0.65 ± 0.09 | 0.887 (6.9) | 0.847 (4.8) | 0.773 (1.9) | +0.24 |
| Sj03 burst r50 -120 | 0.79 ± 0.08 | 0.903 (1.8) | 0.898 (1.6) | 0.774 (0.0) | +0.11 |
| Sj03 burst r50 -200 | 0.79 ± 0.08 | 0.933 (2.8) | 0.918 (2.3) | 0.819 (0.1) | +0.14 |
| L2/3-L5 Sj06 50Hz +10 distal | 0.86 ± 0.09 | 0.845 (0.0) | 0.872 (0.0) | 0.866 (0.0) | -0.02 |
| L2/3-L5 Letzkus 1AP +10 | 0.72 ± 0.03 | 0.751 (1.0) | 0.737 (0.3) | 0.731 (0.1) | +0.03 |
| L2/3-L5 Letzkus 3AP +10 distal | 0.79 ± 0.06 | 0.839 (0.7) | 0.842 (0.8) | 0.790 (0.0) | +0.05 |
| L2/3-L5 Letzkus 3AP -10 proximal | 0.89 ± 0.03 | 0.811 (7.0) | 0.883 (0.1) | 0.860 (1.0) | -0.08 (too much) |
| **sum z^2, LTD controls (11)** | | **45.6** | **24.8** | **18.1** | |
| no change: Sj03 0.1Hz -120 / -200 | 1.05 ± 0.10 | 0.927 / 0.981 (1.5 / 0.5) | 0.918 / 0.938 (1.8 / 1.3) | 1.011 / 1.013 (0.2 / 0.1) | over-LTD |
| no change: Sj01 0.1Hz +10 | 0.97 ± 0.04 | 1.028 (2.1) | 1.024 (1.8) | 1.061 (5.2) | |
| no change: Sj07 pre only / post only | 1.00 ± 0.07 | 0.897 / 1.034 (2.1 / 0.2) | 0.656 / 1.036 (23.6 / 0.3) | 0.903 / 1.055 (1.9 / 0.6) | pre-only LTD |
| no change: Letzkus no post / -500 | 0.98 / 0.99 | 1.043 / 1.054 (1.1 / 1.6) | 0.991 / 1.000 (0.0 / 0.0) | 1.043 / 1.034 (1.1 / 0.8) | |
| blocks: Sj 0.1Hz -10 mGluR / NMDAR | 1.01 / 1.07 | 1.013 / 1.013 (0.0 / 2.0) | 0.978 / 1.013 | 1.019 / 1.013 | |
| blocks: Sj 20Hz -10 mGluR / NMDAR | 1.02 / 1.02 | 1.113 / 1.013 (1.8 / 0.0) | 1.017 / 1.013 | 1.202 / 1.013 | |
| blocks: Markram 10Hz -10 mGluR | 1.07 ± 0.08 | 0.948 (2.3) | 0.964 (1.7) | 1.023 (0.3) | |

L2/3-L2/3 (Zilberter), split1 54-target fit v5_S1_V5cz_s5_39 (chi2_l23l23 111.4):

| target | data ± SEM | pred (z^2) | shortfall |
|---|---|---|---|
| 1AP -10 | 0.56 ± 0.06 | 0.895 (31.1) | +0.33 |
| 1AP +10 | 0.64 ± 0.07 | 0.789 (4.5) | +0.15 |
| 5AP 10Hz +10 | 0.76 ± 0.07 | 0.805 (0.4) | +0.05 |
| 5AP 20Hz -10 | 0.93 ± 0.07 | 0.869 (0.8) | -0.06 |
| train10 50Hz -10 last | 0.72 ± 0.05 | 0.823 (4.2) | +0.10 |
| train10 50Hz -10 last, AM251 (mglu_block lane) | 0.73 ± 0.07 | 0.923 (7.6) | +0.19 |
| train4 50Hz +4 last | 0.76 ± 0.07 | 0.821 (0.7) | +0.06 |
| **sum z^2, LTD rows (7)** | | **49.4** of 111.4 | |
| no change: pre only / post only | 0.98 / 1.03 | 1.010 / 0.970 (0.1 / 2.3) | |
| no change: train10 -4 last | 0.99 ± 0.09 | 0.808 (4.1) | over-LTD |

Notes: in all three fits the mGluR- and NMDAR-block lanes of Sj01 0.1Hz -10 sit at 1.013 (= no rho change), so
every bit of the post-before-pre LTD at 0.1-20 Hz is the eCB step; rho contributes nothing. The eCB-only ceiling is
d_min -0.29 at every synapse, ratio ~0.72; Sj01 10/20 Hz -10, Sj03 -25 and Zilberter 1AP -10 are below it.
Zilberter train10 -10 LTD survives AM251 (CB1 block): it can only be postsynaptic (rho) LTD.

## 2. Per-synapse mechanism (split1 v5c; Zilberter rows at the V5cz_s5_39 params)
diag_ltd.py, base variant (repro of the fit csv preds: 2e-16). The fractions are of synapse-records in the scored records.
eCB = at least one eCB step; down = rho0 >= 0.5 -> rho_f < 0.5; dep = c* > theta_d without pot at some step; W = max W at own
arrivals / theta_eCB. Section types: 2 = basal, 3 = apical (edges afferent_section_type).

| target | n | eCB all (basal / apical) | down | dep | none | W < .01 | W .01-1 | W >= 1 | c* pk / theta_d |
|---|---|---|---|---|---|---|---|---|---|
| Sj01 0.1Hz -10 | 178 | .38 (.45 / .16) | .00 | .09 | .62 | .17 | .45 | .38 | 0.75 |
| Sj01 10Hz -10 | 172 | .41 (.46 / .26) | .03 | .66 | .55 | .10 | .49 | .41 | 1.12 |
| Sj01 20Hz -10 | 172 | .42 (.46 / .29) | .05 (up .08) | .85 | .51 | .10 | .48 | .42 | 1.33 |
| Markram 10Hz -10 | 172 | .41 (.46 / .26) | .03 | .56 | .55 | .10 | .49 | .41 | 1.05 |
| Sj03 0.1Hz -25 | 178 | .36 (.42 / .16) | .00 | .03 | .64 | .15 | .49 | .36 | 0.71 |
| Sj03 -120 (data 1.05) | 178 | .22 | .00 | .00 | .78 | .28 | .50 | .22 | 0.67 |
| Sj03 -200 (data 1.05) | 178 | .07 | .00 | .00 | .93 | .39 | .54 | .07 | 0.67 |
| Sj03 burst r50 -120 | 178 | .31 (.36 / .14) | .00 | .00 | .69 | .21 | .48 | .31 | 0.70 |
| Sj03 burst r50 -200 | 178 | .20 | .00 | .00 | .80 | .30 | .50 | .20 | 0.70 |
| Sj07 pre only (data 1.00) | 171 | .01 | .10 | .58 | .88 | .81 | .18 | .01 | 1.18 |
| L2/3-L5 Letzkus 1AP +10 | 690 | .00 | .26 | .70 | .73 | 1.00 | 0 | 0 | 1.11 |
| L2/3-L5 Letzkus 3AP +10 distal | 228 | .00 | .23 | .79 | .73 | 1.00 | 0 | 0 | 1.17 |
| L2/3-L5 Letzkus 3AP -10 prox | 123 | .39 (.62 / .27) | .14 | .57 | .46 | .35 | .26 | .39 | 1.02 |
| L2/3-L5 Sj06 50Hz +10 distal | 490 | .02 | .32 | .94 | .59 | .84 | .14 | .02 | 1.94 |
| L2/3-L2/3 Z 1AP -10 | 461 | .25 (.16 / .68) | .02 | .07 | .72 | .17 | .58 | .25 | 0.78 |
| L2/3-L2/3 Z 1AP +10 | 461 | .00 | .15 | .51 | .84 | 1.00 | 0 | 0 | 1.01 |
| L2/3-L2/3 Z train10 -10 last | 455 | .38 (.29 / .73) | .09 | .19 | .62 | .04 | .58 | .38 | 0.83 |
| L2/3-L2/3 Z train4 +4 last | 459 | .31 (.21 / .72) | .16 | .71 | .54 | .08 | .61 | .31 | 1.08 |

- The trigger misses most synapses because W stays under theta_eCB (W .01-1: 45-58 %), not because the VDCC current is absent (W < .01: 10-17 % at L5 -10).
  The 3-synapse BCL figure (1 of 3 with a pool) is the typical case, not an exception.
- The basal/apical split goes in opposite directions: L5-L5 apical contacts trigger less (.16 vs .45); L2/3-L2/3 apical contacts trigger more (.68 vs .16).
  So "location" means bAP-VDCC amplitude, not a basal/apical rule.
- The eCB window is too wide and too shallow at the same time: -120 still steps 22 % (data no change), while -10 steps only 38 %.
  Single bAP vs 5-bAP burst at -120 differ only 0.22 vs 0.31, while the data differ 1.05 vs 0.79.
- theta_d / theta_p median: L5-L5 0.65, L2/3-L5 0.64, L2/3-L2/3 0.96.

## 3. Classes, minimal changes, offline chi2 (fixed params, no refit)
chi2 per pathway at the split1 v5c params (l23l23 at V5cz_s5_39), round 1:

| variant | l5 all / LTD rows / rest | l23 all | l23l23 all / LTD / rest |
|---|---|---|---|
| base | 73.2 / 39.0 / 34.2 | 14.7 | 111.4 / 78.2 / 33.3 |
| noE (eCB off) | 173.4 / 130.9 / 42.5 | **11.6** | 126.7 / 107.7 / 19.0 |
| allE (every arrival steps) | 272.7 / 63.4 / 209.3 | 187.8 | 139.9 / 65.8 / 74.1 |
| thE x0.3 / x0.1 / x0.01 | 78.8 / 101.1 / 174.2 (LTD 28.8 / 23.6 / 24.6) | 21.1 / 26.7 / 42.0 | 117.1 / 126.8 / 150.0 |
| Cs (c* copy, tau*) f .2-.8 | 181-212 (LTD 27-45) | 120-157 | 149-190 |
| Cf (c* copy, tau_E1) f .1-.6 | 180-222 (LTD 27-45) | 57-114 | 150-169 |
| d_min -0.40 / -0.50 | **65.4** / 72.7 | 27.2 / 44.0 | 112.8 / 117.6 |
| theta_d x0.9 / x0.8 | 102.2 / 202.8 | 75.7 / 182.7 | 122.8 / 193.9 (LTD 48.4 / 73.9) |

Classes:
1. **eCB trigger fails (L5-L5: Sj01 0.1/10/20 Hz -10, Sj03 -25, burst -120/-200).** At full recruitment these rows are inside 1.3 SEM of the
   0.72 ceiling (allE: z^2 36.1 -> 5.9). Minimal change: a trigger that fires for one bAP at nearly every synapse at -10/-25 but not at -120,
   plus a veto that switches eCB off when a bAP follows the own pre spike (+10 trains). Lowering theta_eCB alone fails: the LTP rows break
   (LTP z^2 34 -> 150), because the (1 - b) weighting does not stop the trigger in +10 trains once the threshold is low. c*-copy triggers
   (Cs/Cf) also fail; c* is too slow (tau* 278 ms) and too large in pairing trains. Round 2 tests the veto (below).
   New parameters: veto window Tv (bracket 15-40 ms from Sj01: 40 Hz -10 = LTP needs a veto at 15 ms, 20 Hz -10 = LTD needs none at 40 ms);
   veto threshold (here theta_eCB, the fitted one; if theta_eCB is refit low, the veto needs its own anchor).
2. **Timing window too wide (Sj03 -120/-200 no change, but burst -120/-200 LTD).** The pool time constant (tau_E1 100 ms) gives W(-120) ~ 0.3 W(-10).
   A burst is only about 1.4x a single bAP in the pool, so a threshold cannot separate them. LTD_LIT point 3: an eCB trace with production ~ VDCC-Ca integral
   and its own lifetime (FAAH-block anchor -400 ms ~ 0.68). This is a new tau, so it needs a user OK.
3. **d_min is not the problem.** dm.40 is the best single change in l5 (-7.8), but l23 gets +12.5 and LTD_LIT anchors -0.29 (ACEA 0.71). Keep it.
4. **rho-LTD needed, theta_d never crossed (Zilberter 1AP -10, train10 -10 incl. AM251).** theta_d/theta_p = 0.96 in L2/3 cells (c_post term
   dominates, a01 ~ a11). theta_d x0.9 improves the Zilberter LTD rows (78 -> 48) but breaks LTP (+41) and L5/L2/3-L5 (+29/+61).
   No uniform c*-threshold move fixes it. A synapse-local option is a depression band on spine VDCC calcium (LTD_LIT 8). Not tested here, and it needs a user decision.
5. **rho-LTD present and correct (Letzkus +10 rows, Sj06 50 Hz +10 distal).** No change needed. LTD_LIT 7 suggests the distal row may be eCB in vivo; that is a flag only.
6. **Over-LTD (Letzkus 3AP -10 proximal 0.81 vs 0.89; Sj07 pre only 0.90 vs 1.00).** The first comes from eCB at 62 % of basal contacts. The second is rho-LTD from pre spikes alone
   (c* 1.18 theta_d at 58 % of synapses): theta_d is too low relative to c_pre (a00). Both should be fixed by the refit of a00 once the trigger change is in.


### Round 2: bAP veto (fixed params = split1 v5c for all three pathways)
Veto: a triggered own arrival takes no eCB step if, within Tv ms after it, (vV) own unweighted spine VDCC influx > theta_eCB,
(vC) own c* influx > c_pre + kappa c_post, or (vP) an own pot event occurs. Triggers: base W, thE x0.1 / x0.01, allE, Cf0.5.

| variant | l5 (LTD / rest) | l23 | l23l23 | 39 targets | 54 targets |
|---|---|---|---|---|---|
| base | 73.2 (39.0 / 34.2) | 14.7 | 210.8 | 87.8 | 298.7 |
| **vV25** | 76.8 (39.1 / 37.7) | **8.8** | 201.8 | **85.6** | 287.4 |
| vP25 | 73.4 | 14.7 | 211.1 | 88.0 | 299.1 |
| thE x0.1 + vV25 | 74.7 (**22.9** / 51.8) | 13.4 | 198.1 | 88.1 | **286.2** |
| thE x0.1 (no veto) | 101.1 (23.6 / 77.5) | 26.7 | 208.3 | 127.8 | 336.1 |
| thE x0.01 + vV25 | 116.1 (24.9 / 91.2) | 23.0 | 203.7 | 139.1 | 342.9 |
| thE x0.01 + vC25 / vP25 | 174.2 / 173.0 | 37.9 / 39.5 | 223.1 | 212 | 435 |
| allE + vV15 / 25 / 35 | 173.9 / 171.1 / 171.1 | 128.3 / 126.5 / 125.3 | 147-150 | 296-302 | 446-449 |
| allE + vC25 / vP25 | 272.7 / 270.3 | 180.4 / 183.4 | 160 | 453 | 613 |
| Cf0.5 + vV25 | 107.5 | 45.0 | 207.7 | 152.5 | 360.2 |
| dm.40 (round 1) | 65.3 | 27.2 | 209.3 | 92.5 | 301.8 |

Readings:
- vV vetoes the eCB in +10 trains. 50 Hz +10 goes 1.393 -> 1.421 at base and 1.157 -> 1.331 at thE x0.1.
  Sj07 step pairing keeps part of its eCB (1.513 -> 1.598; AM251 lane 1.71).
- allE + vV still fails, for two reasons. Synapses with almost no VDCC current trigger but cannot be vetoed (veto threshold
  theta_eCB). And -120/-200 still step.
- c*-based vetoes do nothing within 25 ms (aC25 = allE exactly), so the veto must read the fast VDCC current.
- With the veto threshold tied to theta_eCB, a refit can lower theta_eCB without losing the +10 LTP. That is the refit below.

## 4. Ranked refit candidates (split1 39 targets; gpu_v7_rho.py = v6 kernel + options, gpu_v5/v6 untouched)
1. **V7v: v5c + bAP veto** (SET veto_T 25). Fixed-param gain -2.2 (39 targets) / -11 (54 targets). The refit is expected to lower
   theta_eCB (thE x0.1 + veto: L5 LTD rows 39 -> 23).
   - New parameter: Tv = 25 ms, fixed, not fitted. It is bracketed to 15-40 ms by Sj01 (40 Hz -10 is LTP, so the bAP 15 ms after the pre must veto;
     20 Hz -10 is LTD, so a bAP 40 ms after must not), and is insensitive inside that range (15/25/35 give equal chi2).
   - Veto threshold: tied to theta_eCB, so no new threshold.
   - Mechanism anchor: LTD_LIT 5 (ifenprodil leaves 50 Hz +10 LTP unchanged; AM251 raises step-pairing LTP). Needs the user's OK as a new mechanism.
2. **V7vn: V7v + trigger relative to the own bAP** (SET ecb_ref 1).
   - theta_eCB,i = k_E x P1_i, where P1_i is the synapse's own single-bAP VDCC pool jump (an independent per-synapse
     calibration, like Chindemi's c_post). k_E is free in the theta_eCB slot; the seed is 0.5.
   - Expected anchor for k_E: Sj03 timing, LTD at -25 and none at -100 (Fig 9B), gives exp(-100/tau_E1) < k_E < exp(-25/tau_E1), i.e. 0.37-0.78.
   - Replaces the absolute threshold; no other new parameter. Untested offline.
   - Overlaps A5's E2 (SCALE_E cpost_cai, peak free Ca); here it is the VDCC pool from the records and it carries the veto.
3. Not recommended:
   - d_min -0.40 (anchored -0.29; costs L2/3->L5 +12).
   - theta_d x0.9 / x0.8 (breaks L5 and LTP).
   - c*-copy triggers or vetoes.
   - allE.
4. Next after these: an eCB-lifetime trace (LTD_LIT 3) for burst -120/-200; a VDCC-Ca depression band for Zilberter (user decision).

Jobs (2g MIG, def-emuller; logs/v7_fit_<id>.out; results/v7_<TAG>.json):

| job | tag | what | size |
|---|---|---|---|
| 22156547 | S1_V7off_r | MAXITER 0, options off: must print REPRO OK vs v5_S1_V5c_s6_39 (87.835835) | 50G 0:15 |
| 22156548 | S1_V7v_r | MAXITER 0, veto on at the v5c params: should give about 85.6 (diag vV25: l5 76.84, l23 8.77). Prints REPRO DIFF by design | 50G 0:15 |
| 22156549 | S1_V7vn_r | MAXITER 0, veto + ecb_ref at k_E 0.5 (seed /scratch/dhuruva/ltd_diag/seeds/v5_S1_V5c_s6_39_kE0.5.json); first offline score of candidate 2 (REPRO DIFF by design) | 50G 0:15 |
| 22156550 / 22156551 | S1_V7v_s5 / _s6 | V7v refit, seeded from v5_S1_V5c_s6_39 / unseeded; afterok 22156547 | 50G 0:30 |
| 22156552 / 22156553 | S1_V7vn_s5 / _s6 | V7vn refit, seeded from the kE0.5 json / unseeded; afterok 22156547 | 50G 0:30 |

Rescore size: 50G, from A5's 22155608 (MaxRSS 39.2 GB + 25 %); the 40G suggested was below that.
Diag sizing: run_diag_ltd.sh header (round 2 used 3.99 GB of 4G, at the limit -> 5G next).

### V7 results and the ecb_ref uniformity fix

| tag | chi2 total (BIC) | l5 / l23 | theta_V | theta_eCB | gamma_d / gamma_p |
|---|---|---|---|---|---|
| V7off_r | 87.84, REPRO OK | | | | |
| V7v_r / V7vn_r | 85.61 / 73.82 | | | | |
| V7v_s5 / _s6 | 67.77 (97.08) / 73.41 | 63.71 / 4.06 | 1.010 | 1.352 (absolute; v5c 20.68) | 26.4 / 590 |
| V7vn_s5 / _s6 | 57.56 (86.87) / 60.13; not valid | 48.96 / 8.60 | 1.774 | k_E 0.492 | 21.6 / 564 |

V7vn is not uniform. ecb_ref 1 takes P1 from isolated bAPs in the records. L2/3->L5 has none (3171 of 3171 synapses fall back to P1 = 1), so there theta_eCB is absolute, while at L5 it is k_E x P1 (median 3.92). Its 57.56 is not used. k_E 0.49 is inside the 0.37-0.78 anchor, but only for L5.

Fix, `ecb_ref: 2`: one source for every pathway. P1_i = 1e5 x vdcc_q_post_i (A4 cexp, /scratch/dhuruva/split1/cexp), through fit_v6's loader and SCALE_E. NaN rows get the pooled median over all pathways. ecb_ref 1 code is unchanged. The seed is V7vn_s5 with theta_eCB = 0.4916 x 3.92 = 1.927 (absolute at the median L5 synapse; /scratch/dhuruva/ltd_diag/seeds/v7_V7vn_s5_abs.json). fit_v6 divides it by the pooled median uE.

| job | tag | what | size |
|---|---|---|---|
| 22159502 | S1_V7vn2_r | MAXITER 0 at the V7vn_s5 params, ecb_ref 2 (repro skipped: spec differs) | 47G 0:15 |
| 22159503 / 22159504 | S1_V7vn2_s5 / _s6 | refit seeded from the converted V7vn_s5 / unseeded; afterok 22159502 | 38G 0:30 |
| 22159505 | S1_V7v_s7 | V7v basin check, seed 7, unseeded | 38G 0:30 |
