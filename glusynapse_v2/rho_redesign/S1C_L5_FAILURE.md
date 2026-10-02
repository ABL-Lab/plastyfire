# Why S1C misses the L5->L5 dt-0 / +25 / -10 high-frequency and Sj07 targets (2026-10-02)

Fit: s1C_s.json (G8X shaft licence 0.15 uM in 100 ms + eCB, theta_eCB 0.437, gamma_d 20.8, gamma_p 165, chi2 10.3 on 7 core).
Predictions are taken from s1C_s_val.csv, with s1A_s_val.csv (Chindemi, no licence, no eCB) and s1B_s_val.csv (licence, no eCB) for comparison.
Kernel facts used below come from gpu_v10_rho.py and LTD_DIAG.md:
- **Step size.** With A_eCB = 1, one triggered, unvetoed own arrival sets d = d_min (-0.29). At that synapse this multiplies the ratio by about 0.71.
- **Trigger.** W is the own spine VDCC pool (tau_E1 100 ms), weighted by (1 - b). The b trace counts own arrivals and decays with tau 70 ms. So a bAP that follows an own arrival closely barely fills W. A trigger needs W > k_E P1_i.
- **Veto.** A triggered arrival is vetoed if the own VDCC influx over [t_a, t_a + 25 ms] exceeds the same threshold.

## 1. Data check (PROTOCOLS.md, L5_EXTRA.md, CHINDEMI_DATA.md)
- **Source.** The 20 Hz dt 0 / +25 / -25 rows and the 40 and 50 Hz dt 0 rows are Chindemi's biodata/paired_recordings.csv rows sjs01_11-17. They are Sjostrom 2001 paired L5->L5 recordings: 5 pre + 5 post spikes, 15 sweeps every 10 s, **n = 5 each**. The values are digitised by Chindemi and were not re-read from the 2001 PDF here, so treat them as unverified. They are paired data, so they are allowed as fit targets.
- **40 Hz dt 0.** Its SEM is 0.038 at n = 5 (0.93 +- 0.04), the tightest in the dt-0 set. That is why it alone gives z^2 132.
- **50 Hz dt 0** (0.92 +- 0.035, n 5) was extracted (L5_EXTRA) but is **not in the scored set**: it is absent from every val csv. It is the row that brackets the veto window (section 3), so it should be scored as validation.
- **40 Hz +10 SEM.** Ours is 0.14 (n 11, from the text); Chindemi's csv has 0.113 (n 6).

## 2. What S1C does at high frequency (pred; data in brackets)
| f | +10 | 0 | +25 | -25 | -10 | -10 AM251 lane |
|---|---|---|---|---|---|---|
| 20 Hz | 1.31 (1.31) | 1.30 (0.76) | 1.30 (0.65) | 0.96 (0.70) | 0.94 (0.65) | 1.29 (1.02) |
| 40 Hz | 1.37 (1.53) | 1.37 (0.93) | - | - | 1.22 (1.51) | - |
| 50 Hz | 1.40 (1.57), AM251 lane 1.42 (1.59) | not scored (0.92) | - | - | 1.19 (1.70) | - |

The rho part is timing-blind at >= 20 Hz. At 20 Hz every dt gives 1.29-1.31 without eCB, including the AM251 lane of -10. All the timing comes from the eCB step, and the eCB step is placed wrongly.

## 3. Per-target causes
**(a) 20 Hz dt 0, 20 Hz +25, 40 Hz dt 0: LTP predicted, LTD or no change measured. STRUCTURAL (eCB veto window and (1 - b) weighting).**
- Every own arrival is followed by a bAP 0 ms (dt 0) or 25 ms (+25, and 40 Hz dt 0) later, so the 25 ms veto window covers each arrival.
- Every bAP lies within 25 ms of an own arrival, so (1 - b) keeps W near zero.
- eCB therefore never steps, and the prediction equals the rho-only +10 value. Evidence: pred(dt 0) = pred(+10) to 0.01 at both 20 and 40 Hz.
- No theta_eCB value fixes this, because the veto uses the same threshold as the trigger.
- With full eCB on top of the current rho, the predictions would be 0.71 x 1.30 = 0.92 at 20 Hz dt 0 and +25 (z +1.6 / +1.5) and 0.71 x 1.37 = 0.97 at 40 Hz dt 0 (z +1.1). So the eCB repair alone brings these rows within 1.6 SEM.
- **Edge artefact.** In steady state, 20 Hz +25 and 20 Hz -25 are the same train (ISI 50 ms). The data agree (0.65 / 0.70), but the model gives 1.30 / 0.96. The difference comes only from which spike starts and ends the train.

**(b) 50 Hz -10 (1.19 vs 1.70) and 40 Hz -10 (1.22 vs 1.51). STRUCTURAL, two parts.**
- *eCB leak in the -10 order.*
  - In S1C, -10/+10 = 1.19/1.40 at 50 Hz and 1.22/1.37 at 40 Hz. In S1B (no eCB) the two orders are nearly equal: 1.38/1.40 and 1.34/1.41.
  - So eCB removes 11-15 % of the ratio, and only in the -10 order. The data show no eCB effect in either order (the 50 Hz +10 AM251 arm is 1.59 vs 1.57).
  - Most likely cause: the last pre spike of each train has a bAP 10 ms before it (a trigger) and none after it (no veto). With A_eCB = 1, one such arrival in 15 trains is enough for full d_min at that synapse. The Ae02 job tests this.
- *rho-LTP cap*, as in (c). Even without eCB the prediction is 1.40 (z -1.6).

**(c) LTP amplitude at >= 20 Hz is capped at 1.30-1.42.**
- Affected rows (pred vs data): 40 Hz +10 1.37 vs 1.53, 40Hz_5ms 1.35 vs 1.50, 50 Hz +10 1.40 vs 1.57, Sj07 AM251 lane 1.52 vs 2.13.
- S1B (licence, no eCB) shows the same cap: 1.41 / 1.46 / 1.40 / 1.54. S1A (no licence) reaches 1.66 / 1.77 / 1.70 / 1.82.
- Markram +5 is flat in S1C from 20 to 40 Hz (1.31 / 1.37 / 1.35), while the data rise (1.37 / 1.42 / 1.50).
- Reading (hypothesis): the G1 licence (peak shaft dCa within 100 ms > 0.15 uM) is close to binary per synapse once one bAP has passed. Higher frequency then recruits no new synapses, and the licensed ones saturate at rho = 1.
- If this is right, it is STRUCTURAL under the anchored 0.15 uM peak licence and not fixable by gamma_p. The gp2 and lic0 jobs decide it.

**(d) Sj07 step pairing: AM251 lane 1.52 vs 2.13, control 1.39 vs 1.62. Mixed.**
- *LTP ceiling.* This is the cap from (c): S1B, with no eCB, gives 1.54. STRUCTURAL.
- *NO is off.* The data NO-block arm (1.36) implies about x1.19 from NO. Even if that factor were also in the AM251 arm, rho alone would need about 2.13 / 1.19 = 1.79, so NO closes only half the gap. Parameter-limited: stage 2C adds NO.
- *eCB depth.*
  - Data control/AM251 = 0.76; model 1.39/1.52 = 0.92.
  - The veto sums the own VDCC influx over 25 ms, and the 200 ms depolarising step supplies sustained VDCC influx even between bAPs. Ongoing pre arrivals also keep (1 - b) low, so most triggers are vetoed or never form.
  - This is the opposite of Sjostrom 2004 dLTD, where pre spikes plus a subthreshold step give CB1-dependent LTD. STRUCTURAL: the veto should read a bAP-like transient, not integrated VDCC charge.

**(e) 20 Hz -10 (0.94 vs 0.65) and its AM251 lane (1.29 vs 1.02). Partly PARAMETER-limited.**
- c* integrates the whole 250 ms train, so lags {-10, +40} reach theta_p like lags {+10, -40}. Even full eCB gives at best 0.71 x 1.29 = 0.92.
- The AM251 lane directly measures rho at 20 Hz without eCB. Fits that included it reached 1.11-1.20 (LTD_DIAG section 1), at a cost to the LTP rows.
- **Data conflict, not fixable.** Sj01 10 Hz -10 (0.57) and Markram 10 Hz -10 (0.79, fitted 0.79) are the same protocol in the model, an irreducible z^2 of about 3.7 (LTD_DIAG 6).

## 4. Are the targets mutually consistent?
Each protocol reduces to the lags each pre spike sees: pre to the next post (+) and the previous post to pre (-).

| protocol | lags | data | protocol | lags | data |
|---|---|---|---|---|---|
| 20 Hz +10 | +10, -40 | LTP 1.31 | 40 Hz -10 | +15, -10 | LTP 1.51 |
| 20 Hz 0 | 0, +-50 | LTD 0.76 | 40 Hz 0 | 0, +-25 | none 0.93 |
| 20 Hz +-25 | +25, -25 | LTD 0.65/0.70 | 50 Hz +-10 | +10, -10 | LTP 1.57/1.70 |
| 20 Hz -10 | +40, -10 | LTD 0.65 | 50 Hz 0 | 0, +-20 | none 0.92 |
| 40 Hz +10 | +10, -15 | LTP 1.53 | 100 Hz 0 | 0, +-10 | LTP 1.25 (dropped, biased subset) |

**Under any rule monotone in spine Ca or c* (any thresholds), 40 and 50 Hz dt 0 contradict 40 and 50 Hz +-10.**
- Let f(lag) be the Ca a pre spike gets from the bAP at that lag. 40 Hz dt 0 receives f(0) + f(+25) + f(-25), against f(+15) + f(-10) at 40 Hz -10. For dt 0 to have less Ca, f(0) must be below about 0.15 f(+15).
- At 50 Hz the same comparison needs f(0) < f(+10) - f(+20), about 0.13 f(+10), with NMDA tau_d 70 ms.
- The model has f(0) close to f(+10): the dt-0 predictions equal the +10 ones. Real NMDARs do not reach that suppression either, since they rise in a few ms. So no Ca-amplitude rule can give these signs.
- **20 Hz +10 LTP against +25 LTD is not a strict contradiction.** f(+25)/f(+10) is about 0.8, so a theta_p in between would do, but it is a knife edge, and the same c* then gives LTP at 20 Hz -10 (AM251 lane).

**With a timing veto, the set is self-consistent.** LTP appears exactly when a pre spike is followed by a post spike within about +10 to +15 ms. A post spike in the 25-50 ms before a pre spike gives eCB-LTD unless a +10..+15 pairing vetoes it. Lags 0, +20 and +25 must not veto. This brackets the veto window:
- **Start t0 in (0, 10) ms:** the coincident bAP at dt 0 must not veto, and the +10 bAP must.
- **End Tv in [15, 20) ms:** +15 must veto (40 Hz -10), and +20 must not (50 Hz dt 0).

This tightens LTD_DIAG's 15-40 ms bracket. The only true conflict is Sj01 vs Markram 10 Hz -10.

## 5. Parameter-limited or structural, and what to add to the core
| target | class | add to core? |
|---|---|---|
| 40 Hz dt 0 (z^2 132) | structural: veto window | yes, after the veto change (tight SEM, sets t0); before it, it only wrecks the fit |
| 20 Hz dt 0, 20 Hz +25 | structural: same cause | no: redundant with 40 Hz dt 0 + 20 Hz -10 (validation) |
| 20 Hz -25 | structural: edge artefact | no: the same train as +25 |
| 50 Hz dt 0 | not scored yet | score as validation (brackets Tv < 20) |
| 50 Hz -10 (core), 40 Hz -10 | structural: end-of-train eCB leak + LTP cap | 40 Hz -10 no (consistent with 50 Hz -10, validation) |
| 40 Hz +10 (core), 40Hz_5ms, 50 Hz +10 | structural if gp2 confirms saturation (licence cap); parameter-limited (gamma_p) otherwise | no new rows: 40 Hz +10 is already core |
| Sj07 AM251 lane | cap + NO (stage 2C) + veto reads sustained Ca | after the licence/veto fix; adding it now pulls gamma_p and thresholds against Markram |
| 20 Hz -10 AM251 lane | partly parameter-limited (theta_p vs the 20 Hz c* integral) | **yes, now**: the only direct rho-only timing constraint at 20 Hz |
| Sj01 10 Hz -10, 0.1 Hz dt 0 | conflict with Markram / already fine (z 0.9) | no |

## 6. Smallest rule change (uniform, synapse-local, no new free parameter for steps 1-2)
1. **Veto window (t0, Tv] instead of [0, 25].**
   - The veto still reads the own spine VDCC influx, but only between t_a + t0 and t_a + Tv. Fixed, bracketed by Sjostrom 2001 as above: t0 = 3 ms (bAP after glutamate binding; any value in (0, 10) works), Tv = 17 ms. These two constants replace veto_T 25.
   - Read it as a peak (transient) instead of a 25 ms integral, so that a depolarising step does not veto (fixes the Sj07 eCB depth; Sjostrom 2004 dLTD).
   - It needs a new kernel file (gpu_v11_rho.py, copy of v10, options veto_t0 and veto_peak). Later it needs a new mod for the BCL. v10 is not edited.
2. **Drop the (1 - b) weighting (v5_mode 1).** With a correct veto it is redundant, and it blocks the dt-0 / +25 triggers.
3. **Only if Ae02 confirms the end-of-train leak: graded step, A_eCB free** (k 7 to 8; unanchored, PARAM_ANCHORS lists A_mglu as "none"). Expected value 0.1-0.3, since tLTD saturates over 50-100 pairings (Sjostrom 2003).
4. **Only if gp2 shows saturation: the low-pass shaft licence G2** (gate_src 2, gate_win 0; GATE_ALT 2). Its plateau grows with AP frequency (Helmchen 1996), so 40-50 Hz recruits more synapses. Its threshold is set by the single-bAP calibration, with no new free parameter. It already exists in v10, so it can be tested at once.
- No change is proposed to rho or the c* thresholds. The fit stays at a00-a11, gamma_d, gamma_p and theta_eCB (7), plus A_eCB if needed.

## 7. Diagnostic jobs (submitted; fixed-parameter L5-only rescores of s1C_s.json, maxiter 0, all 40 L5 targets)
Script: /scratch/dhuruva/s1c_diag/run_s1c_diag.sh. Outputs: /scratch/dhuruva/s1c_diag/<tag>.csv. Logs: plastyfire/logs/s1c_diag_{ecb,ltp}_<id>.out (grep "=== |chi2 L5|v10 counts").
Sizing: 10G, 0:15, 3g MIG, 1 CPU, from the measured stage-1 L5-only fits (1.5-1.7 min, 6.5-7.8 GB). Run seff on both jobs afterwards.

**22289864 (ecb):** base, vT15 (veto_T 15), vT15w1 (+ v5_mode 1), Ae02 (A_eCB 0.2).
- base must reproduce the L5 rows of s1C_s_val.csv.
- vT15:
  - 20 Hz +25 falls (no veto at +25) but 40 Hz dt 0 stays at 1.37: the coincident bAP vetoes, so t0 is needed.
  - Both fall: Tv alone is enough.
- vT15w1: 20 Hz dt 0 and +25 move towards 0.92 while the +10 rows hold, so step 2 is confirmed.
- Ae02: 40 and 50 Hz -10 rise to the +10 values while 0.1 Hz -10 and Markram -10 stay at or below 0.80, so the leak is confirmed (step 3).

**22289865 (ltp):** gp2 (gamma_p x2, seed /scratch/dhuruva/s1c_diag/s1C_gp2.json), lic0 (lic_src 0, licence off), lic0gp2.
- gp2: 40/50 Hz +10 and the Sj07 AM251 lane move by less than 0.03, so LTP is saturated: the licence cap, structural, step 4. If they rise, LTP is gamma_p-limited: parameter-limited, and the fitter traded it against Markram 10 Hz +10.
- lic0 shows how much of the cap is the licence.
