# Why S1C (fitted on 7 L5 shape targets) fails on the L2/3 pathways (2026-10-02)

This is a diagnosis only. No rule, kernel or fit was changed.

Inputs:
- s1C_s and s1C_u (fit + `_val*` csvs), s1A_s, and r1E_s, which is the same G8X + eCB rule fitted on all 65 targets.
- The rescore log `logs/s1C_s_val_22266779.out`, for the v8 calib lines and the v10 counters.
- The prior diagnoses listed in the brief.

Correction to the brief: `letzkus_3ap_200hz_dt-10ms@proximal` is predicted at **1.27** (data 0.89 +- 0.03, z 12.7, z^2 161). The value 1.71 is `3ap +10 proximal`, whose data are 1.30.

## 0. Bottom line

1. **S1C has no rho-LTD.**
   - theta_p - theta_d = (2.31 - 2.28) c_pre + (1.27 - 0.67) c_post = **0.04 c_pre + 0.60 c_post**, so the depression band is closed wherever c_post is small.
   - gamma_d is 20.8, at the bottom of its 20-250 box, against gamma_p 165.
   - Every LTD in S1C is therefore the eCB step. That step is post-before-pre only (the (1 - b) weight plus the 25 ms bAP veto), and its ceiling is 1 + d_min = 0.71.
   - The L5 shape targets do not need rho-LTD: L5 LTD is CB1-dependent (LTD_LIT 1). So the 7-target fit had no reason to keep a band.
   - Every L2/3 pre-before-post LTD row needs rho-LTD (LTD_LIT 6, 8), so they all fail.
2. **The 7 L5 targets leave the band undetermined.**
   - s1C_u has the same rule and the same 7 targets, with an equally good L5 shape (Markram +10 1.22, -10 0.80; 20 Hz 1.37).
   - Its a00/a10 are 1.95/2.50, i.e. a band of 0.55 c_pre, and its L2/3 validation scores are 108 / 196 (against 355 / 336 for s1C_s).
   - So part of the failure is only an unconstrained direction, and a refit can fix it.
3. **Part of it is structural.** r1E_s fitted this exact rule with every L2/3 target in the objective. It still misses:
   - Zilberter 1AP +10 (0.97), 1AP -10 (0.74), train +5 first (1.27) and train -10 under AM251 (1.11);
   - Letzkus 1AP +10 (0.85, z 4.4).
   
   It also lost Markram 10 Hz +10 (0.96 against 1.20), 20 Hz +10 (1.05) and 50 Hz -10 (0.77). The band that gives distal L2/3 LTD also takes away L5 LTP.
4. **Verdict.**
   - L2/3->L5: adding the targets fixes about half of the chi2 (3AP +/-10 proximal) and only part of the +10 distal LTD rows.
   - L2/3->L2/3: adding the targets is **not enough**. The single-pairing LTD and the AM251-resistant train LTD need a mechanism the rule does not have, or a user decision on uniformity (MECH_NECESSITY 6.3).
   - Test fits 22289873/74 (section 4) check whether the 7 L5 shape targets survive once the L2/3 rows are in the fit.

## 1. Mechanism evidence (S1C_s rescore log, per pathway)

| quantity | L5->L5 | L2/3->L5 | L2/3->L2/3 |
|---|---|---|---|
| own single-bAP VDCC charge uE, median (q10-q90) | 109 (23-382) | **17** (9-97) | 47 (10-168) |
| single-bAP shaft dCa, median | 0.062 uM (proximal 0.53) | no isolated bAP in the records; **3AP 200 Hz record peak median 1.2 nM**, q90 1.93 uM | 0.115 uM (proximal 0.23); trains up to ~16 uM (GATE_IMPL) |
| rows with licensed pot / blocked pot | 1216 / 509 of 5176 | **250 / 260 of 2688** | 2116 / 1519 of 6716 |
| rows with an eCB step / a vetoed trigger | 1368 / 1009 | **42 / 131** | 1966 / 1719 |
| null baseline (frozen rho, nmdar_block) | 1.013 | **1.034** | 1.03 |

Reading the table:
- **Licence (theta_G 0.15 uM).**
  - At L2/3->L5 it is bimodal. Most Letzkus synapses never open it, even with a 3AP burst; a proximal minority opens it and potentiates hard (3AP +10 proximal 1.71).
  - At L2/3->L2/3 it is open in every train, so it gives no protection there.
- **eCB.** It is practically absent at L2/3->L5 (1.6 % of rows step).
- **Readout offset.** The null baseline of 1.034 is a readout / rho0 offset shared by S1A and S1C. It costs about 1 SEM on Letzkus 1AP (SEM 0.03), so it is not the main cause.

## 2. Failing targets: cause and refit outlook

The prediction columns are s1C_s / s1C_u / r1E_s. r1E_s is the same rule with these targets fitted (but 40 L5 targets).

| target | data | s1C_s | s1C_u | r1E | likely cause (evidence) | refit fixes? |
|---|---|---|---|---|---|---|
| Letzkus 3AP -10 proximal | 0.89 | 1.27 | 0.93 | 0.885 | Closed band: licensed proximal synapses cross theta_p in both orders. Pre lands on the 3rd AP, so the veto fires and the eCB (42 rows) cannot help | **yes** |
| Letzkus 3AP +10 proximal | 1.30 | 1.71 | 1.64 | 1.25 | All the LTP comes from the licensed minority at gamma_p/gamma_d = 8 | **yes** |
| Letzkus 1AP +10 (distal-weighted) | 0.72 | 1.04 | 0.93 | 0.85 | Needs rho-LTD (pre comes before the only bAP, so no eCB). theta_d (2.28 c_pre) sits above the single-pairing c*, so nothing happens: 1.04 = the null. At distal sites the bAP adds ~1 nM of shaft Ca (uE 6x below L5), so c* there is close to the EPSP-alone value (nopost 0.98) | partial (best 0.85) |
| Letzkus 3AP +10 distal | 0.79 | 1.14 | 1.09 | 0.96 | Same cause. Also, 13 % of the "distal" synapses are basal and a quarter of the pairs sit below 281 um, so some of them get licensed LTP (LETZKUS_LOCATION) | partial |
| S&H 50 Hz +10 distal / all | 0.86 / 1.06 | 1.23 / 1.32 | 1.20 / 1.29 | 1.10 / 1.21 | The data point to CB1 LTD (EPSP-alone LTD is blocked by AM251, LTD_LIT 7). The 25 ms veto turns eCB off in +10 trains, which L5 50 Hz +10 needs | partial; conflicts with the veto |
| Zilberter 1AP +10 | 0.64 | 1.01 | 0.96 | 0.97 | CB1-independent post LTD at a single pairing, in the same local state as L5 0.1 Hz +10 (0.97, no change; L23L23_DIAG 5) | **no** |
| Zilberter 1AP -10 | 0.56 | 0.74 | 0.74 | 0.74 | Pure eCB at the 0.71 ceiling (d_min -0.29, anchored). Going lower needs rho-LTD that L5 -10 under AM251 (1.01) forbids | **no** |
| 5ap 10 Hz +10 / 20 Hz +10 | 0.76 / 1.07 | 1.36 / 1.37 | 1.15 / 1.22 | 0.90 / 1.04 | 5 bAPs with an open licence push c* over theta_p. With no band, pre-post LTD cannot appear | 20 Hz yes; 10 Hz partial |
| train10 post only | 1.03 | 1.27 | 1.19 | 1.09 | **The bAP train alone crosses theta_p** (10 summed c_post against a11 = 1.27) while the licence is open: rho-LTP with no pre spike at all | mostly |
| train10 +5 first | 0.97 | 1.41 | 1.36 | 1.27 | Order inversion: the largest c* with an open licence, while the data show no change. No c*-threshold rule separates it from +4 last (L23L23_DIAG 5) | **no** |
| train10 -10 last / AM251 | 0.72 / 0.73 | 0.94 / 1.30 | 0.88 / 1.22 | 0.81 / 1.11 | rho gives train LTP and the eCB pulls it down, so the AM251 lane exposes the rho-LTP. The data show postsynaptic, CB1- and NMDAR-independent, mGluR-dependent LTD (LTD_LIT 8) | **no** (uniformity conflict) |
| train4 / train8 / train10 +4 last | 0.76 / 1.15 / 1.49 | 1.13 / 1.39 / 1.41 | 0.96 / 1.34 / 1.36 | 0.90 / 1.13 / 1.19 | The AP-count threshold sits too low: 4 APs already potentiate. Lowering it with a refit also lowers train10 +4 | partial (a trade-off) |

## 3. The causes in the brief, one by one

1. **Licence or eCB never engaging at L2/3.**
   - The licence: closed at most Letzkus synapses (burst shaft median 1.2 nM against 0.15 uM) and open in all Zilberter trains. Where it is closed, c* stays under theta_d anyway, so the closed licence is not the reason the LTD is missing. The open licence at L2/3->L2/3 is the reason post-only and +5-first trains potentiate.
   - The eCB: 42 of 2688 rows at L2/3->L5. This matters for 3AP -10 proximal only, and a refit fixes that through the band.
2. **c_pre / c_post scaling.** This is the main fit-level cause.
   - With a00 ~ a10 the band depends only on 0.60 c_post. It is therefore narrowest exactly where c_post is smallest, at distal L2/3->L5 sites.
   - With a11 small, a 4-10 bAP train crosses theta_p at the large-c_post Zilberter contacts.
   - Every C fit has a01 ~ 0 or below 1 (s1C_s 0.67, s1C_u 0.011, r1E 0.007), so no fit uses the c_post term to protect L5 single pairings.
3. **bAP attenuation at distal Letzkus sites.**
   - Confirmed: uE median 17 against 109 at L5, and a burst shaft median of 1.2 nM.
   - The emodel under-drives bAP Ca beyond ~130 um by about 3x (DECISIONS T33). So at distal sites 1AP +10 is locally almost the same as the EPSP alone (nopost 0.98), and no threshold on c* can give 0.72 against 0.98 there.
   - Part of the Letzkus +10 LTD is therefore an emodel limit (ion_fitter).
4. **rho0 / initial state.** The 1.034 frozen baseline is a constant offset: z +1.1 on Letzkus 1AP, about +0.5 elsewhere. It is not rule-specific and is minor.
5. **Doublet-filtered keep lists.**
   - L2/3->L5 keeps 102 of 120 pairs.
   - The split2 L5TTPC often fails 200 Hz (201 of 300 burst tasks failed, DECISIONS T39). The 3AP targets therefore rest on the cells that do follow 200 Hz, probably the more excitable, more often licensed ones. That would bias 3AP toward LTP.
   - The size of this bias is not known. It cannot explain the single-AP failures.
6. **eCB LTD needing post-pre timing.**
   - Yes, by design: the (1 - b) weight and the veto, which L5 50 Hz +10 and sj07 need.
   - So every L2/3 pre-before-post LTD row has to come from rho-LTD, and S1C switched rho-LTD off.
   - The S&H distal row is the one case where the data suggest CB1 LTD inside a +10 train. That conflicts with the veto.

## 4. Can refitting with these targets fix it?

**L2/3->L5: largely yes. L2/3->L2/3: no.**
- The refit must open the band (a10 - a00 ~ 0.5-1 c_pre, as in s1C_u / r1E).
- The open question is whether the band can stay open without the L5 LTP loss seen in r1E. That loss may also come from the 33 near-1 L5 targets r1E carried.

Test fits submitted (S1C rule, 7 free, 7 L5 shape + 8 L2/3->L5 + 15 L2/3->L2/3 = 30 targets; Egger and Letzkus -10 distal dropped):
- Script: `rho_redesign/run_s1c_l23fit.sh s|u`.
- **22289873** `s1CL_s`: seeded from s1C_s + s1C_u, seed 5. **22289874** `s1CL_u`: unseeded, seed 6.
- Sizing: 91G 0:45 on 3g.40gb, from r1E_s 22235501 (22:29, 72.9 GB).
- Read `/scratch/dhuruva/s1c_l23fit/s1CL_{s,u}{.csv,_l23.csv,_l23l23.csv}` and `logs/s1CL_{s,u}_<id>.out` (grep `chi2 L5|chi2 l23|v10 counts`).

How to read the result:
- **Markram +/-10 within 1 SEM, and L2/3->L5 under ~40.** The L2/3->L5 part is a pure parameter problem; put those rows into the stage-3 core.
- **Markram +10 below ~1.1.** The distal-LTD band and L5 LTP conflict inside this rule (r1E repeated), and section 5 applies.
- **L2/3->L2/3 rows marked "no" in section 2.** Expect them to stay out (r1E). If they do, they are confirmed as structural.

Run `seff` on both jobs afterwards and record the result.

## 5. Smallest rule changes, if needed (all synapse-local, uniform, new names)

1. **No new mechanism for L2/3->L5 proximal.** Refit with an open band (section 4).
2. **Distal pre-post LTD (Letzkus 1AP / 3AP +10, S&H distal).** The local cue the data use is the weak bAP: LTD gets deeper with distance (r -0.69). Two routes:
   - (a) The emodel side: a stronger distal bAP Ca (ion_fitter, new emodel name). This needs no rule change.
   - (b) The rule side, if (a) does not come: let the Chindemi c_post term carry it. a01 > 1 protects synapses with a large c_post (L5, ~95 um), while 1 < a00 < 1 + s_distal opens the band where c_post ~ 0. This is a 0-parameter re-parameterisation, but its margin rests on the distal NMDA-bAP supralinearity, which the under-driven emodel makes small. Test it by seeding a refit at a01 ~ 2, a00 ~ 1.2 before adding anything.
3. **L2/3->L2/3 post LTD (train -10 under AM251, 5ap 10 Hz, train4) and the train order (post-only, +5 first).** The data name the mechanism: postsynaptic group I/II mGluR, CB1- and NMDAR-independent, set by an intermediate L-type Ca level (D890 and BAPTA, LTD_LIT 8). The smallest synapse-local form is a new arm, **"mGluR-gated VDCC depression" (vdep)**:
   - rho depresses at gamma_d while the synapse's own (1 - b)-unweighted spine VDCC pool (the existing c_VDCC, tau_E1 100 ms) lies above theta_vdep x P1_i (own single-bAP pool, as ecb_ref 2), but only within an own-glutamate window (b > 0.5, the own pre spike as the mGluR proxy).
   - It adds one parameter, theta_vdep, anchored to Zilberter's D890 Ca ratio (0.37: LTP turns into LTD).
   - It reads only own Ca and own pre spikes. It catches train -10 / AM251 and train4, because the L2/3 pool there is ~3x the L5 single-bAP pool (MECH_NECESSITY 4: 1.24 vs 0.36 theta).
   - It **cannot** fix Zilberter 1AP +/-10. Their local state equals the L5 single pairings, which AM251 shows must not change. These two rows remain a uniformity conflict: accept about 20 chi2 or let the user decide (MECH_NECESSITY 6.3).
   - The +5-first / post-only LTP also needs potentiation to read own-glutamate x bAP coincidence (Gb, L23L23_DIAG 3). Under the current "failed licence = depression" convention Gb was fatal to L5. A neutral failed-Gb form is untested.
   - Test vdep only after 22289873/74 show what refitting alone leaves.
