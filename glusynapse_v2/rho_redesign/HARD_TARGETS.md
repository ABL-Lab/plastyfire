# Hard targets for GluSynapse_v2: which ones, and why (2026-10-02)

**Question.** Which plasticity targets does no rule variant fit so far, which ones are hit only when other targets are given up, and what is the cause in each case?

**Short answer.**
- **Missed by every model.** Three targets are missed by all 27 rescored models (|z| >= 2): `sjostrom_40hz_dt0ms` (best 1.05 against 0.93 +- 0.04), Zilberter train10 +5 first (best 1.10 against 0.97) and Zilberter train10 -10 under AM251 (best 0.88 against 0.73).
- **Nearly always missed.** Two more targets are hit by a single model: `sjostrom_20hz_dt0ms` and `sjostrom_20hz_dt-25ms`.
- **The main trade-off.** Everything else can be hit by some model, but never all at once. There is one dominant axis: "L5 high-frequency and Markram LTP" against "LTD at 20 Hz and on L2/3".
  - 50 Hz -10 (1.70) is never hit together with 20 Hz -10 (0.65), its AM251 arm, Letzkus 3AP distal, or any Zilberter LTD row.
  - Markram +10 is never hit together with Zilberter post-only or train4.
- **By class:**
  - **Structural:** the eCB veto timing at 20-40 Hz, and the order and coincidence rows in Zilberter.
  - **Trade-off or parameter-limited:** 50 Hz -10, 20 Hz -10 and its AM251 arm, Markram +10, Letzkus 3AP proximal, Zilberter 5ap 10 Hz and train4.
  - **Input-limited (emodel):** Letzkus 1AP +10 and 3AP +10 distal.
  - **Data conflict under a uniform synapse-local rule:** Zilberter train -10 under AM251, Zilberter 1AP +/-10 (soft), and Sjostrom 10 Hz -10 against Markram 10 Hz -10.

## Data and method
- **Models.** 27 validation rescores, each scoring all 65 targets: s1A/B/C _s/_u, s1C_w, s2C/E/N _s/_u, s3C/N _s/_u, s3V_s, s3W_s, s4N_s/_u, s4V_s, s5N_s, g2t1/g2t2 _s/_u.
  - The files are RS/<model>_val.csv, _val_l23.csv and _val_l23l23.csv.
  - No r1*/r2* val csvs exist; the R1/R2 fits used the rejected 65-target objective.
- **Definitions.**
  - "Hit" means |z| < 2.
  - "Best pred" is the prediction with the smallest |z| over the 27 models.
  - "Joint hit" means a single model hits both targets.
- **Trade-off index.** This is the across-model correlation of the predictions, signed by the direction each target needs to move. Its sign is negative when moving one target towards its data moves the other away. It is computed over the 25 models without s1A, whose predictions sit near maximal LTP everywhere and would dominate.
- **Fitted sets.** s1*/g2t1 were fitted on 7 L5 targets; s2*/g2t2 on 11 L5; s3* on 20 cores; s4* on 35; s5N_s on all 64. s1C_w is the S1C rule fitted on 7 L5 + 23 L2/3 targets with L2/3 weight x0.25 (DECISIONS 2026-10-02).
- **Compute.** awk on the csvs only. Working tables are in the session scratchpad and are not kept.

**Per-pathway validation chi2, L5 / L2/3->L5 / L2/3->L2/3 (40 / 9 / 16 targets), with Markram +10 and -10 (data 1.20, 0.79):**

| model | L5 | L2/3->L5 | L2/3->L2/3 | Markram +10 | Markram -10 | note |
|---|---|---|---|---|---|---|
| s5N_s | 148 | 89 | 119 | 1.02 | 0.72 | best L5 |
| s3C_s | 202 | 59 | 158 | | | |
| s3W_s | 237 | 44 | 146 | 1.17 | 0.79 | |
| s3C_u | 291 | 30 | 114 | 1.14 | | |
| s4N_s | 313 | 24 | 107 | 1.00 | | |
| s4N_u | 566 | 54 | 72 | 1.07 | | best L2/3->L2/3 |
| s1C_w | 337 | 16 | 139 | 1.19 | 0.79 | |
| s1C_s | 268 | 354 | 336 | | | |
| s1A_s | 918 | 564 | 729 | | | |

- Only s1C_w and s3W_s keep Markram within 1 SEM while fitting L2/3.
- So L2/3->L5 can be fitted with Markram held (s1C_w 16, s1B_s 17). L2/3->L2/3 cannot: the best Markram-holding score is 139, against 72 when Markram is given up.

## Hard-target table
Hits are counted out of the 27 models. Data are mean +- SEM. Classes: S = structural, P = parameter-limited or trade-off, I = input-limited (emodel), D = data conflict under a uniform synapse-local rule.

| target | pathway | data | best pred (model); hits | class | reason (short) |
|---|---|---|---|---|---|
| sjostrom_40hz_dt0ms | L5 | 0.93 +- 0.04 | 1.05 (s3W_s); 0. s5N_s fitted it: 1.09 | S | Coincident bAP vetoes eCB; any Ca-amplitude rule contradicts 40/50 Hz +-10 |
| sjostrom_20hz_dt0ms | L5 | 0.76 +- 0.10 | 0.94 (s3W_s); 1. s5N fitted: 1.04 | S | Same veto / (1-b) cause |
| sjostrom_20hz_dt+25 / -25 | L5 | 0.65 +- 0.18 / 0.70 +- 0.07 | 0.85 (s3W_s) / 0.78 (s5N_s); 4 / 1 | S | +25 inside the 25 ms veto; -25 is the same train (edge artefact) |
| sjostrom_50hz_dt-10ms | L5 | 1.70 +- 0.19 | 1.61 (s1A_u, L2/3 not fitted). Best with L2/3 fitted 1.34 (s4N_u); 6. s5N fitted: 0.79 | P+S | Never jointly hit with 20 Hz -10, its AM251 arm or the L2/3 LTD rows. Not gamma_p- or licence-limited. eCB leaks at the end of the train |
| sjostrom_20hz_dt-10ms, control / AM251 | L5 | 0.65 +- 0.09 / 1.02 +- 0.07 | 0.70 / 0.99 (s5N_s); 3 / 7, all with Markram +10 lost | P | Knife edge: c* over the train gives lags {+40,-10} about the same as {+10,-40} |
| 10Hz_10ms (Markram) | L5 | 1.20 +- 0.06 | 1.20 in every L5-only fit. With L2/3 fitted: 1.19 (s1C_w), 1.17 (s3W_s); 4N/5N 1.00-1.07 | P | Global LTP/LTD axis against the L2/3 LTD rows (r 0.88-0.90) |
| sjostrom07_step200ms_pair, control / AM251 | L5 | 1.62 +- 0.07 / 2.13 +- 0.22 | 1.63 (s5N_s) / 1.97 (s4V_s); 8 / 11 | P+S | Licence cap (licence off gives 1.57); needs NO. The veto reads the step's sustained Ca. v11 3V/3W fall to 1.03-1.26 |
| sjostrom_10hz_dt-10ms | L5 | 0.57 +- 0.11 | 0.69 (s4V_s); 18 | D | Same protocol as Markram -10 (0.79 +- 0.026) |
| sjostrom_0.1hz_dt-10ms, nmdar_block | L5 | 1.07 +- 0.04 | 1.01 in all 27; 27 (z -1.5) | lane | The frozen-rho lane cannot leave the 1.01 baseline. Minor |
| letzkus_1ap_dt+10ms | L2/3->L5 | 0.72 +- 0.03 | 0.69-0.71 (s3C_u, s4N_u, s4V_s); 4. s5N fitted: 0.97 | I | Distal bAP 2-13 mV in the model vs 45-65 mV in the data. The target is pooled, the model sample distal-weighted. Never jointly hit with 3AP distal or Zilberter 1AP |
| letzkus_3ap +10 @distal | L2/3->L5 | 0.79 +- 0.06 | 0.77 (s1C_w); 6. Models that hit 1AP give 0.48-0.53 | I | Model LTD comes from burst Ca at near/basal "distal" sites, not from a weak distal bAP |
| letzkus_3ap +10 @proximal | L2/3->L5 | 1.30 +- 0.10 | 1.35 (s3C_u); 5. 4N/4V 0.94-0.96 | P | Trade-off with Markram +10 (r 0.90) and with the distal/1AP LTD level |
| zilberter train10 -10_last, AM251 | L2/3->L2/3 | 0.73 +- 0.07 | 0.88 (s3C_s); 0 | D | CB1-independent post LTD. No synapse-local cue separates it from the L5 AM251 nulls (C_MVD) |
| zilberter train10 +5 (first AP) | L2/3->L2/3 | 0.97 +- 0.06 | 1.10 (s4N_u); 0 | S | Order inversion against +4 last (1.49); a c*-threshold rule cannot separate them |
| zilberter_1ap_dt-10ms | L2/3->L2/3 | 0.56 +- 0.06 | 0.63 (s3W_s), 0.67 (s1C_w); 2 | D (soft) | Below the eCB ceiling 0.71; the extra rho-LTD conflicts with the L5 AM251 nulls |
| zilberter_1ap_dt+10ms | L2/3->L2/3 | 0.64 +- 0.07 | 0.62-0.73 (s3V_s, s3W_s, s1C_w, s1B_u); 4 | D (soft) | Same local state as L5 0.1 Hz +10 (0.97). Only hit with side effects |
| zilberter_5ap_10hz_dt+10ms | L2/3->L2/3 | 0.76 +- 0.07 | 0.78 (s3C_s, s4N_u); 8 | P/D | Same protocol class as L5 10 Hz +10 LTP. Trade-off with Markram +10 (r 0.88) |
| zilberter train10 post_only | L2/3->L2/3 | 1.03 +- 0.04 | 1.02-1.04 (s4N_s, s5N_s); 3. None holds Markram | S | The bAP train alone crosses theta_p. Never jointly hit with Markram +10 |
| zilberter train4 +4_last | L2/3->L2/3 | 0.76 +- 0.07 | 0.76 (s3V_s); 4 | P | AP-count threshold. Never jointly hit with Markram +10, train10 -10 or 1AP -10 |

**Not hard any more:**
- Letzkus 3AP -10 proximal (0.89; s3N_u 0.89, 8 hits) became parameter-limited and was fixed once the band opened (S1C_L23_FAILURE section 2).
- S&H 50 Hz +10 distal (0.86) is hit by 14 models. Most of the stage-3 fits get it through rho-LTD, whereas the paper's mechanism is CB1 (DISTAL_BAP_LIT).

## Groups by mechanism

### A. eCB veto window: 20 Hz dt 0 / +25 / -25, 40 Hz dt 0. Structural
- **Cause.** In S1C every own arrival has a bAP 0-25 ms after it, so the 25 ms veto blocks the eCB step. The (1 - b) weighting also keeps the trigger pool near zero. The prediction then equals the rho-only +10 value (S1C_L5_FAILURE section 3a).
- **No amplitude rule can fix it.** Under any rule monotone in spine Ca or c*, 40/50 Hz dt 0 contradicts 40/50 Hz +-10 (S1C_L5_FAILURE section 4). A timing veto (t0, Tv] makes the set self-consistent.
- **The fits confirm this.**
  - s5N_s fitted both dt-0 rows and still misses them: 1.09 (z +4.2) and 1.04 (z +2.8).
  - The v11 variants 3V/3W (veto 2-17 ms, unweighted trigger) move 40 Hz dt 0 to 1.05-1.16 and 20 Hz +25 to 0.85.
  - The trade-off index is about -0.7 between 40 Hz dt 0 and the 20/40 Hz +10 and 50 Hz -10 rows.
- **The veto fix moves the error rather than removing it.** s3W_s loses 50 Hz -10 (0.93), the step pair (1.03) and Sj07 pre-only (0.72 against 1.00, z -3.9); s3V_s also loses 0.1 Hz +10 (0.81). DECISIONS (L5 diagnostics) shows the same: veto 15 + v5_mode 1 drops 50 Hz -10 to 0.96 and the S07 pair to 1.08.
  - Inferred: the unweighted trigger fires on the last pre spike of each -10 train and during the depolarising step, where nothing vetoes it. This is the end-of-train leak named in S1C_L5_FAILURE 3b.
  - With A_eCB = 1, one such event gives the full d_min.
- **The 20 Hz -25 miss is an edge artefact.** In steady state it is the same train as +25 (S1C_L5_FAILURE 3a).
- **Bracket change (new).**
  - S1C_L5_FAILURE bounds Tv < 20 ms using 50 Hz dt 0. TARGET_AUDIT section 2 found that 50 Hz dt 0 is the same measurement as 40 Hz dt 0 and confirmed the drop.
  - So the only remaining upper bound is the +25 bAP of 40 Hz dt 0 (Tv < 25). Tv = 17 is still inside the bracket, but the bracket is wider than stated.
  - C_MVD_CALIB section 1 confirms the lower edge (t0 = 2 lies between the dt-0 and +4 current peaks).

### B. High-frequency LTP ceiling: 50 Hz -10, 40 Hz +10, Sj07 step pair (+ AM251). Trade-off plus structural
- **Licence cap.**
  - S1B and S1C cap LTP at 1.30-1.42 at >= 20 Hz, while S1A without the licence reaches 1.59-1.82 (S1C_L5_FAILURE 3c).
  - Licence off lifts the S07 pair to 1.57, so S07 is licence-capped (DECISIONS L5 diagnostics).
  - The G2 low-pass licence does not lift 50 Hz -10 (1.18 against 1.20 for G1; DECISIONS stage-3 entry).
- **50 Hz -10 is a different case.** It stays at about 1.2 under gamma_p x2, licence off and A_eCB 0.2 (DECISIONS L5 diagnostics), so it is not gamma_p- or licence-limited.
  - The csvs show it is a trade-off: index -0.87 with 20 Hz -10, never jointly hit with 20 Hz -10, 20 Hz -10 AM251, 20 Hz +25 or any Zilberter/Letzkus LTD row.
  - s5N_s, which fits all 64 targets, puts it at 0.79 (z -4.8). The only models that hit it are the L5-only or NO fits (s1A, s1B_s, s2C, s4N_u at 1.34).
  - Inferred mechanism: the rho part must stay neutral for 20 Hz -10 (lags +40/-10, AM251 1.02) and potentiate to 1.7 at 50 Hz -10 (lags +10/-10), while 20 Hz +10 (lags +10/-40) stays at 1.31. With c* integrated over the train, this is a knife edge on theta_p (S1C_L5_FAILURE 3e, section 4).
  - On top of this comes the end-of-train eCB leak in the -10 order (S1C_L5_FAILURE 3b).
- **Sj07 eCB depth (AM251 / control 0.76 in the data, 0.92 in the model) is structural.** The veto integrates VDCC influx, which the step supplies (S1C_L5_FAILURE 3d). NO models (2C, 4N, 4V, 5N) reach 1.63-1.74 on the pair, so that part is parameter-limited once NO is in.
- **Target value.** The 40 Hz +10 core value is wrong: it pools +-10, n 11. Use 1.54 +- 0.113 (n 6) (TARGET_AUDIT section 1). This makes the row slightly harder, not easier.

### C. Markram 10 Hz +10 against the L2/3 LTD rows. Parameter trade-off
- Every L5-only fit hits Markram. Every fit that also fits the L2/3 targets at full weight loses it: 3C_s 1.02, 3N 1.04, 4N 1.00-1.07, 4V 1.04, 5N 1.02 (DECISIONS stage-3, 4N/5N entries).
- Trade-off indices: -0.90 with Letzkus 3AP +10 proximal, -0.88 with Zilberter 5ap 10 Hz, -0.83 with Zilberter post-only, -0.79 with Zilberter train -10 AM251.
- It is one global axis. The L2/3 rows need less LTP for the same pairing counts, and the uniform thresholds respond by lowering L5 LTP.
- It is not absolute:
  - s1C_w holds Markram (1.19/0.79) with L2/3->L5 at 16;
  - s3W_s holds it (1.17/0.79) with Zilberter 1AP +/-10 hit.
  - Both pay on the L2/3->L2/3 trains (139-146) and on L5 high-frequency LTP.
- The weight-sweep fits show the frontier (DECISIONS: WL23 0.25 / 1 / 4). The Markram-weight-8 refits 22300522-27 will place it for the 4N/5N rules.

### D. Postsynaptic, CB1-independent LTD on L2/3->L2/3. Structural or data conflict
- **Zilberter train10 -10 AM251.** Missed by all 27; best 0.88, mean z +6.8.
  - The data say this LTD is postsynaptic, mGluR-dependent and CB1/NMDAR-independent (S1C_L23_FAILURE 2, 5.3).
  - The MVD cue search found no synapse-local quantity whose margin exceeds 0.21 (pool, low-Ca gate, mGluR low-pass, VDCC/Ca ratio, shaft Ca). L5 AM251 rows have equal or larger VDCC pools (C_MVD_CALIB section 2).
  - Under the current inputs this is a data conflict for a uniform synapse-local rule, and MVD stays off (DECISIONS MVD cue search).
- **Zilberter 1AP +/-10.**
  - The eCB ceiling is 0.71 (d_min -0.29), and the L5 AM251 nulls forbid extra rho-LTD at the same local state (S1C_L23_FAILURE 2).
  - It is soft, not absolute: s3W_s (0.69/0.63) and s1C_w (0.73/0.67) hit both while keeping 0.1 Hz -10 AM251 at 1.01. They pay elsewhere: s3W_s Sj07 pre-only 0.72 and 20 Hz -10 AM251 1.22; s1C_w train rows.
  - Neither ever hits Letzkus 1AP, Zilberter 5ap 10 Hz or post-only at the same time.
- **Train +5 first (0.97) against +4 last (1.49), and post-only (1.03). Structural.**
  - The +5-first train has the largest c* with an open licence, and the bAP train alone crosses theta_p (S1C_L23_FAILURE 2).
  - Telling them apart needs potentiation that reads own glutamate x bAP coincidence (Gb), which is untested in a neutral form.
  - Post-only is hit only by 4N_s/4V_s/5N_s, all of which lose Markram.
- **Train4, 5ap 10 Hz, train10 -10 control. Parameter trade-offs.**
  - These are the AP-count threshold placed for L2/3 against the L5 frequency axis.
  - Each is hit by 4-8 models, never together with Markram +10 (train4) or 50 Hz -10 (train -10).
- **Selection caveat.** Zilberter discarded low-release-probability connections; our 120 pairs are unselected (TARGET_AUDIT section 3). This affects every Zilberter row and is not modelled. No doublet bias here (120/120 kept; DOUBLET_BIAS).

### E. Distal bAP input on L2/3->L5: Letzkus 1AP +10, 3AP +10 distal. Input-limited
- **Model vs data.**
  - The model's distal bAP is 2-13 mV against 45-65 mV in the data.
  - Its distal/proximal single-bAP Ca ratio is 3e-4 to 3e-3 against an inferred 0.1-0.3.
  - split2 fails bap_620 in 30/30 cells (DISTAL_BAP_LIT, table).
- **Consequence.** Distal 1AP +10 is locally the EPSP-alone state (nopost 0.98). Letzkus' own reading is a moderate-NMDA-Ca band LTD that needs the weak bAP (DISTAL_BAP_LIT, Letzkus 2006 p.10427). No threshold choice can separate the two (S1C_L23_FAILURE 3.3).
- **What the csvs show.**
  - The models that hit 1AP +10 (0.69-0.71: s3C_u, s4N_u, s4V_s) over-depress 3AP distal (0.48-0.53).
  - s1C_w hits 3AP distal (0.77) and misses 1AP (0.65; z -2.3).
  - So the 1AP/3AP-distal/3AP-proximal pattern cannot be placed with a near-absent distal bAP.
  - The model's LTD comes from the near or basal part of the "distal" sample: 13 % of it is basal and a quarter lies below 281 um (S1C_L23_FAILURE 2).
- **Scoring and sample effects.**
  - Letzkus 1AP is pooled over locations, but our 120-pair sample weights distal 2:1 (TARGET_AUDIT section 3).
  - The 3AP rows rest on the 102 pairs that pass the doublet filter, and on cells that follow 200 Hz (DOUBLET_BIAS; DECISIONS T39). Bias size unknown.

### F. Other
- **Sjostrom 2001 10 Hz -10 against Markram 10 Hz -10 (0.57 +- 0.11 against 0.79 +- 0.026).** The same protocol in the model; the two values differ by 1.9 combined SEM. This gives an irreducible z^2 of about 3.7 (S1C_L5_FAILURE 3e). Data conflict; leave it as validation.
- **0.1 Hz -10 under NMDAR block (1.07).** Every model gives 1.01, because the lane freezes rho. Only a lane or readout change could move it. Not a rule problem.

## Direct data-conflict check
1. **Same protocol, different papers, same pathway (L5).**
   - Markram against Sjostrom at 10 Hz -10: conflict, 1.9 sigma.
   - Markram 10 Hz +10 / +5 (1.20 / 1.20) against Sj 10 Hz +10 (1.16 +- 0.09): consistent.
   - Markram 20/30/40 Hz +5 against Sj 20/40 Hz +10: consistent within 0.06.
   - Sj 20 Hz +25 against -25 (same steady-state train): consistent (0.65/0.70).
2. **Same pairing class across pathways, opposite sign** (uniform rule: the difference must come from local Ca):
   - (a) Single pairing +10: L5 Sj 0.1 Hz 0.97, Zilberter 0.64, Letzkus 0.72 pooled.
     - Letzkus can be explained by location, once the distal bAP exists (E).
     - Zilberter cannot, under current inputs (S1C_L23_FAILURE 2: same local state).
   - (b) Five pairings at 10 Hz, +10: L5 Markram 1.20 / Sj 1.16 against Zilberter 5ap 0.76.
     - Inferred from target ids only; repetition counts and inter-train intervals were not compared here.
     - Zilberter L2/3 contacts have larger VDCC pools and shaft Ca than L5 (C_MVD_CALIB section 2, cue e). So a Ca-monotone rule predicts the opposite order. This is the core of the axis in C.
   - (c) Post-pre train LTD under AM251: L5 0.1/20 Hz -10 AM251 1.01/1.02 and Markram -10 AM251 1.07, against Zilberter train -10 AM251 0.73. There is no separating local cue (C_MVD_CALIB). This is the hardest conflict.
3. **Within one dataset, rule-dependent.**
   - Sj 40 Hz dt 0 against 40/50 Hz +-10 contradict each other under any Ca-amplitude rule, but not under a timing veto (A).
   - Zilberter +5 first against +4 last contradict each other under any c*-threshold rule, but not under a coincidence-reading potentiation (D).
   - These are structural limits of the rule, not data conflicts.

## Consequences for the model (within the project rules)
- **Do not raise weights on the class-S or class-D rows.** Weighting moves error at constant total (STAGE3_DESIGN section 5.6). The csvs agree: s5N_s fitted 40 Hz dt 0, 50 Hz -10, Letzkus 1AP and Zilberter AM251, and still missed each by 3-8 SEM.
- **Report as structural misses with z:** 40/20 Hz dt 0, 20 Hz +/-25, Zilberter +5 first and post-only. The only fix is rule structure:
  - a timing-read eCB veto that does not leak at the train end;
  - a coincidence-reading potentiation.
  - Both need new kernel and mod names, and no new free parameter beyond the anchored edges.
- **Report as data conflicts for the user:**
  - Zilberter train -10 AM251;
  - Zilberter 1AP +/-10, which are soft and hit only at a cost elsewhere;
  - Sj 10 Hz -10.
  - The options are to accept the chi2 or demote them to validation. MVD cannot be anchored to a separating cue in the current Ca inputs.
- **Report as input-limited:** Letzkus 1AP +10 and 3AP +10 distal. The fix is an apical-bAP emodel variant (ion_fitter, new emodel name), not a rule parameter (DISTAL_BAP_LIT verdict). Until then, score 1AP on a 50/50 rise-time subsample, or treat both rows as validation.
- **Parameter trade-offs:** Markram +10, 50 Hz -10, 20 Hz -10 (+AM251), Letzkus 3AP proximal, Zilberter 5ap 10 Hz and train4.
  - These are the frontier the Chindemi parameters choose on.
  - Show the frontier, e.g. the Markram-weight-8 runs against WL23 0.25 (s1C_w), as fig9 plots, instead of searching for a single optimum.
  - 50 Hz -10 should be scored with its trade-off partner 20 Hz -10 AM251 visible.
- **Before adopting the v11 veto (3V/3W),** check the Sj07 pre-only (0.65-0.72) and 0.1 Hz +10 (3V 0.81) side effects. They suggest that the unweighted trigger fires without a pairing.

## Open points and compute that would settle them (nothing submitted)
1. **Is the s3W_s loss of 50 Hz -10 / S07 pair / pre-only caused by the unweighted trigger (v5_mode 1) or by the veto read?**
   - Run fixed-parameter L5 rescores of s3W_s.json with v5_mode 2, with A_eCB 0.2 (diagnostic only, not a model), and as is.
   - Size: one GPU 3g.40gb slice, 1 CPU, 47G, 0:15. Measured basis: L5-all 3 rescores, 7:31, 37.6 GB (DECISIONS v11 tests).
2. **Tv upper bound.** With 50 Hz dt 0 confirmed as a duplicate (TARGET_AUDIT), check whether Sjostrom 2001 Fig 7 has any other row with a bAP at +18 to +24 ms. Otherwise state the bracket as Tv in [15, 25).
3. **Zilberter 5ap 10 Hz against Markram 10 Hz +10.** Compare the protocols (pairs per train, trains, interval) from the papers before calling (b) a data conflict.
4. **Doublet effect on Letzkus 3AP distal.** Per-bin counts (DOUBLET_BIAS section 2, job 22291773) are still to be filled in.
5. **Pending refits.** The Markram-weight-8 refits (22300522-27) will show whether 4N/5N can hold Markram without losing the L2/3 hits listed above.
