# Stage 3 fit design: 20 core targets on all three pathways (2026-10-02)

Start: S1C (G8X shaft-Ca licence + eCB; free a00 a01 a10 a11 gamma_d gamma_p theta_eCB), chi2 10.3 over the 7 L5 shape
targets. Validation: L5 258 over 40, L2/3->L5 355 over 9, L2/3->L2/3 336 over 16.
User asks: fit L2/3->L2/3 too, put the worst misses into the fit, at most about 20 meaningful targets, a real dynamic
range, and Markram 10 Hz +/-10 ms must stay hit.
Why L2/3->L2/3 was not fitted until now: CORE_TARGETS.md held both L2/3 pathways back for stage 3, to fix the L5 rule
first. That stage is now, and L2/3->L2/3 enters with 6 targets.

## 0. Bug in the stage-2 validation of the NO models (fix before using S2C)

The 2C/2N SET contains `"A_NO": 190.0, "d_NO_max": 0.5`, and gpu_v10_rho.py `_pack10` starts a free v10 parameter
that is in SET at the SET value, overriding the seed json. So `2C|2N val` (`--maxiter 0`) scored 190 / 0.5, not the
fitted 649.4 / 0.405 (see s2C_s_val.json). The accepted REPRO DIFF hid this. S2C_s fit vs val csv: Markram -10
0.796 vs 0.752, +10 1.193 vs 1.146, step pair 1.673 vs 1.720 (S2C_u, S2N_s likewise; S1C, S2E agree exactly).
**The S2C/S2N validation totals (354/106/343 and the rest) are of a different parameter point.** Fix in the val branch:
`[ "$VAL" = val ] && SET=$(echo "$SET" | sed -E 's/, "A_NO": [0-9.]+, "d_NO_max": [0-9.]+//')`, rerun `2C s val`,
and check that the L5 cores match s2C_s.csv to 1e-6. A fit seeded from an NO json needs the same: no A_NO/d_NO_max in SET.

## 1. S1C validation failures, ranked by z^2

D = distinct timing/frequency/location/mechanism, R = redundant with a core; N = non-null, 0 = null or drug block near 1.

| z^2 | path | target | data +- SEM | S1C | info / effect | decision |
|---|---|---|---|---|---|---|
| 160 | 23-5 | letzkus_3ap_200hz_dt-10ms@proximal | 0.89 +- 0.03 | 1.27 | D post-pre burst / N (small) | core, SEM floor |
| 132 | L5 | sjostrom_40hz_dt0ms | 0.93 +- 0.038 | 1.37 | D / 0, tiny SEM | validation |
| 115 | 23-5 | letzkus_1ap_dt+10ms | 0.72 +- 0.03 | 1.04 | D / N | core, SEM floor |
| 74 | 23-23 | zilberter_5ap_10hz_dt+10ms | 0.76 +- 0.07 | 1.36 | D 5 pairs at 10 Hz / N | core |
| 67 | 23-23 | zilberter_train10_50hz_dt-10ms_last AM251 | 0.73 +- 0.07 | 1.30 | D CB1-independent LTD / N | core |
| 54 | 23-23 | zilberter_train10_50hz_dt+5ms | 0.97 +- 0.06 | 1.41 | D / 0 | validation, watch |
| 36 | 23-23 | zilberter_train10_50hz_post_only | 1.03 +- 0.04 | 1.27 | D bAP train alone / 0 | validation, watch |
| 35 | 23-5 | letzkus_3ap_200hz_dt+10ms@distal | 0.79 +- 0.06 | 1.14 | D location / N | core |
| 29 | L5 | sjostrom_20hz_dt0ms | 0.76 +- 0.10 | 1.30 | D dt 0 / N | validation, 1st L5 swap-in |
| 29 | 23-23 | zilberter_1ap_dt+10ms | 0.64 +- 0.07 | 1.01 | D / N | core |
| 29 | 23-23 | zilberter_train4_50hz_dt+4ms_last | 0.76 +- 0.07 | 1.13 | D AP-count threshold / N | validation, 1st L2/3 swap-in |
| 19 | 23-23 | zilberter_train10_50hz_dt-10ms_last | 0.72 +- 0.05 | 0.94 | D / N | core |
| 17 | 23-5 | sjostrom_50hz_dt+10ms@distal (S&H 2006) | 0.86 +- 0.09 | 1.23 | R Letzkus distal / N weak | validation |
| 17 | 23-5 | letzkus_3ap_200hz_dt+10ms@proximal | 1.30 +- 0.10 | 1.71 | D / N | core |
| 15 | L5 | sjostrom_20hz_dt-10ms AM251 | 1.02 +- 0.07 | 1.29 | D eCB-dependent L5 LTD / 0 drug | core, w 0.5 |
| 14 | L5 | sjostrom_20hz_dt-25ms | 0.70 +- 0.07 | 0.96 | R 20 Hz -10 / N | validation |
| 13 | L5 | sjostrom_20hz_dt+25ms | 0.65 +- 0.18 | 1.30 | D / N, large SEM | validation |
| 11 | L5 | sjostrom07_step200ms_pair | 1.62 +- 0.07 | 1.39 | D depolarising step / N | core (stage 2) |
| 10 | L5 | sjostrom_20hz_dt-10ms | 0.65 +- 0.09 | 0.94 | D 20 vs 50 Hz crossover / N | core |
| 9.5 | 23-23 | zilberter_1ap_dt-10ms | 0.56 +- 0.06 | 0.74 | D / N | core |

Below z^2 10, all validation: train8 (1.15, R), sj07 step AM251 (2.13, R), 5ap 20 Hz +10 (1.07, null), S&H 50 Hz +10 on
L2/3->L5 (1.06, null), burst r50 -200 (0.79), 10 Hz -10 (0.57). 50 Hz -10 (1.70) is already core and has never been hit.
Back to validation (still checked): `sjostrom_50hz_dt+10ms|nmdar_block` (nmdar_block freezes rho: 1.01 in every fit,
no information), `sjostrom_0.1hz_dt-10ms|mglu_block` (passes unfitted, z -0.1), `sjostrom07_step200ms_pair|no_block`
(equals the control without NO).
**Main finding: both L2/3 pathways fail on sign, not on flattening** (S1C predicts LTP for 6 LTD targets). On the 20
cores (L5 / L2/3->L5 / L2/3->L2/3) a flat model (all 1.0) scores 221 / 122 / 155 = 498, S1C 46 / 327 / 198 = 571,
S2E 428, S2C_s 415 (wrong A_NO). S1C is worse than flat on both L2/3 pathways.

## 2. Core list: 20 targets (18 non-null, 2 nulls)

L5->L5 (10): the 7 shape targets plus 3.
1. `10Hz_10ms|control` 1.20: Markram LTP, must-hit.
2. `10Hz_-10ms|control` 0.79: Markram LTD, must-hit.
3. `sjostrom_0.1hz_dt+10ms|control` 0.97: no LTP at low frequency (the shape null from stage 1).
4. `sjostrom_0.1hz_dt-10ms|control` 0.69: low-frequency timing LTD.
5. `sjostrom_20hz_dt+10ms|control` 1.31: LTP rises with frequency.
6. `sjostrom_40hz_dt+10ms|control` 1.53: LTP rises with frequency.
7. `sjostrom_50hz_dt-10ms|control` 1.70: high-frequency LTP in either order.
8. `sjostrom_20hz_dt-10ms|control` 0.65: post-pre is still LTD at 20 Hz; with 7, this sets the crossover frequency.
9. `sjostrom_20hz_dt-10ms|mglu_block` 1.02: that LTD is eCB-only, so post rho must stay silent (S1C 1.29). L5 half of the eCB test.
10. `sjostrom07_step200ms_pair|control` 1.62: LTP with a depolarising step (stage-2 core; pulls towards LTP).

L2/3->L5, Letzkus 2006 (4):
11. `letzkus_1ap_dt+10ms|control` 0.72: single-AP pre-post gives LTD at distal inputs.
12. `letzkus_3ap_200hz_dt+10ms@proximal|control` 1.30: burst LTP, proximal.
13. `letzkus_3ap_200hz_dt+10ms@distal|control` 0.79: the same burst at distal inputs gives LTD (location).
14. `letzkus_3ap_200hz_dt-10ms@proximal|control` 0.89: post-pre burst, mild LTD; the worst S1C miss (1.27).

L2/3->L2/3, Zilberter 2009 (6):
15. `zilberter_1ap_dt+10ms|control` 0.64: single-pairing pre-post LTD (anti-Hebbian).
16. `zilberter_1ap_dt-10ms|control` 0.56: single-pairing post-pre LTD (symmetric).
17. `zilberter_train10_50hz_dt+4ms_last|control` 1.49: train-LTP.
18. `zilberter_train10_50hz_dt-10ms_last|control` 0.72: train-LTD.
19. `zilberter_train10_50hz_dt-10ms_last|mglu_block` 0.73: train-LTD survives AM251; post rho must make it (L2/3 half of the eCB test).
20. `zilberter_5ap_10hz_dt+10ms|control` 0.76: 5 pairings at 10 Hz give LTD (S1C 1.36).

3N only (21): add `sjostrom07_step200ms_pair|no_block` 1.36. A mechanism's target enters with that mechanism; without
it, A_NO and d_NO_max are constrained only by the step pair.
Swap-in order if a slot opens: `sjostrom_20hz_dt0ms`, `zilberter_train4_50hz_dt+4ms_last`,
`zilberter_train10_50hz_post_only`. Egger stays validation, and the DROPT pair (Letzkus 3AP -10 distal) stays dropped.

**Keeping the nulls from dominating** (2 of 20 are nulls). The only option is fit_v6 `--weights` (json of z^2 weights,
DE objective only, keys `target|cond`, `NAME/target|cond`, `ltd`/`ltp`; reported chi2 unweighted). gpu_v10_rho.py and
fit_v6.py have no SEM floor, per-group normalisation or dynamic-range penalty, so per-target weights emulate them:
- SEM floor 0.05 as w = (SEM/0.05)^2: 0.1 Hz +10 w 0.64, Letzkus 1AP w 0.36, Letzkus 3AP -10 proximal w 0.36.
- The null drug arm (20 Hz -10 AM251) gets w 0.5.
- Markram is exempt from the floor. Markram +10 gets w 2 (effective SEM 0.044), because the 65-target fit lost it;
  Markram -10 already has the tightest SEM (0.026).
- No per-pathway normalisation in the primary runs: the 10/4/6 split is a design choice, and normalising would
  multiply the tiny-SEM Letzkus rows by 1.7.
- Dynamic range is monitored, not penalised (awk on the csvs): the OLS slope of pred on data over the 20 cores
  (flat = 0; S1C 0.45), and sign agreement on the 18 non-nulls (S1C 12 of 18).

`WEIGHTS='{"l5/10Hz_10ms|control": 2.0, "l5/sjostrom_0.1hz_dt+10ms|control": 0.64, "l5/sjostrom_20hz_dt-10ms|mglu_block": 0.5, "l23/letzkus_1ap_dt+10ms|control": 0.36, "l23/letzkus_3ap_200hz_dt-10ms@proximal|control": 0.36, "l23l23/zilberter_1ap_dt+10ms|control": 0.5, "l23l23/zilberter_1ap_dt-10ms|control": 0.5}'` (the two Zilberter 0.5 weights were added 2026-10-02, orchestrator decision)

## 3. Fit procedure

All runs cover the 3 pathways; the queue runs one GPU fit at a time.

| run | rule | free parameters | seeding |
|---|---|---|---|
| **3C_s** (primary) | S1C | a00 a01 a10 a11 gamma_d gamma_p theta_eCB (7) | members s1C_s.json and s2E_s.json (same rule, refit with the step pair); `--seed-fits` takes a comma list, one member per json |
| **3C_u** | S1C | same 7 | seed 6, no json (basin check) |
| **3N_s** | S1C + NO mode 3 | the 7 + A_NO [32, 1123] and d_NO_max [0.4, 1.0] (9; boxes from C-NO2 / Padamsey 2017; tau_NO 6.7 ms fixed, Hall & Garthwaite) | members s2C_s.json (A_NO 649; d_NO_max 0.405, at its bound) and s1C_s.json (A_NO 0, so it starts at 190 / 0.5) |

**Seed: I recommend S1C for the primary run.** It is the user's choice, and it is the only candidate whose validation
is verified. S2C_s hits Markram in its own fit csv (0.796 / 1.193), but its validation was scored at the wrong A_NO
(section 0), so its L2/3->L5 gain (106 vs 355) does not count until the corrected rescore. S2C_s seeds the NO variant.
maxiter is 300 at popsize 8, as in R1/R2 (27 min measured on 3 pathways). If the s and u runs differ by more than 10%,
resume the better one to maxiter 500 (91G 1:15:00).

run_stage_fit.sh changes, for the technician (described here, not edited):
1. **Cases.** `3C) SET="{$G8X}"; FREE=theta_eCB; SEEDJ=$RS/s1C_s.json,$RS/s2E_s.json ;;` and
   `3N) SET="{$G8X, \"no_mode\": 3, \"tau_NO\": 6.7}"; FREE=theta_eCB,A_NO,d_NO_max; SEEDJ=$RS/s2C_s.json,$RS/s1C_s.json ;;`.
   A_NO and d_NO_max stay out of the 3N SET (section 0). Update the `*)` usage message.
2. **Cores.** `3*) NAME=s$M; NT=20; KEEP2=",sjostrom_20hz_dt-10ms|control,sjostrom_20hz_dt-10ms|mglu_block,sjostrom07_step200ms_pair|control"`;
   for 3N, NT=21 and add `,sjostrom07_step200ms_pair|no_block`.
3. **Per-model drop lists.** An unprefixed `--drop-targets` entry applies to every model (fit_v6 `drop_for`), and
   `sjostrom_50hz_dt+10ms|control` exists in both l5 and l23. So the existing awk on r1D_s.csv prints `l5/<target>`
   (harmless in stages 1-2, where l5 is the only model). The same awk prints `l23/...` with KEEP23 (cores 11-14) on
   s1C_s_val_l23.csv, and `l23l23/...` with KEEP2323 (cores 15-20) on s1C_s_val_l23l23.csv, which also drops Egger.
   Both csvs have the 6 columns the awk checks.
4. **Fit-branch arguments for 3\*.** Add the val branch's `--l23-dirs $L23D`,
   `--extra "l23:paired_l23l5:$L23D:$L23_BASIS_DIR:$GEOM_L23:$K235"` and
   `--extra "l23l23:paired_l23l23,paired_l23l23_egger:$Z:$W/basis_l23l23::$K2323"`. These are the specs of
   run_r1_fit.sh / run_r2_fit.sh, where the 6th-field keep list replaces --joint. Also add
   `--drop-targets "$DROPT,$DROP,$DROP23,$DROP2323" --weights "$WEIGHTS"`. Seeding stays `--seed 5 --seed-fits $SEEDJ
   --seed-set '{}'`, or `--seed 6` for u runs, with no --x0.
5. **Val branch:** only the section 0 SET strip (stays unweighted, all targets). **Header:** sizing below; grep adds
   `model l23|model l23l23|weights`; the log must show 4 / 6 targets for l23 / l23l23 and `over 20 targets` (21 for 3N).

Commands (from plastyfire/, after the edits; sizes are measured):
```
S=glusynapse_v2/sjob.sh; F=glusynapse_v2/rho_redesign/run_stage_fit.sh
BF='MEASURED R1 3-pathway fits 22235493-502 max 72.9 GB 27:23'; BV='MEASURED stage-2 rescores, 84G 0:15'
$S -g -n s3C_s -m 91G -t 00:45:00 -b "$BF" -- "bash $F 3C s"
$S -g -n s3C_u -m 91G -t 00:45:00 -b "$BF" -- "bash $F 3C u"
$S -g -n s3N_s -m 91G -t 00:45:00 -b "$BF" -- "bash $F 3N s"
$S -g -n s3C_s_val -m 84G -t 00:15:00 -b "$BV" -d afterok:<s3C_s id> -- "bash $F 3C s val"   # same for 3C u, 3N s
$S -g -n s2C_s_val2 -m 84G -t 00:15:00 -b "$BV" -- "bash $F 2C s val"                        # corrected S2C rescore
```
After each job, run seff and add a MEASURED line. The fits drop 45 of the 65 targets, so MaxRSS may fall below 72.9 GB;
resize from the seff.

## 4. Acceptance (unweighted chi2 from the fit csv; the rest from the val rescore)

**Hard** (a run that fails any of these is rejected): Markram +10 and -10 both within 1 SEM; at least 16 of the 18
non-null cores on the correct side of 1; OLS slope of pred on data over the cores at least 0.6; prediction range at
least 0.7 (the data span 0.56-1.70).

| pathway (cores) | chi2 limit | conditions |
|---|---|---|
| L5 | <= 25 over 10 | 0.1 Hz +10 < 1.05; 0.1 Hz -10 and 20 Hz -10 < 0.85; 20 and 40 Hz +10 within 1.5 SEM; 20 Hz -10 AM251 within 1.5 SEM; 50 Hz -10 within 2 SEM (soft: this rule has never passed 1.4) |
| L2/3->L5 | <= 25 over 4 | 1AP +10 < 0.85; 3AP distal < 0.95; 3AP proximal > 1.15; 3AP -10 proximal < 1.0 |
| L2/3->L2/3 | <= 30 over 6 | both single pairings < 0.85; train +4 > 1.30; train -10 (both arms) < 0.85; 5ap 10 Hz < 0.9 |

**Validation** (not fitted): L5 at most 258 (the S1C value) over 40, with the stage-2 mechanism rows within 2 SEM; each
L2/3 pathway below its flat-model chi2; Zilberter post_only and train +5 both < 1.15.

## 5. Decision tree

1. **Section 0 first.** Rescore S2C_s correctly. If its corrected L5 validation is still worse than S1C's, or it misses
   Markram, 3N stays a secondary run.
2. **3C_s passes and 3C_u agrees within 10%.** 3C becomes the stage-3 rule. Compare with 3N by BIC on the 20 shared
   cores (k 7 vs 9); NO stays only if it buys more than 2 ln 20 = 6 chi2. Then the val rescore, and the BCL last.
3. **L5 passes, L2/3 fails on sign** (LTP where the data show LTD): the uniform rule cannot place theta_d/p for both
   cell types from the synapse's own c_pre/c_post. In order:
   a. Check the inputs first: are L2/3 EPSP and bAP Ca in range? Beyond 130 um, bAP Ca is about 3x too low (T32). If
      that is the cause, it goes to the emodel/ion_fitter side, not to a new fit parameter.
   b. If 3N fixes L2/3->L5 without breaking L5, take 3N.
   c. Sensitivity runs, not a model: refit the best seed with the L2/3 weights x0.25 and x4 (2 x 91G 0:45), and show
      the user the L5/L2/3 trade-off as fig9-style plots. Pathway-specific parameters are not an option.
   d. If both Zilberter train-LTD arms still fail: the data say this LTD is postsynaptic, mGluR-dependent, and CB1- and
      NMDAR-independent (L23L23_TARGETS.md), but the model's post LTD is NMDA-Ca driven. Propose one anchored mechanism
      (postsynaptic group-I mGluR LTD on own pre spikes, anchored to Zilberter's CPCCOEt arm, not a target) and ask the
      user first.
4. **L2/3 passes, L5 shape degrades** (Markram lost, or slope < 0.6): reject. Rerun with Markram w 4 and the L2/3
   weights x0.5; if L5 still fails, go to 3c and give the user the choice.
5. **Both fail and DE barely moved** (nfev at the cap, or s and u differ): resume to maxiter 500 before any model change.
6. **A core sits at a structural ceiling** (50 Hz -10 about 1.2 and step pair about 1.41 without NO): report it as a
   miss with its z, and do not raise its weight. Weighting moves error between targets at constant total
   (LTD_WEIGHTED.md).
