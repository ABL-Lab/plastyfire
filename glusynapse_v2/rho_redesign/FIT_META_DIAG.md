# How we fit GluSynapse_v2: procedure, objective and pipeline problems (2026-10-02)

Question: apart from which targets are hard, what is wrong with the way the stage-3/4N/5N fits are set up, run, scored
and compared? Sources read: run_stage3_fit.sh, run_stage_fit.sh, gpu_v10_rho.py, gpu_v11_rho.py, and the fitter they
call (fit_v6.py, with fit_v3.rules_row, fit_v2.rules, gpu_v4_rho chi2/_ratios); logs s3C*, s3N*, s4N*, s5N*, v4N*,
v5N*; result csv/json of s1C_s s3C_s s3C_u s3N_s s3N_u s4N_s s4N_u s5N_s s1C_w; the design docs listed in the brief.
No compute was run. All numbers are from those files; lines marked "inferred" are my reading.

## Short answer

The largest problems are in the scorecard and the optimiser, not in the targets.
1. The comparison metric misleads. "L5 validation chi2" mixes in-sample and held-out targets, and one target decides it.
   For 5N, the L5 "validation" 148.03 is its in-sample fit chi2 (all 40 L5 targets were fitted). sjostrom_40hz_dt0ms
   (SEM 0.038) alone gives 156 of 4N_s's 313 and 420 of 4N_u's 566. Without that one row, every model scores 119-175.
   The model ranked best (5N_s) is also the flattest (slope 0.42, LTP under-predicted by 0.30 on average).
2. Seeding makes the DE converge early in a worse basin. 3C_s ended at objective 110.4 against 98.4 unseeded, and 4N_s
   at 144.8 against 133.0. At generation 1 only 3-16 of 56-72 members are admissible, and two of them are the seeds.
3. Three non-Chindemi parameters are free (theta_eCB, A_NO, d_NO_max). They end at their box edges or outside their own
   anchor. gamma_d sits on its floor and gamma_p near its ceiling. Those are constrained optima, not estimates.
4. The objective is a weighted chi2 over a binary readout. It is piecewise constant, so it is flat and degenerate.
   Structurally unreachable rows with tight SEMs steer the shared parameters. Weights move error around but cannot
   enforce "Markram must hit".
5. The pair sampling error (22 L5 pairs) and the hidden admissibility rules (L5-only, and they change with the target
   set) are not in the objective. The offline kernel has been checked against a live simulation on one L5 pair only,
   with the S1C rule (no NO arm, no L2/3 pathway).

## 1. Objective problems (ranked)

1a. The scorecard is in-sample, dominated by one target, and rewards flattening. (Highest impact.)
- The val branch rescores all targets, the fitted cores included (run_stage3_fit.sh:52-57). For 5N every L5 target is
  fitted (l.71, NT=40), so "5N_s L5 val 148" is the fit value: log s5Ns_22296409 "chi2 L5 148.03 over 40", repeated
  exactly in v5Ns_22296419:131 ("repro l5 ... diff 0.0"). 4N fits 11 of the 40 L5 targets, so the two numbers are not
  comparable, and DECISIONS.md l.4 compares them directly.
- Per-target z^2 from the *_val.csv files:

  | model | L5 val | 40 Hz dt0 share | L5 val without it |
  |---|---|---|---|
  | 4N_s | 313.4 | 155.8 | 157.6 |
  | 4N_u | 565.7 | 420.4 (pred 1.717 vs 0.929 +- 0.038) | 145.3 |
  | 5N_s | 148.0 | 17.6 | 130.5 |

  For the others, L5 val without this row is 119-175: 3N_s 119.3, 3N_u 119.9, s1C_s 135.8, 3C_s 142.8, 3C_u 144.9,
  s1C_w 175.1. So the reported "4N_s 313 vs 4N_u 566 basin difference" (brief) is almost all one row with an unfitted
  tight SEM. On the other 39 L5 rows, the unseeded 4N_u is better.
- Six L5 rows with SEM < 0.045 carry 34-75% of the L5 val chi2 (awk over the val csvs; 4N_u 75%, 5N_s 19%). Over the
  64 5N targets, the ten tightest carry 49% of the total precision (sum of 1/SEM^2). SEMs range 0.026-0.31, a 140x
  spread in weight.
- Flattening is not penalised. OLS slope of pred on data over all 65 val targets: 5N_s 0.42 (the lowest), s1C_s 0.43,
  3C_s 0.45, 3N 0.56-0.57, 4N_s 0.55, 4N_u 0.64 (the highest). Mean (pred - data) on the 20 LTP rows (> 1.1): 5N_s
  -0.30, 4N_u -0.11. STAGE3_DESIGN.md section 4 sets "slope >= 0.6" as a hard pass, but nobody computes it. The model
  with the best chi2 fails it, and the one with the worst L5 chi2 is the only one that passes.
- In 5N, 25 of the 64 targets are null-ish (|mean - 1| <= 0.1). They hold 46% of the unweighted precision (7461 of
  16142). This is the "pull to 1" the user rejected for R1D (DECISIONS.md l.37). The 0.05 SEM floor (l.74-75) eases it
  but does not remove it.
- The validation set is reused to choose among more than 20 variants (DECISIONS.md l.4-17), so it is no longer
  independent (inferred: selection bias).

1b. Structurally unreachable rows steer the shared parameters, and weights cannot enforce a must-hit.
- Weighted z^2 is a soft trade. Markram +10 at w 2 is 1 of 20-64 terms, and every multi-pathway fit lost it: 3C_s
  1.017, 3N 1.038, 4N_s 0.997, 4N_u 1.070, 5N_s 1.021 (data 1.201 +- 0.063). STAGE3_DESIGN.md section 5.6 and
  LTD_WEIGHTED.md already say that weighting moves error at constant total. The m8 refits (DECISIONS.md l.4) will
  trade again.
- Rows that no uniform rule here reaches still enter at full weight: 40/20 Hz dt0, the Zilberter AM251 train LTD and
  1AP +-10 (DECISIONS.md l.9: no synapse-local cue separates them), and the LTP ceiling (50 Hz -10, the step pair). Their
  residuals pull the shared a/gamma (section 4). Example: 5N fitted 40 Hz dt0 down to 1.09, but 40 Hz +10 fell to
  1.197, 50 Hz -10 to 0.79 and Markram to 1.02 (s5N_s_val.csv).

1c. The objective silently ignores NaN predictions and frozen lanes. Inferred priority: low-medium.
- In both the weighted wrapper (fit_v6.py:248) and the class chi2 (gpu_v4_rho.py:335), a target whose prediction is
  NaN adds 0, not a penalty. No NaN appears in the current csvs. v11 does report P1 NaN rates of 10-15% (DECISIONS.md
  l.15), so a v11 fit could drop a target without anyone seeing it.
- Frozen-rho lanes are constant in all 9 models: 0.1 Hz -10 nmdar_block 1.011 (z^2 2.2), 20 Hz -10 and 50 Hz +10
  nmdar_block 1.013, Letzkus nmdar_block 1.034. These 4 rows add a fixed chi2 and carry no information, yet they count
  in n for AIC/BIC (fit_v6.py:529).
- The baseline of a no-change prediction is 1.005-1.034, not 1.0, because of the Jensen term
  (1 + min((bs/bm)^2, 0.25)) in gpu_v4_rho.py:321. For SEM 0.04 nulls this sets a floor of z 0.3-0.9 (inferred).

1d. Signs are not in the objective. Sign agreement on the 49 non-null val targets is 36-43 for every model. Several
LTD rows predicted as LTP (Zilberter train +5, AM251 train -10) cost no more than a same-sign miss of the same size.

## 2. Optimiser problems

- Settings. scipy differential_evolution, strategy best1bin (default), popsize 8 (fit_v6.py:327 default; the stage
  scripts never set it), so 56 members for 7 free parameters and 72 for 9. maxiter 300 (run_stage3_fit.sh:82),
  tol 1e-6, polish False (fit_v6.py:492-494). Seeded runs draw a uniform random population and put the seeds in
  members 1 and 2. Unseeded runs use a Latin hypercube (fit_v6.py:468-476).
- Most of the initial population is penalised. fit_v3.rules_row (l.24-35) returns 1e3-1e4 when theta_d > theta_p at
  any L5 synapse, or when fewer than 30% / 15% of L5 synapse-records cross theta_d / theta_p. Generation-1 admissible
  counts: 3C_s 3 of 56, 3C_u 4, 3N_s 12 of 72, 3N_u 13, 4N_s 9, 4N_u 16, 5N_s 7. 5N_s did not improve for its first 4
  generations (500.18). With best1bin, the seeds are the base vector for nearly every trial at the start (inferred).
- Seeding hurts. The seeds are optima of other objectives (s1C: 7 L5 targets; s2C: 11), so they start in the wrong
  basin. s3C_s converged at generation 243 with f 110.45 (s3C_s_22291778.out:545; json nfev 244). The unseeded 3C_u
  reached 98.39 on the same objective. Likewise 4N_s 144.84 against 4N_u 133.04. 3N_s (99.08) and 3N_u (100.89)
  agree, and their parameters agree to within 7% (section 4). In 2 of 3 pairs, the CLAUDE.md rule "seeded beats
  unseeded" did not hold in stage 3.
- No run except 3C_s converged; they all stopped at maxiter, still improving. Last 50 generations: 3N_s 99.48 to
  99.08, 3N_u 101.08 to 100.89, 4N_s 144.89 to 144.84, 5N_s 284.08 to 283.59 (ckpt lines in the logs). "converged"
  never appears in a log. The json nfev (301) counts vectorised calls, i.e. generations, not evaluations.
- A generation costs 1.2-4.8 s (ckpt lines), so a fit costs 5-17 min (json "minutes"). A bigger search is cheap.
- One start per seed mode with one random seed is not a basin census. When s and u differ, the CLAUDE.md fix is to
  resume the better run, but that does not test multiplicity.

## 3. Model-evaluation problems

- Binary readout. rho_sigma defaults to 0 (rho_v4.py:29), so each synapse counts as potentiated only if rho >= 0.5
  (gpu_v4_rho.py step / ndtr branch, ~l.313-321). Each target prediction is a mean over about 22 L5 pair-records of
  EPSP ratios (gpu_v4_rho.py:330-333; L5 has 195 records and 1511 synapses over 9 protocols, about 7.7 synapses per
  record). The objective is therefore a staircase in parameter space. Inferred: one synapse flipping moves an L5
  target by roughly 0.005-0.02, i.e. 0.1-0.5 SEM on the tight rows. Flat steps give exact ties, which makes tol-based
  "convergence" a population collapse onto one step, not an optimum, and widens the degenerate sets (section 4).
- Pair sampling. L5 uses 22 pairs (subset24 intersected with l5_keep_pairs: 24 to 22, all logs "l5 used 22").
  L2/3->L5 uses 102 and L2/3->L2/3 120. The model prediction is deterministic, but it is a mean over a fixed draw of
  22 cells. Its draw-to-draw SE is not in the chi2, and the fit can tune to these 22 pairs. If the between-pair SD of
  the L5 ratio is about 0.2-0.3, the SE is 0.04-0.06, the same size as the tight L5 SEMs (0.026-0.04) (inferred,
  unmeasured; check B below).
- Hidden priors that depend on the target set. The admissibility rules use L5 synapses only (the G0 peak,
  fit_v6.py:206 and :95, fit_v2.py:24 MIN_ACTIVE 0.30 / MIN_POT 0.15). The peak array comes only from the L5 protocols
  of the targets kept (fit_v6.py:89-95). So 3C (9 L5 protocols, 195 records) and 5N (31 protocols, 670 records) search
  different admissible regions. The rules are not uniform across pathways and are not anchored in data.
- Offline vs BCL. The kernel integrates rho on frozen prefire Ca traces and reads out through a linear EPSP basis.
  Drug lanes are edits of the control run: nmdar_block freezes rho, and AM251/APV reuse the control VDCC trace
  (DECISIONS.md l.18). The one BCL check is the V9 pilot: 1 L5 pair, S1C rule, 39/40 targets, live = offline 19.66,
  20 Hz -25 0.913 vs 0.972 (DECISIONS.md l.7). The NO arm (no_mode 3) in 3N/4N/5N and both L2/3 pathways
  (vseg-rs) have no live check.
- Input consistency (emodel side). In every log, cexp and npz disagree for the L2/3 pathways: c_pre rel diff median
  0.67% / max 56% (L2/3->L5), 1.8% / 144% (L2/3->L2/3); c_post max 97% (L2/3->L5). The thresholds use npz c_pre/c_post,
  but the eCB scale uE uses cexp vdcc_q_post, so one synapse is described by two different measurements. None of the
  635 L2/3->L5 synapses has an isolated single-bAP shaft reference ("0 of 635", s5Ns log l.53). The distal bAP is
  under-driven, 2-13 mV vs 45-65 mV (DECISIONS.md l.22). The uniform a/gamma absorb any input error, which is one route
  for the L5 vs L2/3 trade-off (inferred).
- Doublets. The keep lists drop 2 of 24 L5 pairs and 18 of 120 L2/3->L5 pairs, the latter possibly from the
  strongest-Ca distal sites (DOUBLET_BIAS.md sections 2-3; its per-bin table was still pending). The fits still use
  these pair-level lists (run_stage3_fit.sh:29-31). Check that this matches the latest doublet decision ("doublet =
  one firing", PLASTYFIRE_DOUBLET_TOLERANT=2). If the lists predate that decision, the fitted pair set is out of date.

## 4. Identifiability (fitted json values; boxes: a 0-5, gamma_d 20-250, gamma_p 100-600, A_NO 32-1123, d_NO_max 0.4-1)

| run | a00 | a01 | a10 | a11 | gamma_d | gamma_p | theta_eCB | A_NO | d_NO_max |
|---|---|---|---|---|---|---|---|---|---|
| s1C_s | 2.28 | 0.67 | 2.31 | 1.27 | 20.8 | 165 | 0.437 | - | - |
| s1C_u | 1.95 | 0.011 | 2.50 | 0.97 | 24.0 | 585 | 0.442 | - | - |
| s3C_s | 0.94 | 2.98 | 0.96 | 4.95 | 20.0 | 589 | 0.146 | - | - |
| s3C_u | 0.89 | 2.78 | 1.05 | 3.21 | 97.2 | 599 | 0.263 | - | - |
| s3N_s | 0.88 | 3.35 | 0.95 | 4.87 | 28.7 | 480 | 0.525 | 1104 | 0.401 |
| s3N_u | 0.82 | 3.30 | 0.95 | 4.86 | 24.6 | 489 | 0.474 | 1114 | 0.402 |
| s4N_s | 0.80 | 3.02 | 1.59 | 3.29 | 20.0 | 581 | 0.526 | 1123 | 0.400 |
| s4N_u | 0.82 | 3.46 | 0.84 | 4.06 | 159 | 448 | 0.430 | 1122 | 0.400 |
| s5N_s | 1.37 | 1.83 | 2.68 | 1.69 | 24.6 | 303 | 0.174 | 33.4 | 0.76 |

- Well constrained within stage 3: a00 (0.80-0.94) and mostly a10. That fits PARAM_REDUCTION.md section 2, where
  a00 - a10 was the stiffest direction.
- Sloppy: gamma_d (20-159, x8, on its floor in 3 of 7 stage-3 runs), a11 (3.2-4.95, at the a box edge in 3C_s),
  theta_eCB (x3.6), and a01. gamma_d/gamma_p is 0.034 in 3C_s/4N_s against the literature 0.13-0.54
  (PARAM_REDUCTION.md section 1).
- Parameters on bounds: A_NO at its ceiling in 3N/4N (1104-1123) and at its floor in 5N (33). d_NO_max at its floor
  0.400 in every 3N/4N run. So the NO arm behaves as a switch that saturates at once (or is off), and the box is the
  active constraint. theta_eCB ends at 0.146-0.263 in 3C_s/3C_u/5N_s, outside the 0.37-0.78 bracket that EXP_PARAM.md
  section 4 gives as its anchor.
- Adding the L2/3 rows moves the c_post coefficients a lot: a01 goes from 0.01-0.67 (stage 1) to 2.8-3.5, and a11
  from 0.97-1.27 to 3.2-4.95. The pathway trade-off lives in the (a01, a11, gamma_d) directions (inferred).
- Rule check: only a00..a11, gamma_d and gamma_p may be free. theta_eCB, A_NO and d_NO_max are free in every stage-3
  model (run_stage3_fit.sh:36-40). theta_V 3.575 is a value fitted earlier and carried in SET. It is documented as
  unused once gate_theta > 0 (gpu_v10_rho.py doc, gate section); confirm it is inert in v5_mode 2.

## 5. Targets that may conflict under any uniform synapse-local rule (flags only; the other agent maps them)
- 40 Hz dt0 0.93 +- 0.038 and 20 Hz dt0 0.76 vs 40 Hz +-10 1.53/1.51 and 20 Hz +10 1.31. Every fit gives dt0 >= +10.
  This needs a timing-specific veto (the v11 window), not a Ca level.
- Zilberter 1AP +10 LTD 0.64 vs L5 0.1 Hz +10 null 0.97; Zilberter AM251 train LTD vs the L5 no-fire rows. No
  separating cue was found (DECISIONS.md l.9, l.13).
- 0.1 Hz -120/-200 null (1.05) vs burst r50 -120/-200 LTD (0.79).
- LTP ceiling: 50 Hz -10 (1.70) and the step pair (1.62) were never reached together with Markram.
- Not a conflict: s1C_w gets both Zilberter 1AP +10 (0.727) and Markram +10 (1.193), and 4N_u gets 2/5 Hz null with
  Markram +5 at 1.205. Those pairs are reachable.

## 6. Ranked fixes to the procedure

1. Fix the scorecard. Cost: minutes (awk on the csvs). Gain: high. Allowed. Report held-out-only chi2 per pathway and
   never mix in fitted rows. Report the 40 Hz dt0 row on its own line. Add robust metrics: median |z|, number within
   2 SEM, sign agreement, OLS slope and LTP deficit. Apply the STAGE3_DESIGN.md section 4 hard passes as gates before
   anyone looks at chi2. Freeze a holdout list that is never used to pick variants.
2. Multi-start DE with admissible initialisation. Cost: 3-5x GPU per fit (5-17 min becomes about 0.5-1.5 h). Gain:
   high. Allowed. Rejection-sample the initial population inside rules_row = 0. Use popsize 15 and maxiter 600. Run at
   least 4 starts: 2 seeded, plus 2 or more unseeded with seeds 6-9. Accept a basin only if 2 starts agree within
   delta f <= 2. Select by the fitted objective, not by validation. Consider rand1bin for the first ~100 generations.
3. Fix theta_eCB, A_NO and d_NO_max at anchored values (theta_eCB inside 0.37-0.78; NO size from the Sjostrom 2007
   1.62/1.36 split). This leaves 6 free (Chindemi only). Cost: one refit round. Gain: rule compliance and fewer
   degenerate or bound-hitting directions. Allowed and required by the rules. Keeping them free needs a user decision.
4. Encode "must hit" as a constraint. Use a hinge penalty for Markram +-10 (zero inside +-1 SEM, steep outside)
   instead of w 2/8, and add a slope/sign term only if the user wants it. Cost: a small change in a new fitter file.
   Gain: the user's hard pass is enforced, and the price elsewhere shows up explicitly. Allowed for Markram (it encodes
   a user rule). The sign/slope term needs a user decision.
5. Take structurally unreachable rows out of the objective until a mechanism exists that could reach them (dt0 rows
   before v11; the Zilberter AM251/1AP rows). Keep them as validation. Cost: none. Gain: stops irreducible residuals
   from steering a01/a11/gamma_d. Needs a user decision, because it changes the core list.
6. Replace the ad hoc SEM floor with a measured model-discrepancy term: SEM_eff^2 = SEM^2 + SE_pair^2, with SE_pair
   from the pair jackknife (check B). Cost: about 3 GPU-h. Gain: principled weights, and tight rows stop dominating.
   Allowed (independent calibration, not a fitted value).
7. Smoothed search. Anneal rho_sigma > 0 (the existing ndtr readout) during the DE, then rescore at rho_sigma 0. Cost:
   moderate (a fitter option in a new file). Gain: fewer plateaus and meaningful convergence. Allowed (procedure only;
   the model stays binary).
8. Make the admissibility rules pathway-uniform (all models' synapses, a fixed protocol set), or drop them and log
   whether they bind at the optimum. Cost: small. Allowed (log in DECISIONS.md).
9. Input hygiene before the next round. Reconcile cexp vs npz c_pre/c_post for the L2/3 pathways. Confirm the doublet
   keep lists match the latest doublet decision. Run a BCL pilot on 1 L2/3->L5 pair, 1 L2/3->L2/3 pair and 1 L5 pair
   with the NO arm. Cost: small to moderate. Allowed. Emodel fixes go to ion_fitter, under a new name.
10. More L5 pairs than the 24 of subset24 (prefire on more n120 pairs). Cost: heavy (>30 CPU-h, must be logged).
   Gain: lower L5 sampling error. Needs a user decision.

## 7. Checks that need compute (specified, not submitted; restate the Slurm sizing and no-login-node rules in any brief)
A. DE-population identifiability. Load the *_ckpt.npz of 3C_s/u, 3N_s/u, 4N_s/u and 5N_s (x, f). Unpack them to
   physical units, as fit_v6.unpack5 does (log params in log). For each run, count the members within delta f 2/4/10
   of the best, and give the per-parameter min-max and the correlation matrix. Basis: PARAM_REDUCTION job 22133878,
   39 MB, < 1 s. Size: 1 CPU, 256M, 0:15.
B. Pair jackknife for SE_pair. This is an L5-only rescore (no --extra), the 40 L5 targets, --maxiter 0, --seed-fits on
   each of s3N_s, s4N_u and s5N_s, with --pairs = the 22 pairs minus one, giving 22 x 3 = 66 rescores. Per target,
   compute SE_pair = sqrt((n-1)/n * sum (p_i - p_mean)^2) and the ratio SE_pair/SEM. Basis: "L5-all 3 rescores 7:31,
   37.6 GB" (DECISIONS.md l.15), so about 2.5 min each. Size: 3 jobs (one per json, 22 rescores serially), 3g.40gb MIG
   slice, 1 CPU, 47G, 1:30.
C. Basin census for 4N (after fixes 2-3). Six starts at popsize 15 and maxiter 500. Basis: 4N fits 23-26 min, max
   71.5 GB at popsize 8 / maxiter 300, so about 3.2x the work. Size: 3g.40gb, 1 CPU, 90G, 2:15 each. That is about
   13.5 GPU-h in total, which is not a heavy CPU run.
D. Admissibility binding. At each fitted json, print the rules_row margins: the fraction of L5 records with
   theta_d < peak and with theta_p < peak, and min(theta_p - theta_d). Add this to any rescore of A/B.

## Open points
- Whether the m8 refits (22300522-27) hold Markram at w 8, and at what cost in L2/3. That is a trade, not a fix (1b).
- The DOUBLET_BIAS.md per-bin table (job 22291773) is still unfilled.
- The audit fix for the sjostrom_40hz_dt+10ms target (1.54 +- 0.113) is deferred (DECISIONS.md l.21). All current
  fits use the double-counted 1.53 +- 0.14.
