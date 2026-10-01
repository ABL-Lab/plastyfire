# GluSynapse_v2 fit targets: difficulty, discrimination, cost (2026-09-30)

Fits used. Current rule (cur): td4 s1-s3, td3 s1, ljp25 td4 s2/s3, v4 td4g s2. A0: v4 A0 s2/s3 and A0g s3/s4. B: v4 B s2/s3.
Joint (J): v4 A0 joint s3/s4 (logs 22114141/2). Rule A (s2/s3) is counted only in med and spread. D (L5 chi2 ~313, rule broken) is counted only in "+D".
L2/3 rows: the transfer evals (l23_*.csv) for cur, A0 and B; for J they were fitted.
- **Per-family |z|:** the median over that family's fits.
- **med:** the median |z| over all non-D fits.
- **spread:** (max - min pred over the non-D fits) / SEM. "+D" adds the two D fits.
- **cost:** the protocol's share of synapse-steps (GPU kernel work; out/target_cost.csv, job 22115351). It is given as % of the L5-only fit (805 M steps) / % of the joint fit (1463 M = 805 L5 + 658 L2/3). Extra conditions on a protocol (mglu_block, no_block: the same kernel thread carries their dpre; nmdar_block: readout only) add no GPU work, so they are marked "+0".
- **Records:** 24 records (191 synapses) per L5 protocol (Sj07 pair and pre_only: 23). Trace length per synapse: 0.1 Hz ~300-338k steps, burst r50 361-376k, Sj07 200k, Markram 168k, 10-50 Hz 95-114k. L2/3: Letzkus 401k steps with 110/55/58 records, 50 Hz 95k steps with 111 records.

## 1. All 38 targets

| # | target | cond | data | cur | A0 | B | J | med | spread / +D | cost % | tag |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 10Hz_-10ms (Markram) | ctrl | 0.79+-0.03 | 1.3 | 0.7 | 1.4 | 1.0 | 1.1 | 2.0 / 2.3 | 4.0 / 2.2 | core |
| 2 | 10Hz_5ms | ctrl | 1.20+-0.06 | 0.3 | 0.9 | 0.4 | 2.0 | 0.4 | 2.7 / 2.7 | 4.0 / 2.2 | validation |
| 3 | 10Hz_10ms | ctrl | 1.20+-0.06 | 1.0 | 0.6 | 1.4 | 3.1 | 1.1 | 2.6 / 2.6 | 4.0 / 2.2 | core |
| 4 | sj 0.1hz -10 | mglu_block | 1.01+-0.04 | 0.1 | 0.2 | 0.2 | 0.1 | 0.2 | 1.2 / 1.2 | +0 | core |
| 5 | sj 20hz -10 | mglu_block | 1.02+-0.07 | 1.6 | 0.9 | 0.8 | 0.1 | 1.2 | 2.9 / 2.9 | +0 | core |
| 6 | sj 0.1hz -10 | nmdar_block | 1.07+-0.04 | 1.4 | 1.4 | 1.4 | 1.4 | 1.4 | 0 / 0 | +0 | redundant |
| 7 | sj 20hz -10 | nmdar_block | 1.02+-0.07 | 0.1 | 0.1 | 0.1 | 0.1 | 0.1 | 0 / 0 | +0 | redundant |
| 8 | sj 50hz +10 | mglu_block | 1.59+-0.19 | 0.1 | 0.1 | 0.1 | 1.0 | 0.1 | 1.4 / 2.4 | +0 | core |
| 9 | sj 50hz +10 | nmdar_block | 1.04+-0.05 | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 | 0 / 0 | +0 | redundant |
| 10 | sj 0.1hz -25 | ctrl | 0.65+-0.09 | 1.5 | 1.3 | 1.6 | 1.8 | 1.5 | 1.0 / 2.9 | 7.2 / 4.0 | validation |
| 11 | sj 0.1hz -120 | ctrl | 1.05+-0.10 | 0.4 | 0.4 | 0.7 | 0.8 | 0.4 | 1.9 / 1.9 | 7.7 / 4.2 | validation |
| 12 | sj 0.1hz -200 | ctrl | 1.05+-0.10 | 0.4 | 0.4 | 0.4 | 0.4 | 0.4 | 1.9 / 1.9 | 8.0 / 4.4 | redundant (same data as 11) |
| 13 | burst5x20 r50 -120 | ctrl | 0.79+-0.08 | 0.1 | 0.3 | 0.1 | 0.5 | 0.3 | 2.4 / 3.0 | 8.6 / 4.7 | redundant (same data as 14) |
| 14 | burst5x20 r50 -200 | ctrl | 0.79+-0.08 | 0.4 | 0.4 | 0.5 | 0.5 | 0.4 | 3.0 / 3.7 | 8.9 / 4.9 | core |
| 15 | sj 0.1hz +10 | ctrl | 0.97+-0.04 | 0.8 | 2.5 | 1.5 | 0.4 | 1.1 | 3.8 / 3.8 | 7.2 / 3.9 | core |
| 16 | sj 0.1hz -10 | ctrl | 0.69+-0.07 | 1.5 | 1.3 | 1.7 | 4.6 | 1.5 | 3.6 / 3.6 | 7.2 / 3.9 | core |
| 17 | sj 10hz +10 | ctrl | 1.16+-0.09 | 0.5 | 0.3 | 0.6 | 1.4 | 0.5 | 2.2 / 2.2 | 2.7 / 1.5 | validation |
| 18 | sj 10hz -10 | ctrl | 0.57+-0.11 | 1.7 | 1.8 | 1.5 | 1.4 | 1.7 | 1.4 / 1.4 | 2.7 / 1.5 | validation |
| 19 | sj 20hz +10 | ctrl | 1.31+-0.14 | 0.1 | 0.3 | 0.2 | 0.7 | 0.2 | 1.4 / 1.8 | 2.4 / 1.3 | validation |
| 20 | sj 20hz -10 | ctrl | 0.65+-0.09 | 2.6 | 1.9 | 2.1 | 2.2 | 2.2 | 1.4 / 1.4 | 2.4 / 1.3 | core |
| 21 | sj 40hz +10 | ctrl | 1.53+-0.14 | 0.2 | 0.2 | 0.2 | 0.8 | 0.2 | 1.4 / 1.7 | 2.3 / 1.3 | validation |
| 22 | sj 40hz -10 | ctrl | 1.51+-0.31 | 1.4 | 1.4 | 1.6 | 1.4 | 1.4 | 0.7 / 1.2 | 2.3 / 1.3 | validation |
| 23 | sj 50hz +10 | ctrl | 1.57+-0.26 | 0.1 | 0.1 | 0.1 | 0.6 | 0.1 | 1.0 / 1.7 | 2.3 / 1.2 | core |
| 24 | sj 50hz -10 | ctrl | 1.70+-0.19 | 3.0 | 3.0 | 3.4 | 3.2 | 3.1 | 1.1 / 3.4 | 2.3 / 1.2 | core |
| 25 | Sj07 step200 pair | ctrl | 1.62+-0.07 | 0.5 | 0.5 | 0.4 | 1.2 | 0.5 | 2.8 / 11.6 | 4.6 / 2.5 | core |
| 26 | Sj07 step200 pair | mglu_block | 2.13+-0.22 | 0.4 | 0.3 | 1.0 | 0.3 | 0.4 | 1.9 / 6.9 | +0 | core |
| 27 | Sj07 step200 pair | no_block | 1.36+-0.07 | 0.4 | 1.1 | 0.1 | 1.2 | 0.6 | 2.5 / 12.4 | +0 | core |
| 28 | Sj07 pre_only | ctrl | 1.00+-0.07 | 1.1 | 1.0 | 0.7 | 1.6 | 0.9 | 3.3 / 4.3 | 4.6 / 2.5 | core |
| 29 | Sj07 post_only | ctrl | 1.00+-0.07 | 0.4 | 0.6 | 0.5 | 0.1 | 0.6 | 1.7 / 1.7 | 4.8 / 2.6 | validation |
| 30 | L23 Letzkus 1AP +10 | ctrl | 0.72+-0.03 | 14.5 | 11.0 | 11.3 | 4.6 | 11.0 | 10.7 / 10.7 | - / 20.1 | core |
| 31 | L23 3AP200 +10 distal | ctrl | 0.79+-0.06 | 9.9 | 6.0 | 8.3 | 3.4 | 6.3 | 7.4 / 9.7 | - / 9.9 | validation |
| 32 | L23 3AP200 +10 distal | nmdar_block | 1.03+-0.07 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0 / 0 | +0 | redundant |
| 33 | L23 3AP200 +10 proximal | ctrl | 1.30+-0.10 | 0.3 | 1.3 | 1.1 | 0.9 | 1.0 | 3.0 / 9.5 | (with 31) | validation |
| 34 | L23 3AP200 -10 distal | ctrl | 1.42+-0.09 | 6.1 | 5.6 | 5.3 | 4.6 | 5.5 | 2.1 / 2.1 | - / 10.2 | core (structural miss) |
| 35 | L23 3AP200 -10 distal | nmdar_block | 0.97+-0.06 | 1.1 | 1.1 | 1.1 | 1.1 | 1.1 | 0 / 0 | +0 | redundant |
| 36 | L23 3AP200 -10 proximal | ctrl | 0.89+-0.03 | 5.6 | 1.2 | 7.4 | 0.8 | 2.7 | 9.2 / 9.5 | (with 34) | core |
| 37 | L23 sj 50hz +10 | ctrl | 1.06+-0.09 | 6.0 | 4.2 | 6.3 | 2.5 | 5.3 | 5.1 / 6.4 | - / 4.8 | core |
| 38 | L23 sj 50hz +10 distal | ctrl | 0.86+-0.09 | 8.8 | 5.4 | 8.3 | 2.8 | 6.4 | 6.6 / 6.7 | (with 37) | core |

Reading the table:
- **Hard targets (med >= 1.5):** 50 Hz -10, 20 Hz -10, 10 Hz -10, 0.1 Hz -10, 0.1 Hz -25, and every L2/3 control target except 3AP +10 proximal.
- **Discriminating targets (spread >= 2.5):** all L2/3 control targets, 0.1 Hz +10 and -10, Sj07 pre_only and pair, burst -200, Markram 10 Hz +5/+10, and 20 Hz -10 mglu_block.
- **Hard but flat:** 50 Hz -10, 20 Hz -10, 10 Hz -10, 0.1 Hz -25 and 40 Hz -10 are missed by every non-D rule. They measure rule structure rather than separate parameter sets.
- **nmdar_block targets (6, 7, 9, 32, 35):** the prediction is a function of rho0 only, because rho = rho0 and dpre = 0. It is identical in every fit (L5 1.013; L2/3 1.032 and 1.036), so these 5 targets add a constant 2.33 (L5) + 1.21 (L2/3) to chi2. They can be dropped from the objective with no change to the optimum (`--conditions control,mglu_block,post_nmdar,no_block`); only the reported chi2 shifts. Dropping them saves no GPU work, because their protocols are kept for the control targets.
- **Duplicate data:** 0.1 Hz -120/-200 share the same data (1.05+-0.10), and so do burst r50 -120/-200 (0.79+-0.08). In each pair, one target is enough.

## 2. Core set: 10 L5 protocols (15 informative targets) + 3 L2/3 protocols (5 targets)

**L5:**
- 0.1 Hz +10 (15)
- 0.1 Hz -10 (16, 4 [+6])
- 20 Hz -10 (20, 5 [+7])
- 50 Hz -10 (24)
- 50 Hz +10 (23, 8 [+9])
- Markram 10 Hz -10 and +10 (1, 3)
- burst r50 -200 (14)
- Sj07 pair (25, 26, 27)
- Sj07 pre_only (28)

**L2/3 (joint only):**
- Letzkus 1AP +10 (30)
- 3AP200 -10 distal and proximal (34, 36 [+35])
- 50 Hz +10 and +10 distal (37, 38)

**Coverage:**
- **Pathways:**
  - eCB-LTD: 16, 20, 1, 14.
  - NO-LTP: 25, 27, with 28 as the pre-alone control.
  - Post rho LTP: 24, 23, 3, 26.
  - Post rho LTD: 30 (dpre = 0, so it is pure rho LTD), with 15 and 5 as rho must-not-move constraints.
- **Pharmacology arms:** mglu_block 4, 5, 8, 26 constrain A_mglu; no_block 27 constrains A_NO; nmdar_block 6, 9, 35 are free to keep, but carry no information.
- **Frequency dependence:** -10 at 0.1/20/50 Hz (+ Markram 10 Hz); +10 at 0.1/10/50 Hz.
- **Timing:** ±10 at 0.1, 10 and 50 Hz, plus -200 in bursts.
- **Sjostrom 2007:** 25-28.
- **L2/3 distal targets:** 34, 38.

| parameter | role | constraining core targets | dropped targets that also constrain it |
|---|---|---|---|
| a00, a01 | theta_d = a00 c_pre + a01 c_post (post LTD threshold) | 30, 36, 15, 5, 28 | 29, 11, 33 |
| a10, a11 | theta_p (post LTP threshold) | 24, 23, 3, 26, 15, 37, 38 | 17, 19, 21, 22, 2, 31, 33 |
| gamma_p | rho LTP rate | 24, 23, 3, 25, 37 | 17, 19, 21 |
| gamma_d | rho LTD rate | 30, 36, 15 | 31, 33 |
| A_mglu | eCB-LTD strength (mglu_block arms) | 4, 5, 8, 26 vs 16, 20, 25 | - |
| theta_Te, tau_T | T integrator of post impulses (eCB timing window) | 16, 15, 14, 1, 20 | 10, 11, 12, 13 |
| theta_Tg | T gate on eCB-LTD (burst vs single pairing) | 14, 16, 20, 24 | 11, 18, 22 |
| dpre_min | eCB-LTD floor (depth) | 1, 16, 20 | 18, 10 |
| A_NO | NO-LTP strength | 27 vs 25 | - |
| theta_NOi (= theta_Ti), tau_NO | NO drive N (pre VDCC) | 25, 27, 28, 23, 37 | 29, 21 |
| tau_Z, theta_Z | presynaptic count Z (frequency of NO-LTP) | 25, 28, 23, 3 | 17, 19, 21, 29 |

## 3. Speed-up and risk

- **L5-only fit:** the core keeps 381.5 M of 805.2 M synapse-steps (47%), about 2.1x less kernel work. The ~2 min load also shrinks (10 of 21 protocols, 238 of 502 npz). Expected time is about 7 min instead of 12, and GPU trace memory drops by about half.
- **Joint fit:** the core keeps 895 of 1463 M (61%), about 1.6x. L2/3 Letzkus 1AP alone is 20% of the joint fit (110 records x 401k steps). Moving 3AP +10 to core as well would give 71%, about 1.4x.
- **Caveat:** steps measure kernel work only. The CPU readout and DE overhead do not shrink in proportion.

**Dropped targets most likely to break**, in order. Check them after every core fit.
1. **0.1 Hz -120/-200 (11, 12):** no LTD at long single-pairing delays. The core has burst -200 LTD but nothing that forbids single-pairing LTD at -200, so tau_T or theta_Tg can drift. Rule A already broke -120 (z 1.3).
2. **10 Hz -10 Sjostrom (18):** the deepest LTD, already z 1.7. Only Markram 10 Hz -10 stands in for it.
3. **40 Hz -10 (22):** LTD to LTP transition.
4. **0.1 Hz -25 (10):** window width (D: z 4.0).
5. **Markram 10 Hz +5 and Sjostrom 10 Hz +10 (2, 17):** they rise to z 2.0 and 1.4 under joint pressure.
6. **L2/3 3AP +10 distal/proximal (31, 33):** spread 7.4 and 3.0; must be checked in every joint core fit.

Low risk: 20 and 40 Hz +10 (19, 21; bracketed by 3 and 23), post_only (29), burst -120 (13).

## 4. Running a core fit

fit_v4.py has no per-target filter. What it supports now:
- `--groups` selects whole groups (paired_l5, sjostrom07).
- `--conditions` drops a condition everywhere, e.g. nmdar_block (no effect on the optimum).
- `--dirs` drops whole dirs. Dropping sj03 (from the default run_fit_v4.sh DIRS list) removes targets 10-12 and saves 23% of L5 steps, but none of those targets are core.
- `--pairs` subsamples records. It is a separate speed lever with more noise.

A protocol-level core needs one small change (not implemented):
1. In `build()`, after `T = {...}` is built from `load_targets`, add `if keep: T = {k: v for k, v in T.items() if k[0].split("@")[0] in keep}`. `protos` and the npz loading follow from T, so only the core records load.
2. Pass `keep` through a new `--protos` argument (comma list, applied to both the L5 and the `--joint` L2/3 build), plus an env `PROTOS` in run_fit_v4.sh.

Core list for `--protos`: `sjostrom_0.1hz_dt+10ms,sjostrom_0.1hz_dt-10ms,sjostrom_20hz_dt-10ms,sjostrom_50hz_dt-10ms,sjostrom_50hz_dt+10ms,10Hz_-10ms,10Hz_10ms,sjostrom_burst5x20hz_r50_dt-200ms,sjostrom07_step200ms_pair,sjostrom07_step200ms_pre_only,letzkus_1ap_dt+10ms,letzkus_3ap_200hz_dt-10ms`. sjostrom_50hz_dt+10ms is shared by L5 and L2/3 (the name matches in both builds).

Caveat: `peak` (the rules_row admissibility penalty) is built from the loaded records only, so the core fit's admissible set is slightly looser than the full one.

Validation after each core fit: re-score on all targets with the existing zero-generation path. Run run_fit_v4.sh with `MAXITER=0 SEEDFITS=<core>.json` (plus `JOINT=1` for the L2/3 rows). With `--maxiter 0`, fit_v4 evaluates `seeds[0]` and writes the full 29-target (+9) table. That also re-checks admissibility on all records. Size it from the full fit's MaxRSS (24-26 GB), with a few minutes of wall time.
