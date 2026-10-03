# s4wide diagnosis (2026-10-03)

Scope: 32 s4wide kimchi runs present (w_4N_u12, w_5N_s11, w_5N_u12, w_5N_u13 still missing; their vulcan attempts timed out)
plus the 16 rorqual census fits (f3Wdah, f4Vdah, f4Ndah; there are no f3V/f3N/f5N census fits). Sources: ledger `.results`, phase-A jsons in
/scratch/dhuruva/stage4 and /scratch/dhuruva/glusyn_runs, census `*_val.csv`, select_*.csv. Pure jq/awk, no jobs submitted.
Notation: A = phase-A de_fun (rho_sigma 0.1, smoothed readout); B = final de_fun (phase B, binary readout). Both include the hinge.

## A. Did the DE run? Yes. The summary fields describe phase B only

- run_stage4_fit.sh fits in two calls: phase A (450 generations max at rho_sigma 0.1), then phase B (`--resume` from phase A's final
  population, 150 generations max, binary readout). The summary `nfev`, `minutes` and `de_fun` come from the phase-B json. Only
  `de_fun_phaseA` comes from phase A.
- With `vectorized=True`, scipy 1.17 counts one nfev per generation call. So nfev = generations + 1. "nfev 7, 1.18 min" for w_3N_s11
  means phase B ran 6 generations. The 3525 s `elapsed_fit_s` is mostly phase A on the L40S.
- Phase A did the work. Census phase A ran 271-450 generations: 13 of 16 stopped on the tol 1e-6 test (the population collapsed) and 3 hit
  the 450 cap (f4Ndah_u8, f4Vdah_s6, f4Vdah_u11). The kimchi logs on rorqual agree: w_3V_u15 converged at generation 354, w_5N_u15 at 429.
- Phase B is close to a no-op. It starts from a collapsed population on a piecewise-constant binary objective, so all energies match after
  1-20 generations and the tol test stops it. Only w_4V_s12 ran all 150 (nfev 151). Its phase A comes from the V3 basin, which hit the cap.

| model | runs | phase-B nfev | B - A | fit wall, H100 (s) | fit wall, L40S (s) |
|---|---|---|---|---|---|
| 3W | 6 | 2-21 | -3.2 to +1.0 | 601-960 | 2445-3117 |
| 3V | 6 | 7-56 | -1.3 to 0.0 | 755-772 | 2393-2675 |
| 3N | 6 | 2-19 | -6.3 to +12.6 | 893-1045 | 3398-4098 |
| 4V | 7 | 9-151 | -8.0 to +37.4 (census f4Vdah_s5: +94) | 1693-2435 | timeout (1:30) |
| 4N | 5 | 2-43 | -14.7 to +8.0 | 1788-2048 | timeout |
| 5N | 3 | 2-13 | -11.9 to +27.0 | 2226-2696 | timeout |

- **Seeded vs unseeded.** No seeded run kept its seed: no final parameter set is near a seed json. s* and u* starts fall into the same
  basins: 3N s12 = u13; 4N s11 = census u7; 4V s11 = census u8/u9 and w u14. There is no systematic s/u difference. The single best start
  is seeded for 3N (s11), 3V (s11) and 5N (s12), and unseeded for 3W, 4V and 4N.
- **No optimiser bug.** There are two bookkeeping issues:
  1. The summary should carry phase-A nfev and minutes.
  2. basin_s4wide.csv ranks by B with df <= 2, so binary jitter inside one basin hides convergence (see B). Its header line is also
     sorted to the bottom.

## B. Basin structure (A is smooth and consistent within a basin; B jitters)

`*` marks a parameter on its bound: gamma_d 20 (lower), gamma_p 100 / 600, a = 0. a = (a00, a01, a10, a11); a2x/a3x are tied to them.

| model | basin: a; gamma_d; gamma_p | A | B | starts |
|---|---|---|---|---|
| 3V | (0.79, 2.37, 0.91, 3.54); 20*; 481-486 | 115.71 | 114.42-115.67 | 6/6 |
| 3W | W1 (0.79-0.84, 1.88-2.37, 0.86-0.91, 2.95-3.50); 20-23*; 264-549 | 127.71-131.66 | 126.91-129.66 | 5/10 |
| 3W | W2 (1.17-1.19, 0.52, 1.34, 1.79-1.83); 41-45; 508-595 | 137.36-137.44 | 137.67-138.66 | 4/10 |
| 3W | w_u15 (0.81, 2.16, 1.22, 2.62); 22; 412 | 137.44 | 138.32 | 1/10 (same A as W2) |
| 3N | s11 (1.12, 2.47, 1.33, 2.63); 151; 485 | 152.79 | **146.52** | 1/6 |
| 3N | u15 (2.60, 0*, 2.69, 0.36); 128; 100* | 157.89 | 154.57 | 1/6 |
| 3N | (1.95-2.08, 0-0.29, 3.45-3.47, 0.20-0.25); 20*; 432-446 | 160.69-160.71 | 155.79-158.89 | 3/6 |
| 3N | u12 (2.49, 0*, 2.60, 0.51); 249; 436 | 162.75 | 175.32 | 1/6 |
| 4V | V1 (0.97-0.98, 0.36, 2.39-2.45, 0.41-0.51); 20-22*; 463-600 | 198.78-198.96 | **190.54-193.46** | 6/14 |
| 4V | V4 (1.85, 0*, 1.94, 0.67); 91; 100-108* | 204.46 | 241.86 / 298.66 | 2/14 |
| 4V | V2 (0.90, 2.62, 1.09, 3.71); 20*; 533-600* | 206.14-206.17 | 215.36-217.03 | 2/14 |
| 4V | u7 (0.80, 2.79, 0.80, 4.21); 20*; 550 | 208.84 | 224.72 | 1/14 |
| 4V | V3 (1.62, 0*, 1.80, 0.75); 104-106; 155-157 | 210.03-210.04 | 227.34-240.79 | 3/14 |
| 4N | N2 f4Ndah_u8 (1.80, 0*, 1.86, 0.73); 99; 107 | **159.59** | 170.37 | 1/9 |
| 4N | w_u14 (1.72, 0*, 1.76, 0.92); 209; 520 | 169.33 | 177.31 | 1/9 |
| 4N | N1 (0.94, 0.52, 2.37, 0.51); 20*; 443-493 | 173.37-173.39 | **161.41-162.52** | 3/9 |
| 4N | N3 (0.87, 2.67, 1.16, 3.56); 20*; 600* | 193.27 | 188.86 | 2/9 |
| 4N | N4 (1.06, 2.50, 1.36, 2.52); 123-128; 429-447 | 206.70-206.80 | 192.06-192.65 | 2/9 |
| 5N | s12 (2.37, 0*, 2.50, 0.29); 159; 100* | **360.14** | **387.14** | 1/3 |
| 5N | (1.68-1.72, 2.03, 2.03-2.06, 1.99-2.04); 152-184; 504-532 | 402.36-403.38 | 391.48-391.56 | 2/3 |

- **Converged:** 3V (6/6), 3W (best basin W1 in 5/10) and 4V (V1 in 6/14). The basin census reports 1 start within 2 for 4V only
  because B jitters by 3 inside V1. **Not converged:** 3N (best basin 1/6) and 5N (3 of 6 runs; the two basins swap order between A and B).
  4N is ambiguous: A ranks N2 first by 14 and B ranks N1 first by 8.
- **Binary jitter.** In one basin, parameters that differ at the 1e-3 level give B values up to 3 apart (4V V1 190.5-193.5; 3N
  155.8-158.9), Markram +10 moving by 0.02 (V1 1.133-1.152), and up to 94 apart on the V4 corner. Identical B down to 1e-13 also occurs for
  different parameters: 3N s12/u13 have gamma_p 432 vs 438. The objective is piecewise constant (22 L5 pairs), so B cannot rank
  neighbouring points.
- **Smoothing reorders basins.** The phase-A optimum is not the binary optimum, and phase B cannot move between basins (see A). Which basin
  "wins" therefore depends on the readout: 4N N2 vs N1, and the V4 corner (2nd best by A, worst by B).
- **Degeneracy.** 3W W2 and w_3W_u15 have the same A (137.44) but a00 0.81 vs 1.17, a01 2.16 vs 0.52 and gamma_d 22 vs 41. That is an
  a00/a01 ridge (c_pre and c_post co-vary). gamma_p is weakly set: 264-549 inside 3W W1, 463-600 inside 4V V1.
- **Bounds.** gamma_d sits on its literature floor of 20 in the best basin of 3V, 3W, 4V and 4N. The corner basins put a01 at 0 and gamma_p
  at 100.

## C. Why no run passes the L5 slope gate

- **Noise ceiling.** The data SEMs dilute the OLS slope: var(data) over the 40 rows is 0.119 and mean SEM^2 is 0.014. A model that
  predicts the true values exactly would score about 0.88 on average, not 1.
- **The same rows in every model.** Write 1 - slope = sum_i (x_i - xbar)(x_i - p_i) / Sxx, with x = data and p = pred. Mean over the 31
  runs analysed (the 4V_u12 and census vals came in after this pass): slope 0.47. The 10 rows below make up 0.41 of the 0.53 deficit. Each
  row's share has the same sign and similar size in all six models.

| L5 row | data | pred (best runs) | share of deficit | in fit core? |
|---|---|---|---|---|
| S07 step pair, mGluR block | 2.13 | 1.30-2.07 | 0.099 | no |
| Sj 50 Hz -10 | 1.70 | 0.96-1.24 | 0.076 | yes |
| Sj 20 Hz +25 | 0.65 | 0.90-1.32 | 0.044 | no |
| Sj 40 Hz -10 | 1.51 | 0.93-1.26 | 0.034 | no |
| Sj 20 Hz dt 0 | 0.76 | 1.08-1.46 | 0.033 | no |
| Sj 20 Hz -10 | 0.65 | 0.88-1.06 | 0.031 | no |
| Sj 20 Hz -25 | 0.70 | 0.90-1.12 | 0.026 | no |
| Sj 50 Hz +10 | 1.57 | 1.16-1.41 | 0.026 | no |
| Sj 10 Hz -10 (demoted) | 0.57 | 0.71-0.82 | 0.022 | no |
| S07 step pair, control | 1.62 | 1.13-1.97 | 0.021 | yes |

- **Structural, not a search problem.** The data show a narrow LTP window at 20 Hz (only +10 potentiates; 0, +25, -10 and -25 depress)
  and LTP in either order at 40-50 Hz. Every rule variant does the reverse: it potentiates 20 Hz at 0 and +25 and holds 40-50 Hz -10 near
  1.0-1.25. Across runs, pred(20 Hz +25) and pred(50 Hz -10) correlate at r = 0.93, and pred(20 Hz 0) with pred(50 Hz -10) at r = 0.77.
  A run that lifts 50 Hz -10 also lifts the false 20 Hz LTP.
  - 3W/3V hold 20 Hz 0/+25 at 0.9-1.1, but then 50 Hz -10 is 0.96 and the S07 pair 1.13.
  - 3N/4N reach 1.2 at 50 Hz -10, with 20 Hz 0 at 1.25-1.46.
  - 5N fits all 40 L5 rows and still scores 0.35-0.52, so the gate is not missed because rows are left out of the fit.
  - EMODEL_LIMITS.md already lists these rows as rule limits that a perfect emodel would not fix.
- **Trade-off with Markram +10.** Over all 48 runs, r(Markram +10 pred, slope) = -0.68. The front:

| run | slope | Markram +10 (z) | Markram -10 (z) | gate |
|---|---|---|---|---|
| f4Vdah_u7 | 0.695 | 1.104 (-1.56) | 0.784 | MISS, PASS |
| f4Vdah_s9 / w_4V_u12 (V2) | 0.652 / 0.650 | 1.112 (-1.43) | 0.779 / 0.766 | MISS, PASS |
| f4Ndah_u7 / w_4N_s11 (N3) | 0.625 | 1.119 (-1.32) | 0.769 | MISS, PASS |
| w_3N_s11 (best HIT) | 0.550 | 1.145 (-0.91) | 0.767 (-0.96) | HIT, low |
| f4Vdah_u8/u9 (V1) | 0.537-0.539 | 1.147 (-0.87) | 0.768 | HIT, low |
| 3W, all HITs | 0.35-0.41 | 1.17-1.20 | 0.77-0.80 | HIT, low |

  - Every run with slope >= 0.6 has Markram +10 at 1.10-1.12. No HIT exceeds 0.55.
  - The high-slope runs gain slope mainly through the S07 step-pair rows: mGluR block 2.11 vs 1.67 in the HIT runs, and control 1.94
    against data of 1.62 (an overshoot of z +4.6). They also gain through the 20-40 Hz 5 ms rows (+0.10-0.13). They do not gain at
    20 Hz 0/+25 (1.25/1.08, still false LTP).
  - The S07 mGluR-block row alone carries 22 % of the leverage (x = 2.13). Without it, the best slope is 0.62 (f4Vdah_u7), every other run
    is <= 0.56, and w_3N_s11 drops to 0.52.

## D. Markram +10 at 1.12-1.15: the hinge binds, and the gate sits on its edge

- **How the hinge works.** It is lam 100 x max(|z| - 1, 0)^2 on both Markram rows. Its gradient is zero at |z| = 1, so a fitted point
  ends just outside the band whenever the other targets pull. +10 at |z| = 1 is 1.139, and the data mean is 1.201.
- **Where it binds** (hinge > 0 at the final point): 3V 6/6, 4V 5/7, 4N 4/5, 3N 4/6, 5N 3/3, 3W 0/6.
  - 3V misses at z -1.019 (penalty 0.034) and 4V V1 at -1.02 to -1.09. These misses are within the binary jitter (see B), so HIT/MISS
    inside one basin is a coin flip: V1 is HIT for u8/u9 and MISS for u10/u13/u14/s11.
  - 5N misses at z -1.25 to -1.29 (penalty 6-11). There the other 40 L5 rows pull hard.
- **Implied pull of the rest of the objective** at the final point (2 lam x overshoot): about 4 chi2 per sigma of +10 in 3V and 3N s12;
  17-19 in 4N u15 and 4V u13; 50-64 in 5N, 4N s11 and 3N u12/u14. 3W is unconstrained and sits at 1.17-1.20.
- **Cost of hitting Markram, from basin pairs:**
  - 4N: N2 (HIT 1.187) vs N1 (MISS 1.133): chi2_eff/n 2.69 vs 2.66 (about +1 over 32 core rows), B +8, slope 0.42 vs 0.54.
  - 4V: V1 HIT (u8) vs V1 MISS (u10): chi2_eff/n 3.02 vs 2.98.
  - So hitting the band costs about 1-8 chi2 in 3V/4V/4N, but about 0.1 of slope where it means a basin change. In 5N it costs tens of
    chi2.
  - Moving to the data mean (1.20) only happens in 3W, whose chi2_eff/n is 3.8-4.0 against 2.7-3.0.

## E. Cluster / GPU

- **Same answers on L40S and H100.** No same-start pair finished on two GPU types (the vulcan 4V/4N/5N attempts timed out). Within one
  basin the values agree across clusters:
  - 3V: B 115.66639258288839 on vulcan L40S (u12), fir 3g.40gb (s12, u13) and rorqual 3g.40gb (u15), identical to 14 digits, with A
    115.71 everywhere.
  - w_4N_s11 (fir) and f4Ndah_u7 (rorqual): B 188.8632693203795, identical.
  - The 3V B outliers (114.42, 114.85) are binary jitter, not GPU: vulcan u12 equals fir.
- **Timing.** L40S full GPU is 3.2-3.5x slower than an H100 3g.40gb slice (3V 2393-2675 vs 755-772 s). That is why all six vulcan
  4V/4N/5N attempts hit 1:30.
- **Memory.** Peak RSS is 41-68 GB with no GPU or cluster pattern (max 67.9 GB, w_5N_u15 rorqual). Val takes 348-463 s on H100.

## F. Bottom line and next steps

- **Search is done for 3V, 3W, 4V and 4N.** Their best basins come back from several starts and from both seeded and unseeded starts.
  3N and 5N need more starts only if they stay in contention.
- **The two gates exclude each other for all six models.** That is structural:
  - Markram HIT needs +10 >= 1.139 and leaves slope <= 0.55.
  - Slope >= 0.6 only comes with +10 <= 1.12 and an S07 overshoot.
  - The slope deficit comes from the 20-50 Hz timing rows (false LTP at 20 Hz 0/+25/-25, missing LTP at 40-50 Hz -10). They move
    together in every variant.
- **The current lead is noise.** Under the pinned rule, w_3N_s11 (slope 0.55, chi2_eff/n 3.59, its basin found by 1 of 6 starts) leads
  f4Vdah_u8 (0.54, 3.02) by 0.01 of slope, which is below the binary jitter.
- **No fitter bug.** Two procedure issues remain: phase B is a near no-op, and the A/B readouts rank basins differently.

**Next steps, ranked by expected value:**
1. **User decision on the gates (no compute).** As it stands no model can pass, so the gate picks noise. Options:
   - accept the best-slope HIT with a tie band (for example 0.05 of slope);
   - use an SEM-weighted or leverage-robust slope (S07 mGluR block holds 22 % of the leverage);
   - set the slope bar against the 0.88 noise ceiling.
   With a 0.05 tie band, chi2_eff/n decides, which gives f4Vdah_u8 (V1) over w_3N_s11.
2. **Hinge margin refit (6-8 fits, about 0.5 h each on 3g.40gb).** Use k = 0.8 on both Markram rows, seeded from 4N N1, 4V V1, 3V and
   3N s11. This would turn the boundary misses into robust HITs at a small chi2 cost (pull 4-19 per sigma). N1 has the best chi2_eff/n of
   all runs (2.66) but misses at 1.133. It does not fix the slope.
3. **Fitter procedure (cheap).**
   - Choose the final point by binary rescoring (`--maxiter 0`) of the distinct members of every start's phase-A ckpt population, or
     re-spread the phase-B population.
   - Count basins by A.
   - Add phase-A nfev/minutes to the summary.
   This removes the A/B basin reordering (4N) and the jitter-driven HIT/MISS flips.
4. **Rule research (synapse-local, Ca only).** Find a local Ca feature that separates 50 Hz -10 (LTP) from 20 Hz 0/+25 (LTD). The
   existing variants (veto window, NO, eCB low-pass) all move these rows together (r = 0.93). More DE starts will not help here.
