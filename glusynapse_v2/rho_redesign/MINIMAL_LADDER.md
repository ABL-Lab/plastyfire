# Minimal extension of Chindemi 2022: nested-model ladder (joint fit)

Goal (user): "minimal extension of chindemi's work with less number of params, that fits all pathways."
Rule (2026-10-01): free = a00, a01, a10, a11, gamma_d, gamma_p only; every added parameter needs an experimental anchor.

Common setup: run_fit_v4.sh, JOINT=1, FITGAMMA=1, default L5 dirs/groups (no sj04), DROPT Letzkus 3AP -10 distal
(control + nmdar_block) -> 39 targets (30 L5 + 9 L2/3->L5). Seeded runs: SEEDFITS v4_C1Ajn_s5.json (SEED 5, SEEDSET {}:
the seed has every free slot); basin check: unseeded (SEED 6; pop[0] a's from x0 delta-cooker as in all v4 fits).
2g MIG, 1 CPU, 50G, 0:45 (basis 22132568/22133630: 39.4-39.5 GB, 20-26 min); check 0:15 (22133806: 2:42, 39.1 GB).
AIC = chi2 + 2k, BIC = chi2 + k ln 39 (ln 39 = 3.664).

## How eCB and NO are switched off (exactly)

- FIX (--fix-params) A_mglu = 0 and A_NO = 0. Their DE slots are log10 amplitudes, so the fixed slot is -inf and
  fit_v2.unpack gives 10**-inf = 0.0 exactly.
- In the kernel (gpu_v3 / gpu_v4_rho, same lines): am = 0 gives fm = 1, d = dmin + (d - dmin) = 0 exactly; r_no = 0 gives
  ex = 1, d = dmax - (dmax - d) = 0 exactly. So dpre = dpre0 = 0 in every lane (control, mglu_block, no_block).
- The 8 eCB/NO filters that run_fit_v4.sh always frees (theta_Te, tau_T, theta_Tg, dpre_min, theta_NOi, tau_NO, tau_Z,
  theta_Z) are inert once both amplitudes are 0. They are fixed via FIX (at the C1Ajn_s5 values, tau_Z > tau_NO so the
  rules penalty is not hit) only to take them out of the DE vector. They are not counted in k and are not a reduction.
- rho_gamma = 1 via SET (Chindemi cq = c_post, bitwise); vamp_mode 0 (M0) or 1 (M1).
- Check (job 22135018, MAXITER=0, M0 at the C1Ajn_s5 a's and gammas): "off" holds if, in results/v4_LM0chk.csv,
  pred(mglu_block) == pred(control) for 10Hz_-10ms, sjostrom_0.1hz_dt-10ms, sjostrom_20hz_dt-10ms,
  sjostrom_50hz_dt+10ms, sjostrom07_step200ms_pair, and pred(no_block) == pred(control) for sjostrom07_step200ms_pair
  (in C1Ajn_s5 these differ, e.g. step200 pair 1.674 / 2.137 / 1.237). The 4 fits run afterok on it.

FIX used for M0 and M1:
`{"A_mglu": 0, "A_NO": 0, "theta_Te": 0.2245, "tau_T": 30.44, "theta_Tg": 6.18, "dpre_min": -0.245, "theta_NOi": 0.000188, "tau_NO": 7.0, "tau_Z": 48.74, "theta_Z": 1.75}`

## Rungs

| rung | model | free (k) | added constants and their anchor |
|---|---|---|---|
| M0 | pure Chindemi (no eCB, no NO, no gate, rho_gamma 1) | a00 a01 a10 a11 gd gp (6) | none |
| M1 | M0 + spine VDCC-Ca potentiation gate (C1, vamp_mode 1) | + theta_V (7) | theta_V free: requirement anchor only (Ni2+-sensitive Ca needed for LTP: Kampa 2006, Letzkus 2006), no absolute value. tau_E1 = 100 ms fixed (gate memory): its current anchor is the eCB tLTD timescale, not a gate measurement -> needs an anchor. i_scale = units (absorbed by theta_V). |
| M1b | M1 + rho_gamma free | 8 | rho_gamma: no anchor (Chindemi / GB use 1) -> needs the user's OK |
| M2 | M1 + eCB-LTD at PARAM_ANCHORS values | 7 + unanchored eCB params | waits for PARAM_ANCHORS.md |
| M3 | M2 + NO-LTP at PARAM_ANCHORS values | | waits for PARAM_ANCHORS.md |

## Jobs and results

| job | tag | rung | seed | chi2 L5 / L2/3 / total | k | AIC | BIC | worst 5 (z) |
|---|---|---|---|---|---|---|---|---|
| 22135018 | LM0chk | M0 check, MAXITER=0 | C1Ajn_s5 | pending | 6 | | | |
| 22135019 | LM0_s5 | M0 | C1Ajn_s5 | pending | 6 | | | |
| 22135021 | LM0_s6 | M0 | unseeded | pending | 6 | | | |
| 22135020 | LM1_s5 | M1 | C1Ajn_s5 | pending | 7 | | | |
| 22135022 | LM1_s6 | M1 | unseeded | pending | 7 | | | |
| ref | C1Ajn_s5 | current rule | | 46.89 / 5.08 / 51.98 | 18 | 87.98 | 117.92 | |

A rung beats the current rule when chi2 < 75.98 (AIC) / 95.94 (BIC) for M0, 73.98 / 92.27 for M1, 71.98 / 88.61 for M1b.

## Three pathways (L5 30 + L2/3->L5 9 + L2/3->L2/3 Zilberter 15 = 54 targets)

fit_v4n.py (FITPY), JOINT=1 + EXTRA `l23l23:paired_l23l23:glusynapse_v2/extracted/zilberter_l23l23_delta-prefire-vseg-rs:basis_results_edges_zilberter_l23l23_delta_rs`
(the spec of pilot 22134473: 1527 records, 15 targets). Same FIX / SET / DROPT as above. fit_v4n.py got --fix-params
(port of fit_v4.py, default {} unchanged), checked by LMzrepro (2-pathway, 4 params fixed at the C1Ajn_s5 values: the
log must say REPRO OK). Seed C1Ajn_s5 only (LM0_s5 / LM1_s5 had not run at submission). 3g MIG, 1 CPU, 79G (63.3 GB of
22134473 + 25%), 1:45 (26 min x 2533 / 1006 records = 66 min, +50%); checks 0:15. Fits run afterok on both checks.
BIC = chi2 + k ln 54 (3.989). Reference C1Ajn_s5 3-pathway (22134473): 46.89 + 5.08 + 110.27 = 162.24, k 18, AIC 198.24,
BIC 234.04. Beats it when chi2 < 186.24 (AIC) / 210.11 (BIC) for M0, 184.24 / 206.12 for M1.

| job | tag | rung | seed | chi2 L5 / L2/3->L5 / L2/3->L2/3 / total | k | AIC | BIC | worst 5 (z) |
|---|---|---|---|---|---|---|---|---|
| 22135203 | LMzrepro | fit_v4n --fix-params repro (2 pathways) | C1Ajn_s5 | must equal 46.891 / 5.084 | 18 | | | |
| 22135204 | LM0zchk | M0 off-check, MAXITER=0 | C1Ajn_s5 | pending | 6 | | | |
| 22135205 | LM0z_s5 | M0 | C1Ajn_s5 | pending | 6 | | | |
| 22135207 | LM0z_s6 | M0 | unseeded | pending | 6 | | | |
| 22135206 | LM1z_s5 | M1 | C1Ajn_s5 | pending | 7 | | | |
| 22135208 | LM1z_s6 | M1 | unseeded | pending | 7 | | | |

## Anchored eCB / NO rungs (39 targets; PARAM_ANCHORS.md)

All: M1 base (vamp_mode 1, theta_V free, rho_gamma 1), fit_v4.py, 2g 50G 0:45, seed C1Ajn_s5 (SEED 5); M2a also unseeded (SEED 6).
eCB anchors: dpre_min = -0.29 (E, Sjostrom 2003 ACEA ceiling). Threshold collapse = PARAM_ANCHORS section 3 option B:
theta_Te = 0 (structural: the offset is a reparameterisation of the threshold) and tau_T = 50 ms (E upper bound,
Heinbockel 2005 at 32-37 C), so the one free eCB threshold is theta_Tg (= S* x tau_T in the quasi-steady reading).
v4 has no step-on-S kernel, so this is the closest exact emulation; the T stage stays a 50 ms low-pass (lag).
Risk: with tau_T 50 a single event's T peaks ~69 ms after the bAP, so the -100/-120 no-LTD targets may depress.
NO anchors (M3): tau_NO = 6.7 ms (E, Hall & Garthwaite 2006), theta_NOi = 6e-6 nA = 6 fA (C-NO1, nNOS EC50 0.2-0.3 uM),
tau_Z = 10 ms (E lower bound, Padamsey 2017), theta_Z = 0 (E hippocampus, one spike suffices). NO off in M2: A_NO = 0.
k = number of fitted parameters. AIC = chi2 + 2k, BIC = chi2 + k ln 39.

| job | tag | rung | free (k) | fixed by anchor | flagged unanchored (fitted) | beats C1Ajn_s5 if chi2 < AIC / BIC | result |
|---|---|---|---|---|---|---|---|
| 22135509 | LM2a_s5 | M2a: M1 + eCB, A_mglu = 1 (saturation test) | 6 + theta_V, theta_Tg (8) | dpre_min, theta_Te, tau_T; A_mglu = 1 | theta_V (C-V null), theta_Tg (S*) | 71.98 / 88.61 | pending |
| 22135510 | LM2a_s6 | M2a unseeded | 8 | same | same | same | pending |
| 22135511 | LM2b_s5 | M2b: M2a, A_mglu free | 9 | dpre_min, theta_Te, tau_T | theta_V, theta_Tg, A_mglu | 69.98 / 84.95 | pending |
| 22135512 | LM3a_s5 | M3a: M2a + NO | 6 + theta_V, theta_Tg, A_NO (9) | M2a + tau_NO, theta_NOi, tau_Z, theta_Z = 0 | theta_V, theta_Tg, A_NO (C-NO2 not done) | 69.98 / 84.95 | pending |
| 22135513 | LM3b_s5 | M3b: M3a, theta_Z free | 10 | M3a minus theta_Z | theta_V, theta_Tg, A_NO, theta_Z | 67.98 / 81.28 | pending |

Fixed constants with weak anchors in every gated rung: tau_E1 100 ms (L-M bracket), K (C-K insensitivity check not run),
dpre_max 1.0. theta_NOi's mapping from bulk to nNOS Ca is L confidence (PARAM_ANCHORS section 4: high risk).
FIX strings: M2a `{"A_mglu": 1, "dpre_min": -0.29, "theta_Te": 0, "tau_T": 50, "A_NO": 0, "theta_NOi": 0.000188, "tau_NO": 7.0, "tau_Z": 48.74, "theta_Z": 1.75}`
(NO entries inert); M2b = M2a without A_mglu; M3a `{"A_mglu": 1, "dpre_min": -0.29, "theta_Te": 0, "tau_T": 50, "theta_NOi": 6e-6, "tau_NO": 6.7, "tau_Z": 10, "theta_Z": 0}`; M3b = M3a without theta_Z.

## Results M0/M1, 39 targets (filled in by the coordinator)
| rung | job | k | chi2 L5 | chi2 L2/3->L5 | total | AIC | BIC (n 39) |
|---|---|---|---|---|---|---|---|
| M0 seeded | LM0_s5 22135019 | 6 | 151.01 | 29.01 | 180.02 | 192.0 | 202.0 |
| M0 unseeded | LM0_s6 22135021 | 6 | 148.69 | 15.36 | **164.05** | 176.1 | 186.0 |
| M1 seeded | LM1_s5 22135020 | 7 | 154.87 | 2.38 | **157.25** | 171.3 | 182.9 |
| M1 unseeded | LM1_s6 22135022 | 7 | 157.33 | 5.50 | 162.83 | 176.8 | 188.5 |
| C1Ajn_s5 (ref) | | 18 | 46.89 | 5.08 | 51.98 | 88.0 | 117.9 |

The worst targets in both rungs are the Sjostrom post-before-pre LTD rows: 0.1 Hz -10/-25, 10 and 20 Hz -10, and the 5x20 Hz bursts at -120/-200. The model predicts about 1.0 there, against data 0.57-0.79. The Sj07 step-pair NO-block arm is also off (1.36 vs 1.57-1.62). The Ca-threshold rule alone cannot produce this LTD, so eCB presynaptic LTD is necessary. Next: M2a/M2b/M3a/M3b (22135509-13) and v5/v5b/v5c.
Sizing: fits 7:46-12:57 and 27.8-36.0 GB, so 50G stays and 0:30 is enough next time.

## Results M0/M1, 54 targets (3 pathways; coordinator)
| rung | job | k | L5 | L2/3->L5 | L2/3->L2/3 | total | AIC | BIC (ln 54) |
|---|---|---|---|---|---|---|---|---|
| M0 s5 | 22135205 | 6 | 157.35 | 23.77 | 137.46 | **318.58** | 330.6 | 342.5 |
| M0 s6 | 22135207 | 6 | 172.04 | 17.76 | 134.57 | 324.36 | 336.4 | 348.3 |
| M1 s5 | 22135206 | 7 | 169.58 | 10.71 | 91.96 | **272.25** | 286.2 | 300.2 |
| M1 s6 | 22135208 | 7 | 165.72 | 21.91 | 87.64 | 275.28 | 289.3 | 303.2 |
| C1Ajn_s5 (not refit) | | 18 | 46.89 | 5.08 | 110.27 | 162.24 | 198.2 | 234.0 |
The gate helps L2/3->L2/3 even without eCB (about 135 -> 88-92, better than C1Ajn's 110), so the eCB in C1Ajn hurts L2/3->L2/3. L5 needs eCB in any case.
seff: 11:29-14:18, MaxRSS 41.9-63.1 GB (79G ok), so 0:30 next time.

## Results M2/M3 (anchored), 39 targets (coordinator)
| rung | job | k | L5 | L2/3->L5 | total | AIC | BIC (n 39) |
|---|---|---|---|---|---|---|---|
| M2a s5 | 22135509 | 8 | 97.02 | 5.87 | **102.89** | 118.9 | 132.2 |
| M2a s6 (unseeded) | 22135510 | 8 | 106.16 | 4.40 | 110.56 | 126.6 | 139.9 |
| M2b s5 (+A_mglu) | 22135511 | 9 | 99.80 | 9.47 | 109.27 | 127.3 | 142.2 |
| M3a s5 (+NO) | 22135512 | 9 | 100.07 | 7.01 | 107.08 | 125.1 | 140.0 |
| M3b s5 (+NO, θZ) | 22135513 | 10 | 90.76 | 11.48 | 102.24 | 122.2 | 138.9 |
| C1Ajn_s5 (ref) | | 18 | 46.89 | 5.08 | 51.98 | 88.0 | 117.9 |

eCB cuts M1's 157 to 103 with 1 extra parameter, so it is the key mechanism. NO and a free A_mglu add nothing (M2b and M3a score worse than M2a, which means DE has not converged at 150 iterations).
M2a's worst targets: sjostrom 0.1 Hz -120/-200 predicted 0.73 against data 1.05. This was the predicted cost of the τ_T 50 ms pass-through emulation: one event's eCB acts at long post-pre intervals. Also bad: Sj07 step pair (1.29 vs 1.62), 50 Hz -10 (1.03 vs 1.70) and 0.1 Hz +10 (1.11 vs 0.97). The -120/-200 miss needs a proper timing-selective eCB trigger, which v5b/E2 provide.
seff: 13:28-22:44, 27-39 GB.
