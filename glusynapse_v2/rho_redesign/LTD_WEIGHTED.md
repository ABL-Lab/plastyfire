# LTD-weighted refits of v5c (delta-split1, 39 targets, 14 LTD = data mean < 1.0, 25 LTP)

Weighting: z^2 x w on LTD targets in the DE objective only. All chi2 below are UNWEIGHTED. Jobs 22154961-5, tags S1_V5cW*, results/v5_<tag>.json.
Rescore of v5_S1_V5c_s6_39 (no refit): REPRO OK 1e-9, 87.8358 = LTD 50.415 + LTP 37.421 (so baseline LTD already carries 57% of chi2).

| w | seed | chi2 L5 | chi2 L23 | total | chi2_LTD | chi2_LTP | BIC |
|---|---|---|---|---|---|---|---|
| 1 | V5c_s6_39 (unseeded) | 73.17 | 14.66 | 87.84 | 50.42 | 37.42 | 117.15 |
| 2 | s5 (seeded) | 72.33 | 15.38 | 87.71 | 47.38 | 40.33 | 117.02 |
| 2 | s6 (unseeded) | 71.32 | 16.79 | 88.11 | 45.55 | 42.56 | 117.42 |
| 4 | s5 (seeded) | 69.28 | 19.33 | 88.61 | 39.13 | 49.48 | 117.92 |
| 4 | s6 (unseeded) | 86.49 | 18.33 | 104.81 | 39.24 | 65.57 | 134.12 (worse basin) |

## Fitted parameters
| fit | a00 | a01 | a10 | a11 | theta_V | theta_eCB | gamma_d | gamma_p |
|---|---|---|---|---|---|---|---|---|
| w1 s6_39 | 1.066 | 4.455 | 1.657 | 4.674 | 1.03 | 20.7 | 34.8 | 507.0 |
| W2 s5 | 1.064 | 4.622 | 1.659 | 4.670 | 1.02 | 13.1 | 34.6 | 496.6 |
| W2 s6 | 1.187 | 2.741 | 1.234 | 4.807 | 2.94 | 11.1 | 97.4 | 516.6 |
| W4 s5 | 1.206 | 2.949 | 1.248 | 4.565 | 3.30 | 11.1 | 147.5 | 598.1 |
| W4 s6 | 1.111 | 3.753 | 1.924 | 4.053 | 0.92 | 13.1 | 36.4 | 599.4 |

## LTD targets: data vs unweighted fit (w1) vs weighted fits (pred)
L5 (rows 1-8) then L2/3 (rows 9-14); SEM = data SEM.
| target | data | SEM | w1 | W2 s5 | W2 s6 | W4 s5 | W4 s6 |
|---|---|---|---|---|---|---|---|
| markram 10Hz -10ms | 0.79 | 0.03 | 0.82 | 0.80 | 0.82 | 0.81 | 0.80 |
| sj 0.1Hz -25ms | 0.65 | 0.09 | 0.89 | 0.88 | 0.87 | 0.87 | 0.88 |
| sj burst r50 -120ms | 0.79 | 0.08 | 0.90 | 0.88 | 0.88 | 0.88 | 0.88 |
| sj burst r50 -200ms | 0.79 | 0.08 | 0.93 | 0.91 | 0.90 | 0.90 | 0.91 |
| sj 0.1Hz +10ms | 0.97 | 0.04 | 1.03 | 1.03 | 1.09 | 1.05 | 1.03 |
| sj 0.1Hz -10ms | 0.69 | 0.07 | 0.88 | 0.87 | 0.87 | 0.87 | 0.87 |
| sj 10Hz -10ms | 0.57 | 0.11 | 0.86 | 0.84 | 0.86 | 0.84 | 0.80 |
| sj 20Hz -10ms | 0.65 | 0.09 | 0.94 | 0.93 | 0.85 | 0.81 | 0.85 |
| sj 50Hz +10ms distal (L2/3) | 0.86 | 0.09 | 0.84 | 0.84 | 0.76 | 0.71 | 0.87 |
| letzkus 1AP +10ms | 0.72 | 0.03 | 0.75 | 0.75 | 0.75 | 0.75 | 0.74 |
| letzkus 3AP 200Hz +10ms distal | 0.79 | 0.06 | 0.84 | 0.84 | 0.82 | 0.82 | 0.85 |
| letzkus 3AP 200Hz -10ms proximal | 0.89 | 0.03 | 0.81 | 0.80 | 0.84 | 0.83 | 0.81 |
| letzkus nopost | 0.98 | 0.06 | 1.04 | 1.04 | 1.04 | 1.04 | 1.04 |
| letzkus 3AP -500ms | 0.99 | 0.05 | 1.05 | 1.05 | 1.05 | 1.05 | 1.05 |

## Conclusion
Weighting does not fix LTD: it moves error from LTD to LTP at about constant total. Seeded fits go 87.8 -> 87.7 (w2) -> 88.6 (w4) in unweighted total while chi2_LTD falls 50.4 -> 47.4 -> 39.1 and chi2_LTP rises 37.4 -> 40.3 -> 49.5. The weighted objective buys 11 LTD chi2 at a cost of 12 LTP chi2 (BIC +0.3 at w4). The unseeded w4 basin is simply worse (105). So the model structurally cannot do both with this rule and these params; the unweighted fit is already near the best compromise, not an LTD-starved local optimum.

Which LTD rows weighting can reach: only the L5 rows that sit at 0.80-0.90 (10 Hz / 20 Hz -10 ms, 50 Hz distal L2/3 and Letzkus 3AP distal) which move by 0.02-0.15, and those moves are at the expense of LTP rows (W2 s6 and W4 s5 trade the 50 Hz +10 ms distal L2/3 row from 0.84 to 0.76/0.71 against 0.86 data, i.e. overshoot). Cannot reach: the deep Sjostrom LTD rows (0.1 Hz -25/-10 ms, 10 Hz -10 ms, 20 Hz -10 ms; data 0.57-0.69, all fits 0.80-0.94, the 4x weight does not move 0.1 Hz at all, 0.88 -> 0.87) and the Sjostrom burst r50 rows. Those depress by about 0.1-0.3 too little in every fit, so the rule lacks a strong low-frequency / post-before-pre LTD arm that is independent of the LTP arm. Rows already fine (Letzkus 1AP, nopost, -500 ms) stay put. One miss goes the other way: Letzkus 3AP -10 ms proximal (0.89 data, 0.81 predicted) is over-depressed, so LTD is not uniformly short.

Measured (seff): rescore 22154961 4:43 elapsed, 37.25 GB of 40G, 52% CPU; fits 11:43-13:15, 26.2-37.2 GB of 50G, 95-99% CPU.
