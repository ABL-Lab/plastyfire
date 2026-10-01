# BCL validation: G8X_G1a_u6 (G1 shaft gate, live GluSynapseV8)

Fit: rho_redesign/results/v7_G8X_G1a_u6.json (gpu_v8_rho gate_src 2, gate_win 100 ms, theta_G 0.15 uM fixed,
t_exact 1, v5_mode 2, veto_T 25, ecb_ref 2; gd 35.16, gp 370.68, k_E 0.4958; chi2 292.21, BIC 321.4, k 7).
Live: GluSynapseV8 gate_src 2 (mod/GluSynapseV8.mod, mod_build_v8), live_v8.py; offline: compare_prefire_v8.py (CPU port
of the gpu_v8 G1 kernel with the v7x exact-arrival corrections). Outputs: /scratch/dhuruva/bcl_valid_v8_G8X_G1a_u6/.

## Pre-checks
- compile + smoke_v8 22193266 (47 s, 281 MB): gate_src 0 == V7 bit for bit, licence timeline, licence replaces Vg.
- G1a pilot 22193267 (2:15, 1.44 GB): L5 192879-186028 Sj01 50/20 Hz +10: rho disagreements 0/12, licence openings
  live 39 / offline 39, live chi2 = offline 36.11.

## Jobs (submit_valid_v8.sh)
- gate smoke 22195075 (1800M 0:15) -> shards L5 22195076-78 (20G 0:45), L2/3->L5 22195079-80 (23G 0:35), L2/3->L2/3 22195081-88 (16G 0:30), 32 CPU -> compare 22195089 (20G 0:30). Sizing in the script header and DECISIONS.md.

## Results (compare 22197386, 27:58, 21.8 GB)
| pathway | targets | BCL live | offline | PASS | rho disagree | licence openings live / offline |
|---|---|---|---|---|---|---|
| L5 | 40 | 154.94 | 150.66 | 40/40 | 1 / 7044 | 112109 / 112242 |
| L2/3->L5 | 9 (+2 val) | 18.53 | 17.46 | 11/11 | 6 / 3650 | 50595 / 50834 |
| L2/3->L2/3 | 16 (+1 val) | 124.65 | 124.08 | 17/17 | 4 / 7769 | 163817 / 163746 |
| total | 65 | 298.12 | 292.21 | 68/68 | 11 / 18463 | 326521 / 326822 |

Verdict: PASS. The live mod reproduces the offline shaft-gate fit; the W2_N3 dt0 failures are gone (exact-arrival timing).
