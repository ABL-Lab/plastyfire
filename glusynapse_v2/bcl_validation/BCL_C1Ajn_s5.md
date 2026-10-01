# BCL prefire validation of fit v4_C1Ajn_s5

**Request (user):** "run prefire BCL on C1Ajn_s5 across all targets and see if the fitted values match what you are seeing".

**Status:** pilot submitted 2026-10-01 (jobs 22135321 live, 22135323 compare). The full run is pending the pilot's seff.

## What is compared
- **Fit:** `rho_redesign/results/v4_C1Ajn_s5.json`, χ² 51.97 over 39 targets: 30 L5→L5 (paired_l5 + sjostrom07) and 9 L2/3→L5 (paired_l23l5). The rule is C1 (vamp_mode 1, θ_V 5.50), rho_gamma 0.52, γd 51.7, γp 158.4, plus eCB (t_drive 4) and NO.
- **Live:** bluecellulab with `mod/GluSynapseV4.mod` (`live_v4.py --from-fit`). Each record's prefire config is rerun with cooker expression, so the induction reproduces the traces the fit was made on, while the v4 rule runs live with the fit's globals and per-synapse θ from the record's c_pre/c_post. The rule is continuous (`V4_BIN=0`), as it would run in a simulation.
- **Readout:** the final live ρ and dpre of each synapse go through the pair's EPSP basis (the fit's readout). This gives the live ratio per record, which is then averaged with the fit's pair selections per target.
- **Conditions:** these are live switches. mglu_block sets A_mglu 0. no_block sets A_NO 0. nmdar_block sets gmax_NMDA 0 for the whole prefire, plus A_mglu = A_NO = 0. Offline, nmdar_block assumes ρ = ρ0 and dpre 0; live, it tests that assumption.
- **Records:** 716 L5 (the fit's 24 pairs) and 579 L2/3 (111 all_protocols pairs; located nmdar targets on their distal pairs only), 1295 in total. This includes the 2 dropped Letzkus 3AP −10 distal targets, which are validation only.
- **Tolerance per target:** |BCL live − offline (same records)| ≤ max(0.02, 0.25 SEM).

## Why this mode
In equiv 22131696 (fit C1Ajd_s5), the continuous prefire runs had 0 ρ disagreements in 7/7 records, and the ratio was within 0.039 of the offline value. The offset comes from dpre +0.01 to +0.1, because the live rule rectifies the instantaneous VDCC current, while the offline rule rectifies the 0.25 ms bin mean. Bin mode (`V4_BIN=0.25`) is 6–18× slower per task, so the full set would cost about 350–1100 CPU·h. Its equivalence run (22134590) was still going when this was submitted.

## Verdict (full run 22135837 + compare 22135838, 2026-10-01)
**PASS 41/41.** All 1295 records ran (716 L5, 579 L2/3->L5). χ² over 39 fitted targets: BCL live 52.40, offline on the same records 51.98, fit 51.98 (max |offline - offline_fit| 0.0000).
- ρ disagreements: L5 2 / 5702 synapses, L2/3 15 / 3767. |d_basis| mean 0.0011 (L5) / 0.0081 (L2/3).
- Largest per-target live-offline gap: sjostrom07_step200ms_pair control +0.013 (dpre_maxdiff 1.12 at one synapse, eCB step timing), within tol 0.02.
- Figures: results/C1Ajn_s5_full.png/.pdf; tables results/C1Ajn_s5_full_{targets,records}.csv.
- seff: live 1:58:16, 16.99 GB of 17G, 89% CPU (32 workers); compare 4:36, 17.99 GB of 18G (both at the limit: 22G next time).
