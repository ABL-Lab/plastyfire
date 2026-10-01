# BCL prefire validation of fit v5_V5c_s7

**Status:** submitted 2026-10-01. Compile 22138635, then pilot 22138636 (afterok), then compare 22138637 (afterok). The full run is pending the pilot's seff.

## What is compared
- **Fit:** `rho_redesign/results/v5_V5c_s7.json`, χ² 84.69 over 39 targets (30 L5→L5 + 9 L2/3→L5), k 8. The rule is v5c (`v5_mode` 2): Chindemi ρ (ρ-γ 1, γd 168.1, γp 258.4), C1 gate θ_V 0.245, and eCB on the (1−b)-weighted pool W with θ_eCB 6.91, giving d → d_min −0.29 (A_eCB 1). There is no NO arm. Design: `rho_redesign/V5_DESIGN.md` §4.
- **Live:** bluecellulab with `mod/GluSynapseV5.mod` (new SUFFIX GluSynapseV5, built into `mod_build_v5/` by `compile_v5.sh`). The script is `live_v5.py --from-fit`, and the pipeline is the v4 one (BCL_C1Ajn_s5.md): prefire configs, cooker expression and a continuous rule.
- **Offline:** `compare_prefire_v5.py`. `v5_rec` is a CPU port of the gpu_v5_rho kernel, applied to the same extracted records.
- **Conditions (live):**
  - mglu_block: A_eCB 0.
  - nmdar_block: gmax_NMDA 0 and A_eCB 0. Offline, ρ is frozen and d = 0.
  - no_block: the same as control.
- **Deliberate discretisation difference:** the kernel snaps each own arrival to the first 0.25 ms grid sample at or after it, and steps W, V and b on the sampled −ica_VDCC. Live, the eCB check and the b jump happen at the exact arrival time, and the integration is continuous (CVODE).
- **Pilot:** the 38 tasks of run_valid_C1Ajn_s5.sh STAGE=pilot. That is L5 180351-198084 (30 records) and L2/3 10149-186264 (8 records). Outputs go to `/scratch/dhuruva/bcl_valid_V5c_s7/`.
- **Pass:** rho disagreements ≈ 0 and |d_basis| ≤ ~0.01. The v4 C1Ajn_s5 run gave 0/258 and ≤ 0.0094. The per-target tolerance is max(0.02, 0.25 SEM).

## Pilot result
Pending (log `logs/valid_V5c_s7_cmp_22138637.out`).

## Full run
22139631 (32 workers, 1295 records: 716 L5, 579 L2/3->L5): 1:45:37, MaxRSS 21.99 GB of 22G, 93% CPU. Compare 22139632: 2:44, 21.04 GB of 22G.

## Verdict (2026-10-01)
**PASS 39/41.** All 1295 records ran. χ² over 39 fitted targets: BCL live 89.45, offline on the same records 84.69 = fit (max |offline - offline_fit| 0.0000).
- ρ disagreements: L5 4 / 5702 synapses, L2/3 14 / 3767; |d_basis| mean 0.0036 (L5) / 0.0054 (L2/3).
- Fails: sjostrom_10hz_dt-10ms control (live 0.809 vs offline 0.776, d +0.032, tol 0.028; 1 ρ flip) and letzkus_3ap_200hz_dt-10ms@distal nmdar_block (validation-only target, d +0.028; 7 flips).
- Cause: the live eCB step (dpre -0.29, all-or-none) fires at exact arrival times while the offline kernel uses the 0.25 ms grid; at synapses where W sits at θ_eCB, one step more or less (dpre_maxdiff 0.29) moves a synapse across ρ* (bistable flip). v4 C1 (graded eCB) had 41/41 with 17 flips. The v5c gap (+4.8 χ²) is the price of the all-or-none eCB step near threshold, not a model/port error.
- Figures: results/V5c_s7_full.png/.pdf; tables results/V5c_s7_full_{targets,records}.csv.
