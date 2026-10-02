# glusynapse_v2 orchestrator state (overwrite at each milestone; history lives in DECISIONS.md)

Updated: 2026-10-02

## In flight
- STAGED FIT (user 2026-10-02, picked S1C): plan rho_redesign/IMPROVE_PLAN.md. Running: 22289873/74 S1C rule + 23 L2/3 targets (s/u, out /scratch/dhuruva/s1c_l23fit); 22289864/65 L5 diagnostics (/scratch/dhuruva/s1c_diag); 22290236-39 corrected S2C/S2N val. Watchers b4ad990ck, bz28venwe. Next: stage 3 (20 cores, STAGE3_DESIGN.md) seeded from S1C; plots only via plot_stage_fit.py into figs/stage_fits/.
- PENDING after the stage-3/W/G2 jobs finish: edit ebner/ebner_targets.csv row 8 (40 Hz +10) to 1.54, 0.113, n 6, source Fig 1D; check whether the kernel reads it at runtime or it is baked into the extracted dirs.
- 10 AGENTS (2026-10-02): K gpu_v11 kernel+MVD+veto window (SPEC_V11.md, test_v11.sh); S3 run_stage3_fit.sh 3C/3N fits; A ANCHORS_V11.md; D DISTAL_BAP_LIT.md; T TARGET_AUDIT.md; G run_g2_test.sh; W run_s1c_l23w.sh weight sensitivity; B V9 BCL pilot on S1C; F DOUBLET_BIAS.md; R plan review. Next after K: GluSynapseV10.mod (MVD + veto window) mod agent, then stage-3 v11 models.
- T31 DONE (review): basal spine bAP dV <60 um 77 +/- 17 mV vs Chindemi 80-100. Emodel OK; Ca gap = steep VDCC V dependence. Next: T32.
  Script cexp/bap_vpk.py (job 22218182).
- T32 scan (L5 CEXP_SABATINI): arrays 22218430/32/34/36, merges 22218431/33/35/37; out /scratch/dhuruva/split1/cexp_sabatini/t32/<cfg>. Grid ljp/gca: 25/0.0321, 12/0.0744, 12/0.110, 0/0.160.
- T32 LOCKED (review): ljp 25 / gca 0.0321, spine/delta_ljp25_g0321.json; lit cexp/BAP_CA_DISTANCE_LIT.md. Rebinned check job 22222133.
- T33 step 1+1b on og-delta CANCELLED (wrong emodel; 22222192-97, 22222649-55). Redo on delta-split1, all 7 sets incl. l5extra + l23l23extra (technician brief 2026-10-02).
- T33 step 1 DONE on delta-split1 (spine/run_ljp25g0321.sh, ~155 CPU-h; K 5.65e-07): caches 22223209-11, prefire 22223212-19, K 22223220, extracts 22223221-23 -> /scratch/dhuruva/split1/extracted/*_delta-split1-ljp25g0321-prefire-{vca,vseg-rs}. Reuses split1 workdirs, find calibration, EPSP bases (VDCC current enters V, ~6% of EPSP Ca; E1 also reused). Watch 22223218 (Zilberter split, OOM risk). Next: T33 step 2 refit.
- T35 G9 test (rho_redesign/submit_t35_g9.sh, gpu_v9_rho.py): (resubmitted, abs GLUSYN_GLOBALS) pilot 22229198 + check 22229199, arrays 22229200-02, merges 22229203-05, check 22229206; GPU G9_r0 22229207 (must = T33_G8X_r0), G9A_r 22229208, G9A_s 22229209, G9A_u 22229210, robustness 22229211-14 -> results/v7_T35_*.json. Adopt G9A if chi2 <= T33_G8X + 4.2.
- USER: shaft-Ca licence chosen (T35), W2X jobs cancelled. T33 step 2 (rho_redesign/submit_t33_refit.sh): cexp re-measure new VDCC 22224026-31 -> /scratch/dhuruva/split1/cexp_ljp25g0321 (eCB pool scale_E reads cexp); fits afterok extracts+merges: W2X v7x r0 22224032, s7 22224033, n3 22224034, u 22224035; G8X v8 r0 22224036, s 22224037, u 22224038 -> rho_redesign/results/v7_T33_*.json. Refs W2_X_s7 278.0/311.4, W2_N3 276.4/309.8, G8X_G1a_u6 292.2/321.4.
- 7-h plan v2 (2026-10-02, user revised): split2 only. R0 split2 inputs (chain + bases + cexp) + kernel gpu_v9 + anchors research; R1-R5 five streams (A licence cpre/cpost spine Ca, B licence shaft rel. own bAP, C eCB on shaft Ca, D NO presyn LTP, E baselines/integrator), 25 models x (seeded + unseeded); last hour prefire BCL of the winner with GluSynapseV9. See DECISIONS 2026-10-02.
- R0 split2 chain s2g0321 (spine/submit_s2g0321.sh, resubmitted 22230193-22230235 after 1-min limit bug, ~390 CPU-h): extracts /scratch/dhuruva/s2g0321/extracted/*_delta-split2-ljp25g0321-*, bases /scratch/dhuruva/s2g0321/basis_{l5l5,l23l5,l23l23}, CEXP_DIR /scratch/dhuruva/s2g0321/cexp, K /scratch/dhuruva/s2g0321/k.txt. Kernel gpu_v9 + ANCHORS_V9.md agents running.
- Kernel is rho_redesign/gpu_v10_rho.py (gpu_v9 = T35-G9 gate_rel kernel of an earlier session; its split1 jobs 22229204-214 left alone, owner unknown). Tests 22230184-88 (a must equal 292.206). Anchors rho_redesign/ANCHORS_V9.md: spine licence n_L(cpre+cpost) n_L [1.8,2.4] a20=a21, 100 ms low-pass; shaft kappa [3.0,4.6] (Kampa 2006); eCB k_E rel [0.65,1.0], tau_e [25,125]; NO tau 6.7 ms fixed, gated by LTP licence, d_NO_max >= 0.4. Pending kernel additions: low-pass licence option, NO gated by licence (no_mode 3), tau_E1 freeable.
- T39 audit 22230444-47 (afterok caches): /scratch/dhuruva/s2g0321/doublets/{l5,l23l5,l23l23}_{cpost_spikes.csv,exclude_pairs.txt,keep_pairs.txt}, rejections.csv. Re-run audit_fin after the prefire jobs to add guardrail failures. R1 fits exclude pairs as follows. L5: --pairs = subset24 intersected with l5_keep. L2/3 to L5: replace --joint with --extra l23:paired_l23l5:<dirs>:$L23_BASIS_DIR:$GEOM_L23:<l23l5_keep>. L2/3 to L2/3: add a 6th field to --extra l23l23 pointing at <l23l23_keep>. The Letzkus yaml already uses 1 Hz with 100 reps.
- Phase-2 kernel tests 22231539-45 (lic_lp/tau_L, no_mode 3 as a rate dd_NO/dt = A_NO*N, tau_E1 free, l23 tier-2 bAPs for lic5/ecb1). split2 chain: bases and caches done; sims sj03/sj07/extra, prefire and extracts are still waiting on queue priority (rrg fairshare 0.37). Watchers: b1jlgtlz1 (chain), b9xf33mcq (audit), b4dn1yx6x (tests2).
- R1 (T41) resubmitted 22235493-502 (first submit 22233493-502 cancelled by a dep) (r1{A..E}_{s,u}, run_r1_fit.sh, 85G 0:45, afterok on the ext/cexp/audit/kernel tests): A lic4 + a20 tied (8 free), B lic5 + kappa tied (8), C ecb_src1 with k_E [0.65,1] (7), D no_mode 3 (9), E G8X baseline on split2 (7). New default-off SET option lic_tie in gpu_v10. Results: rho_redesign/r1*.json, logs/r1X_*.out.
- T39 keep lists /scratch/dhuruva/s2g0321/doublets/*_keep_pairs.txt (exclude c_post doublets only; L5 22, L2/3->L5 102, L2/3->L2/3 120); the strict whole-pair lists are *_strict.txt. Split2 L2/3->L5 Letzkus 3AP bursts fail in about half the cells (2 of 3 spikes at 200 Hz).
- R1 done: winner D (G8X + NO mode 3) 297.6, BIC 335.2 (r1D_s.json); baseline E 329.6. R2 = 5 NO3 variants (technician); V9.mod being written (G8X + NO3).
- V9 BCL: rho_redesign/submit_valid_v9.sh [-n dry] [-p pilot] FIT TAG (full run 185 CPU-h, 15 jobs); pilot 22243712/13 on r1D_s.
- USER 2026-10-02: Egger 1999 validation only (drop it from fits from R3 on). Companion updated to R1D and pushed (4f06d89).
- USER 2026-10-02: 65-target fit rejected (flattens to 1, Markram 10 Hz LTP missed). Staged core-target fitting (CORE_TARGETS.md, 18 core, the rest validation). Stage 1 (7 L5 STDP-shape targets; S1A Chindemi, S1B + shaft licence, S1C R1E rule) being submitted.
- Known limit: distal basal (>130 um) bAP Ca ~3x below Kampa/Gordon in every VDCC setting (emodel bAP); possible ion_fitter item.
- Other agents: none. Team: see DECISIONS 2026-10-02 (glusyn-* agents, sjob.sh, guard hook).

## Best fits (65 targets, 3 pathways; results in rho_redesign/results/)
| fit | licence | chi2 / BIC | k | prefire BCL (live) | full-protocol BCL | neurodamus |
|---|---|---|---|---|---|---|
| v7_W2_N3 | VDCC charge theta_V | 276.4 / 309.8 | 8 | PASS 66/68, live 286.8 | no | no |
| v7_W2_X_s7 | VDCC charge, exact timing | 278.0 / 311.4 | 8 | no | no | no |
| v7_G8X_G1a_u6 | shaft Ca >= 0.15 uM | 292.2 / 321.4 | 7 | PASS 68/68, live 298.1 | no | no |
Formalism: rho_redesign/EXP_PARAM.md sections 4 and 4b. BCL notes: bcl_validation/BCL_W2_N3.md, BCL_G8X_G1a_u6.md.

## VDCC / spine Ca findings (cexp/)
- Synapse mod = Chindemi exactly (CHINDEMI_SPINE_CA.md). Chindemi VDCC literature-only, never fitted, no Ba->Ca shift (CHINDEMI_VDCC_FIT.md); thesis has no numbers (CHINDEMI_THESIS_VDCC.md).
- Sabatini/Chindemi Supp Fig 2 replication (SABATINI_VALID.md): synaptic Ca passes; bAP Ca ljp0 0.79 uM FAIL, E1 1.22 PASS (E1 gca tuned on this target).
- delta-split1 (measure_cexp) basal bAP only ~3-20 mV below Chindemi at <60 um (T31); falls to 48 mV at 100-150 um.
- User plan (2026-10-02): vshift approved as calibration, not free fit: (1) shift from Ba->Ca (25 mV), (2) gca to Sabatini 1.7 uM + Koester, (3) refit 6 Chindemi params. Tasks T32, T33.

## Waiting on user (tsk blocked)
T36 spine neck E2 + Zilberter LTD band, T37 delete ckpt.npz. Also: archive T10/T11/T12?

## Next
T31 result -> T32 calibration -> T33 re-prefire + refit (log heavy runs first) -> T34 full-protocol BCL -> T38 tex + push.
