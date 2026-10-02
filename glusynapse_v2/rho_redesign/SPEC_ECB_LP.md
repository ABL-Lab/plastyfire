# SPEC_ECB_LP: eCB trigger from a 10 ms low-pass of the own VDCC Ca (gpu_v11_rho.py, ecb_lp_mode; 2026-10-02)

User-approved idea (EBNER_HOW section 5 item 2): replace the 25 ms bAP veto with Ebner's pre-LTD read-out, translated
to Ca. The default is ecb_lp_mode 0, which runs v11 operation for operation (the new branches are compile-time dead).

## Rule (per synapse i; notation as SPEC_V11)
- **L**, the own eCB low-pass: L' = -L / tau_lp + x / i_scale. x is the same input as the existing trigger pool: s (the
  own spine VDCC current) with v5_mode 1, or (1 - b) s with v5_mode 2. Synapse-local, uniform, with no voltage involved.
  L is in pool units, so one isolated bAP lifts L by about P1_i and L(Delta ms later) ~ P1_i e^(-Delta / tau_lp).
- **Trigger** at each own arrival t_a (exact with t_exact 1, the same v7x correction as V / W): L(t_a) > theta_eCB x P1_i.
  This is the existing thE, so no new threshold is added and theta_eCB stays the only free eCB parameter. For a single bAP
  Delta ms before t_a, the trigger fires while Delta < tau_lp ln(1 / theta_eCB). That gives 13 ms at theta 0.26 (s3C_u),
  10 ms at 0.37 (s3W_s) and 27 ms at 0.07.
- **ecb_lp_mode 1** (Ebner, no veto): every coincident or post-before-pre arrival triggers.
- **ecb_lp_mode 2**: the LP trigger plus the configured veto. Use it with the anchored v11 window (veto_t0 2, veto_Tv 17).
  The v7 25 ms veto is the one that kills dt 0 and +25.
- **ecb_lp_tau**: 10 ms by default (anchored below), not free. It is SET only, for a sensitivity check.
- The eCB step itself is unchanged (A_eCB 1, d_min -0.29).

## Why mode 2 and the (1 - b) weighting are needed (inferred, to be checked by the test)
- **The first pre spike of a -10 train looks the same from behind at 10 Hz and at 50 Hz.** In both, one bAP arrives
  ~8.5 ms before t_a (lag L lands at L - d_syn + d_bAP, d_syn 1-2 ms, ANCHORS_V11 section 2).
  - With A_eCB 1, one trigger gives the full d_min. A pure read-back (mode 1) then puts 50 Hz -10 at rho-ratio x 0.71,
    while the data say 1.70.
  - Ebner avoids this with a graded per-spike LTD that the LTP arms outweigh. Our step is all-or-none.
  - Only a forward look separates them: the next bAP at t_a + 8.5 ms in 50 Hz -10. The narrow window (2, 17] keeps that
    look but excludes the dt 0 bAP (t_a - 1 ms) and the +25 bAP (t_a + 23.5 ms). This is what makes 20/40 Hz dt 0 and
    20 Hz +-25 LTD-capable.
- **The end-of-train leak** (HARD_TARGETS A: the last -10 pre has no bAP after it, so nothing vetoes it) is removed by the
  weighting, not by the window.
  - With v5_mode 2, Bg exceeds 1 after the first arrival at 40/50 Hz (tau_d_NMDA 70 ms). Later bAPs then add nothing to
    L, and the last pre reads L ~ 0.
  - At 20 Hz, Bg is still below 1 for the first bAPs, so 20 Hz -10, dt 0 and -25 keep their trigger.
  - Unweighted (v5_mode 1, as in 3V/3W), the leak stays, so W1 / W2 are expected to lose 50 Hz -10 and S07 pre-only again.
- **Predicted L(t_a) / P1, one bAP:**

  | Rows | bAP position | L(t_a) / P1 |
  |---|---|---|
  | dt 0 | t_a - 1 ms, partly before t_a | ~0.5 |
  | -10 | 8.5 ms back | ~0.43 |
  | 40 Hz +10 | 13.5 ms back | ~0.26 |
  | +-25, 0.1 Hz -25 | 23.5 ms back | ~0.095 |
  | 20 Hz +10 | 41.5 ms back | ~0.016 |
  | Markram +10 | 91.5 ms back | ~0 |

  So 20 Hz +-25 and 0.1 Hz -25 need theta_eCB below ~0.09. The seeds' theta (0.26, 0.37) cannot reach them, which is why
  the test adds C2t at theta 0.07. Refitting theta_eCB (with the Chindemi parameters) is the real test.
- **Expected loss:** sjostrom_burst5x20hz_r50_dt-120/-200 (0.79). A 10 ms LP cannot see a bAP 120 ms back (e^-12); only
  the 100 ms pool can. Sjostrom 0.1 Hz -120/-200 (1.05, no LTD) is consistent with the LP.

## Anchors
| Quantity | Value | Source |
|---|---|---|
| tau_lp | 10 ms | Ebner, Clopath, Jedlicka & Cuntz 2019 Cell Rep 29:4295, doi 10.1016/j.celrep.2019.11.068, ModelDB 251493 (syn_4p.mod pre-LTD: rectified u low-pass, tau 10 ms, read by a delta pulse at each pre spike). Model-derived: Ebner set it against Sjostrom 2001 and Nevian 2006. |
| post-before-pre window, L5-L5 | LTD at -10 (0.69) and -25 (0.65), none at -120/-200 (1.05), 0.1 Hz | Sjostrom, Turrigiano & Nelson 2003 Neuron 39:641, doi 10.1016/s0896-6273(03)00476-8. Window width = tau_lp ln(1/theta_eCB), so -25 fixes theta_eCB <~ 0.09 at tau 10; the window is set by eCB handling (FAAH block widens it), not by Ca decay. |
| Ca must precede mGluR activation; LTD Ca threshold ~1/2 the LTP one; coincidence excess decays tau ~29 ms | | Nevian & Sakmann 2006 J Neurosci 26:11001, doi 10.1523/JNEUROSCI.1749-06.2006. Supports a read-out at t_a of Ca that came before it (L), and the (1 - b) order weighting. |
| eCB tLTD window -100 to +25 ms (L4->L2/3) | | Bender, Bender, Brasier & Feldman 2006 J Neurosci 26:4166, doi 10.1523/JNEUROSCI.0176-06.2006. Wider than tau 10; reached only with a small theta_eCB. Upper side (+25) agrees with the (2, 17] veto letting +25 through. |
| veto window (2, 17] ms (mode 2) | | ANCHORS_V11 section 2 (Sjostrom 2001 lag table, Bender +5, Zilberter +4). |

No new free parameter. theta_eCB (existing) and the Chindemi parameters stay the only free ones.

## Mod implementation (later GluSynapseV10.mod, new name)
- Add STATE L_eCB: L_eCB' = -L_eCB / tau_lp + w * (-ica_VDCC) / i_scale, with w = 1 - min(Bg, 1) under v5_mode 2.
- In NET_RECEIVE (own pre): if L_eCB > theta_eCB * P1, then mode 1 applies the dpre step at once; mode 2 starts the
  SPEC_V11 section 4.3 pending/window mechanism with this trigger in place of the c_VDCC one.

## Test (test_ecb_lp.sh, outputs /scratch/dhuruva/ecb_lp_test)
- REPRO11: v11 core "a" against s1C_s.csv (1e-6).
- L5-all fixed-parameter rescores:
  - C0 / C1 / C2 / C2t: s3C_u, weighted (v5_mode 2). Mode off / 1 / 2 / 2 at theta 0.07.
  - W0 / W1 / W2: s3W_s, unweighted (v5_mode 1). Mode off / 1 / 2.
- Report: focus rows (20/40 Hz dt 0, 20 Hz +-25, Markram +-10, 50 Hz -10, S07) and the total L5 chi2 per run.
