# SPEC_V11: eCB veto window, peak veto, MVD arm (gpu_v11_rho.py, 2026-10-02)

gpu_v11_rho.py is a copy of gpu_v10_rho.py (v10 is not edited). With every v11 option at its default, the kernel runs the v10
operations in the same order, so predictions are bit-identical. test_v11.sh `core` checks this against s1C_s.csv (tolerance 1e-6).
NO (no_mode 1/2/3), the licences (lic_src) and the shaft eCB (ecb_src) are unchanged. All drives are synapse-local.

## Notation (per synapse i, kernel grid step k of length h_k ms)
- s_k: own spine VDCC current (step mean of -ica_VDCC). V: unweighted own VDCC pool, V' = -V/tau_E1 + s/i_scale (the existing c_VDCC).
- b = min(Bg, 1): own glutamate-bound state. Bg += 1 per own arrival; Bg' = -Bg/tau_d_NMDA (70 ms). Same as the v5b weight.
- P1_i = 1e5 x cexp vdcc_q_post_i (ecb_ref 2): the own single-bAP pool jump. Missing rows get the pooled median (fit_v6 rule).
- dep = [c* > theta_d]. pot = [c* > theta_p] x licence. Existing rho' = -rho(1-rho)(rho*-rho) + pot gamma_p (1-rho) - dep (1-pot) gamma_d rho.

## 1. Veto window and peak veto (spine eCB, ecb_src 0 only)
| option | default | meaning |
|---|---|---|
| veto_t0 | 0 ms | window start after the exact arrival t_a (t_exact 1) |
| veto_Tv | = veto_T | window end. The window is (t_a + veto_t0, t_a + veto_Tv] |
| veto_peak | 0 | 0: integral; 1: transient peak |
| veto_peak_k | 0.5 (placeholder; free-able, box 0.1-3) | peak threshold in single-bAP peak units |
- Any non-default value moves the veto sum into the kernel at the candidate's tau_E1. Each grid step is weighted by its overlap with the window, including the part of step k-1 after t_a.
- Integral (veto_peak 0): Q = sum_j tau_E1 (1 - e^(-h_j/tau_E1))/i_scale x s_j x overlap_j/h_j, and the trigger is vetoed if Q > theta_eCB P1_i (v7 Q, new window).
- Peak (veto_peak 1): Q = max s_j over the steps overlapping the window, and the trigger is vetoed if Q > veto_peak_k x I1_i.
  - I1_i is its own reference in peak units, not k_E P1_i. It is the own single-bAP peak VDCC current: the median over isolated bAPs of max s over [x-1, x+8] ms.
  - Isolation follows v7 ecb_ref 1: no other post spike in [x-100, x+20]. Tier 1 also requires no own arrival in [x-150, x+8]; tier 2 allows arrivals and is used only where tier 1 is empty; the rest get the pooled median.
  - It comes from the records, because cexp has no peak column. The log prints the I1/P1 ratio per pathway. Adding a cexp `vdcc_ipk_post` column would give one calibration source with P1.
- Anchored values (ANCHORS_V11 2): **veto_t0 2, veto_Tv 17, veto_peak 1, veto_peak_k 0.5**. Pass them in SET; the kernel defaults stay at v10 for repro. t_a here = pre spike + NetCon delay (model_v2.arrival_index), i.e. glutamate onset at the synapse, so a protocol lag L puts the bAP at t_a + L - delay + d_bAP. t0 2 is confirmed once C_MVD_CALIB lists the dt-0 and +4 bAP peak lags (dt 0 must fall before t_a + 2, +4 after it).

## 2. Unweighted eCB trigger
Already available: **v5_mode 1**. In every kernel (v5..v11), v5_mode only switches the trigger between W, the (1-b)-weighted pool (mode 2), and V (mode 1). Nothing else depends on it, so no new option was added.

## 3. MVD, "mGluR-gated VDCC depression" (mvd_mode 1; needs ecb_ref 2)
- M_i(t) = [lo_i < V <= theta_MVD_hi P1_i] x [b > b_c]: own VDCC pool above the lower edge (no upper edge by default), own glutamate window open. b = min(Bg, 1) is already normalised (jump 1 per arrival, b_peak = 1); b_c = exp(-mvd_W/tau_d_NMDA) = 0.49 at mvd_W 50 ms (Marcaggi 2009), mvd_b >= 0 overrides it.
- mvd_pot 0 (default): D = max(dep, M), and rho' = ... - D (1-pot) gamma_d rho. MVD is silent during licensed potentiation, exactly like dep.
- mvd_pot 1: rho' = ... - max(dep (1-pot), M) gamma_d rho. MVD also competes with potentiation.
- MVD adds to the depression *indicator*, so gamma_d enters once (no double count). A failed licence still counts as dep, as in v8.
| option | default | meaning |
|---|---|---|
| mvd_mode | 0 | 1 = MVD on |
| theta_MVD_lo | 0.5 (ANCHORS_V11 midpoint; 0.37 is NOT an MVD edge) | lower edge x R_i (mvd_ref 0/1); free-able, box 0.2-1.0 |
| theta_MVD_abs | 100 (placeholder) | lower edge for mvd_ref 2, absolute c_VDCC pool units (1e5 x charge, P1 units); box 1-1e4 log; value from C_MVD_CALIB |
| theta_MVD_hi | 1e30 (= no upper edge) | upper edge, x R_i; free-able, box 0.5-10 (give a SET start value if freed) |
| mvd_ref | 0 | 0 = P1_i; 1 = P1_i sum_{k<n} e^(-k isi/tau_E1) (not recommended: per-protocol, ANCHORS_V11); **2 = absolute theta_MVD_abs** (recommended if C-MVD shows a gap) |
| mvd_ref_n, mvd_ref_isi | 10, 20 ms (placeholders) | train used by mvd_ref 1 |
| mvd_W, mvd_b | 50 ms, -1 | glutamate window; b_c = exp(-mvd_W/tau_d_NMDA); mvd_b >= 0 sets b_c directly |
| mvd_pot | 0 | see above |
- **Edges (ANCHORS_V11 1).** 0.37 is the D890 train-Ca ratio, which constrains the LTP licence (theta_VDCC in (3.0, 4.4] P1), not MVD. MVD has no upper edge (APV train and train-LTD are LTD at the full pool), so theta_MVD_hi stays 1e30. In P1_i units no uniform lower edge exists (Zilberter LTD at <= 1-1.8 P1, L5 AM251 no LTD at <= 2.3 P1); hence mvd_ref 2.
- **mvd_ref 1.** A per-synapse pool *measured* in the control protocol is not available for every synapse: only the L2/3->L2/3 records contain Zilberter trains, and L5 or L2/3->L5 synapses have none.
  - A self-normalised reference (the row's own peak in the same protocol) would be non-causal and cannot go into a mod.
  - mvd_ref 1 therefore uses the linear train sum of the own P1_i, which is synapse-local and causal. It is a rescaling of P1_i by G(n, isi, tau_E1), so it matters only when tau_E1 is free or when lo/hi are anchored as train ratios.
- theta_MVD_abs comes from C_MVD_CALIB.md (job 22293154). theta_MVD_lo 0.5 is the ANCHORS_V11 midpoint (sensitivity 0.2 / 0.9). If C-MVD shows no gap, MVD cannot be uniform with this emodel (user decision).
- Counter 6 (v10_counts 1) counts MVD-depressing steps: steps where M raised the depression drive.

## 4. Implementation in a later GluSynapseV10.mod (new name; for the full-network BCL run)
1. STATE c_VDCC (unweighted, tau_E1, i_scale), Bg (decay tau_d_NMDA; NET_RECEIVE own pre: Bg = Bg + 1). PARAMETER/RANGE P1, I1,
   theta_MVD_abs, theta_MVD_lo, theta_MVD_hi, b_c, veto_t0, veto_Tv, veto_peak, veto_peak_k. P1 and I1 are set per synapse from the same
   calibration as the kernel (cexp vdcc_q_post x 1e5; record-based I1, or a future cexp vdcc_ipk_post).
2. MVD in the rho DERIVATIVE: mvd = (c_VDCC > lo)*(c_VDCC <= theta_MVD_hi*P1)*(Bg > b_c), lo = theta_MVD_abs (mvd_ref 2) or theta_MVD_lo*P1 (0);
   D = max(depress, mvd); rho' = (-rho(1-rho)(rho*-rho) + pot gamma_p (1-rho) - D (1-pot) gamma_d rho)/(1e3 tau_ind).
   Keep the kernel's (1-pot) convention and the licence masking of pot.
3. Veto window:
   - At an own pre arrival: if c_VDCC > theta_eCB P1 (v5_mode 1 trigger), set pending = 1 and Qv = 0. Then net_send(veto_t0, 2), which opens the window (win = 1), and net_send(veto_Tv, 3), which closes it.
   - At flag 3: if Qv <= threshold, apply dpre = dmin + (dpre - dmin)(1 - A_eCB). Then win = 0 and pending = 0.
   - While win: in integral mode Qv' = -ica_VDCC/i_scale (O(h/tau_E1) from the kernel sum); in peak mode, in BREAKPOINT, Qv = max(Qv, -ica_VDCC) against veto_peak_k I1.
   - The kernel applies the step at t_a using future influx. The mod applies it veto_Tv later. That is harmless while the pre ISI is above veto_Tv (rates < 59 Hz at 17 ms). For faster pre trains, keep two alternating accumulators.
4. Validate the mod against the kernel on single records first: Markram 10 Hz +/-10, Sj07 step, Zilberter train10 -10 AM251.

## 5. Tests (test_v11.sh; outputs /scratch/dhuruva/v11_test)
- `core` (7 L5 core targets): a = defaults (REPRO11 vs s1C_s.csv); b = mvd_mode 1 (lo 0.5 P1, b_c 0.49 when it runs); c = veto_t0 3, Tv 17, peak 1 (submitted before the t0 2 anchor).
- C-MVD (c_mvd_calib.py, run_c_mvd.sh, job 22293154): per-row glutamate-window pool, abs and P1 units -> C_MVD_CALIB.md.
- `l5all` (40 L5 targets, pilot size): a40 / b40 / c40. a40 is compared with s1C_s_val.csv, where small differences may come from the cross-pathway
  pooled-median fills. b40 shows whether MVD breaks the L5 cores (20 Hz -10 AM251, 40/50 Hz LTP, Markram).
- Both modes print the P1 NaN fraction per pathway table (cexp), and per model the `v11 calib` P1 / I1 quantiles, filled rows and the counters.
