# LTP-licence alternatives: data inventory and implementation (A14, 2026-10-01; design in GATE_ALT.md)
## Inventory (split1 npz, extract_v2.py; key check job 22176661: 13 s, 179 MB)
- Every split1 npz (L5 vca, L2/3->L5 and L2/3->L2/3 vseg-rs) holds per synapse: effcai, **shaft_cai**, vdcc (-ica_VDCC step mean),
  cev/vev events, pre/post spikes. extract_v2 refuses traces without shaft_cai. No ica_NMDA, spine cai_CR or v traces are stored;
  spine [Ca] above rest is recovered exactly from effcai (batch_v2.cacr_from_effcai). **No re-prefire needed, so no pilot was run.**
- (1) shaft_cai = `cai` of the segment carrying the synapse (GluSynapse `shaft_cai = cai`), point-sampled on the effcai grid:
  0.25 ms everywhere in L5 Markram, 0.25 ms in spike windows and 5 ms outside them in the windowed L2/3 records. Rest = 65 nM
  (cad cainf) before the first stimulus. Markram 10 Hz -10 peaks at 0.56-1.03 uM; Letzkus 1AP +10 (L2/3->L5 tuft) peaks at 0.065-0.066 uM,
  so one bAP adds about 1 nM there. Zilberter peaks range 0.07-16 uM.
- (2) plastyfire/GluSynapse.mod and GluSynapseV4/V5/V7 only `USEION ca READ cai`, with no WRITE ica. Synaptic NMDA/VDCC Ca stays in the private
  spine pool cai_CR, so **shaft Ca is not NMDA-contaminated**. It is the emodel's dendritic Ca_HVA2/Ca_LVAst influx through `cad`.
- (3) delta-split1 L5TPC: `cad` (modified Kampa CaDynamics) in basal and apical with the same parameters (gamma_0 0.24, kE 62,
  cainf 65 nM, diam/4 shell). Decay tau = (1+kE) diam / (4 gamma_0), about 66 ms per um of diameter, so trunk, oblique and tuft
  differ only through diameter. Ca_HVA2/Ca_LVAst are uniform within basal (2.39e-3 / 2.39e-3) and within apical (3.19e-3 / 1.75e-3), with
  no distal hotspot and trunk = tuft. Soma and axon use CaDynamics_DC0 (decay ~282 ms). L2TPC/L3TPC have the same layout.
- (4) Spearman of single-bAP shaft dCa against cexp vdcc_q_post, per pathway and pooled: printed as `v8 calib` lines in the G1r015 / G5r050 logs, together with the
  single-bAP peak by cexp loc, the proximal (<= 100 um) median against 0.15-0.30 uM, and the 3AP 200 Hz record peaks.
## Implementation: gpu_v8_rho.py (subclass of GPUModelV7; v7, v6, fit_v6 not modified)
- `gate_src` in SET: 0 = Vg (pure v7 path, bit-for-bit); 1 = spine [Ca] - min_ca_CR (uM, from effcai); 2 = shaft dCa above
  the synapse's pre-stimulus rest (uM; on the GPU as uint16 nM, +6.9 GB, 31 -> 38 GB of 3g.40gb); 3 = c* itself (uM s = mM ms).
- `gate_win` > 0 (src 1/2): peak-in-window licence `[max x over [t - win, t] > theta_G]` (GATE_ALT G1, win 100 ms). `gate_win 0`:
  100 ms low-pass `dG/dt = (x - G)/tau_E1`, licence `[G > theta_G]` (G2). A failed licence counts as depression, as in V7vn2.
- Threshold theta_G: one natural-unit value for every pathway. `gate_theta` > 0 fixes it (anchored, k 8 -> 7). Otherwise
  theta_G = kappa x theta_V (DE slot), with kappa = Q50 / 3.575. Q50 is the median over synapses (pooled over the 3 models) of each synapse's gate
  quantity at the moments its Vg crosses W2_N3's theta_V, so a W2_N3 seed starts at the converted median-synapse threshold.
  `gate_kappa` pins kappa. The json gets a `v8` block {gate_src, gate_win, kappa, Q50, theta_G, unit}.
- For src != 0, fit_v6's REPRO DIFF at MAXITER 0 is expected and caught (exit 0). run_joint_w2.sh takes env KERNEL (default v7).
## Jobs (submit_gate_v8.sh; 65-target W2 env; outputs results/v7_G8_*.json)
- R0 22177165: gate_src 0 MAXITER 0 on W2_N3 -> must print `REPRO OK` (78G 0:15).
- Rescores at W2_N3 params: G1r015 22177166 (G1, theta 0.15 uM) and G5r050 22177168 (G5, 0.5 uM), both 87G 0:20, with calibration.
- Refits, afterok 22177166, seeded s5 (W2_N3) / unseeded u6: G1 (src 2, win 100) 22177169/70, G2 (src 2 low-pass) 22177171/72,
  S1 (spine, win 100) 22177173/74, C3 (c*) 22177175/76. G1a 22177177 = G1 anchored at 0.15 uM, k 7.
  Sizes: src 2 87G 0:40, S1 78G 0:35, C3 78G 0:30. Run seff on each: the src-2 size is estimated (measured 61.6 GB + 8 GB of traces), not measured.
## Results and v7x wave (2026-10-01)
- G8: R0 REPRO OK. G1_s5 285.41 (theta_G 0.169 uM); G1a 290.41 (k 7, 0.15 uM); G2 302.7; S1 304-314; C3 306-307; G5 419. W2_N3 = 276.39.
  Calibration: Spearman 0.86 vs vdcc_q_post; proximal single-bAP shaft dCa 0.53 uM (L5), 0.23 uM (L23L23). src-2 MaxRSS 68.0 GB.
- GPUModelV8 now inherits gpu_v7x_rho.GPUModelV7X. SET t_exact 1 applies the v7x exact-arrival trigger and veto corrections in the v8 kernel;
  t_exact 0 keeps the v7 arithmetic. submit_gate_v8x.sh (87G): G1a_r0 22192847 (t_exact 0 MAXITER 0 on G8_G1a_s5, log must print REPRO OK),
  G8X_G1a_s5 22192848 / G8X_G1a_u6 22192849 (0.15 uM fixed, k 7), G8X_G1_s5 22192850 (theta free, seeded G8_G1_s5); fits 0:35.
