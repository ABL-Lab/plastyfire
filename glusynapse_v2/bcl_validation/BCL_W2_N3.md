# Prefire BCL validation of v7_W2_N3 (V7vn2 joint fit, delta-split1; 65 targets, 3 pathways)

Live GluSynapseV7 (live_v7.py: bAP veto 25 ms, ecb_ref 2 theta_eCB,i = k_E x 1e5 vdcc_q_post) on every record behind the 65 targets (L5->L5 40 incl. paired_l5_extra, L2/3->L5 9, L2/3->L2/3 16 Zilberter + Egger), Banerjee 2014 (validation only) and the 2 dropped Letzkus targets, against the CPU port of the gpu_v7_rho kernel (compare_prefire_v7.py). Fit: k_E 0.4869, theta_V 3.575, gamma_d 47.18, gamma_p 343.86. Outputs: /scratch/dhuruva/bcl_valid_v7_W2_N3/ (cmp_full_{records,targets}.csv, cmp_full.png/.pdf); log logs/valid_v7_cmp_22186907.out.

## Verdict (2026-10-01)
**PASS 66/68**, 3469/3469 records ok. chi2 over 65 fitted targets: live 286.77, offline on the same records 276.39 = fit csv (max |offline - offline_fit| 0).
- L5->L5 (908 records): live 137.90 vs offline 128.39 (40 targets); rho disagreements 1 / 7044; eCB steps 53574 / 53789, vetoes 58555 / 60534.
- L2/3->L5 (564): 19.25 vs 19.12 (9); rho 4 / 3650; eCB 7041 / 7084, vetoes 11878 / 11933.
- L2/3->L2/3 (1963): 129.62 vs 128.88 (16); rho 3 / 7769; eCB 143281 / 142808, vetoes 66240 / 66094.
- nvpend 0 and uE equal live / offline on every synapse.
- Banerjee 2014 (validation only): 0.730 vs data 0.77 +/- 0.07, match.
- FAIL: L5 sjostrom_20hz_dt0ms, live 1.073 vs offline 1.027 (d +0.046, tol 0.025), and L5 sjostrom_40hz_dt0ms, live 1.130 vs offline 1.109 (d +0.021, tol 0.020). sjostrom_0.1hz_dt0ms passes only just (d +0.034, tol 0.035). Together the three account for most of the +9.5 L5 chi2 gap.

## dt 0 diagnosis
- Only dt 0 is affected. Triggered arrivals (steps + vetoes, live / offline): 0.1 Hz 700 / 1400, 20 Hz 2105 / 2624, 40 Hz 2242 / 2628. Every other L5 protocol agrees within about 3%: dt +/-10, the 5 ms trains.
- So the trigger (W > theta_eCB,i at the arrival) is the main difference, not only the veto window. Offline triggers more, so it makes more steps and gives lower ratios.
- Mechanism: the kernel reads W at the first 0.25 ms sample at or after t_a, so W includes the VDCC input of [t_a, t_k] at full weight. Live reads W(t_a) exactly, and after the arrival bglu = 1 sets the W input to 0. At dt 0 the bAP's VDCC current starts inside that sub-bin, so W sits at the threshold edge. The veto Q has the same sub-bin offset, at the start of its window.
- Test: diag_dt0_v7.py, job 22190200 (0:42, 2.60 GB). Triggered arrivals, live / V0 kernel / V1 exact t_a / V2 sample before t_a: 0.1 Hz 700 / 1400 / 627 / 550; 20 Hz 2105 / 2624 / 2242 / 2132; 40 Hz 2242 / 2628 / 2365 / 2319. Synapses matching live under V1: 176/178, 137/172, 141/178. The hypothesis is partly confirmed: V1 removes most of the gap, and a residual remains. The nearest post spike is 0.6-4 ms before t_a, so the bAP reaches the synapse before the arrival, and W rises steeply inside the 0.25 ms bin, where left-point sampling is crude.
- dt 0 is ill-conditioned. The trigger compares W(t_a) with theta_eCB,i while the bAP's VDCC influx is rising at the synapse, so sub-ms differences in bAP timing at the synapse (dendritic distance, AP latency, jitter) switch the eCB step on or off. This is a physical property of the rule, not only a numerical one: the dt 0 predictions are sensitive to the bAP arrival time, and live / offline agreement there is limited by the sampling.

## Fix (gpu_v7x_rho.py, gpu_v7_rho.py unchanged)
SET t_exact 1 (with veto_T > 0): per arrival entry f = (t_k - t_a)/h_j (host, veto_q order). Trigger W(t_a) = W(t_k) - f x the bin-j W input, and veto Q + f x the bin-j unweighted input, i.e. W just before t_a and the window [t_a, t_a + T] (= V1). Run with KERNEL=gpu_v7x_rho.py run_joint_w2.sh, using the W2_N2 env.
- 22190861 W2_X_r: MAXITER 0 at the W2_N3 params, 78G 0:15. It ends in REPRO DIFF (exit 1) by design. Compare its chi2 with live 286.77 / offline 276.39 (per pathway: live L5 137.90, L23 19.25, L23L23 129.62).
- 22190862 W2_X_s5: refit seeded from W2_N3. 22190863 W2_X_s7: unseeded. Both 78G 0:30, independent of the rescore.

## Runs
- Gate 22174635. Shards (32 CPU, 84-95% eff): L5 24-29 min, 14.5-15.4 GB; L2/3->L5 20-22 min, 17.8-18.5 GB; L2/3->L2/3 13-16 min, 11.5-12.3 GB. Next: 20G / 23G / 16G.
- Compare 22174649 timed out at 0:15 (15.7 GB). 22186907: 12:26, 15.75 GB of 25G. Next: 20G, 0:30.
