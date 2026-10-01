# Why C1Ajn_s5 scores chi2 110.27 on L2/3 -> L2/3 (Zilberter 2009)

Diagnosis only. Nothing in the model, fitters or edges was changed.

**Inputs**
- Score: pilot 22134473 (fit_v4n, 3 pathways, MAXITER 2; the DE never left the seed). Files: `results/v4_C1Ajz_pilot{.json,_l23l23.csv}`, `logs/v4_fit_22134473.out`.
- Rule: C1 (vamp_mode 1, theta_V 5.50), t_drive 4, tau_E1 100 ms, A_mglu 0.308, theta_Tg 6.18, tau_T 30.4 ms, dpre_min -0.245, theta_Z 1.75, tau_Z 48.7 ms, rho_gamma 0.52, gamma_d 51.7, gamma_p 158.4.

**Script**
- `rho_redesign/diag_l23l23.py` + `run_diag_l23l23.sh`.
- It replays the kernel record by record on the CPU, in the same step order as `diag_dltd.py`.
- In the same pass it evaluates parameter-free rule variants at the fixed C1Ajn_s5 parameters (section 3).

**Run 1: 22135726 / 22135727 / 22135728 (l23l23 / l5 / l23)**
- The replay reproduces every fit-csv prediction: max |diff| is 3e-16 (l23l23), 3e-12 (l5) and 1.5e-9 (l23).
- Each job then crashed in the chi2 print (a `D.sem` bug), after the target table had been printed. So the target-level numbers below are final.
- The per-protocol local-quantity tables were not written.

**Run 2: 22135934 / 22135935 / 22135936** (fixed script, plus the latched-gate variant GbL)
- It writes `/scratch/dhuruva/l23l23_diag/<model>_{syn,rec,targets,proto}.csv`. Logs are `logs/diag_l23l23_<id>.out`.
- Results are in section 5. The v5c follow-up (Run 3) is in sections 6–7.

## 0. The three worst targets: model and data were swapped in the brief

| target | data | model | z |
|---|---|---|---|
| train10 +4 ms last (train-LTP) | **1.49** ± 0.12 | 0.89 | -5.0 |
| train10 -10 ms last, AM251 (mglu_block) | **0.73** ± 0.07 | 1.06 | +4.7 |
| 1ap -10 ms | **0.56** ± 0.06 | 0.76 | +3.4 |

The model gives **too little** LTP at train +4, and **no** LTD under AM251 where the data keep it.
- In the model, train-LTD is 100% eCB (pre-side dpre), so AM251 removes it.
- In the data it is AM251-insensitive and postsynaptic (PPR and CV unchanged), so it has to come from rho.

## 1. Decomposition (replay, l23l23, 15 targets)

**NO is off in every Zilberter protocol.**
- NO needs Z > theta_Z = 1.75, where Z counts the synapse's own arrivals and decays with tau_Z = 48.7 ms.
- One pre spike gives Z = 1; five at 20 Hz settle at 1.56.
- Consequently rho_only == mglu_block exactly. The prediction is rho × eCB with no third term.

**eCB is saturated.**
- ecb_only = 0.765 in every train-with-arrival, 1ap -10 and 5ap 20 Hz -10 protocol.
- That is the floor 1 + dpre_min in the readout: d0 = dpre_min at practically every synapse.

| protocol | data (AM251) | model | rho part (rho_only) | eCB part (ecb_only) |
|---|---|---|---|---|
| 1ap +10 | 0.64 | 0.83 | **0.83** | 1.01 |
| 1ap -10 | 0.56 | 0.76 | **0.98** | 0.79 |
| 5ap 10 Hz +10 | 0.76 | 0.85 | 0.85 | 1.01 |
| 5ap 20 Hz +10 | 1.07 | 1.03 | 1.03 | 1.01 |
| 5ap 20 Hz -10 | 0.93 | 0.78 | 1.03 | 0.77 |
| train4 +4 last | 0.76 | 0.69 | 0.91 | 0.77 |
| train8 +4 last | 1.15 | 0.85 | 1.13 | 0.77 |
| train10 +4 last | 1.49 (1.73) | 0.89 | **1.18** | 0.77 |
| train10 -4 last | 0.99 | 0.80 | 1.05 | 0.77 |
| train10 -10 last | 0.72 (0.73) | 0.80 | **1.06** | 0.77 |
| train10 +5 first | 0.97 | 1.16 | **1.16** | 1.01 |
| train10 post only | 1.03 | 1.06 | 1.06 | 1.01 |

**No eCB at all (control targets scored with rho_only) gives chi2 157.3 rather than 110.3.** Removing eCB from L2/3 does not help by itself. The failure is in the **post rule**:

1. **No rho-LTD where the data put postsynaptic LTD.** 1ap -10 gives 0.98 (data 0.56) and train -10 gives 1.06 (0.73). 1ap +10 is only 0.83 (0.64).
2. **Train rho is set by the AP count, not by pre timing.**
   - By AP count: 0.91 → 1.13 → 1.18 for 4 → 8 → 10 APs.
   - Every 10-AP timing, including the train alone and pre-before-the-first-AP, gives 1.05–1.18.
   - The data need: +4 last 1.5–1.7, +5 first 1.0, -4 1.0, -10 0.73.
3. **The eCB saturates on any bAP before the own arrival.** Zilberter's pharmacology says L2/3 -> L2/3 has no CB1 component, while L5 -10 LTD is CB1-dependent (Zilberter L5 10 Hz -10 under AM251 is 1.07; Sjostrom 0.1 Hz -10 under AM251 is 1.01).

## 2. Mechanism in synapse-local terms

### Pre side: eCB (t_drive 4)

- Each own VDCC event (K_ca crossing, one per bAP) feeds S1 with weight 1 − b. Here b is the synapse's own glutamate-bound NMDA state: 1 at arrival, decaying with tau_d 70 ms.
- A train before the pre spike gives 9–10 events at weight 1. S1 builds up and T ≫ theta_Tg, so each arrival moves d0 31% of the way to dpre_min, and 40 pairings reach the floor.
- One bAP 10 ms before the arrival is already enough (1ap -10: 0.79).
- In this rule eCB means only "a bAP before my glutamate". **No local quantity in the current state separates L5 single -10 (eCB needed) from L2/3 train -10 (no eCB)**:
  - E1 (no eCB at a high c_VDCC pool) also removes the L5 single -10 eCB (0.1 Hz -10: 0.80 → 0.91). So **one bAP already lifts c_VDCC above theta_V within 10 ms**.
  - E2 (eCB only from events into an unprimed pool) cuts train eCB by about half and keeps every single-bAP eCB (table in section 3).

### Post side: rho (C1)

These are hypotheses for Run 2's `*_proto.csv` to confirm (Vpk, V_arr, t_tp, t_pot, t_conv, P_up).
- **Gate.** c_VDCC (the own VDCC current through tau_E1 = 100 ms) crosses theta_V after a single bAP, as the E1 result shows. In any multi-bAP protocol the gate is open for the whole late part of the pairing, whatever the pre timing.
- **Threshold.** effcai (tau 278 ms) integrates the train's VDCC Ca. The own NMDA Ca of one pre spike is a small increment on top of it. Time above theta_p with the gate open is therefore set by the train. That explains the flat 1.05–1.18 and the absence of depression at -10: no "c* > theta_p with the gate closed" time is left for C1 to turn into depression.
- **What is missing.** The rule has no synapse-local signal for "a bAP arrived while my own glutamate was bound", which is what separates +4 last from -4 / -10. (+5 first also has such a bAP, at the start of the train.)
- **LTP ceiling** (P_up in Run 2). Even the most permissive variant gives at most 1.25 at train +4 under AM251 (C2), against data 1.73. If P_up(train10 +4) is already about 1 for rho0 = 0 synapses, 1.73 is out of reach of the readout: only 46% of synapses start at rho0 = 0 (edges check), and the Use_p / Use_d gap caps the gain. Run 2 decides this.

### Q2: train10 +4, too little LTP

It is not over-potentiation. It is two things:
- **Saturated eCB LTD that should not be there:** ×0.765, while the data give ×0.86 within 1.2 SEM.
- **A train-dominated rho rule that tops out at 1.18:** the pre spike's coincidence with the last bAP adds almost nothing on top of the train's c* and c_VDCC.

Ca saturation of the threshold is not the cause: the time above theta_p is large in every train.

### Q3: train10 -10 under AM251

The model's LTD there is 100% eCB, because rho gives 1.06. To move it to rho, the local quantity has to be **high at -10 and low at +4 / L5 +10 pairings**. The candidate is the **b-weighted own VDCC event count (Gb): bAPs that arrive while my own glutamate is bound**, the mirror of the eCB weight 1 − b.

- **It points the right way at -10.** Gb gives train -10 under AM251 0.84 (data 0.73, base 1.06).
- **As a plain "pot needs Gb > 0.5" gate under C1, it is fatal.** C1 converts every c* > theta_p moment outside the gate into depression: during the train before the pre spike, and in the c* tail after Gb decays (tau_E1 100 ms vs effcai 278 ms).
  - L2/3: train +4 drops to 0.49 and train8 to 0.50.
  - L5: sj07 step pair drops to 0.47 (data 1.62) and Markram 10 Hz -10 to 0.54 (0.79).
  - L2/3 -> L5: Letzkus 3AP proximal drops to 0.42 (1.30).
- The latched form GbL is in Run 2. It opens at Gb > 0.5 while c* > theta_p and stays open until c* falls below theta_p, which removes the tail artefact but not the pre-arrival train part.

## 3. Uniform remedies tested at fixed C1Ajn_s5 parameters (Run 1)

All variants below are synapse-local and add no free parameter. They reuse theta_V, tau_E1, tau_d_NMDA and the 0.5 count threshold of the existing (D) vgate. chi2 is per pathway: l5 (30 targets) / l23 = L2/3 -> L5 (9) / l23l23 (15). Values are computed by hand from the printed target tables; Run 2 prints them.

| variant | change | l5 | l23 | l23l23 | total | verdict |
|---|---|---|---|---|---|---|
| base (C1) | – | 46.9 | 5.1 | 110.3 | 162.2 | |
| **E2** | own VDCC events feed eCB S1 only while c_VDCC <= theta_V | ~50.1 (r50 bursts -120 / -200: 0.77 / 0.82 → 0.89 / 0.91, data 0.79) | 5.1 | **98.1** (train +4 1.04, train8 1.00, -4 0.91; -10 0.92 worse) | **~153** | only variant that helps without refit; partial |
| E1 | no eCB step at an arrival while c_VDCC > theta_V | worse: 0.1 Hz -10 0.91, 10 Hz -10 0.91, -25 0.88 | – | 109.1 (1ap -10 0.87 worse) | worse | reject: one bAP primes the pool, so L5 single -10 loses eCB |
| C2 | a failed pot step is neutral | worse: 0.1 Hz +10 1.11, 20 Hz -10 1.02 | 1ap +10 0.97 (data 0.72) | train +5 1.41, 1ap +10 0.95 | worse | reject (as in earlier fits) |
| Gb / Gbp | pot needs Gb > 0.5 | sj07 pair 0.47, 10 Hz -10 0.54 | 3AP proximal 0.42 | ~227 (train +4 0.49) | far worse | reject in this form; -10 AM251 moves the right way |
| C2+Gbp | | 10 Hz -10 0.65, sj07 0.88 | 3AP proximal 0.56 | train +4 0.80, -10 0.65 / AM251 0.86 | worse | reject |
| GbL, GbL+E2 | latched Gb (see Q3) | Run 2 | Run 2 | Run 2 | | pending |

### What this means for the minimal change

1. **eCB.** E2 is the smallest uniform change that moves the pre side toward the data: about -12 on L2/3 and +3 on L5. It keeps every single-bAP eCB and costs only the L5 5x20 Hz r50 bursts. The biological case is weak: 2-AG synthesis rises with Ca, so treat E2 as a phenomenological "eCB from unprimed bAPs" gate.
2. **Post rule.** None of the fixed-parameter post variants helps. The missing piece is rho-LTD at single pairings and at train -10, together with timing-selective rho-LTP. Both need a change in what the potentiation gate reads: the own glutamate × bAP coincidence (Gb), without C1 converting non-coincident train Ca into depression. Run 2's GbL tests the gentlest form.
   - If GbL also fails, the next form to test is **a gate on potentiation only, with the failed-gate time made neutral** (C2 semantics for the Gb gate but C1 for the theta_V gate). That needs one more kernel flag but no parameter.
3. **Refit.** For any variant that passes Run 2:
   - add it as a kernel flag with no new parameter;
   - refit only a00..a11, gamma_d and gamma_p, seeded from C1Ajn_s5 plus one unseeded seed, on all 3 pathways (79G, 1:45 on 3g, from 22134473);
   - E2 alone is worth refitting first (pre side only, no L5 risk beyond the r50 bursts).
4. **Fit-impact estimate.**
   - E2 + refit: about -10 to -20 total.
   - A working coincidence gate: up to about -50 on L2/3. The AM251 pair (27.4), train +4 (25.0), train8 (13.7) and 1ap -10 (11.6) are the reachable part.
   - The train +4 AM251 target at 1.73 may stay out of reach of the readout (P_up check, section 5).

Not proposed:
- A pathway-specific CB1 (no eCB on L2/3 PC terminals) would fit the AM251 data trivially, but it is against the uniform rule.
- Raising K_ca or theta_Tg uniformly would remove the L5 single -10 eCB.

## 4. Model limits seen here

- APV (Z20/Z21) stays out: nmdar_block freezes rho. A test of the NMDAR-independence of train-LTD would need effcai without the NMDA Ca, which is not extracted.
- Expression locus: rho and dpre both act on Use, so the PPR split (LTP pre, LTD post) cannot discriminate between them.

## 5. Run 2 results (22135934 l23l23, 22135935 l5, 22135936 l23; done, 3:49 / 1:59 / 2:39, ~2.1 GB)

chi2 per variant at fixed C1Ajn_s5 parameters (base reproduces 110.27 / 46.89 / 5.08):

| variant | l23l23 | l5 | l23 | total |
|---|---|---|---|---|
| base | 110.27 | 46.89 | 5.08 | 162.2 |
| **E2** | **98.06** | **50.13** | **5.08** | **153.3** |
| E1 | 109.07 | 139.6 | 38.5 | 287 |
| C2 | 164.1 | 88.9 | 149.2 | 402 |
| GbL | 112.99 | 461.7 | 309.4 | 884 |
| GbL+E2 | 84.71 | 463.9 | 309.4 | 858 |
| Gb / Gbp / C2+Gbp (+E2) | 135–227 | 305–655 | 280–434 | > 700 |

- GbL+E2 is the best L2/3→L2/3 score of any rule (84.7), but it kills L5 LTP. E2 is the only variant that improves the total; V5_DESIGN §5 has adopted it as v5_mode 3.
- **The local quantities at single pairings are nearly identical in L2/3→L2/3 and L5.**
  - 1ap −10 (L2/3) vs Sjöström 0.1 Hz −10 (L5): pairing peak c* 0.024 mM in both, V at own arrival 0.36 θV(5.5) in both, 1 own VDCC event in both, and t_dep = 0 in both.
  - The only local difference is the c_pre-normalised pairing Ca: 0.52 c_pre (L2/3) vs 0.40 c_pre (L5). At 1ap +10 it is 0.61 vs 0.50.
- **10-AP trains invert any rule that is monotone in total spine Ca.**
  - train +5 first: c* peak 1.36 c_pre, the largest of all, and the data show no change (0.97).
  - train +4 last: c* peak 0.79 c_pre, and the data need LTP (1.49).
  - The local quantity that separates them is the order, i.e. how much VDCC Ca arrives before own glutamate. That is W (VDCC Ca weighted by 1 − b), not c*.
- **The V pool:**
  - L2/3 trains reach Vpk ≈ 10 θV(5.5) units, vs 1.4–1.9 in the L5 10 / 20 / 50 Hz protocols.
  - But the sj07 step pairing reaches Vpk 13, so a V-level LTD rule also hits sj07.

## 6. v5c (gpu_v5_rho v5_mode 2): which targets drive the L2/3→L2/3 misfit

**V5cz_s5** (3-pathway fit, `results/v5_V5cz_s5_l23l23.csv`, chi2 104.9) per target:

| target | data | pred | z |
|---|---|---|---|
| **train10 +5 first** | 0.97 | 1.355 | **+6.42** |
| **train10 +4 last** | 1.49 | 0.970 | **−4.34** |
| **train10 −4 last** | 0.99 | 0.674 | **−3.52** |
| **1ap −10** | 0.56 | 0.730 | **+2.83** |
| +4 last, AM251 | 1.73 | 1.113 | −2.57 |
| train8 +4 | 1.15 | 0.947 | −2.53 |
| 5ap 10 Hz +10 | 0.76 | 0.927 | +2.38 |
| −10 last, AM251 | 0.73 | 0.831 | +1.44 |
| others | | | \|z\| < 1 |

- The top three are all train-position errors. Together they make up 41 + 19 + 12 = 72 of the 105 chi2. LTP lands on the wrong pairing: +5 first gets LTP and +4 last does not. −4 also gets eCB LTD that should not be there.
- 1ap −10 and 5ap 10 Hz +10 are single or low-frequency pairings that need more LTD than the eCB ceiling, 1 + d_min = 0.71, can give. So they need ρ-LTD at a Ca level where L5 single pairings must not change.
- The −10 AM251 target needs postsynaptic LTD (data 0.73, model 0.83 with eCB blocked).
- **V5c_s7** (L5 + L2/3→L5 fit) has no L2/3→L2/3 csv. Jobs 22138665–7 score it on all three pathways.

**Run 3: the v5c replay plus 13 fixed-parameter variants, 6 jobs**
- Jobs: V5c_s7 → 22138665 (l23l23), 22138666 (l5), 22138667 (l23). V5cz_s5 → 22138668 (l23l23), 22138669 (l5), 22138670 (l23).
- Script: `diag_l23l23_v5.py` + `run_diag_l23l23_v5.sh` (1 CPU, 2G, 0:15).
- Output: `/scratch/dhuruva/l23l23_diag_v5/<fit>_<model>_{syn,rec,targets,proto}.csv`. Logs: `logs/diag_l23v5_<id>.out`.
- Every variant is uniform and synapse-local, reads only own Ca quantities, and uses the fitted parameters (no refit):

| variant | change | what it tests |
|---|---|---|
| E2 | v5_mode 3: W fed only by excursions that start with the pool unprimed | the adopted fix |
| noE | eCB off | post-only bound |
| DW / DV | Chindemi depression also while W (or V) > θ_eCB | **post, AM251-insensitive LTD from the eCB trigger state** |
| PW / PWo | ρ → 0 at an own arrival with W > θ_eCB, with / without the eCB step | the same, as a one-shot post LTD |
| GW / GWE | LTP gate reads W > θ_V (or > θ_eCB) instead of V > θ_V | **order-sensitive LTP: fixes +5 first vs +4 last?** |
| td.8 / .6 / .5 | θ_d × 0.8 / 0.6 / 0.5 (an a00/a01 move, scored unrefit) | single-pairing ρ-LTD from the 0.52 vs 0.40 c_pre margin |
| E2+DW | combination | |

**What to read from Run 3**
1. **Repro:** the V5cz_s5 logs must print `repro max |base - fit_pred|` < 1e-6 for all three models. The V5c_s7 l5 / l23 logs must do the same against `v5_V5c_s7{,_l23}.csv`.
2. **Scoring:** sum the `chi2 per variant` lines over the three models. A candidate counts only if l5 + l23 stays within about 10 of base.
3. **Reachability (`*_targets.csv`):**
   - The `reach` column flags targets outside [all_min (post_min under AM251), ltp_max].
   - `ltp_max` at +4 AM251 (1.73) and at +4 last (1.49) shows whether those values are readout-reachable at all.
4. **Local state (`*_proto.csv`), L2/3 1ap ±10, 5ap, train10 vs L5 sjostrom −10 and Letzkus 1ap:** check W_arr, VE_arr, trig (fraction of arrivals that fire eCB), Bg_arr, Vpk, Wpk, t_dep and pk_tp.
5. **Absolute spine Ca (µM above rest), per protocol:** check bAP_ref, EPSP_ref, bAP_EPSP and supra at 1ap +10 / −10, compared with the literature in section 7.

## 7. Is Zilberter reachable on this L2/3 emodel? (argument; numbers are confirmed by Run 3 section 5)

**Literature (according to PubMed metadata; local full text of Nevian & Sakmann 2006 in `nevian/`)**
- **L2/3 basal spines** (Nevian & Sakmann 2006, J Neurosci 26:11001, [doi:10.1523/JNEUROSCI.1749-06.2006](https://doi.org/10.1523/JNEUROSCI.1749-06.2006)):
  - Peak spine Ca rises linearly with the number of APs in the burst (r² > 0.98), for both orders.
  - EPSP → 3 APs at +10 is supralinear (1.8 ± 0.1). AP-before-EPSP sums linearly.
  - LTP and LTD protocols reach the same peak Ca (Fig 8). The direction is set by an mGluR → PLC → eCB sequence detector, "VDCC before mGluR".
- **L5 basal spines** (Koester & Sakmann 1998, PNAS 95:9596, [doi:10.1073/pnas.95.16.9596](https://doi.org/10.1073/pnas.95.16.9596)):
  - The spine Ca from a single AP is comparable to that from a single EPSP, up to 80 µm.
  - EPSP → AP is supralinear and AP → EPSP is sublinear.
- **Koester & Sakmann 2000** (J Physiol 529:625, [doi:10.1111/j.1469-7793.2000.00625.x](https://doi.org/10.1111/j.1469-7793.2000.00625.x)) is about **axonal boutons** of L2/3 PCs (about 500 nM residual free Ca per AP), not spines. It anchors only presynaptic Ca, and the model has no presynaptic Ca variable.
- In both cell types the bAP Ca is therefore of the same order as the unitary EPSP Ca in the active spine.

**The model**
- c_post / c_pre is 0.011 (L2/3) and 0.005 (L5) in effcai units. The bAP spine Ca is about 1 % of the EPSP Ca. An earlier audit (EMODEL_CONSTRAINTS_FROM_PLASTICITY.md) found the L5 pooled spine cai_CR per bAP at about 0.1 µM, against 1.1 µM in the literature.
- The L2/3 hocs share every dendritic density with L5 (L23_PATHWAYS.md). Their bAP Ca has not been audited; Run 3 bAP_ref / bAP_EPSP measures it.
- Consequence: adding 9 APs barely raises c*. The train10 −10 peak is 0.66 c_pre vs 0.52 at 1ap −10, where Nevian finds a linear rise with AP number. So **no rule that reads c* can tell a 10-AP train from a single AP**. Only the VDCC pools V and W carry the AP count and the order.

**Verdict (pending Run 3)**
- (i) **Train-position targets** (+5 first, +4 last, −4 last; about 72 chi2): not reachable by any c*-only rule, because c* inverts their order. They are reachable in principle with an order-sensitive VDCC quantity (W) on the LTP side, which is what GW / GWE test. If GW also breaks L5 LTP (Markram / Sjöström +10 have little VDCC Ca before own glutamate), the data are not reachable with the present bAP:EPSP spine-Ca balance.
- (ii) **Single-pairing LTD** (1ap ±10 at 0.56 / 0.64, below the 0.71 eCB ceiling): needs ρ-LTD at the same local state where L5 must not change (§5). The only lever is the 0.52 vs 0.40 c_pre margin (td scan). It is narrow and would be set by per-synapse spread, so it is not robust.
- (iii) **Post LTD in train −10 under AM251**: DW / DV give it from the VDCC state. Run 3 shows whether sj07 (Vpk 13) survives.

**Emodel constraint if (i) / (ii) fail**
- The L2/3 basal spine bAP Ca per AP must be of the same order as the unitary EPSP Ca (bAP:EPSP about 0.5–1, as in Koester & Sakmann 1998 and Nevian & Sakmann 2006).
- Spine Ca must grow about linearly with AP number in a burst.
- EPSP → AP at +10 must be supralinear (about 1.8) and AP → EPSP linear.
- In this emodel that means a larger bAP-driven spine-head Ca entry at basal L2/3 synapses. It would come as a new, renamed channel or emodel (new-names rule), calibrated on those Ca observables, not on plasticity.
- Only then can a uniform c*-threshold rule separate 1 AP from 10 APs. The train-position effects would still need the order-sensitive VDCC term.
