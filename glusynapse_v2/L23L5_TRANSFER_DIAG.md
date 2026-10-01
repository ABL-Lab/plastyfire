# L2/3→L5 transfer failure: diagnosis (2026-09-30, no refit, no model edits)

**Setup.** Fit `reduced_gpu_subset_pl5r50_td2_s2` (and td2_s1, which gives the same picture within 0.05). Records: `ebner_l23l5_delta-prefire` (111 `all_protocols` pairs) and `ebner_delta-prefire` (all 111 L5→L5 pairs plus the 24-pair fit subset). Protocol `sjostrom_50hz_dt+10ms` is identical in both yamls except the pulse width (5 vs 3 ms). My per-record ratios reproduce `l23l5_pl5r50_td2_s2_rise_records.csv` exactly (max diff 0) and the fit csv's 24-pair values (1.582, 0.966). With dpre = 0 the readout is exact: `R_post = (1 + Σ_i Δb_i·g_i)(1 + CV0²)`, where Δb_i = binary flip and g_i = δ_i / EPSP0. That lets each flip be attributed to its synapse. Scripts: `diag_l23l5/diag.py`, `agg2.py`, `diag2.py` (login node, about 2 min, no jax). Job 22074079: 15.2 GB MaxRSS of 20 GB, 1:30, 76% CPU.

## 1. The "L2/3 gets more LTP than L5" comparison is a subset artefact
| sjostrom_50hz_dt+10ms | model R | post-rule R_post | pre part | data |
|---|---|---|---|---|
| L5→L5, 24 fit pairs | 1.582 | 1.486 | +0.097 | 1.57 ± 0.26 |
| L5→L5, all 111 pairs | **1.847** | 1.631 | +0.216 | — |
| L2/3→L5, 111 pairs | **1.845** | 2.059 | −0.215 | 1.06 ± 0.09 |

At population level the model is **pathway-blind**: 1.847 vs 1.845. The rho0 mean is similar too (0.53 vs 0.51). What differs is the composition. In L2/3 the LTP is all postsynaptic, and the pre rule is already pushing it down. In L5 the pre rule supplies about a quarter of the LTP.

## 2. Where the L2/3 post-rule LTP comes from: synapses with no bAP Ca
c_post is bimodal. At the resting floor (below 1e-3 mM, mode about 1e-4) sit **72% of L2/3→L5 synapses vs 52% of L5→L5**. This matches LETZKUS_CA: the bAP is about 2 mV beyond 450 µm. At these synapses θ = a·c_pre + b·c_post collapses to **θ_p ≈ 1.10·c_pre and θ_d ≈ 0.80·c_pre**. The pairing peak is 2.2·c_pre (50 Hz) or 1.27·c_pre (1 AP), so the synapse's own NMDA Ca crosses θ_p without any bAP. The resulting LTP does not depend on post timing.

| net EPSP gain Σ Δb·g per pair | floor synapses | bAP-reached | total |
|---|---|---|---|
| L2/3→L5, 50 Hz | **0.600** (80% of floor rho0 = 0 flip up) | 0.362 | 0.961 |
| L5→L5 all, 50 Hz | 0.133 (59% flip up, 33% of the rho0 = 1 synapses flip down) | 0.477 | 0.609 |
| L2/3→L5, Letzkus 1AP +10 | **0.320** | 0.081 | 0.401 |

## 3. Decomposition of the L2/3 50 Hz miss (model 1.845 vs data 1.06, Δ = 0.785)
- **(a) Synapse population × θ scaling: −0.63 (80%).** Freezing the post rule at floor synapses gives 1.213. S&H distal drops from 1.79 to **1.00** (data 0.86 ± 0.09) and proximal drops from 1.99 to 1.74 (derived data about 1.21).
- **(c) EPSP mapping: −0.11 (14%).** The gain per flip is 0.289 in L2/3 vs 0.239 in L5 (×1.21; all-potentiated ceiling 2.18 vs 1.88; with Use_p/Use_d 10.6 vs 8.4 and unitary EPSP 0.39 vs 2.0 mV). The CV correction is 1.049 vs 1.013. This is also why the nmdar_block rows are 1.03. With both (a) and (c) removed the prediction is 1.10.
- **(b) Protocol encoding: about 0.** For Letzkus 1AP, truncating to 15/50/100 pairings gives 1.43/1.51/1.48. Emulating 0.1 Hz by adding 9 s of rho relaxation per pairing changes it by ≤ 0.01. The rule saturates within about 15 pairings, and effcai (τ 278 ms) does not sum at 1 Hz.
- **(d) The rule proper: the pre rule is not the cause** (dpre −0.07 in L2/3 vs +0.24 in L5). The degenerate floor thresholds are part of the rule's θ normalisation. They matter only because of the population in (a).
- **Residual: +0.04,** within the SEM.

**Letzkus 1AP +10 (1.48 vs 0.72):** gating the floor synapses gives 1.14 and adding (c) gives 1.08. The remaining −0.36 is a missing +10 ms single-AP LTD mechanism. That is a cell and rule limit, already predicted in L23_PATHWAYS §4, and no gate fixes it.

## Verdict: **(a)**
The cause is the L2/3→L5 synapse population: 72% sit at the c_post floor, where the c_pre/c_post threshold scaling makes their post rule respond to pre-only Ca. The EPSP mapping (c) adds a second-order 14%. Protocol (b) and the pre rule are cleared.

## Smallest fix within the rules
**Post-rule bAP gate.** For synapses with c_post at the resting floor (< 1e-3 mM, the gap in the bimodal c_post), set θ_d = θ_p = ∞, so rho is frozen and the pre rule is kept. This adds no new parameters. If the cutoff is freed it is one parameter, fitted on L5 paired data only.
- **Offline estimate before refit:**
  - L2/3 all: 1.21 (z 8.7 → 1.7).
  - L2/3 distal: 1.00 (z 10.3 → 1.6).
  - L5 24-pair set: 50 Hz 1.58 → 1.50, 0.1 Hz +10 0.97 → 1.03.
- **Cutoff sensitivity:** 3e-4 to 3e-3 gives 1.29 to 1.16 on L2/3 50 Hz.
- **Next steps:** add a gate flag to thetas (batch_v2/jax_v2, owner's call), refit paired_l5 on the 24 pairs, then re-run `run_eval_l23l5.sh` unchanged. Nothing is fitted on L2/3 data. The mapping part (c) stays as it is, because it follows from the circuit's Use/EPSP values.
- **Caveat:** the floor fraction is inflated by the cell's too-steep bAP attenuation (LETZKUS_CA). The gate encodes "no bAP Ca, no Hebbian post signal", which agrees with S&H 2006, where distal inputs get no LTP.
