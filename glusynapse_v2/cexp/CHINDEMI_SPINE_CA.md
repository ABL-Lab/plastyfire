# Spine Ca and spine VDCC in Chindemi et al. 2022 vs our model (A15, 2026-10-01)

Paper: Chindemi et al. 2022 Nat Commun ([doi](https://doi.org/10.1038/s41467-022-30214-w), PMC9160074), full text via PubMed. Section names are from Methods; the PMC text has no equation bodies, so formulas come from the mod. Mods: Chindemi `DEES_cell_packages/SSCx_og_mods/GluSynapse.mod` (og), `plastyfire/plastyfire/GluSynapse.mod`, ours `glusynapse_v2/mod/GluSynapseV5.mod` / `GluSynapseV7.mod` (V7 line numbers below). Supplementary PDF not read (not in PMC text). Earlier audits reused: `cell_audit/NMDA_VDCC_CALIBRATION.md`, `cell_audit/ljp25_validation/VALIDATION.md`, `cell_audit/CA_VS_DATA.md`, `cell_audit/VDCC_DECODE.md`.

## (a) Chindemi's equations (point process, no neck, no diffusion to dendrite; v = dendritic segment v)
- Ca ODE (Methods "Postsynaptic calcium dynamics"; og 325-326, V7 604-605): `dcai_CR/dt = -(ica_NMDA + ica_VDCC)·gamma_ca_CR / (2·F·volume_CR) - (cai_CR - min_ca_CR)/tau_ca_CR`. gamma_ca_CR 0.04, tau_ca_CR 12 ms, min_ca_CR 70 nM.
- NMDA Ca (Methods "NMDAR"; og 294-295, V7 526-527): a fixed fraction of g_NMDA, not GHK. `ica_NMDA = Pf·g_NMDA·(v - 40 mV)`, `Pf = 4cao/(4cao + 120/1.38)·0.6` = 0.0505 at cao_CR 2 mM. Erev 40 mV and Pf(V→-inf) are from Schneggenburger 1993 ([doi](https://doi.org/10.1016/0896-6273(93)90277-x)). Mg block `1/(1+exp(-0.072 v)·mg/2.552)`, refit to Vargas-Caballero & Robinson 2003 ([doi](https://doi.org/10.1152/jn.01038.2002)).
- VDCC (Methods "VDCC"; og 297-300, 320-323; V7 529-532, 574-577): R-type, `ica_VDCC = 1e-3·gbar_abs·m²h·(v - E_Ca)`, with E_Ca Nernst from cai_CR and cao_CR. `m_inf = 1/(1+exp((vhm-ljp-v)/km))`, vhm -5.9 mV, km 9.5; `h_inf` vhh -39 mV, kh -9.2; tau_m 1 ms, tau_h 27 ms; ljp_VDCC 0. `gbar_abs = gca_bar_VDCC·4π·(3·volume_CR/4π)^(2/3)`, sphere surface, gca_bar_VDCC 0.0744 nS/µm². At 0.087 µm³ that is 0.071 nS, about 19 channels.
- Spine volume (Methods "Dendritic spines", "Correlation of synaptic parameters"): per synapse (RANGE volume_CR in the edges). L5 basal reconstructions give 0.087 ± 0.088 µm³, lognormal μ -2.8, σ 0.87. The volume is mapped linearly from synaptic conductance (ratio of means, L5 basal, applied to all pathways), so gbar_VDCC ∝ V^(2/3) and is tied to g_NMDA.
- In og and ours, gca_bar_VDCC, ljp_VDCC, gamma, tau and rest are GLOBAL (one value per run). Only volume_CR is per synapse.

## (b) How each was set: none were fitted to plasticity; only the 8 thresholds, tau_effca and gamma_d/p were fitted
| quantity | source / method |
|---|---|
| Pf, Erev 40 mV | GHK fractional current from Schneggenburger 1993. alpha 0.6 tuned so the fixed-fraction current matches GHK (self-consistency, Supp Fig). |
| Mg block 0.072 / 2.552 | refit by inspection to the neocortical steady-state unblock of Vargas-Caballero & Robinson 2003 (LJP-corrected) |
| VDCC gating | Magee & Johnston 1995 cell-attached dendritic patches in 110 mM BaCl2 ([doi](https://doi.org/10.1113/jphysiol.1995.sp020862)); tau_h corrected to 34 °C; tau_m chosen to fit the rise. No Ba→Ca correction (ljp 0). |
| VDCC density | 20 channels/µm² (Sabatini & Svoboda 2000, [doi](https://doi.org/10.1038/35046076)) × 3.72 pS unitary (Bartol 2015, mod comment) = 0.0744 nS/µm² |
| gamma 0.04, tau 12 ms, rest 70 nM | Sabatini, Oertner & Svoboda 2002 ([doi](https://doi.org/10.1016/s0896-6273(02)00573-1)): kappa_E ≈ 24, tau 14 (12-20) ms, rest 70 ± 29 nM |
| volume | Benavides-Piccione P14 L5 basal spines (8423 protrusions); correlation with g from Harris & Stevens and Arellano |
| **validation only** (Results, Supp Fig 2) | vs Sabatini 2002 (CA1). Synaptic: 0.67 ± 0.44 µM model vs 0.7 ± 0.4. bAP: 1.4 ± 0.6 vs 1.7 ± 0.6, **L5 basal, path < 60 µm only**, because "beyond 100 µm [bAPs are] insufficient to evoke VDCC calcium". No APV-residual, Koester, Nevian or Kampa target was used. The text says "isolated presynaptic spikes primarily trigger NMDAR" Ca and "isolated postsynaptic spikes only activate VDCCs". |

## (c) Chindemi vs our current model (delta-split1 prefire: `prefire_simulation_config.json` conditions.mechanisms sets only cao_CR, tau_effca, gamma_d/p, init_depleted)
| item | Chindemi (paper / og mod) | ours (V5/V7 mod + runtime) | diff |
|---|---|---|---|
| Ca ODE, gamma, tau, rest | as (a) | identical (V7 234-236, 604-605) | none |
| NMDA Ca: Pf·0.6, Erev 40, Mg 0.072/2.552, E_NMDA -3, tau 0.29/70 ms | as (a) | identical (V7 208-213, 522-527) | none |
| VDCC gating, ljp_VDCC | M&J, ljp 0 | identical, ljp 0 (not overridden) | none |
| gca_bar_VDCC (GLOBAL) | 0.0744 | 0.0744 (not overridden, cannot be 0 per synapse) | none |
| VDCC written into cai_CR, and into i | yes (og 300, 302) | yes (V7 532, 536, 604) | none |
| calcium_current_flag | absent | RANGE, default 1, never set in any .py | none in effect |
| volume_CR | 0.087 mean, from g | edges per synapse: L5→L5 median 0.153 (IQR 0.105-0.211), 68 synapses | **×1.76 → VDCC [Ca] ×0.83, NMDA [Ca] about unchanged (V ∝ g)** |
| emodel | cADpyr232 (2015, not on disk) | og-delta / delta-split1 | **bAP at < 60 µm peaks ~0 mV here (0.86 µM vs 1.4)** |
| c_pre definition | all RRP vesicles at once; effcai units (tau_effca fitted) | same (cexp "cpre" = full RRP, effcai tau 278 ms) | same, but **not the same quantity as a single-release EPSP** |

## (d) Why ~0 % VDCC share and a low bAP/EPSP ratio (hypothesis, with evidence)
1. **This is Chindemi's model working as designed. It is not a parameter error.** Every Ca parameter and the runtime globals equal the paper's. gca_bar is not 0, ljp is not missing from the code (it is 0 by design), and VDCC Ca does enter cai_CR.
2. **No spine-head voltage.** The VDCC sees dendritic v. A unitary EPSP depolarises the dendrite by a few mV, and with m² gating m_inf² is 4e-6 at -65 mV and 3e-5 at -55 mV, against 0.42 at 0 mV. So EPSP VDCC Ca is about 1e-3 of NMDA (CEXP: vdcc_q/nmda_q ~1e-3), and APV leaves the 70 nM floor. To reach a 10-30 % APV-resistant share, the spine head has to depolarise (neck resistance; Bloodgood 2009 [doi](https://doi.org/10.1371/journal.pbio.1000190), Grunditz 2008 [doi](https://doi.org/10.1523/JNEUROSCI.2702-08.2008)). That is outside Chindemi's formulation.
3. **The bAP side is the Ba-shifted gating plus location.** M&J gating is from 110 mM Ba²⁺ with no surface-charge correction. That shifts activation positive by about 20-30 mV (`cell_audit/VDCC_DECODE.md`). With ljp 0, a bAP that peaks at -27 to -44 mV (60-150 µm and beyond) gives m² of 0.01 or less. Chindemi calibrated only at < 60 µm. Most of our basal synapses lie beyond 60 µm (60 of 68 in the audit), so the pooled bAP Ca is 0.08 µM.
4. **The ratio is partly a definition issue.** The CEXP cpost/cpre (0.033 L5 basal, 0.14 L2/3→L5 basal) divides by a full-RRP EPSP in effcai units. tau_effca 278 ms weights the 70 ms NMDA transient far more than the 12 ms VDCC one. Koester compares peak ΔF for one AP against one EPSP. In free Ca, L5 basal is 0.21/1.6 = 0.13, and against a single release (0.52 µM given release) the pooled ratio is about 0.15, about 2.8 at < 60 µm and about 0.5 at 60-150 µm.
5. The volume being 1.76× the default contributes only ×0.83 to VDCC [Ca], so it is minor.

## (e) Changes for the implementing agent (synapse level only; emodel untouched)
**E1 (no mod change, only Chindemi's own GLOBALs): adopt `spine/delta_ljp25.json`, i.e. `ljp_VDCC_GluSynapse 25`, `gca_bar_VDCC_GluSynapse 0.02306` (0.31 × 0.0744).**
- How to apply it: add both to `conditions.mechanisms.GluSynapse` in the prefire configs, or pass them through pairrunner `--glusyn-globals`. For V5/V7 runs, use the `_GluSynapseV5/_V7` suffixes. They are GLOBAL, so they are uniform across pathways.
- Anchoring: the Ba→Ca shift (Smith 1993 [doi](https://doi.org/10.1085/jgp.101.5.767), Zamponi & Snutch 1996 [doi](https://doi.org/10.1007/BF02207290); VDCC_DECODE.md) sets ljp. gca is then recalibrated to Chindemi's own target, Sabatini's bAP Ca at basal < 60 µm. Neither is fitted to plasticity data.
- Validated in VALIDATION.md (22090528): bAP 1.50 ± 0.32 µM at < 60 µm (target 1.4-1.7), **0.99 µM at 60-150 µm** (delta 0.07; Koester & Sakmann 1998 [doi](https://doi.org/10.1073/pnas.95.16.9596) "comparable to EPSP up to 80 µm"), EPSP VDCC share 0 → 6 %, +10 supralinear, +10 > -10 kept.
- Known costs:
  - spine/shaft Ca ratio too high (CA_VS_DATA.md);
  - L5-only refit was not better (43.4 vs 39.5, DECISIONS 2026-09-30 evening, under td4, not under v7);
  - every VDCC-referenced quantity (c_VDCC pool, theta_V/theta_eCB, uE_GB = 1e5·vdcc_q_post, K_ca_GB) shifts, so it needs a new prefire, a new cexp and a refit.
- Recheck after the change: rerun `cexp/run_measure_cexp.sh` with the globals. Expected: cpost/cpre basal up about 10× (free Ca); NMDA share stays about 0.94-0.99.
- Variant: `ljp_VDCC 5` (gca unchanged) matches Chindemi's < 60 µm number only (1.36 µM). It leaves 60-150 µm low, so it does not reach Koester.

**E2 (only if the 10-30 % APV-resistant EPSP share is required; beyond Chindemi; ask the user first): new mod `GluSynapseV8.mod`, copied from V7 with SUFFIX GluSynapseV8.**
- Add a spine-head voltage for the VDCC and the Mg block only, `v_sp = v - R_neck·i_syn` (i_syn = i_AMPA + i_NMDA + ica_VDCC, inward negative). The membrane current stays at the segment.
- R_neck is a GLOBAL anchored from the literature (~500 MΩ, Harnett 2012 [doi](https://doi.org/10.1038/nature11554), CA1 apical; no L5 basal value checked), not fitted. Default 0 means V7 exactly.
- It is synapse-local and adds no emodel channel. It raises both the EPSP VDCC share and the NMDA nonlinearity, so recheck Mg0/cpre and Sabatini 0.7 µM synaptic Ca.
- Never edit V5/V7 in place.
