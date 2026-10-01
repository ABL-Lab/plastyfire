# GluSynapse NMDA / VDCC spine-Ca calibration vs Chindemi et al. 2022 (2026-09-30)

**Question:** why is synaptic Ca 92–98% NMDA, and why is spine Ca per bAP only 0.08 µM (pooled median) when Sabatini 2002 reports ~1.1–1.7 µM?

**Short answer:** mostly (c). The mod matches the paper exactly. The comparison mixed pools: Chindemi calibrated the bAP against basal spines < 60 µm from the soma only. In that pool our cell gives 0.86 µM against Chindemi's 1.4 ± 0.6 µM. Explanation (b) accounts for a small remaining factor of 0.6×.

Figure: `figs/nmda_vdcc_calibration.png`. Scripts and data: `nmda_vdcc/` (summarize.py, clamp_vdcc.py, plot.py, run.sh, rows_all.json). The sims are diag_burst runs `results/diag_burst_*_{sp1,cal1}.json` (8 L5→L5 pairs, 68 synapses, job 22066817/21, ~0.6 CPU·h).

## 1. Calibration table

| Parameter | Chindemi 2022 (Methods) | `SSCx_og_mods/GluSynapse.mod` | Runtime (prefire config / edges) |
|---|---|---|---|
| NMDA Ca fraction | s = α·Pf(V→−∞), α = 0.6, Pf from Schneggenburger; E = 40 mV (eqs 21–23) | `Pf = 4cao/(4cao + 120/1.38)·0.6` = 0.0505 at 2 mM; E = 40 mV | cao_CR = 2.0 (config `conditions.mechanisms`) |
| VDCC type / gating | R-type, m²h; Vhm −5.9, km 9.5; Vhh −39, kh −9.2; τm 1, τh 27 ms (eqs 26–32) | identical; ljp_VDCC = 0 | not overridden |
| VDCC density | 0.0744 nS/µm² on a spherical head (eq 28) | gca_bar_VDCC = 0.0744 | not overridden |
| Spine volume | lognormal μ −2.8, σ 0.87 (0.087 ± 0.088 µm³), scaled with synaptic conductance | default 0.087 | edges `volume_CR`: for these 68 L5→L5 synapses, median 0.153 (IQR 0.105–0.211) µm³ |
| η (free fraction) | 0.04 | gamma_ca_CR = 0.04 | — |
| τCa, [Ca]rest | 12 ms, 70 nM | 12 ms, 70e-6 mM | — |
| [Ca]o | 2 mM in vitro | cao_CR = 2.0 | 2.0 |
| **Targets** | synaptic 0.67 ± 0.44 µM (Sabatini 0.7 ± 0.4); **bAP 1.4 ± 0.6 µM (Sabatini 1.7 ± 0.6), basal only, path < 60 µm**, "beyond 100 µm the bAP is insufficient to evoke VDCC Ca" | | |

The plastyfire and v2 mods differ from SSCx_og only in `calcium_current_flag` and the v2 presynaptic terms; the Ca equations are the same. Nothing in pairrunner sets VDCC or Ca globals, except `--glusyn-globals` when it is passed.

## 2. Fig 2-type validation in our cell (og-delta; 8 pairs, 68 synapses: 50 basal, 18 apical)

Spine Ca is `cai_CR` peak minus rest. VDCC share is the VDCC fraction of spine Ca charge (∝ integrated Ca, which is what effcai sees).

| Case | All 68: median [IQR] µM | Basal < 60 µm (n = 8) | Basal 60–150 (n = 33) | VDCC share (all / < 60) | Reference |
|---|---|---|---|---|---|
| EPSP only | 0.22 [0.0002–0.56]; mean 0.33 (63% release); **0.52 ± 0.34 given release** | 0.31 ± 0.27 | 0.13 | 0.00 / 0.00 | 0.67 ± 0.44 (Sabatini 0.7 ± 0.4) |
| bAP only | 0.08 [0.007–0.8] | **0.86 ± 0.28** | 0.07 | 1 / 1 | **1.4 ± 0.6** (Sabatini 1.7 ± 0.6) |
| pre→post +10 ms | 1.17 [0.26–1.9] | 1.88 ± 0.72 | 0.66 | 0.13 / 0.19 | Fig 1e: supralinear |
| post→pre −10 ms | 0.92 [0.25–1.3] | 1.05 ± 0.22 | 0.66 | 0.16 / 0.24 | Fig 1e: ≈ linear sum |

Peak local v during the bAP: +0.1 mV at basal < 60 µm, −27 mV at 60–150 µm, −44 mV beyond 150 µm, −17 mV apical (Fig B). Synaptic (NMDA) Ca given release is within 1 SD of Chindemi and Sabatini; the lower unconditional mean reflects failures (init_depleted, one spike). The bAP pairing ordering (+10 > −10) is reproduced.

## 3. Diagnosis

- **(a) Mod mis-parameterised: no.** Every Ca parameter equals the paper's value (table above), and there is no runtime override. The larger L5→L5 spine volume (0.153 vs 0.087 µm³) lowers VDCC Ca by only (0.153/0.087)^(−1/3) = 0.83×. It is also the paper's own conductance-to-volume prescription.
- **Mod-alone clamp (Fig C):** a single compartment clamped to a bAP waveform, using the compiled mod. Ca per bAP is set by peak V and width alone. It rises e-fold per ~6–7 mV: at half-width 2.1 ms it is 0.02 / 0.10 / 0.38 / 1.04 / 2.0 µM at −30 / −20 / −10 / 0 / +10 mV. The cell's synapses lie on this curve, so the mod works as designed. Chindemi's 1.4 µM needs about +4 mV peak, or about 2.5 ms half-width at 0 mV, at the spine.
- **(b) Cell bAP weaker than Chindemi's: only a small factor.** At < 60 µm our bAP reaches ~0 mV, which gives 0.86 µM, or 0.6× Chindemi (0.9 SD below). Chindemi's 2015 L5TTPC emodel (cADpyr232) is not on disk; no `*232*` exists in DEES_cell_packages or plastyfire. The closest available relative is the SSCx og emodel with passive basal dendrites (no NaTg, Ka, Kv3 or BK). Applied to the same synapses, it gives an even weaker proximal bAP: 0.78 ± 0.31 µM at < 60 µm and 0.05 µM pooled. So og-delta's basal Na did not cause a regression. It already gives slightly more.
- **(c) Pool mismatch plus m² steepness: the main cause.** 60 of 68 synapses lie beyond 60 µm, where the bAP peaks at −27 to −44 mV, and there R-type m²(V) ≈ 0.005–0.05. This is exactly the attenuation Chindemi describe ("beyond ~100 µm insufficient to evoke VDCC Ca"). The 0.08 µM pooled median and the small VDCC share are therefore properties of Chindemi's validated model at these synapse locations. They are not a calibration error. The effcai integration is fine: NMDA Ca lasts ~5× longer than the ~12 ms VDCC transient, so even at < 60 µm a −10 ms pairing is only 24% VDCC by charge. The "2–8%" figure was burst-protocol effcai from syn1; for single pairings the pooled share is 13–16%.

## 4. Recommended fix

Keep the mod and the cell. If the aim is to match Chindemi's own calibration target exactly, add one runtime GLOBAL: **`ljp_VDCC_GluSynapse = 5` (mV)**. This shifts R-type activation and inactivation 5 mV negative. Tested in diag_burst (cal1):

| Variant | bAP Ca, basal < 60 µm | bAP Ca, pooled median | VDCC share +10 / −10 (pooled) |
|---|---|---|---|
| og (mod as is) | 0.86 ± 0.28 | 0.08 | 0.13 / 0.16 |
| **s5 (ljp +5)** | **1.36 ± 0.35** (Chindemi 1.4 ± 0.6) | 0.20 | **0.21 / 0.25** |
| g2 (gbar ×2) | 1.63 ± 0.53 (Sabatini 1.7 ± 0.6) | 0.16 | 0.23 / 0.27 |
| s10 | 2.0 ± 0.43 | 0.45 | 0.31 / 0.36 |
| s10g3 (delta-sv) | 5.5 ± 1.2 (3× too high) | 1.26 | 0.55 / 0.61 |

- **Why s5 over g2:** the Magee & Johnston 1995 gating comes from cell-attached recordings with a high-Ba²⁺ pipette and no LJP correction. Both effects shift apparent activation positive by several mV, which is the reason `ljp_VDCC` exists in the mod. +5 mV is the smallest shift that restores the Chindemi < 60 µm number, and it leaves NMDA Ca unchanged (EPSP share 0). g2 is the density-based alternative, if a denser-than-20/µm² channel count is preferred; it gives similar shares.
- **What it does not do:** the pooled per-bAP median stays ~0.2 µM, and the VDCC share stays ~20–25%. This is inherent to R-type m² kinetics at the attenuated bAP of 60–200 µm synapses. Pushing the pooled median to 1.1 µM (s10g3) overshoots the proximal calibration 3×. It would be fitting to plasticity data, not to Ca data.
- **Alternative: accept the current state.** It is within 1 SD of Chindemi at the calibrated locations. Changing the cell's basal bAP was already shown not to help (BASAL_NAKV.md), and it costs spike fidelity.
