# Spine VDCC Ba->Ca correction (ljp_VDCC 25 mV) vs the Chindemi 2022 / Sabatini 2002 validation (2026-09-30)

Setup: same 8 pairs, 68 synapses, protocols and diag_burst as NMDA_VDCC_CALIBRATION.md s2. Job 22090528 (8 x 1 CPU, 6 min,
<= 0.99 GB); summary 22090852 (`summary_22090852.txt`). Spine Ca = cai_CR peak - rest (uM). The recommended variant is
**ljp25_g031: ljp_VDCC 25, gca_bar_VDCC 0.02306 (0.31 x 0.0744)**.

| Check (Chindemi Fig 2 / Sabatini) | Target | current delta | ljp25_g031 | holds? |
|---|---|---|---|---|
| bAP spine Ca, basal < 60 um | 1.4 +- 0.6 (Sabatini 1.7 +- 0.6) | 0.86 +- 0.28 | **1.50 +- 0.32** | yes (better) |
| EPSP spine Ca, all (incl. failures) | 0.67 +- 0.44 (given release) | mean 0.33 | mean 0.37 | unchanged; VDCC share of EPSP Ca 0 -> 6% |
| pre->post +10 supralinear vs EPSP + bAP (basal < 60) | supralinear | 1.79 vs 0.23 + 0.83 | 2.41 vs 0.24 + 1.43 | yes |
| post->pre -10 ~ linear sum (basal < 60) | ~ linear | 1.09 vs 1.06 | 1.57 vs 1.67 | yes |
| +10 > -10 ordering | yes | 1.79 > 1.09 | 2.41 > 1.57 | yes |
| "beyond ~100 um the bAP is insufficient" (Chindemi text) | model statement | 60-150 um: 0.07; >= 150: 0.004 | 60-150 um: **0.99**; >= 150: 0.14 | **no** |

- The one deviation is Chindemi's qualitative statement. It is a property of their model, which used the same uncorrected
  Ba gating, not a measurement. Data point the other way: AP-evoked spine Ca is measurable at 50-150 um (Nevian & Sakmann
  2006, L2/3 basal) and is present in L5 basal dendrites (Kampa & Stuart 2006).
- Consequence: the VDCC share of pairing Ca (effcai input) rises from 13-16% to 34-40%. effcai changes everywhere, so the
  post rule must be refit on a new prefire. g025 / g040 bracket the proximal target (1.23 / 1.90 uM).
