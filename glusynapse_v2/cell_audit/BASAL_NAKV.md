# Basal Na / Kv audit, og-delta L5TTPC (2026-09-29)

Question: are basal Na and Kv too low for the fitted L5→L5 protocols to leave distinct synaptic signals? Code and data are in `cell_audit/basal_nakv/`: `grid.py` (diag_burst with pre trains), `analyze.py`, `plot.py`, `contrasts.csv`, `same_syn.csv`, and `out/*.json`. Figure: `figs/basal_nakv.png`. Grid job 22065638: 78 tasks, about 6.5 CPU·h.

## 1. Channels (S/cm², from `SSCx-AAD-delta-emodels/cADpyr_L5TPC.hoc`)
| | soma | basal (og-delta) | apical | antic-delta basal |
|---|---|---|---|---|
| NaTg | 0.290 | **0.003** | 0.0229 | 0.0143 |
| SKv3_1 (Kv3) | 0.388 | **0.00182** | 0.00182 | 0.00182 |
| K_Pst / K_Tst | 0.100 / 0.013 | not inserted | not inserted | — |
| Ka_kampa | — | 0.002 | 0.0228 | 0.025 |
| kBK | — | 0.0427 | 0.00356 | 0.0427 |
| Ca_HVA2 / Ca_LVAst | 8.4e-4 / 2.2e-4 | 2.39e-3 / 2.39e-3 | **2.39e-3 / 2.39e-3** | same as og |

- The only delayed-rectifier-type K channel in the basal dendrites is SKv3_1, at 1/213 of the somatic density.
- **The apical Ca claim is confirmed.** Apical Ca_HVA2 and Ca_LVAst are both 0.0023897659621395684, identical to 17 digits, so they are one tied value, not two fitted ones.
  - The tie is identical in beta, gamma, delta and antic-delta. It started in beta: alpha had 5.37e-3 / 2.81e-3, and SSCx-og had 3.5e-6 / 9.9e-4.
  - The og-delta basal block copies the same value (8× gamma basal).

## 2. Basal bAP compared with data
bAP amplitude at the synapse (peak − rest, one somatic AP, median; Fig. A):

| Variant | 0–50 µm | 50–100 µm | 100–150 µm | 150–200 µm | 200–300 µm |
|---|---|---|---|---|---|
| og-delta | 78 | 56 | 35 | 34 | 15 |
| antic-delta | 83 | 65 | 41 | 42 | 17 |

- The 9-cell `screen_bd` results (s1) agree with the og-delta row: 58, 39 and 28 mV at 50–100, 100–150 and ≥200 µm.
- The only numeric target encoded (`score_bd` bap200) is ≤ 30 mV at ≥ 200 µm, from Nevian 2007 ([10.1038/nn1826](https://doi.org/10.1038/nn1826)) and Kampa. og-delta meets it.
- Nevian 2007 describes the bAP as "significantly attenuated", which is consistent with og-delta.
- Antic 2003 ([10.1113/jphysiol.2002.033746](https://doi.org/10.1113/jphysiol.2002.033746), VSD, relative units only) reports "very little amplitude modulation" within 200 µm. og-delta, which halves by about 100 µm, attenuates more steeply than that.
- **Verdict:** the og-delta basal bAP sits on the attenuating edge of the data, not outside it. antic-delta and Na×3 add only 7–10 mV at 50–150 µm.

## 3–4. Grid on the fitted protocols
Setup:
- 6 pairs and 38 basal synapses. Variants: basal Na ×{1, 3, 5} × SKv3_1 ×{1, 3} × Ka {og, antic}, plus antic-delta.
- Protocols:
  - Sjöström 2001: 5 pairings at 0.1 (a single pairing), 10, 20 and 50 Hz, ±10 ms. Markram 10 Hz is the same as the 10 Hz row.
  - Sjöström 2003: single post AP at −25 and −120 ms, and a 5-AP 20 Hz burst ending 120 ms before pre.
- Only spike-matched runs are used.
- Cells show the median per-synapse **effcai peak ratio** for the same synapses, og → variant (n synapses in brackets).
- The paired fraction A > B for og is 0.94–1.00 for every contrast except 0.1 Hz ±10 (0.71) and −10 50 vs 20 Hz (0.76). So the direction is consistent, but the effect size is small.

| Contrast (target) | og | Kv ×3 | Na ×3 | Na ×5 + Ka antic | antic-delta |
|---|---|---|---|---|---|
| +10 vs −10, 0.1 Hz (0.97 / 0.69) | 1.06 | 1.07 (38) | 1.13→1.15 (21) | 1.13→1.24 (11) | 1.14→1.21 (20) |
| +10 vs −10, 10 Hz (1.16 / 0.57) | 1.15 | 1.15 | 1.14→1.18 | 1.14→1.18 | 1.14→1.17 |
| +10 vs −10, 20 Hz (1.31 / 0.65) | 1.09 | 1.09 | 1.08→1.07 | 1.09→1.11 | 1.09→1.10 |
| +10 vs −10, 50 Hz (1.57 / 1.70) | 1.09 | 1.09 | 1.05→1.10 | 1.05→1.07 | 1.09→1.10 |
| −10, 50 vs 20 Hz (1.70 / 0.65) | 1.13 | 1.12 | 1.21→1.09 | 1.21→1.19 | 1.21→1.19 |
| −10, 50 vs 0.1 Hz | 3.13 | 3.14 | 2.28→2.38 | 2.19→2.40 | 2.70→2.53 |
| −25 vs −120 (0.65 / 1.0) | 1.05 | 1.05 | 1.05→1.05 | 1.05→1.05 | 1.05→1.05 |
| burst vs single −120 (0.79 / 1.0) | 1.37 | 1.36 | 1.95→1.68 | 1.00→1.00 (5) | 1.00→1.00 (5) |

**Firing.** These are the fractions of protocol runs with the intended somatic AP count, followed by f-I spikes relative to og (0.5 nA, 600 ms).

| Variant | AP count as intended | f-I vs og |
|---|---|---|
| og | 97% | 1.00 |
| Kv ×3 | 95% | 1.00 |
| Na ×3 | 55% | 0.98 |
| Na ×5 | 41% | 0.84 |
| Na ×5 + Ka antic | 61% | 0.83 |
| antic-delta | 68% | 0.83 |

- Na ×3 and ×5 with og Ka make every pulse a 2–6 AP burst, and the basal bAP becomes a regenerative spike: 70–90 mV at 150–300 µm.
- Antic Ka blocks the single 3 ms pulse in 3/6 cells (0 APs).
- The NMDA share of synaptic Ca (sj20 +10) is 95% for og, 88% for antic and 57% for Na ×5 (spiking).

## Verdict
- **Raising basal Kv does nothing.** SKv3_1 ×3 changes every contrast by ≤ 0.01. Basal SKv3_1 is already negligible, and repolarisation there is set by Ka and BK.
- **Raising basal Na (with the Ka that keeps the soma firing) gives only a small increase in timing separation.**
  - Single pairing +10/−10 goes from 1.13 to 1.24, and 10 Hz from 1.14 to 1.18.
  - 20 Hz and 50 Hz ±10 stay at 1.07–1.11.
  - −25 vs −120 stays at 1.05.
  - Burst vs single does not improve.
  - It costs spike fidelity: 30–60% of runs lose the intended AP count, and f-I drops by 17%.
- In og-delta, frequency (50 vs 0.1 Hz: ×3) and burst vs single (×1.37, VDCC ×4) are already strongly separated.
- The targets ask for opposite signs on effcai differences of 5–15%:
  - 20 Hz +10 is LTP 1.31 and −10 is LTD 0.65, at a 1.09 effcai ratio;
  - −25 ms is LTD and −120 ms is unchanged, at a 1.05 ratio;
  - 50 Hz −10 flips to LTP from 20 Hz −10 at a 1.13 ratio.
- No basal Na/Kv/Ka setting within the physiological range widens these gaps.
- **The rule is the limiter, not the cell.** Timing sign has to come from a non-Ca signal (pre-before-post trace, mGluR/eCB), consistent with the delta-bd and sp1 conclusions in DECISIONS.md.
