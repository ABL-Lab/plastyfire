# bAP-evoked Ca vs path distance in the L5TTPC delta emodel (spine and shaft), 2026-09-30

## Verdict
- **Comparison with Nevian 2007 is qualitative only.** Its full text is not on PMC and the figures are inaccessible, so nothing is digitised (only the abstract and text-level bounds). Nevian & Sakmann 2006 (the pdf) gives no Ca-vs-distance data, only ΔG/R in spines at 50–150 µm, with no µM scale.
- **bAP amplitude is within the text bounds.** The bAP is "significantly attenuated" in the data. The bound encoded in BASAL_NAKV.md (Nevian 2007 / Kampa) is ≤ 30 mV at ≥ 200 µm, and the model gives 26–27 mV at 160–250 µm.
- **Ca falls far more steeply than V.** Because of a Ca-vs-V e-fold of about 9 mV, a 55% loss of bAP amplitude by 100 µm becomes a 50–100× loss of Ca:
  - Basal spine Ca per bAP is 1.1 µM at < 40 µm, 0.10 at 40–80 µm and 0.007–0.017 at 80–160 µm. That is 1–2% of the proximal value at 100 µm.
  - Chindemi/Sabatini targets 1.4 ± 0.6 µM (< 60 µm) and says the bAP is ineffective beyond ~100 µm. The model is consistent with that, but it leaves essentially no bAP Ca at the 50–150 µm basal spines where Nevian & Sakmann measured AP-evoked spine transients.
- **Spine Ca follows the local peak V, not an independent shaft signal.**
  - Spine and shaft Ca have the same exponential dependence on the local V peak (e-fold 8–10 mV, both).
  - Spearman: spine~V 1.00, spine~shaft 0.98 (basal, 1 AP). The log-log slope of spine vs shaft is 0.99 basal and 0.8 apical.
  - Spine/shaft is ≈ 1 (0.7–2) in basal dendrites and 4–5 on the distal apical trunk.
- **og vs antic:** within the cell-mix noise. The one clear difference is in cell 181015-184976: antic-delta does not fire on the 1 AP pulse.
- **The L2/3 tree is too short to test the 50–150 µm range.**

## Method
- **Script and cells:** `bap_ca_map/bap_ca_map.py`. It uses the real post cell (bluecellulab, delta circuit, SK_E2 off, as in diag_burst.py). The cells are L5 181015-184976 and 182339-200396, and an L2/3 PC (gid 143065, Ebner L23PC_L23PC_basal pilot config). Variants are `delta` (= og-delta) and `antic_delta`, with the diag_burst VARIANTS applied to the basal Na/Ka.
- **Sites:** a GluSynapse (volume_CR 0.153 µm³, no NetCon) every ~20 µm path distance:
  - basal: the 6 terminal paths covering the most new segments (L5: 64 sites and 130 in total for the two cells);
  - apical: the trunk plus 3 further paths (L5: 121 and 223 sites in total).
- **Stimuli:** somatic pulses from each protocol's own config, 5 ms wide. 1 AP = nevian_1ap (1.33 / 1.89 nA); 3 AP 50 Hz = nevian_3ap_50hz; 5 AP 20 Hz uses the nevian_3ap_20hz amplitude for five pulses.
- **Two sims per (cell, variant, protocol):** with spines (local V, spine `cai_CR`, `ica_VDCC`) and without (the unperturbed shaft `cai`). The spine's VDCC current adds to the segment's ica, hence the separate no-spine run.
- **Reported values:** peak minus rest over the burst window. V is the local peak above rest.
- **Excluded runs:** cases where the soma fired fewer APs than pulses (the pair's own pulse, the same as in the induction runs). The excluded runs are:
  - 181015 delta 3AP50 (2 of 3 APs);
  - 181015 antic 1AP (0 of 1) and 3AP50 (2 of 3);
  - 182339 antic 5AP20 (4 of 5).
- **Caveats:** the absolute columns therefore mix different cells per protocol, so only the per-site ratios (burst / 1 AP) are cell-matched. The L5 pulses are strong (1.3–1.9 nA).
- **Files:** `out_*.json` (raw), `plot_summ.py`, `tables.md` (all bins, both variants), `run.sbatch`, `figs/bap_ca_{L5_basal,L5_apical,L23_basal}.png` (top: absolute values with the Chindemi/Sabatini point and the ≤ 30 mV bound; bottom: normalised to the most proximal bin, plus spine/shaft).

## L5 basal, og-delta (median over sites; µM; 3AP50/1AP and 5AP20/1AP are per-site ratios of the spine Ca peak)
| bin µm | n | bAP mV | shaft 1AP | spine 1AP | spine / prox | 3AP50/1AP | 5AP20/1AP | spine/shaft | Data |
|---|---|---|---|---|---|---|---|---|---|
| 0–40 | 21 | 85 | 0.49 | 1.12 | 1 | 1.20 | 1.03 | 2.1 | Chindemi/Sabatini spine 1.4 ± 0.6 µM at < 60 µm (text) |
| 40–80 | 23 | 56 | 0.16 | 0.10 | 0.09 | 1.29 | 1.01 | 0.85 | no numbers |
| 80–120 | 23 | 43 | 0.021 | 0.017 | 0.015 | 1.39 | 1.00 | 0.87 | NS06 spines at 50–150 µm respond to APs [text], no numbers |
| 120–160 | 22 | 35 | 0.0076 | 0.0069 | 0.006 | 1.48 | 0.88 | 0.82 | idem |
| 160–250 | 37 | 26 | 0.0032 | 0.0026 | 0.002 | 1.56 | 1.00 | 1.0 | Nevian 2007 ≤ 30 mV at ≥ 200 µm (bound) |
| 250–400 | 4 | 11 | 0.0003 | 0.0004 | 0.0004 | 1.58 | 1.59 | 1.4 | idem |

bAP relative to proximal: 1 / 0.65 / 0.50 / 0.41 / 0.31 / 0.13. antic-delta (1 AP from cell 182339 only): bAP 69 / 54 / 44 / 39 / 27 / 13 mV; spine 1 AP 0.55 / 0.07 / 0.022 / 0.010 / 0.0024 / 0.0004 µM. Same shape and same burst ratios (1.2–1.6).

## L5 apical (trunk + obliques), og-delta
| bin µm | n | bAP mV | shaft 1AP | spine 1AP | 3AP50/1AP | 5AP20/1AP | spine/shaft | Data |
|---|---|---|---|---|---|---|---|---|
| 0–100 | 13 | 73 | 0.18 | 0.74 | 1.63 | 1.04 | 4.0 | none digitised |
| 100–250 | 22 | 50 | 0.060 | 0.071 | 2.2 | 1.13 | 2.4 | idem |
| 250–450 | 23 | 23 | 0.00037 | 0.0020 | 1.66 | 1.25 | 4.4 | idem |
| 450–700 | 46 | 8.9 | 0.00008 | 0.00044 | 1.33 | 0.99 | 5.1 | idem |
| 700–1000 | 85 | 3.6 | 0.00006 | 0.00020 | 1.25 | 1.10 | 4.1 | idem |

- **Ca falls faster than V:** the bAP falls to 5% of the proximal value by 700–1000 µm, while spine Ca falls to 3e-4 of proximal at 450–700 µm.
- **bAP amplitude:** the literature cited in LETZKUS_CA.md (Stuart 1997; Larkum 1999) gives a few tens of mV at 450–650 µm. The model's 9 mV at 450–700 µm is several-fold too low, matching the earlier letzkus_ca result.
- **Regenerative events:** at some 100–250 µm apical sites a regenerative Ca event occurs. The 3AP50 median spine Ca is 4.2 µM (og-delta), and antic-delta 182339 has a 1 AP shaft median of 8.8 µM. These bins are dominated by a few sites and are not reliable medians.
- **Burst gain:** 3AP50/1AP spine Ca is 1.3–2.2 and never rises with distance. Beyond ~450 µm the spine Ca is at leak level, as in LETZKUS_CA.md.

## L2/3 basal (one PC, gid 143065, pulse 0.5 nA / 5 ms from the pipette config)
| bin µm | n | bAP mV | shaft 1AP | spine 1AP | 3AP50/1AP | spine/shaft |
|---|---|---|---|---|---|---|
| 0–40 | 11 | 100 (antic 103) | 1.9 | 2.0 | 1.0 | 1.0 |
| 40–80 | 6 | 83 (antic 90) | 3.2 | 1.65 | 1.0 | 0.5 |

- The built tree reaches only 65 µm in path distance, so the 50–150 µm basal-spine range of Nevian & Sakmann 2006 is not reachable. This is the tree's real extent.
- **Spine Ca:** 1.6–2.0 µM at 0–80 µm, at the Sabatini/Chindemi level (1.4–1.7), and about 16× the L5 basal value at 40–80 µm (bAP 83 vs 56 mV).
- **Bursts add nothing:** 3AP50/1AP and 5AP20/1AP are 1.0, so the Ca peak is set by the first AP, and the cai_CR decay (τ 12 ms) prevents summation.
- **Sample size:** n = 17, so the spine~shaft correlation there (−0.4) is meaningless.

## Deviations from the data and constraints
1. **Distance scale.** The model bAP Ca is essentially gone (< 2% of proximal) beyond ~100 µm on basal dendrites. This agrees with Chindemi's own statement but not with Nevian & Sakmann's spine responses at 50–150 µm. It is the unchecked deviation to verify against the Nevian 2007 figures.
2. **Apical trunk.** The bAP falls to ≈ 9 mV at 450–700 µm. Ca in spines and shaft is 3–4 orders of magnitude below proximal.
3. **The spine has no Ca "filter" of the V decay.** Its Ca is a fixed exponential of local V, so any V error (basal Na/Ka, apical Na) is amplified ~e^(ΔV/9 mV). Ca-vs-distance can only be corrected by fixing V (bAP), or the VDCC voltage dependence (owned by the parallel vdcc_decode audit).
4. **Burst boost.** 3AP50/1AP is only 1.2–1.6 and 5AP20/1AP is ≈ 1.0, because cai_CR decays with τ 12 ms and the pulse gives at most 2 APs in some cells.

## Compute (all 1 CPU)
| Job | Task | Elapsed | MaxRSS | CPU eff | Request |
|---|---|---|---|---|---|
| 22088253 pilot (181015 delta) | – | 2:28 | 1.07 GB | 97% | 1.5G / 30 min |
| 22088355_0..3 (L5 2 cells × 2 variants) | – | 2:24–2:29 | 0.71–0.94 GB | 97–98% | 1.4G / 15 min |
| 22088355_4..5 (L2/3) | – | 0:26–0:27 | 0.33–0.53 GB | 88–89% | 1.4G / 15 min |

Total ≈ 0.15 CPU·h. The requests came from the pilot (MaxRSS 1.12 GB + 25% = 1.4G; elapsed 2.5 min + 50% rounded up to 15 min). The measured usage is below 60% of the request for the L5 tasks and the L2/3 tasks. Next run: L5 ≈ 1.2G, L2/3 ≈ 0.7G. The 15 min limit is the minimum step.
