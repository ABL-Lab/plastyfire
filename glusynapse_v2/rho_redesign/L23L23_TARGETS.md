# L2/3 PC -> L2/3 PC paired-recording targets (2026-10-01)

Sources (papers/): Zilberter et al. 2009 Cereb Cortex 19:2308 (Z) and Hardingham et al. 2007 J Neurophysiol 97:2965 (H).
Ratio = EPSP after / before. Every value below is printed in the paper text unless marked fig-read. Z reports mean +- SEM.
Model pathway: L2/3 PC (L2_TPC / L3_TPC, og-delta cADpyr_L2TPC / cADpyr_L3TPC) -> L2/3 PC in the same circuit as L2/3 -> L5.
Yaml: configs/Zilberter2009_L23PC_L23PC.yaml. Fit block: targets.py PAIRED_L23L23 (group `paired_l23l23`).

## Preparations

| | Zilberter 2009 | Hardingham 2007 |
|---|---|---|
| tissue | rat visual cortex, P14-21, 300 um parasagittal | rat visual cortex, P19-27, 400 um |
| temperature | 32-34 C | 23 C (main set); 11 extra at 35 C |
| ACSF Ca / Mg | 2 / 1 mM | 2.5 / 1 mM |
| pipette | K-gluconate, no Ca buffer (fura-2 100 uM only in imaging cells) | K-gluconate + 10 mM EGTA / 2 mM Ca |
| pairing | 40 pairings every 5 s (0.2 Hz); post APs by current pulses (width not given) | post AP 5 ms after pre AP, 20 pairs at 20 Hz per train, 10 trains at 0.5 Hz, 3 groups 1 min apart = 600 pairs |
| synapses | 5.0 +- 1.0 contacts, 36.5 +- 5.4 um from soma, mostly basal | N = 2.5 +- 0.2 release sites |

## Zilberter 2009: all paired results

dt > 0 = pre before post. "last" = dt to the last post AP. Sim id = configs/Zilberter2009_L23PC_L23PC.yaml.

| # | protocol | condition | mean +- SEM | n | where | sim id / status |
|---|---|---|---|---|---|---|
| Z1 | 1 pre + 1 post, +10 ms | control | 0.64 +- 0.07 | 6 | p.2311, Fig 2B, 3 | zilberter_1ap_dt+10ms: **fit** |
| Z2 | 1 + 1, +2-5 ms | control | 0.77 +- 0.09 | 5 | p.2311 (data not shown) | representable, not simulated (dt range only; would need an assumed dt) |
| Z3 | 1 + 1, -10 ms | control | 0.56 +- 0.06 | 4 | p.2311 (data not shown) | zilberter_1ap_dt-10ms: **fit** |
| Z4 | pre alone, same frequency | control | 0.98 +- 0.08 | 5 | p.2311, Suppl Fig 1B | zilberter_pre_only: **fit** |
| Z5 | 1 + 1, +10, post held at -45.5 +- 1.77 mV (100 ms step) | control | 0.84 +- 0.05 | 8 | p.2311, Fig 2C | no: needs a step + AP post stimulus (simwriter has step OR pulse) |
| Z6 | 1 + 1, +10, + extracellular EPSP 8.3 +- 0.65 mV | control | 0.82 +- 0.06 | 5 | p.2311, Fig 2D | no: extracellular co-stimulation |
| Z7 | 1 + 1, +10, Mg-free ACSF | control | 0.72 +- 0.07 | 5 | p.2311, Fig 2F | no: Mg-free needs new sims (no offline condition) |
| Z8 | 5 + 5 at 10 Hz, +10 | control | 0.76 +- 0.07 | 19 | p.2311, Fig 2G | zilberter_5ap_10hz_dt+10ms: **fit** |
| Z9 | 5 + 5 at 20 Hz, +10 | control | 1.07 +- 0.11 | 6 | p.2311, Fig 2H | zilberter_5ap_20hz_dt+10ms: **fit** |
| Z10 | 5 + 5 at 20 Hz, -10 | control | 0.93 +- 0.07 | 5 | p.2311 (data not shown) | zilberter_5ap_20hz_dt-10ms: **fit** |
| Z11 | 10 post APs at 50 Hz, pre 3-5 ms before the 10th (train-LTP) | control | 1.49 +- 0.12 | 11 | p.2312, Fig 4A,D | zilberter_train10_50hz_dt+4ms_last (Fig 5: 4 ms): **fit** |
| Z11b | train-LTP after a prolonged wait (D890 control) | control | 1.34 +- 0.11 | 5 | p.2315, Suppl Fig 1A | not fitted separately (same protocol; Z11 is the main group) |
| Z12 | pre 5 ms before the 1st AP of the train | control | 0.97 +- 0.06 | 4 | p.2312, Fig 4B | zilberter_train10_50hz_dt+5ms: **fit** |
| Z13 | train alone, no pre | control | 1.03 +- 0.04 | 4 | p.2312, Suppl Fig | zilberter_train10_50hz_post_only: **fit** |
| Z14 | pre 3-5 ms after the 10th AP | control | 0.99 +- 0.09 | 6 | p.2312, Fig 4D | zilberter_train10_50hz_dt-4ms_last: **fit** |
| Z15 | pre 5-12 ms after the 10th AP (train-LTD) | control | 0.72 +- 0.05 | 13 | p.2312, Fig 4C,D | zilberter_train10_50hz_dt-10ms_last (Fig 4C inset 10 ms): **fit** |
| Z16 | 4 post APs, pre 4 ms before the last | control | 0.76 +- 0.07 | 4 | p.2315, Fig 5D | zilberter_train4_50hz_dt+4ms_last: **fit** |
| Z17 | 8 post APs, pre 4 ms before the last | control | 1.15 +- 0.08 | 8 | p.2315, Fig 5D | zilberter_train8_50hz_dt+4ms_last: **fit** |
| Z18 | train-LTP, AM251 2 uM | CB1 block -> mglu_block | 1.73 +- 0.24 | 4 | p.2315, Fig 6C | **fit** |
| Z19 | train-LTD, AM251 2 uM | CB1 block -> mglu_block | 0.73 +- 0.07 | 7 | p.2315, Fig 6C | **fit** |
| Z20 | train-LTP, APV 50 uM | nmdar_block | 0.70 +- 0.07 | 4 | p.2315, Fig 6D | mapped, not fitted: nmdar_block makes the model exactly 1.0 (T25) |
| Z21 | train-LTD, APV 50 uM | nmdar_block | 0.73 +- 0.08 | 7 | p.2315, Fig 6D | mapped, not fitted: same |
| Z22 | train-LTD, CPCCOEt 25 uM + EGLU 50 uM | group I/II mGluR block | 1.02 +- 0.05 | 4 | p.2316-7, Fig 6D | no: postsynaptic mGluR; mglu_block is the presynaptic eCB/CB1 chain and AM251 (Z19) shows it is not that |
| Z23 | train-LTP, mGluR antagonists | same | 1.23 +- 0.12 | 4 | p.2317 (data not shown) | no: same |
| Z24 | train-LTP, D890 200 uM in post pipette (L-type VGCC) | VDCC block | 0.57 +- 0.07 | 5 | p.2315, Fig 5B | no: needs sims with VDCC block (cf. Letzkus Ni); D890 cuts train Ca to 0.37 +- 0.04 (n=4, Fig 5A) |
| Z25 | train-LTP, BAPTA 0.01 / 0.05 / 0.25 mM post | Ca buffer | ~0.97 / ~0.68 / ~0.82, fig-read +- 0.03 (bar on 0.05 ~ +- 0.1); text: 0.01 blocks LTP, 0.05 gives LTD, 0.25 blocks LTD | 3-11 per point | Fig 5E, p.2315 | no: new sims with a buffer |
| Z26 | train-LTP / train-LTD, Mg-free | Mg-free | 2.0-5.0 (range, n=3) / 1.39 +- 0.13 (n=6) | | p.2314 (data not shown) | no: Mg-free sims |
| Z27 | pre near the 8th AP (-5 < dt < 5) | control | no mean printed (Fig 5F points) | | Fig 5F | no |

Also printed, not targets: extracellular L2/3 +10 ms 1.34 +- 0.09 (n=15; gabazine 1.33 +- 0.15 n=5 vs 1.34 +- 0.12 n=10),
extracellular, so excluded. L5 -> L5 pairs: +10 with ext. EPSP 1.34 +- 0.01 (n=5, extracellular); **5 + 5, 10 Hz, -10 ms,
AM251: 1.07 +- 0.08 (n=3, p.2315)**, a paired L5->L5 drug result that could join paired_l5 later (not added here).
Expression locus (validation of the pre/post split, not a ratio target): LTP lowers PPR 1.1 +- 0.04 -> 0.87 +- 0.05 (n=26),
LTD leaves it 0.95 +- 0.04 -> 0.95 +- 0.05 (n=41), CV analysis agrees (Fig 6A,B), so LTP is presynaptic and LTD postsynaptic.

## Hardingham 2007: all paired results (none fitted)

| result | value | n | where |
|---|---|---|---|
| LTP group (>= +20%) | +109 +- 20% (ratio 2.09 +- 0.20) | 16 of 50 | p.2967, Fig 1A, 3 |
| LTD group (<= -20%) | -33 +- 2% (ratio 0.67 +- 0.02) | 20 of 50 | p.2967, 2972, Fig 1A, 5 |
| no change | no mean printed | 14 of 50 | p.2967 |
| 35 C set | 3 LTP / 4 nc / 4 LTD | 11 | p.2967 (data not shown) |
| zero EGTA / zero Ca pipette | 4 LTP / 4 LTD / 4 nc | 12 | p.2972 |
| no pairing, 1 h | 1 +- 12% depression | 10 | p.2972 (data not shown) |
| initial Pr by group | LTP 0.26 +- 0.03, nc 0.47 +- 0.04, LTD 0.69 +- 0.04; after pairing no difference | quantal subset 28 (9/9/10) | p.2969, Fig 2E |
| change vs initial amplitude / failure rate / Pr / PPR | r = -0.42 / +0.68 / -0.66 / +0.37 | 50 / 50 / 28 / 36 | Figs 1D, 1E, 2B, 2A |
| LTP locus | m +91 +- 19%, Pr 0.26 -> 0.44, 82 +- 8% presynaptic | 16 | p.2971, Fig 3C |
| initial amplitude | LTP 176 +- 39 uV, LTD 790 +- 139 uV, all 484 +- 69 uV | | p.2970, 2972 |

Why no fit target: (1) the paper gives the LTP/nc/LTD split, not a population mean; a mean needs the nc-group mean, which
is not printed. (2) 10 mM EGTA in the post pipette and 23 C have no model condition (the zero-EGTA set has only counts).
Use as validation: with our per-synapse U_SE, the predicted per-pair ratio for a 20 Hz, +5 ms, 600-pair protocol should
fall with U_SE (r ~ -0.66 with Pr), with potentiation carried by the presynaptic rule.

## Summary

- **15 fit targets** (13 protocols): 13 control + 2 AM251 (mglu_block), all from the Zilberter text.
- Mapped but not fitted (2): APV under nmdar_block. The data show NMDAR-independent LTD (APV turns train-LTP into LTD and
  leaves train-LTD at 0.73), but nmdar_block freezes rho and zeroes A_mglu/A_NO, so the model gives exactly 1.0. Fitting
  them would add a constant chi2 of ~30. Structural point: L2/3->L2/3 LTD is postsynaptic, mGluR-dependent, CB1- and
  NMDAR-independent; the model's post rule is NMDAR-Ca driven, so this LTD has to come from VDCC/bAP Ca in the model.
- Not representable now (11): step + AP depolarisation, extracellular co-stimulation, Mg-free (3), D890, BAPTA (3),
  postsynaptic mGluR block (2); plus Z2 (dt range) not simulated; all of Hardingham (validation only).
- Biggest test for the shared rule: single pairings give LTD at both +10 and -10 (symmetric, anti-Hebbian), and a
  preceding 50 Hz train switches +4 ms to LTP with an AP-count threshold (4 APs 0.76, 8 APs 1.15, 10 APs 1.49).
