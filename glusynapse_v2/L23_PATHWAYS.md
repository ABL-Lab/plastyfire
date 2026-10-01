# L2/3 pathways: paired-recording data and transfer test (2026-09-30)

Rule: fit targets come from **unitary paired recordings only**. Extracellular stimulation is validation only. No new protocols are proposed here. Only the existing L23PC→L5TTPC workdirs are used.

## 1. Data

### L2/3 PC → L5 TTPC, paired (unitary). Group `paired_l23l5` in `targets.py`, 9 targets

| # | Protocol (existing workdir id) | Δt | Pairing | Location | Ratio ± SEM | n | Source, figure | Text/digitised |
|---|---|---|---|---|---|---|---|---|
| 1 | `letzkus_1ap_dt+10ms`: 1 pre / 1 post AP | +10 | 100× at 1 Hz | pooled | 0.72 ± 0.03 | 15 | Letzkus 2006, Results, Fig 2 | text |
| 2 | `letzkus_3ap_200hz_dt+10ms`: 3 AP @ 200 Hz | +10 (to 1st AP) | 100–200× at 1 Hz | proximal (rise < 2.7 ms) | 1.30 ± 0.10 | 10 | Letzkus 2006, Fig 5 | text |
| 3 | same | +10 | same | distal (rise > 2.7 ms) | 0.79 ± 0.06 | 10 | Letzkus 2006, Fig 5 | text |
| 4 | `letzkus_3ap_200hz_dt-10ms` | −10 (to 1st AP) | same | proximal | 0.89 ± 0.03 | 16 | Letzkus 2006, Fig 5 | text |
| 5 | same | −10 | same | distal | 1.42 ± 0.09 | 7 | Letzkus 2006, Fig 5 | text |
| 6 | #3 with APV (`nmdar_block`) | +10 | same | distal | 1.03 ± 0.07 | fig only | Letzkus 2006, Fig 7A | text |
| 7 | #5 with APV (`nmdar_block`) | −10 | same | distal | 0.97 ± 0.06 | fig only | Letzkus 2006, Fig 7B | text |
| 8 | `sjostrom_50hz_dt+10ms`: 5 + 5 @ 50 Hz | +10 | 15× at 0.1 Hz | all unitary L2/3→L5 | 1.06 ± 0.09 | 19 | Sjöström & Häusser 2006, Results (Fig 1) | text |
| 9 | same | +10 | same | distal (> 200 µm, EPSP < 1 mV) | 0.86 ± 0.09 | 8 | S&H 2006, Fig 4C open blue circles | text |

Paired data that are **not scored**:
- **They would need a new condition or protocol:**
  - Letzkus NiCl₂ at distal inputs: 0.81 ± 0.01 at +10 and 0.99 ± 0.10 at −10. There is no Ni condition.
  - Letzkus no-post-activity 0.98 ± 0.06 (n 6), burst offset +500 ms 0.99 ± 0.05 (n 5), and 1-AP unpaired/offset 0.98 ± 0.04 (n 11).
  - Letzkus Fig 3: AP plus a dendritic pipette.
- **Derived only:** S&H proximal unitary ≈ 1.21 (n ≈ 11, no SEM), computed as (19·1.06 − 8·0.86)/11.
- **Excluded because extracellular stimulation is involved:**
  - S&H "boosted" distal unitary 1.26 ± 0.07 (n 4), boosted by extracellular stimulation.
  - S&H pooled proximal 1.37 and distal 0.80, which include extracellular experiments.
  - S&H Fig 5B unpaired distal EPSPs 0.71 ± 0.06 and AM251 0.99 ± 0.03. The preparation is not stated, and Fig 5A is extracellular.

### L2/3 PC → L2/3 PC, paired: **0 usable numbers**. Group `paired_l23l23` is empty
- **Zilberter et al. 2009**, Cereb Cortex 19:2308. Unitary L2/3 pyramid pairs, rat.
  - Low-frequency pairing gives only LTD, for either spike order.
  - Pairing becomes Hebbian after a preceding postsynaptic AP train.
  - Above 30 Hz there is LTP regardless of timing.
  - Numbers are paywalled, so I need the PDF.
- **Hardingham et al. 2007**, J Neurophysiol 97:2965 ([doi](https://doi.org/10.1152/jn.01352.2006)). L2/3 pairs; 600 paired APs in bursts of 20 at 20 Hz, post +5 ms.
  - A third of connections show LTP, a third LTD and a third no change. The sign follows initial Pr.
  - There is no population mean in the abstract. I need the PDF.
- Out of scope: Egger 1999 (L4→L4), and Banerjee, Rodríguez-Moreno & Paulsen, and Koester & Johnston (L4→L2/3).

### Extracellular stimulation: validation only
- **Nevian & Sakmann 2006:** L2/3 PC, basal pipette. 16 csv rows. Only pipette-mode pilot workdirs exist (3 posts, `Ebner2019_L23PC_L23PC_basal`).
- **Froemke & Dan 2002 and Froemke, Poo & Dan 2005:** L2/3 V1, afferent stimulation.
- **Feldman 2000 and Bender et al. 2006:** L4→L2/3.
- **Banerjee et al. 2014:** L4 vertical and L2/3 horizontal inputs.

I did not read these papers. The preparation is from the abstracts and web summaries; PubMed was disconnected.

### Disagreements with `ebner_targets.csv`
- The csv row `sjostrom_50hz_dt+10ms` with location "all" (1.57) is **L5→L5**. Never score it on L2/3→L5 pairs. The L2/3→L5 all-pairs value is 1.06 ± 0.09 (n 19), which is missing from the csv. It is added in `PAIRED_L23L5_EXTRA`.
- S&H distal: the csv already has 0.86 ± 0.09 (n 8). Ebner's 0.79 is a power-law fit value, not data.
- Letzkus: all rows match the text. Two points to note:
  - Fig 2D's rise-time correlation uses n = 20, while the pooled mean uses n = 15.
  - Ebner ran these protocols at 0.1 Hz, but the paper and our yaml use 1 Hz.

## 2. Transfer test design (L2/3→L5 only)
- **Workdirs:** all 4 protocols above exist on 120 pairs, 466 workdirs, `n120/seed20262009/Ebner2019_L23PC_L5TTPC`. The EPSP basis is done (`basis_results_edges_ebner_l23l5_delta`, 120 pairs). Pair splits are in `ebner/pair_geometry_L23PC_L5TTPC.csv`:
  - `letzkus_distal`: 80/120 pairs, 72 of them in `all_protocols`;
  - `letzkus_distal_median`;
  - `sh_distal`: 87 pairs, 79 in `all_protocols`;
  - `all_protocols`: 111 pairs.
- **The brief is out of date: the prefire is already running.** Job 22064522 (other session, `run_prefire_l23l5.sh`, 32 workers, 195G) was at 200/466 when I checked (about 01:10), with an ETA of about 25 min, at 100–260 s per workdir. Extraction 22064523 (262G) is `afterok` on it and writes `extracted/ebner_l23l5_delta-prefire`.
- **Problem for the owner session:** 53 of the first 200 workdirs failed, almost all Letzkus 3AP@200 Hz runs. The error is the spike-fidelity guard: "296 post spikes, expected 300 (99%)". The pool exits 0, so extraction still runs, but about 25% of burst records will be missing, which biases rows 2–7 toward cells that fire reliably. The fix is theirs, for example a tolerance on the guard.
- **Scoring:** `glusynapse_v2/eval_l23l5.py` (new). It takes the fit json unchanged and records `record_ratios` once, then selects pairs per target:
  - Letzkus rows use the rise-time split (`--split rise|median`);
  - the S&H distal row uses `sh_distal`;
  - pairs are restricted to `all_protocols`.

  It is needed because `BatchV2.objective` without a loc map would score `@proximal` targets on all pairs, and these targets need two different loc criteria. The sbatch for it was not submitted: auto-mode denied it.
- **Command, once the extract is done:** run it for td2_s1 and v22_s2, each with `--split rise` and `--split median`:
  ```
  ANALYTICAL_BASIS_DIR=basis_results_edges_ebner_l23l5_delta .venv/bin/python glusynapse_v2/eval_l23l5.py \
      --fit glusynapse_v2/results/reduced_gpu_subset_pl5_td2_s1.json --save glusynapse_v2/results/l23l5_td2_s1_rise
  ```
  It runs on CPU with numba and no jax. There is no comparable measured run, so per CLAUDE.md the first run is a pilot. My estimate from v22r_evalall (160 GiB for 4054 records, about 410 records here) is roughly 16 GiB. Size the pilot from that, then set the rest from its MaxRSS.
- **L2/3→L2/3:** transfer cannot be tested. There are no paired numbers and no paired workdirs, and I propose none.

## 3. Eval results
None yet. The extracted L2/3 records do not exist: `extracted/` has no l23 set until 22064523 finishes.

## 4. Predictions (from `cell_audit/LETZKUS_CA.md` and the cache note: c_post at the resting floor for 39% of synapses)
- **Proximal (#2, #4):** probably reachable. The bAP is about 74 mV below 200 µm, like the L5→L5 fit regime. The +10 burst should give LTP. With `t_drive 2`, post-before-pre bursts drive eCB-LTD at −10.
- **Distal Letzkus (#3, #5): not reachable.**
  - Beyond 450 µm the bAP is about 2 mV and 7 of 8 cells make no Ca spike, so both timings stay near 1. Expect z ≈ +3.5 at +10 and ≈ −4.7 at −10.
  - The in-silico split at 2.7 ms also puts 200–450 µm synapses in "distal". There the bAP is about 13 mV, which dilutes both signs.
  - This is a cell limit, not a fit issue.
- **Letzkus 1 AP +10 (#1, 0.72): likely a large miss (z ~ +8 to +9).** In the L5 fits, pre-before-post single APs give no change or LTP (Sjöström 0.1 Hz +10 = 0.97). Nothing local drives LTD at distal inputs.
- **S&H distal (#9, 0.86): possibly reachable with td2 only.** The cell-level post-AP eCB trace gates pre spikes at every synapse, so distal inputs with no post Ca can still get presynaptic LTD. That is also the S&H mechanism (CB1-dependent). v2.2's local drive should stay near 1.
- **APV rows (#6, #7):** 1.0 by construction, so z is +0.4 and −0.5.
- **L2/3 post cells:** the regenerated `cADpyr_L2TPC` and `L3TPC` delta hocs share all dendritic densities with `cADpyr_L5TPC`. Only the soma, axon, passive and Ih parameters differ. That includes uniform apical NaTg and the tie Ca_HVA2 = Ca_LVAst. Their bAP and spine Ca are unaudited, which matters only for future L2/3→L2/3 work.
