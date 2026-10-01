# Chindemi 2022: Fig. 2c, training/test split, in-vitro data inventory, L2/3->L2/3 additions (2026-10-01)

Source: Chindemi et al. 2022 Nat Commun 13:3038, doi:10.1038/s41467-022-30214-w (papers/*.md + pdf). Data/code: Zenodo
doi:10.5281/zenodo.5654788 (analysis.zip 10.3 MB, cell_packages.zip 233.5 MB; v2.1 fetches in-vitro data from the original
authors' sites rather than shipping csvs). Our copy of the in-vitro table: biodata/paired_recordings.csv (42 rows).
"fig-read" = read off the low-resolution pdf figure (about +-0.03), not printed in the paper.

## 1. Fig. 2c: what it shows

- **What it is:** mean EPSP ratio of the model (in silico, blue) against the in-vitro data (orange) for the 5 **training**
  protocols (Table 1, dagger) after optimisation, with n = 100 in-silico pairs per protocol. Welch t-test n.s. for all.
- **Protocols:** they are exactly biodata/selected_paired_recordings.csv.

| Fig 2c row (top to bottom) | csv id | pathway | protocol | in vitro mean+-SEM (n) | in silico | p |
|---|---|---|---|---|---|---|
| L2/3 to L5, 50 Hz, 10 ms | sjh06_02 | L2/3 PC -> L5 TTPC (VC) | 5 pre + 5 post, +10 ms, 15 sweeps / 10 s | 1.06 +- 0.09 (19) | ~1.07 fig-read | 0.78 |
| L5 to L5, 10 Hz, -10 ms | mrk97_08 | L5 -> L5 (SSC) | 5+5, -10 ms, 10 sweeps / 4 s | 0.792 +- 0.026 (6) | ~0.80 fig-read | 0.84 |
| L5 to L5, 10 Hz, 10 ms | mrk97_07 | L5 -> L5 | 5+5, +10 ms | 1.201 +- 0.063 (6) | ~1.18 fig-read | 0.67 |
| L5 to L5, 5 Hz, 5 ms | mrk97_02 | L5 -> L5 | 5+5, +5 ms | 1.016 +- 0.073 (5) | ~1.0 fig-read | 0.63 |
| L5 to L5, 2 Hz, 5 ms | mrk97_01 | L5 -> L5 | 5+5, +5 ms | 0.989 +- 0.040 (2) | ~1.0 fig-read | 0.69 |

- **Optimisation (Methods "Model fitting", Table 2):**
  - **Free parameters (11):** tau* (calcium integrator time constant), the apical thresholds a00, a01 (theta_d) and a10, a11
    (theta_p), the basal thresholds b00 to b11, gamma_d and gamma_p.
  - **Best values:** tau* 278.318 (bounds 150-350, printed "s" in our md copy, but by the c* traces of Fig 2a it is ms,
    unverified), a = 1.127 / 2.456 / 5.236 / 1.782, b = 1.002 / 1.954 / 1.159 / 2.483, gamma_d 101.5, gamma_p 216.2.
  - **Algorithm:** a multi-objective GA (IBEA, ref 65) for 103 generations. A surrogate (boosted-tree) chimera was injected
    at generation 25. 10360 solutions were evaluated, and 83 were "valid" (every error < 1 SEM).
  - **Cost:** the per-protocol error is |R_sim - R_vitro| / SEM_vitro. The best solution minimises the **maximum** error over
    the 5 protocols (minimax), not the sum.
  - **Pairs:** 200 (100 L5-L5 and 100 L2/3-L5). The ratio uses the last 60 EPSPs of a 40 min follow-up over a pre-induction
    baseline.
- **Cpre / Cpost (Methods eq. 36-37):**
  - theta_d = a00 Cpre + a01 Cpost and theta_p = a10 Cpre + a11 Cpost, per synapse. Apical and basal synapses have separate
    coefficients; nothing else differs between them.
  - **Cpre:** the peak of the **calcium integrator c\*** in a short simulation where all synapses of the connection are
    activated at once with the **whole RRP released** (main text: "a single EPSP due to the release of all vesicles in the
    RRP").
  - **Cpost:** the peak of c\* in a simulation of a single induced postsynaptic spike (bAP).
  - **Calcium:** free spine [Ca]i from the NMDAR Ca current plus the VDCC current, with eta = 0.04 free fraction, tau_Ca =
    12 ms and rest 70 nM. c\* integrates [Ca]i - [Ca]rest.
  - **Mg block:** fitted to Vargas-Caballero data at 1 mM Mg. In-vitro simulations use in-vitro [Ca]o, and the low-Ca runs
    recompute Cpre and Cpost at 1.2 mM. The exact in-vitro [Ca]o used for fitting is not printed here (unverified).
- **Validation on the trained pathways (new pairs):**
  - **Fig 3d:** Markram +5 ms at 2-40 Hz (mrk97_01-06; above 10 Hz not trained). p 0.40-0.99.
  - **Fig 3e, f:** Markram 10 Hz at -10, +5 and +10.
  - **Supp A.4:** Sjostrom 2001 (sjs01_*), qualitative only. The model switches LTD to LTP near 20 Hz, against about 40 Hz
    in vitro.
  - **Fig 3j:** sjh06_02 (p 0.245).
  - **Fig 3l:** Sjostrom & Hausser distance dependence. The power-law exponents match (-0.32 +- 0.13 vs -0.41 +- 0.04), but
    that in-vitro set includes extracellular data.
- **L2/3->L2/3 transfer without refit (Fig 4, Supp A.6):**
  - **Egger egg99_01** (20 Hz, +10 ms): in silico ~1.75 fig-read against 1.30 +- 0.08. The difference is **significant**
    (p = 0.002, n = 99, one outlier removed), so the model overshoots LTP.
  - **Banerjee bnr14_01** (0.2 Hz, -15 ms, 100 pairings): in silico ~1.0 fig-read against 0.77 +- 0.07. It is n.s.
    (p = 0.069, n = 98, two outliers removed), but the model gives **no LTD**.
  - **Zilberter zlb09_*** (Supp A.6): **not reproduced**, with too much LTP. It matched only after excluding connections
    with low release probability at a cutoff the authors call too large to be reasonable.
  - **L4 PC -> L2/3 (Fig 5, Rodriguez-Moreno):** the LTD protocol matched, and MK-801 was emulated as gamma_d = 0.
  - **L4SS -> L4SS (Supp A.7):** not reproduced (the LTD is NMDAR-independent).
  - **Net:** Chindemi's single parameter set is LTP-biased on L2/3->L2/3, the same direction as our transfer misses.

## 2. In-vitro inventory

| file | what | rows |
|---|---|---|
| biodata/paired_recordings.csv | Chindemi's full in-vitro table (Table 1). Columns: protocol_id, pre/post mtype, area, animal, age, frequency_train, dt_train, n_stimuli_train, n_sweep, period_sweep, mean/sem_epsp_ratio, sample_size, reference, notes. Values digitised or from authors' sites | 42: mrk97 x8, egg99 x5, sjs01 x17, rmp08 x3, bnr14 x2, zlb09 x5, sjh06 x2 |
| biodata/selected_paired_recordings.csv | Chindemi training set = Fig 2c (mrk97_01, 02, 07, 08, sjh06_02) | 5 |
| biodata/README.md | source list + copyright note (data only for reproducing Chindemi figures), points to Zenodo 5654788 | - |
| biodata/recipe_{andras,mod,dhuruva}.csv | **Not plasticity data.** Synapse-physiology recipe per pre x post mtype (18 x 18 excitatory = 324): gsyn+-SD, nrrp, dtc, U, D, F (+SD), NMDA/AMPA (gsynSRSF), Hill coefficient, synapse_type, spine volume, correlations (u-gsyn 0.9, gsyn-nrrp 0.9, spinevol-nrrp 0.92, spinevol-gsyn 0.88 = Chindemi eq. 46-50) and distribution families. **andras** vs **mod**: only L5_TPC:A/B -> L5_TPC:A/B D 365.6 -> 671 ms. **dhuruva**: those 4 L5 rows re-parameterised (gsyn 1.6+-1.1, nrrp 2.9, U 0.5+-0.2, D 672.7, F 15.9). L2/3->L2/3 row is the same in all three (gsyn 1.0+-0.5 nS, nrrp 2.6, U 0.46+-0.26, D 671, F 17, SRSF 0.7, E2_L23PC). Circuit edges match recipe_andras (ebner/audit_l23l5/FINDINGS.md) | 324 each |
| evaluation_results/Chindemi_params*/, data/epsp_ratio_comparison.csv | our older in-silico reproductions (e.g. mrk97_08 sim 1.12 vs 0.79), not in-vitro data | - |

No other Chindemi in-vitro csv in the repo. configs/L23PC_L23PC_STDP.yaml (Egger-labelled 2-20 Hz grid) is a simulation grid, not data.

## 3. Coverage: every Chindemi in-vitro protocol vs our targets

Paired = presynaptic cell patched (fit target allowed). Ours = TARGETS.md / L23L23_TARGETS.md / targets.py.

| csv ids | paper (doi) | pathway | paired? | in our fit? |
|---|---|---|---|---|
| mrk97_03, 07, 08 (10 Hz +5/+10/-10) | Markram 1997, 10.1126/science.275.5297.213 | L5->L5 | yes | yes (10Hz_5ms, 10Hz_10ms, 10Hz_-10ms) |
| mrk97_01, 02 (2, 5 Hz +5) | same | L5->L5 | yes | **no**, although these are Chindemi training rows (cheap null/weak-LTP constraints) |
| mrk97_04-06 (20, 30, 40 Hz +5) | same | L5->L5 | yes | **no** (Chindemi validation, Fig 3d) |
| sjs01_01-10 (0.1-50 Hz +-10) | Sjostrom 2001, 10.1016/S0896-6273(01)00542-6 | L5->L5 | yes | yes (sj * +-10; ours re-digitised in ebner_targets.csv; differences <= 0.03, except the 40 Hz +10 SEM: 0.14 ours (text n=11) vs 0.113 csv (n=6)) |
| sjs01_11-17 (dt 0 at 0.1/20/40/50/100 Hz; 20 Hz +-25) | same | L5->L5 | yes | **no** (ours have 0.1 Hz -25 from Sjostrom 2003 instead) |
| sjh06_01 (L5->L5 50 Hz +10, 1.40+-0.06, n=34) | Sjostrom & Hausser 2006, 10.1016/j.neuron.2006.06.017 | L5->L5 | yes | **no** (our 50 Hz +10 L5 target is the Sjostrom 2001 value) |
| sjh06_02 | same | L2/3->L5 | yes | yes (L23 sj 50hz +10) |
| zlb09_01-05 | Zilberter 2009, 10.1093/cercor/bhn247 | L2/3->L2/3 | yes | yes (5 of our 15 Zilberter targets, same values) |
| egg99_01 | Egger 1999, 10.1038/16026 | L2/3->L2/3 | yes (see note) | added today: group paired_l23l23_egger (sims 22155160-63), not yet in a fit |
| bnr14_01, 02 | Banerjee 2014, 10.1002/phy2.271 | L2/3->L2/3 | **yes, paired** (see note) | listed as VALIDATION_L23L23 "extracellular": **wrong** |
| egg99_02-05 | Egger 1999 | L4SS->L4SS | yes | out of scope (cell type not modelled) |
| rmp08_01-03 | Rodriguez-Moreno & Paulsen 2008, 10.1038/nn.2125 | L4->L2/3 | yes (MK-801 in presynaptic pipette needs a patched pre cell) | out of scope |

Notes on the paired check:
- **Banerjee 2014** (PMC4002250 full text):
  - **Extracellular (validation only):** Figs 1-2 and 3A-B, and 4 (vertical/horizontal extracellular stimulation).
  - **Paired:** Fig 3C uses paired L2/3->L2/3 recordings. 216 pairs were tested, 13 connected and 11 used. Text: "In five
    pairs, the t-LTD protocol induced robust t-LTD (77 +- 7%, n = 5)", and MK-801 in the **presynaptic** pipette gave
    80 +- 6% (n = 6), not different from control.
  - **csv mapping:** bnr14_01/02 are exactly these paired values.
  - **Conditions:** mouse P12-18, **22-24 C**, EPSP slope, 100 pairings at 0.2 Hz, single pre and post spikes.
  - **dt:** -15 ms in Chindemi's csv, read from Fig 3C. The text says the paired protocol is "otherwise identical" to the
    extracellular t-LTD protocol, which is -10 ms. **dt unverified.**
- **Egger 1999:** the abstract (PMID 10570487) describes paired whole-cell recordings, but only of L4SS pairs. Two sources
  place the L2/3 PC row among paired recordings: Chindemi lists it with the "whole-cell paired recordings" datasets, and
  Feldmeyer et al. 2006 (PMC1819447) cites Egger 1999 among studies of L2/3 PC pairs. The full text is paywalled, so the
  L2/3 n, value and conditions come from the csv only (**L2/3 details unverified** against the paper).
- **Extracellular, all validation only:** Banerjee's extracellular horizontal L2/3->L2/3 data (PMC text): +10 ms 1.53 +- 0.05
  (5), -10 ms 0.68 +- 0.08 (5), +50 0.94 +- 0.08 (6), -50 1.09 +- 0.02 (6), -100 0.97 +- 0.03 (6), -200 1.04 +- 0.15 (5),
  D-AP5 1.01 +- 0.02 (4), postsynaptic MK-801 1.11 +- 0.08 (6). Like Nevian, they are all validation only.

## 4. L2/3->L2/3 targets to add

| id (proposed) | data | f | dt | n spikes | sweeps | ratio +- SEM | n | paired | condition / mapping | overlap with our 15 Zilberter |
|---|---|---|---|---|---|---|---|---|---|---|
| egger1999_5ap_20hz_dt+10ms | egg99_01 | 20 Hz | +10 | 5 pre + 5 post | 10 every 10 s | 1.299 +- 0.082 | 12 | yes (L2/3 details unverified) | control (already PAIRED_L23L23_EGGER) | same pattern as Z9 (5+5, 20 Hz, +10: 1.07 +- 0.11, n=6, 40 sweeps at 0.2 Hz). The two disagree (barrel vs visual cortex); Chindemi missed Egger by overshooting |
| banerjee2014_1ap_dt-15ms | bnr14_01 | 0.2 Hz | -15 (unverified, maybe -10) | 1 + 1 | 100 | 0.77 +- 0.07 | 5 | **yes** | control | close to Z3 (1+1, -10, 0.2 Hz, 40 pairings: 0.56 +- 0.06, n=4): same sign, Banerjee shallower |
| banerjee2014_1ap_dt-15ms, pre MK-801 | bnr14_02 | 0.2 Hz | same | 1 + 1 | 100 | 0.80 +- 0.06 | 6 | **yes** | **control, not nmdar_block**: MK-801 in the presynaptic pipette blocks only the presynaptic cell's NMDARs. Our rule has no presynaptic-NMDAR term, so the prediction equals bnr14_01. targets.py maps it to nmdar_block, which would force 1.0 | none; it only checks that there is no presynaptic-NMDAR dependence (automatically satisfied) |

Total: +1 informative protocol beyond Egger (bnr14_01; bnr14_02 duplicates it). The Banerjee prefire pilot is 22155162.
Caveat: Banerjee is mouse at room temperature. The 23 C objection that kept Hardingham's main set out applies here too
(Hardingham was also excluded for EGTA and no population mean), so a low-weight or report-only entry is defensible.

## Findings

1. Fig 2c shows the 5 Chindemi training protocols (= selected_paired_recordings.csv): Markram 2/5 Hz +5 and 10 Hz +-10 on
   L5->L5, and Sjostrom & Hausser 50 Hz +10 on L2/3->L5. The model is within ~1 SEM on each (p 0.63-0.84, n=100). The
   in-silico means are not printed; the fig-read values are within ~0.03.
2. Chindemi's cost is a minimax: the max over protocols of |dR|/SEM, optimised by GA over 11 parameters (tau*, 8
   apical/basal Cpre/Cpost coefficients, gamma_d, gamma_p). Our chi2 (sum of z^2) is a different objective, so best-fit
   comparisons should report both max|z| and chi2.
3. Cpre and Cpost are per-synapse **c\* peaks** (not raw [Ca] peaks): Cpre from a full-RRP single pre activation of all
   the connection's synapses, Cpost from a single postsynaptic AP; NMDAR-Ca + VDCC with eta 0.04 and tau_Ca 12 ms.
4. L2/3->L2/3 was a pure transfer test in Chindemi (no refit). Egger is overshot (~1.75 vs 1.30, p=0.002), Banerjee shows
   no LTD (~1.0 vs 0.77, n.s.), and Zilberter is not reproduced (too much LTP). The original rule is LTP-biased on this
   pathway too.
5. **Correction:** Banerjee bnr14_01/02 are **paired** L2/3->L2/3 recordings (Fig 3C), not extracellular. targets.py
   VALIDATION_L23L23, configs/L23L23extra_split1.yaml and L23L23_EXTRA.md label them extracellular. bnr14_02 must map to
   control, not nmdar_block. Only Banerjee's Figs 1-2/3A-B/4 data are extracellular.
6. Egger egg99_01 is paired per Chindemi and Feldmeyer 2006. The paper's L2/3 numbers could not be checked (paywall).
7. Paired L5->L5 rows in Chindemi's csv that we do not fit: mrk97_01/02 (Chindemi training), mrk97_04-06, sjs01_11-17
   (dt 0, 20 Hz +-25) and sjh06_01. They are cheap additions on existing pairs if the L5 set needs frequency or dt=0 nulls.
8. recipe_*.csv are synapse-physiology recipes (gsyn, nrrp, U/D/F, NMDA/AMPA, spine-volume correlations), not plasticity
   data. The only semantic differences are 4 L5_TPC->L5_TPC rows.
9. Zenodo 5654788 ships code and cell packages. The in-vitro values are equivalent to biodata/paired_recordings.csv, and
   the repo holds no other in-vitro csv.
