# Letzkus 2006 "distal" vs our L2/3 -> L5 "@distal" (2026-10-01)

Source: Letzkus, Kampa & Stuart 2006, **J Neurosci 26:10420-10429** (the brief's "Nat Neurosci" is wrong).
Local PDF: `plastyfire/papers/Letzkus et al. - 2006 - Learning Rules ... Dendritic Synapse Location.pdf`.

## 1. Letzkus' definition
- **Method:** synapse location is inferred from the **somatic uEPSP 10-90% rise time**. Nothing is recorded at the synapse.
  - Calibration, Fig 1F: 5 pairs filled with Alexa 594, 8 putative contacts (close appositions < 2 µm) found by two-photon imaging. Physical distance from the soma correlates linearly with rise time (r = 0.98). Somatic aEPSCs injected at known dendritic sites give the same line (r = 0.99), and so does the passive model.
  - Rise time stays a location measure after h-channel block (ZD7288; supplementary Fig 3).
  - Dendritic recordings (pipette on the apical trunk at 429 ± 13 µm, n = 16) see faster and larger uEPSPs than the soma does.
- **Cut-off:** **rise < 2.7 ms = proximal, > 2.7 ms = distal** (Fig 5; text p. 10425). 2.7 ms is the median of their 191 L2/3 -> L5 pairs (Fig 1C; range 0.7-5.5 ms). The paper gives no µm value. The Fig 1F fit (about 320 µm/ms) puts 2.7 ms at **about 400-450 µm** path distance, so the review's "~400 µm" is a fair reading.
- **Dendrites:** "Contacts were found on proximal oblique dendrites, the main apical trunk and secondary dendrites in the apical tuft". **Basal dendrites are not mentioned.**
  - The Fig 1F contacts run from about 100 µm to about 950 µm. The fast group (1.6-2.5 ms) sits at 100-330 µm; the slow group (4.1-4.2 ms) at 850-950 µm.
  - Their model (Fig 8A) uses an oblique at 200 µm as "proximal" and a **distal tuft site at 800 µm** as "distal".
- **Distance distribution:** only the rise-time histogram (Fig 1C, n = 191) is given, plus the 8 imaged contacts. Each pair had 1-2 contacts, at comparable distances, so inputs were localized.

## 2. Ours
- **Split:** `ebner/pair_geometry.py` sets `letzkus_distal = rise_10_90_ms > 2.7`.
  - The rise time is that of the somatic unitary EPSP in our L5 TTPC delta emodel: GluSynapse, 30 spikes at 1 Hz, averaged.
  - The output is `ebner/pair_geometry_L23PC_L5TTPC.csv`, read by `fit_v4.py` (`GEOM`), `gpu_v4_rho.py:288` and `eval_l23l5.py`. It is the same threshold as Letzkus, but in model milliseconds.
  - `letzkus_distal_median` (rise > our median) is the `--split median` alternative.
  - `sh_distal` (path_mean > 200 µm and EPSP < 1 mV) is for Sjöström & Häusser 2006, not Letzkus.
- **Pairs:** 120, of which 80 are distal (72 with all protocols) and 40 proximal. Synapses come from the circuit's L2/3 PC -> L5 TTPC contacts (edges `dhuruva_modified_edges_l23l5.h5`).
- **Locations** from the existing CSV. In SONATA, "apical" includes obliques, trunk and tuft.

| split | pairs | syn | basal | apical | path_mean q25 / median / q75 / max (µm) | pairs path_mean > 400 | > 600 |
|---|---|---|---|---|---|---|---|
| distal | 80 | 501 | 13% | 87% | 281 / 390 / 509 / 934 | 38 | 11 |
| proximal | 40 | 297 | 30% | 70% | 105 / 170 / 276 / 551 | 2 | 0 |

  - Our model rise times are slower than Letzkus': median 3.22 ms (vs their 2.7), up to 12 ms (vs their 5.5).
  - So a model 2.7 ms falls at a shorter distance than in the slice. Our distal median is 390 µm, and a quarter of distal pairs are below 281 µm. In Letzkus, 2.7 ms corresponds to about 420 µm.
- **Pending (job 22136361):** splitting apical into oblique / trunk / tuft. The output goes to `/scratch/dhuruva/letzkus_loc/`:
  - `summary.txt`;
  - `ours_synapses.csv`;
  - `ours_pairs.csv`;
  - `post_cells.csv`.
  
  Script: `letzkus_loc.py` / `run_letzkus_loc.sh`. It covers:
  - (A) class x split;
  - (B) our rise -> distance line vs Fig 1F;
  - (C) pairs passing apical > X;
  - (D) all circuit L2/3 PC contacts onto the same 120 posts at > 400 µm.

## 3. Verdict (preliminary; refine with job 22136361)
- **Partial mismatch, not a wholesale one.**
  - The review's claim that our distal synapses sit on basal/obliques is too strong. 87% of our distal synapses are apical, with a median pair distance of 390 µm, and 38/80 distal pairs have path_mean > 400 µm.
  - But the same nominal 2.7 ms cut selects synapses closer than Letzkus' cut. About half of our distal pairs sit at 280-400 µm, which in their calibration corresponds to about 2.3-2.7 ms, i.e. proximal.
  - Our distal set also lacks their 800-950 µm tuft tail (only 11 pairs are > 600 µm).
  - A basal share of 13% has no counterpart in Letzkus.
- **Minimal change that matches the preparation:** split by **anatomy, using Letzkus' own calibration**, instead of by model rise time.
  - Distal = all synapses of the pair apical (no basal) with path_mean > 430 µm (2.7 ms on Fig 1F). Proximal = path_mean < 430 µm.
  - Optionally drop pairs with basal contacts from both groups, since Letzkus saw none.
  - It is the same rule for every pathway and needs no new parameter. It replaces `letzkus_distal` in `pair_geometry.py` with a geometric column; neither the fitter nor the kernel changes.
- **Qualifying now:** 29 all_protocols pairs are all-apical with path_mean > 400 µm (9 at > 600 µm). That is enough for a distal group the size of Letzkus' (n = 7-10 per cell). The 430 µm count and the trunk/tuft breakdown come from job 22136361 (summary C).
- **Does the circuit have such contacts?** Yes: our own pairs reach a path_mean of 934 µm and a path_max of 1208 µm, and post apical dendrites reach about 1110 µm on average. The size of the full L2/3 PC pool beyond 400 µm, and how much of it is tuft, comes from job 22136361 (summary D). If C is too small, add pairs from that pool. This needs new workdirs and prefire, i.e. heavy-run sizing.
- **Caveat (from DISTAL_LTP.md):** placing the synapses correctly is necessary but not sufficient for 3AP -10 distal LTP. It also needs a burst-evoked dendritic Ca spike at those sites, which is an emodel matter.

## 4. Alternative geometric split (2026-10-01)
- **File:** `ebner/pair_geometry_L23PC_L5TTPC_geo.csv`, made by `letzkus_geo.py` (job 22136490, `run_letzkus_geo.sh`) from `ours_synapses.csv`. It is a copy of the default geometry csv with `letzkus_distal` replaced; the old column is kept as `letzkus_distal_rise_ms_split`.
- **Rule:** distal = at least half of the pair's synapses on apical dendrites (oblique, trunk, tuft) at path distance > 430 µm (Fig 1F: 2.7 ms ~ 400-450 µm). Proximal = the complement.
- **Not implemented:** excluding proximal pairs with at least half basal synapses beyond 200 µm. Dropping rows would also remove them from the non-Letzkus L2/3 targets that read the same csv. They are flagged in `basal_far_flag` only.
- **Use:** `GEOM_L23=/abs/path/to/..._geo.csv` in the environment of `run_fit_v4.sh` / `run_fit_v5.sh` (read by `fit_v4.py`, `fit_v4n.py` and `fit_v5.py`; unset = default file). The kernels are untouched.
- **Counts (all_protocols pairs):** see `/scratch/dhuruva/letzkus_loc/letzkus_geo_22136490.out`; fill here.
- **Rescore (MAXITER=0, seeded from each json, 39 targets, drop Letzkus 3AP -10 distal):**

| fit | job (old / new) | L5 chi2 | L2/3->L5 chi2 old split | L2/3->L5 chi2 new split | total old / new |
|---|---|---|---|---|---|
| v4 C1Ajn_s5 | 22136500 / 22136501 | 46.89 | | | 51.9757 (repro) / |
| v4 LM1_s5 (M1) | json / 22136502 | 154.87 | 2.38 | | 157.25 / |
| v5 V5_s6 | json / 22136503 | 149.62 | | | 165.98 / |

Fill from the `chi2 L2/3->L5` and `chi2 total` lines of `/project/rrg-emuller/dhuruva/plastyfire/logs/v4_fit_<job>.out` (v5_fit_ for 22136503). The old repro 22136500 must print 51.9757.
