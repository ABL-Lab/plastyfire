# Model spine and shaft Ca vs experimental data (delta = ljp_VDCC 0, and ljp25_g031), 2026-09-30

## Verdict
- **ljp25_g031 is closer overall.** It does better on the quantities that set plasticity at real synapse locations:
  - proximal bAP amplitude;
  - bAP Ca at 50–150 µm (basal distance dependence);
  - the EPSP/bAP ratio.
- **delta's largest miss is bAP Ca beyond ~50 µm.** At 80–120 µm basal it gives 0.017 µM, against Cornelisse's 1.05 µM (range 0.59–7.98) at ~100 µm and Koester & Sakmann's "comparable amplitude up to 80 µm". That is a 35–60× shortfall. ljp25 gives 0.42 µM there, 1.4× below the data minimum.
- **ljp25 has a new, specific deviation: spine/shaft ratio.** Basal spine/shaft is 6–56 beyond 40 µm and apical is 8–200, while the data give 1.2–2.7. The ljp correction shifts only the spine GluSynapse VDCC; the emodel's dendritic Ca channels keep the uncorrected gating, so shaft Ca still falls ~e-fold per 9 mV. delta's spine/shaft ratio (0.8–2.1 basal) is roughly right. If shaft Ca ever becomes a rule drive, it needs the same correction.
- **Both variants are inside the data for:**
  - spine decay τ: 12.5–14 ms against 12–27 ms dye-free;
  - proximal shaft τ: 81 ms against 82 ms dye-free;
  - 3AP50/1AP burst ratio;
  - absolute EPSP Ca given release;
  - −10 ms linearity.
- **Pairing supralinearity (+10 ms):** delta is in the NS06 range (1.43–2.06 against 1.73 ± 0.28). ljp25 is below it (1.28–1.37), because the larger bAP Ca raises the denominator.
- **Both variants fail on high-frequency bursts:**
  - 3AP100 Ca charge is 2.1–2.5× a single AP, against ≈ 3.2× in NS06 (dye-loaded);
  - there is no distal supralinearity at 200 Hz (Kampa & Stuart); the model fires only 2 of 3 APs there and the charge is ≤ linear.
- **Apical Ca vs distance:** no verified dye-free data in the sources, so it is not scored. The model loses 3–4 orders of magnitude by 450 µm (delta) and 2 orders (ljp25).

## Data handling (judgement calls)
- **ΔF/F or ΔG/R vs µM.** Only Cornelisse (Table 1, one-compartment extrapolation to zero dye) and Sabatini (via Chindemi 2022) give dye-free µM. Nevian & Sakmann 2006 (NS06) is ΔG/R with 500 µM OGB-6F (low affinity, so near-linear, but heavily buffered), so it is used only as ratios.
- **Dye-loaded peaks.** With a large added buffer the decay is slow (Cornelisse raw τ is 91 ms at only 100 µM OGB-1), so dye-loaded peak ratios approximate ratios of *total Ca influx*. They are therefore compared with the model's **Ca charge** (∫ ica_VDCC + ica_NMDA; triangles in the figure). Dye-free µM peaks are compared with model `cai_CR` peaks.
- **Dye-free decay.** τ0 = τ·(1+κE)/(1+κE+κD) with Cornelisse's Table 1 κ values:
  - spine: 91.2·20/68 = 26.8 ms;
  - dendrite: 200.9·63/155 = 81.7 ms;
  - their own dye-free multicompartment model (Fig 7A) gives a spine 1/e decay of ~15 ms [digitised].
- **Dye-free burst ratio.** Projecting NS06's constant per-AP increment to dye-free peaks gives 3AP50/1AP = Σ e^(−20k/τ0) = 1.22 (τ0 12 ms) to 1.70 (τ0 27 ms). Cornelisse's dye-free model gives ≈ 1.15 for 50 Hz spines [digitised Fig 7A].
- **Temperature.**
  - Cornelisse 33–35 °C, NS06 32–36 °C.
  - Sabatini 2002 temperature is UNVERIFIED (full text not accessible); its τ = 12 ms is what GluSynapse uses.
  - Koester 1998 and Kampa 2006: temperature UNVERIFIED (abstracts only).
- **Cell type.**
  - Cornelisse: L5, mouse visual cortex, **P6–15** (immature), "secondary dendrites ~100 µm from soma"; Methods say oblique dendrites, so it is basal-or-oblique and is plotted on both panels.
  - NS06: L2/3, rat barrel cortex, P13–15, basal spines 50–150 µm, extracellular stimulation (validation-grade only).
  - Koester, Kampa: rat L5 basal.
  - Sabatini: rat CA1.
  - Model: L5 TTPC og-delta. The model's EPSP and pairing numbers come from 68 L5→L5 synapses (8 pairs).
- **NS06 nonlinearity factor.** The factor is per-spine (Fig 5E bottom). From the panel means it would be 1.23 at +10 and 0.74 at −10 (1 AP), so the plotted per-spine factor (1.73) is the upper reading.

## Table (model: median [IQR] over sites or synapses; L5 cells 181015-184976 and 182339-200396; "in?" = inside the data range)
| Quantity | Data (source, fig, prep, T, dye-corrected?) | delta | ljp25_g031 | in? d / l |
|---|---|---|---|---|
| bAP spine Ca, basal < 60 µm (µM) | 1.7 ± 0.6 (Sabatini 2002 via Chindemi 2022 text; CA1; T UNVERIFIED; yes). Chindemi model target 1.4 ± 0.6 | 0.86 ± 0.28 | 1.50 ± 0.32 | no (−0.8 SD) / yes |
| bAP spine Ca, ~100 µm (µM) | 1.05, range 0.59–7.98 (Cornelisse Tab 1; L5 P6–15 mouse; 33–35 °C; yes, 1-comp); their dye-free model 0.7 (text, Fig 7) | basal 80–120: 0.017 [0.007–0.036]; apical 0–100: 0.74 | basal 0.42 [0.23–0.63]; apical 1.49 | basal no / no (close); apical yes / yes |
| bAP shaft Ca, ~100 µm (µM) | 0.38, range 0.23–1.3 (Cornelisse Tab 1; yes) | basal 80–120: 0.021; apical 0–100: 0.18 | same (shaft unchanged) | no / no |
| Spine/shaft, 1 AP | 2.7 (Tab 1, 0-dye), 1.75 (Fig 7 model); ΔF/F 1.2 [digitised Fig 1E, 100 µM OGB-1] | basal 2.1 (< 40), 0.8–1.0 (40–250); apical 2.4–5 | basal 2.9 (< 40), 6–56 beyond; apical 8–200 | ~yes (close) / no |
| Spine decay τ (ms) | 12 (Sabatini, 0-dye); ~15 [digitised Cornelisse Fig 7A]; 26.8 (Cornelisse Tab 1 κ-corrected); raw 91 ± 13 | 12.5–13.7 (< 120 µm); rises to 18–29 distally (VDCC tail) | 12.9–14.1 | yes / yes |
| Shaft decay τ (ms) | 81.7 (Cornelisse Tab 1 κ-corrected; raw 201 ± 20) | 81 (< 40 µm); 39–51 (40–120) | same | proximal yes / yes; ~100 µm 2× fast |
| Basal distance, spine Ca at 40–80 µm / 0–40 µm | ≈ 1 ("comparable up to 80 µm", Koester & Sakmann 1998 abstract; L5 basal; OGB-type dye; ratio) | 0.09 | 0.67 | no / yes |
| Basal distance, 80–120 µm / 0–40 µm | NS06: 1-AP ΔG/R ≈ 0.75× EPSP at 50–150 µm [digitised Fig 5E]; Kampa 2006: single APs attenuate (abstract) | 0.015 | 0.28 | no / qualitative yes |
| Apical distance | no verified dye-free spine data (Cornelisse point only) | 0.74 → 4e-4 µM at 450–700 | 1.49 → 0.018 | not scored |
| 3AP50 / 1AP, peak (spine) | dye-free 1.22–1.70 (projection above); ≈ 1.15 [digitised Cornelisse Fig 7A] | 1.20 (< 40), 1.39–1.48 (80–160) | 1.04 (< 40), 1.19–1.25 | yes / yes (proximal edge) |
| 3AP50 / 1AP, Ca charge vs dye-loaded peak | 2.3 [digitised NS06 Fig 5E, 0.075/0.033 ΔG/R; L2/3; 32–36 °C] | 2.8 (< 40), 2.1–2.3 | 2.3 (< 40), 1.8–1.9 | yes / borderline |
| 3AP100 / 1AP, Ca charge | 3.2 [digitised NS06 Fig 5F] | 2.4–2.55 | 2.1–2.4 | no / no |
| 200 Hz burst, distal basal | supralinear distally, not proximally (Kampa & Stuart 2006 abstract; critical ~100 Hz) | 2 APs only; charge 2.0 → 1.65 (≤ linear 2) | 1.87 → 1.55 | no / no |
| EPSP spine Ca given release (µM) | 0.7 ± 0.4 (Sabatini via Chindemi; CA1; yes) | < 60: 0.35 ± 0.26; 50–150: 0.59 ± 0.44; all 0.52 | 0.36; 0.71; 0.58 | yes / yes |
| EPSP / 1 AP, peak | 0.41 (Sabatini 0.7/1.7, 0-dye) | 0.28 [0.2–0.57] (< 60) | 0.16 [0.14–0.3] | yes / no |
| EPSP / 1 AP, charge vs dye | 1.36 [digitised NS06 Fig 5E, 50–150 µm]; "comparable" ≈ 1 (Koester, < 80 µm) | 50–150: 22.7 [4.5–210]; < 60: 1.6 | 50–150: 2.0 [0.8–8.2]; < 60: 0.89 | no / yes (IQR) |
| +10 ms nonlinearity (1 AP) | 1.73 ± 0.28 [digitised NS06 Fig 5E]; 3 AP 1.8 ± 0.1 (text, n = 39) | peak 1.51 / 2.06; charge 1.43 / 1.56 (< 60 / 50–150) | peak 1.28 / 1.37; charge 1.33 / 1.33 | yes / no |
| −10 ms nonlinearity (1 AP) | 1.03 ± 0.18 [digitised Fig 5E]; linear (text) | peak 0.85 / 1.21; charge 1.05 / 1.29 | peak 0.88 / 0.94; charge 1.04 / 1.25 | yes / yes |

The model's ljp25 bAP Ca in basal 60–150 µm (diag_burst synapses) is 0.99 µM, against 0.07 for delta. The L2/3 cell (gid 143065, tree ≤ 65 µm) is also computed:
- spine 1 AP: delta 1.65–1.98 µM, ljp25 2.1–3.2;
- spine/shaft: 0.5–1.0 (delta), 1.1–1.3 (ljp25);
- τ: 12.6–13.3 ms;
- 3AP50 peak ratio: 1.0.

## Citations (DOI verified via PubMed metadata)
- Cornelisse et al. 2007 PLoS One 2:e1073, [10.1371/journal.pone.0001073](https://doi.org/10.1371/journal.pone.0001073) (pdf read: Tab 1, Figs 1, 7).
- Nevian & Sakmann 2006 J Neurosci 26:11001, [10.1523/JNEUROSCI.1749-06.2006](https://doi.org/10.1523/JNEUROSCI.1749-06.2006) (pdf read: Methods, Fig 5).
- Sabatini, Oertner & Svoboda 2002 Neuron 33:439, [10.1016/s0896-6273(02)00573-1](https://doi.org/10.1016/s0896-6273(02)00573-1). Abstract only; the numbers are from Chindemi 2022 (local md).
- Kampa & Stuart 2006 J Neurosci 26:7424, [10.1523/JNEUROSCI.3062-05.2006](https://doi.org/10.1523/JNEUROSCI.3062-05.2006). Abstract only; the PMC full text was blocked (captcha / 403).
- Koester & Sakmann 1998 PNAS 95:9596, [10.1073/pnas.95.16.9596](https://doi.org/10.1073/pnas.95.16.9596). Abstract only (same block).
- Nevian et al. 2007 Nat Neurosci 10:206, [10.1038/nn1826](https://doi.org/10.1038/nn1826). Abstract only (bAP "significantly attenuated"; no Ca numbers used).
- UNVERIFIED: the Sabatini temperature and the Kampa/Koester dyes and temperatures; the OGB-6F Kd (~3 µM, from memory) is used only qualitatively.

## Files and method
- **`ca_vs_data/bap_ca_map_glob.py`:** a copy of bap_ca_map.py plus:
  - `--glob` (GluSynapse GLOBALs);
  - the 3ap_100hz and 3ap_200hz protocols;
  - the 1/e decay of cai_CR and shaft cai;
  - 250 ms tails.
- **Run:** `run.sbatch` (job 22092211, 2 L5 cells + the L2/3 cell × {delta, ljp25_g031}, og-delta emodel). The raw output is in `out/bapmap_*.json`.
- **Exclusions:** 181015 fires 2 of 3 APs on 3AP50 and is excluded there. Both L5 cells fire 2 of 3 APs at 200 Hz (Letzkus pulse), so 200 Hz is analysed as 2 APs.
- **Analysis:** `analyze.py` produces `out/summary.txt` (all numbers above) and `figs/ca_vs_data.png` (one panel per quantity). EPSP and pairing numbers come from `ljp25_validation/rows_ljp25.json`, restricted to release in the EPSP, +10 and −10 runs (n = 7 < 60 µm, 18 at 50–150 µm).
- **Digitising:** `render.sbatch` and `crop.py` produce `figs/src/*.png`, the page crops used for the [digitised] values.

## Compute (all 1 CPU, def-emuller)
| Job | What | Elapsed | MaxRSS | CPU eff | Request |
|---|---|---|---|---|---|
| 22092211_0..3 | L5 maps (5 protos) | 4:18–5:02 | 1.11–1.15 GB | ~98% | 1.2G / 15 min (95% of memory; next run 1.45G) |
| 22092211_4..5 | L2/3 maps | 0:25–0:28 | 0.40–0.81 GB | 90–98% | 1.2G (next run 1.0G) |
| 22092241 | pdftoppm render | 0:06 | 43 MB | 83% | 1G (oversized; the crop run used 256M) |
| 22092315 | crop | 0:03 | 56 MB | – | 256M |
| 22092403 | analyze + plot | 0:05 | 143 MB | 40% (5 s job) | 1G, resized to 256M for 22092475 / 22093094 / 22093117 (seff unavailable: slurmdb connection refused) |

Total ≈ 0.3 CPU·h.
