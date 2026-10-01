# Sjostrom 2001 dt 0 / 20 Hz +-25 targets: check (2026-10-01)

Paper: Sjostrom, Turrigiano & Nelson 2001 Neuron 32:1149 ([doi:10.1016/S0896-6273(01)00542-6](https://doi.org/10.1016/S0896-6273(01)00542-6); PubMed 11754844, not in PMC).
Full text read from the author's lab PDF (plasticity.muhc.mcgill.ca/DataPage/PDFs/Sjostrom_Neuron_2001.pdf). The Fig 7D source data come from the lab data page
(DataPage/Sjostrom_2001_Fig_7D/Synch_firing.xlsx: freq 0.1/20/40/100 Hz, mean 91.01/76.08/92.89/124.46 %, SEM 13.91/9.94/3.84/6.36). These match csv sjs01_11/13/15/17 exactly.

## What the paper's protocol is
- The dt 0 rows are real pairings: **Fig 7D "Synchronous firing" (dt = 0 ms)**, also plotted in the Fig 7C timing curves (panels 0.1 Hz, 20 Hz and 40/50 Hz).
  Pre and post APs are evoked simultaneously by 5 ms current steps (0.8-1.5 nA). The random-firing paradigm is a separate set (Fig 8C: Gaussian jitter SD 7 ms,
  mean dt 0, at 0.1/20/35/50 Hz, n 3/3/4/7). The csv does not contain the Fig 8C data.
- Induction: 0.1 Hz = 50 single pairings. 10 Hz and above = 5 pre + 5 post spikes x 15 trains at 0.1 Hz (75 pairings). Spike delays were shifted to set dt.
- n: Fig 7D legend says "n's are between 3 and 5 per point". Fig 7C says 3-7 per point, except 40/50 Hz +10 and -10 (n = 12 and 15). The csv n = 5 is the upper bound for dt 0.
- The paper makes a point of the 40 Hz dt 0 result. Text p 1156: "synchronous firing (dt = 0 ms) produced no LTP at 40 Hz (Figures 7C and 7D)". At 40 Hz each
  post spike is +-25 ms from the neighbouring pre spikes, outside the ~20 ms LTP window. At +-10 the neighbours fall at +15 or +10 ms, which gives LTP. The legend:
  "Synchronous firing at frequencies up to 40 Hz resulted in depression, whereas potentiation was obtained at 100 Hz." So the 0.93 next to 1.57/1.69 is the real finding.
- SEMs: 40 Hz 3.8 % with n 3-5 means SD ~7-9 %. That is tight but plausible, and it matches the author's own file, so it is not a digitisation error.

## Our yaml (/scratch/dhuruva/split1/configs/l5_extra.yaml = configs/L5extra_split1.yaml)
dt 0 = pre spike at the first post AP (dt_ref ap), 5+5 APs, 15 sweeps every 10 s (0.1 Hz: 1+1, 50 sweeps). This is the same paradigm as the paper.
The only difference is the pulse width: 3 ms with amp "find" 0.5-3 nA, against the paper's 5 ms at 0.8-1.5 nA. That is the same as every other sjostrom_* L5 sim, and with
dt referenced to the AP it does not change the timing. It cannot explain z = 11.5.

## Verdict per row
| row | paper | yaml same? | verdict |
|---|---|---|---|
| sjostrom_40hz_dt0ms 0.929+-0.038 (sjs01_15) | Fig 7D 40 Hz point, synchronous, n 3-5 | yes | **keep**: a genuine model failure (too wide or too frequency-driven an LTP window), not a label error |
| sjostrom_50hz_dt0ms 0.923+-0.035 (sjs01_16) | **the same data point**: Fig 7C "40/50 Hz" panel digitised (0.923) = Fig 7D 40 Hz (0.929). No separate 50 Hz dt 0 exists | yes | **drop**: a duplicate that double-counts one measurement (about 132 of the 264 z^2) |
| sjostrom_20hz_dt0ms 0.761+-0.099 (sjs01_13) | Fig 7C/7D 20 Hz synchronous: more LTD than 0.1 Hz (post AP sits in the LTD window of the next pre, -50 ms) | yes | **keep** |
| sjostrom_20hz_dt+25ms 0.652+-0.184 (sjs01_12) | Fig 7C 20 Hz; LTD, because +25 = -25 relative to the next pre spike at 20 Hz | yes | **keep** |
| sjostrom_20hz_dt-25ms 0.705+-0.070 (sjs01_14) | Fig 7C 20 Hz, LTD | yes | **keep** |
| sjostrom_0.1hz_dt0ms 0.910+-0.139 (sjs01_11) | Fig 7C 0.1 Hz / 7D, 50 synchronous pairings, no significant change | yes | **keep** (low weight) |
| (sjs01_17 100 Hz dt 0, 1.245+-0.064) | Fig 7D, synchronous at 100 Hz, LTP | yes | already dropped (10/24 cells); the paper's LTP point is real |

## Notes
- Chindemi 2022 ([doi:10.1038/s41467-022-30214-w](https://doi.org/10.1038/s41467-022-30214-w), PMC9160074, via PubMed) used sjs01 only as held-out validation (Supp Fig, not in the PMC text).
  It reports that its model switches LTD to LTP at a lower frequency than in vitro and that STDP is abolished above 20 Hz. That is the same failure as our 40 Hz dt 0 and 20 Hz dt 0 rows.
- After the duplicate is dropped, 40 Hz dt 0 still gives z ~11. If one SEM dominates, the options are a floor on the SEM or accepting the residual; that is a scoring decision, not a data fix.
- Optional new targets: Fig 8C random firing (0.1/20/35/50 Hz: ~0.77/0.86/1.02/1.36, n 3/3/4/7; values fig-read).
