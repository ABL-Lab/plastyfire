# Spine-Ca validation, Chindemi 2022 Supp Fig 2 replication (2026-10-01)

## Protocol (Chindemi et al. 2022 Methods "Postsynaptic calcium dynamics", via PubMed PMC9160074, doi 10.1038/s41467-022-30214-w; data: Sabatini, Oertner & Svoboda 2002, doi 10.1016/s0896-6273(02)00573-1, not in PMC, methods not read)
1. Pool: synapses of L5TTPC→L5TTPC connections, "following the acceptance criteria in Sabatini et al.": basal dendrites, diameter < 2 µm, 2 ≤ branch order ≤ 4, path length < 150 µm. The bAP pool is further restricted to path < 60 µm (Supp Fig 2f–i). n is not given in the text.
2. Synaptic event: "repeated presynaptic activation at 0.2 Hz" with stochastic release (TM + MVR, RRP ≈ 2). That is the circuit's own Use/Nrrp, one connection driven by its pre cell; it is not full RRP. Readout: mean spine Ca transient per synapse. Target: 0.67 ± 0.44 µM model vs 0.7 ± 0.4 data (mean ± SD).
3. bAP: one somatic AP. Target: 1.4 ± 0.6 µM model vs 1.7 ± 0.6 data.
4. Quantity: free spine [Ca] (cai_CR) peak; baseline subtraction is implied (Sabatini Δ[Ca]) but not stated. The text does not say whether failures are averaged in, so both readouts are reported.

## Implementation (`measure_cexp.py`, opt-in `CEXP_SABATINI=1`; default output unchanged)
- Condition "syn": `CEXP_SAB_NTRIALS` (default 10) pre spikes at 0.2 Hz from 1 s, no rho0/Use_p override, thresholds off.
  - Per trial k: peak cai_CR in [t_k, t_k + 200 ms] minus cai_CR at t_k − 1 ms.
  - `syn_dca_mean` is the mean over all trials. `syn_dca_succ` is the mean over trials with ΔCa > `CEXP_SAB_THR` (0.05 µM). Also `syn_psucc` and `syn_dca_trials`.
- Condition "post": the bAP, as the default cexp. `bap_dca` = peak after 990 ms minus cai_CR at 999 ms.
- Geometry columns: `secname`, `diam_um` (segment), `branch_order` (sections from soma; primary = 1), `dist_um` (path).
- `--merge` prints mean ± SD and n for the accepted pool and for all basal synapses (< 60 µm, and < 150 µm for synaptic), with PASS/FAIL within 1 SD of Sabatini.

## Runs (L5L5 delta-split1, `PATHS["L5"]`; outputs in /scratch/dhuruva/split1/cexp_sabatini/{ljp0,ljp25}/)
| run | globals | pilot (pair 0; 1 CPU, 2G, 0:15) | array | merge/summary |
|---|---|---|---|---|
| ljp0 | current (ljp 0, gca 0.0744) | 22199261 | (after pilot seff) | |
| ljp25 | E1 `spine/delta_ljp25.json` | 22199262 | (after pilot seff) | |

Array, after the pilot (sizing: pilot elapsed per pair × CHUNK + 50 %, MaxRSS + 25 %; the E1 L5 default arrays measured 8:40 and 1.58 GB at CHUNK 12):
`cd /project/rrg-emuller/dhuruva/plastyfire; CEXP_SABATINI=1 PATHNAME=L5 CHUNK=12 CEXP_OUT=/scratch/dhuruva/split1/cexp_sabatini/ljp0 sbatch --array=0-1 --job-name=sab_ljp0 --mem=<M> --time=<T> glusynapse_v2/cexp/run_measure_cexp.sh`. For ljp25, add `GLUSYN_GLOBALS=$PWD/glusynapse_v2/spine/delta_ljp25.json`. The 24 pairs fit in 2 tasks, and pair 0 is skipped because the pilot wrote it.
Summary: the same env with `MERGE=1 sbatch --dependency=afterok:<array> --mem=2G --time=00:15:00 ...`.

## Results (accepted pool; all-basal in brackets)
| run | synaptic all-trial, < 150 µm | synaptic successes, < 150 µm | synaptic < 60 µm | bAP < 60 µm | n syn / n bAP | PASS? |
|---|---|---|---|---|---|---|
| Chindemi model | 0.67 ± 0.44 | | | 1.4 ± 0.6 | | |
| Sabatini 2002 | 0.7 ± 0.4 | | | 1.7 ± 0.6 | | |
| ljp0 | | | | | | |
| ljp25 (E1) | | | | | | |

## Measured results (merges 22199568 ljp0 / 22199569 ljp25; 24 L5->L5 pairs, 191 syn, 101 Sabatini-accepted basal, 10 trials at 0.2 Hz)
| quantity (accepted pool) | ljp0 (current) | ljp25 (E1) | Chindemi model | Sabatini 2002 |
|---|---|---|---|---|
| synaptic, all trials, < 150 um | 0.70 +- 1.17 uM (n 85) | 1.04 +- 2.04 (n 85) | 0.67 +- 0.44 | 0.7 +- 0.4 |
| synaptic, successes, < 150 um | 0.91 +- 1.21 (n 84) | 1.24 +- 2.12 (n 84) | | |
| synaptic, all trials, < 60 um | 0.30 +- 0.14 (n 18) | 0.30 +- 0.14 (n 18) | | |
| **bAP, < 60 um** | **0.79 +- 0.57 (n 18): FAIL** | **1.22 +- 0.46 (n 18): PASS** | 1.4 +- 0.6 | 1.7 +- 0.6 |
| bAP, all basal < 60 um | 0.90 +- 0.55 (n 28) | 1.31 +- 0.44 (n 28) | | |

Verdict: synaptic Ca matches in both settings (mean on Chindemi's 0.67-0.7; our SD is larger, from distal shared-branch synapses). bAP Ca fails with the current globals (0.79 uM, about half of Chindemi's 1.4) and passes with E1 (1.22-1.31 uM). So the current model does NOT reproduce Chindemi's spine-Ca validation; E1 does. Arrays at the 2G limit (2.09 GB) -> 2.7G next time.
