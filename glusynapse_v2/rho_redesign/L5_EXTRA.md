# L5->L5 extra paired targets (split1), 2026-10-01

Source: Chindemi's `biodata/paired_recordings.csv` L5TTPC->L5TTPC rows that no fit simulated (A3: CHINDEMI_DATA.md). 12 new protocols (+1 combined target).
LTD-relevant: Sjostrom 2001 20 Hz +25 / dt 0 / -25 (0.65 / 0.76 / 0.70), 40 and 50 Hz dt 0 (0.93, 0.92), 0.1 Hz dt 0 (0.91).

## Missing-row check (against TARGETS.md, PROTOCOLS.md, targets.py, Ebner csv)
| csv | protocol | existing? | new id | ratio +- sem (n) |
|---|---|---|---|---|
| mrk97_01 | 2 Hz, 5+5, +5 ms, 10 sweeps / 4 s | no (only 10 Hz) | 2Hz_5ms | 0.989 +- 0.040 (2) |
| mrk97_02 | 5 Hz | no | 5Hz_5ms | 1.016 +- 0.073 (5) |
| mrk97_03 | 10 Hz +5 | YES = 10Hz_5ms (1.2038) | - | - |
| mrk97_04/05/06 | 20 / 30 / 40 Hz +5 | no | 20Hz_5ms / 30Hz_5ms / 40Hz_5ms | 1.369 +- 0.091 (11) / 1.424 +- 0.075 (3) / 1.501 +- 0.120 (4) |
| sjs01_11 | 0.1 Hz, 1+1, dt 0, 50 sweeps | no (only +-10, -25, ...) | sjostrom_0.1hz_dt0ms | 0.910 +- 0.139 (5) |
| sjs01_12/13/14 | 20 Hz 5+5, +25 / 0 / -25, 15 sweeps / 10 s | no (only +-10) | sjostrom_20hz_dt+25ms / dt0ms / dt-25ms | 0.652 +- 0.184 / 0.761 +- 0.099 / 0.705 +- 0.070 (5 each) |
| sjs01_15/16/17 | 40 / 50 / 100 Hz, dt 0 | no | sjostrom_40hz_dt0ms / 50hz_dt0ms / 100hz_dt0ms | 0.929 +- 0.038 / 0.923 +- 0.035 / 1.245 +- 0.064 (5 each) |
| sjh06_01 | 50 Hz +10, 15 sweeps | same stimulus as sjostrom_50hz_dt+10ms (Sjostrom 2001: 1.57 +- 0.26, n=6) | none | 1.40 +- 0.06 (34) |

Sjostrom 2001 dt follows the same convention as the existing sjostrom_* ids (dt to the first post AP, dt_ref ap). Markram csv columns (period 4 s, 10 sweeps) = the
10Hz_* protocols (T 4000, nreps 10). The 1-pair 0.1 Hz protocol: the pair 207453-189325 cannot fire exactly one AP on split1 (as in sj03), so 23 of 24 pairs.

## Files
- yaml: /scratch/dhuruva/split1/configs/l5_extra.yaml (copy configs/L5extra_split1.yaml), same 24 pairs (Sjostrom2007_subset24_pairs.csv), delta-split1 circuit, label Sabrina, index_label L5extra_L5TTPC_L5TTPC. Stimulus keys 5ap_{20,40,50}Hz_3ms and 1ap_0.1Hz_3ms are reused from the single_cells pkls; 5ap_{2,5,30,100}Hz_3ms are calibrated by the sim job (per post cell).
- jobs: /scratch/dhuruva/split1/jobs/{sim,pf,ext}_l5_extra.sh: sim 22155650 -> prefire array 0-2 22155651 -> ext 22155652 (afterok chain; sizes in DECISIONS.md).
- extraction: glusynapse_v2/extracted/l5extra_delta-split1-prefire-vca (Markram ids without --window as markram_*, Sjostrom ids with --window); basis /scratch/dhuruva/split1/basis_l5l5 (unchanged).
- targets.py: `PAIRED_L5_EXTRA` (12 targets, group `paired_l5_extra`) and `PAIRED_L5_SJH06` (group `paired_l5_sjh06`: replaces the paired_l5 entry sjostrom_50hz_dt+10ms by the inverse-variance combination of Sjostrom 2001 and Sjostrom & Hausser 2006, ~1.41 +- 0.059, n=40; optional, keeps the target count at 12 new). Existing groups unchanged.

## Use in fit_v5 / split1 refits (39 -> 51 targets)
```
L5DIRS=$X/ebner_delta-split1-prefire-vca,$X/markram_delta-split1-prefire-vca,$X/sj03_delta-split1-prefire-vca,$X/sj03r50_delta-split1-prefire-vca,$X/sj07_delta-split1-prefire-vca,$X/l5extra_delta-split1-prefire-vca
L5GROUPS=paired_l5,sjostrom07,paired_l5_extra          # add ,paired_l5_sjh06 for the combined 50 Hz +10 value
```
(everything else as in run_split1_refits.sh COMMON; ANALYTICAL_BASIS_DIR=/scratch/dhuruva/split1/basis_l5l5). Target keys are protocol ids, so the 10Hz_*/sjostrom_* ids in other dirs do not clash.
Check after the chain: seff of the three jobs; number of records per protocol (24, or 23 for the 0.1 Hz dt 0); 100 Hz calibration failures (warning "no valid stimulus" in the sim log) would drop that protocol for the affected cells.

## Result (2026-10-01) and exact env
Measured: prefire 22155651 tasks 12:31/68.1, 19:12/81.6, 16:39/75.6 GB (88G request; next time 102G); ext 22155652 5:58, 45.0 GB at its 45G cap (next time 57G).
Records (of 24 pairs): 23 for 2/5/30/40Hz_5ms, sjostrom_0.1hz_dt0ms, 20hz +25/-25, 40hz/50hz dt0; 22 for 20Hz_5ms and sjostrom_20hz_dt0ms. Pair 207453-189325 (post 189325) has no valid stimulus for most keys (as in sj03); pair 192879-186028 failed 20Hz_5ms and 20hz_dt0ms.
sjostrom_100hz_dt0ms: no valid 5-AP 100 Hz stimulus for 14 of 24 post cells -> 10 records, a biased subset -> DROPPED from PAIRED_L5_EXTRA (11 extra targets, 39 -> 50). Its records stay in the extracted dir but are not scored.
```
L5DIRS=$X/ebner_delta-split1-prefire-vca,$X/markram_delta-split1-prefire-vca,$X/sj03_delta-split1-prefire-vca,$X/sj03r50_delta-split1-prefire-vca,$X/sj07_delta-split1-prefire-vca,$X/l5extra_delta-split1-prefire-vca
L5GROUPS=paired_l5,sjostrom07,paired_l5_extra      # +,paired_l5_sjh06 for the combined 50 Hz +10 value
ANALYTICAL_BASIS_DIR=/scratch/dhuruva/split1/basis_l5l5
```
