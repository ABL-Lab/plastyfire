# L2/3->L2/3 extra protocols (split1), 2026-10-01

## Mapping of Chindemi's paired_recordings.csv L2/3->L2/3 rows
| csv | paper | protocol | status |
|---|---|---|---|
| zlb09_01..05 | Zilberter 2009 | 1ap +10 (0.64), 1ap -10 (0.56), 5ap 10 Hz +10 (0.76), 5ap 20 Hz +10 (1.07), 5ap 20 Hz -10 (0.93) | already the Zilberter targets zilberter_1ap_dt+10ms / _dt-10ms / _5ap_10hz_dt+10ms / _5ap_20hz_dt+10ms / _5ap_20hz_dt-10ms (same values; 0.2 Hz, 40 sweeps). The 10 other Zilberter targets (pre-only, trains, AM251) are not in the csv |
| egg99_01 | Egger 1999 Nat Neurosci 2:1098 (PMID 10570487, paired whole-cell recordings, rat barrel cortex P14; abstract names L4 spiny stellate, the L2/3 PC row is Chindemi's) | 5 pre + 5 post at 20 Hz, +10 ms, 10 sweeps / 10 s: 1.299 +- 0.082 (n=12) | NEW fit target egger1999_5ap_20hz_dt+10ms (group paired_l23l23_egger) |
| bnr14_01/02 | Banerjee 2014 (Fig 3C), paired L2/3 pairs, mouse, room temperature 22-24 C | 1 pre + 1 post, -15 ms (csv; unverified, text implies -10), 100 sweeps at 0.2 Hz: 0.77 +- 0.07 (n=5); bnr14_02 = same with MK-801 in the PRESYNAPTIC pipette only, 0.80 +- 0.06 (n=6), same model protocol (control) | VALIDATION ONLY (mouse, RT, like Hardingham 23 C): VALIDATION_L23L23 in targets.py (one control entry); protocol banerjee2014_1ap_dt-15ms simulated |

Egger egg99_01 is confirmed a paired recording (A3): fit target.

## Files
- yaml: /scratch/dhuruva/split1/configs/l23l23_extra.yaml (copy: configs/L23L23extra_split1.yaml), label L23L23extra, same 120 pairs as l23l23, single_cells pkls copied from the Zilberter run (same stimulus keys, no recalibration).
- sims: /scratch/dhuruva/split1/refitting_results/fitting/n120/seed20262009/L23L23extra/simulations; job scripts /scratch/dhuruva/split1/jobs/{sim_l23l23x,pf_l23l23x_egg,pf_l23l23x_bnr_pilot,ext_l23l23x_egg}.sh
- extraction: glusynapse_v2/extracted/l23l23extra_delta-split1-prefire-vseg-rs (Egger now; Banerjee after its prefire). Basis: reuse /scratch/dhuruva/split1/basis_l23l23 (per pair, read from any workdir of the pair; pair set unchanged).

## Jobs (all rrg-emuller)
sim 22155160 (4 CPU 3G 0:15) -> pf_egg 22155161 (8 workers 100G 0:30) -> ext_egg 22155163 (4 workers 44G 0:15); Banerjee prefire pilot 22155162 (4 pairs, 4 workers 60G 0:15) afterok sim.
Sizing basis (Zilberter split1): sim 22136808 2.3 GB MaxRSS (calibration dominated, now cached); prefire 22136811 1560 workdirs, 32 workers, 1:16, 310G (limit hit, ~9.7+ GB/worker), 38.8 CPU-h = 1.5 CPU-min/workdir, bio 200 s/workdir; Egger bio 100 s -> ~1 CPU-min x 120 = 2 CPU-h; ext 22136812 1527 records 47:50 4 workers 44G (limit hit) = 32 rec/min -> 120 records ~4 min. Banerjee bio 500 s (2.5x): pilot first for RAM/time, then full 120 pairs with --skip-existing.

## Use in a fit
EXTRA="l23l23:paired_l23l23,paired_l23l23_egger:$X/zilberter_l23l23_delta-split1-prefire-vseg-rs,$X/l23l23extra_delta-split1-prefire-vseg-rs:/scratch/dhuruva/split1/basis_l23l23"
(groups and dirs are comma lists; as in run_split1_refits.sh variable Z). 54 -> 55 targets (one chi2 term, 5 stim; value 1.299 is a ratio of means over 10 sweeps). Banerjee (paired, validation only) is not in any group; score by hand against VALIDATION_L23L23.

## Results and exact use (2026-10-01)
Measured: Egger pf 22155161 9:09 56.9 GB; Egger ext 22155163 1:20 27.1 GB; Banerjee pf 22155971 22:40 80.0 GB of 80G (at cap, next time 100G); Banerjee ext 22155972 4:39 60.0 GB of 60G (at cap, next time 75G). Records: 118 egger1999 + 118 banerjee2014 of 120 pairs each (236 files in the dir, no overlap; 2 pairs per protocol dropped by the prefire spike-count guardrail).

Egger fit target (55 targets):
EXTRA="l23l23:paired_l23l23,paired_l23l23_egger:$X/zilberter_l23l23_delta-split1-prefire-vseg-rs,$X/l23l23extra_delta-split1-prefire-vseg-rs:/scratch/dhuruva/split1/basis_l23l23"

Banerjee validation only (not in any group, never fitted): score a finished fit's params on `banerjee2014_1ap_dt-15ms` with the extracted dir `$X/l23l23extra_delta-split1-prefire-vseg-rs`, basis `/scratch/dhuruva/split1/basis_l23l23`, target `targets.VALIDATION_L23L23[("banerjee2014_1ap_dt-15ms","control")]` = 0.77 +- 0.07 (bnr14_02 0.80 +- 0.06 is the same protocol, compare to both). Env for a score-only run (MAXITER 0, as the S1_*_r rescores): EXTRA="banerjee:paired_l23l23_val:$X/l23l23extra_delta-split1-prefire-vseg-rs:/scratch/dhuruva/split1/basis_l23l23" with MAXITER=0 SEEDFITS=<fit json>. NOTE: group `paired_l23l23_val` does not exist in load_targets yet; add it (t.update(VALIDATION_L23L23)) before using this line.
