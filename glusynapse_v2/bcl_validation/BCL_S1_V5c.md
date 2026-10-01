# Prefire BCL validation of v5c on delta-split1 (fit S1_V5c_s6_39)

Live GluSynapseV5 (bluecellulab prefire runs on the delta-split1 emodel, split1 circuit configs) on every record behind the fit's 39 targets + 2 validation-only, against the offline CPU port of the gpu_v5_rho kernel. Scripts: run_valid_S1_V5c.sh (STAGE pilot / full / shard), run_valid_S1_V5c_cmp.sh (merges shard_*.jsonl).

## Verdict (2026-10-01)
**PASS 41/41.** 1244/1244 records (680 L5, 564 L2/3->L5). χ² over 39 fitted targets: BCL live 91.57, offline on the same records 87.84 = fit.
- L5->L5: live 73.59 vs offline 73.17 (30 targets); per record r 0.9991, ρ disagreements 2 / 5276 synapses.
- L2/3->L5: live 17.97 vs offline 14.66 (9 targets); per record r 0.9899, ρ disagreements 5 / 3650. Largest per-target gaps: letzkus 3AP 200 Hz +10 @proximal +0.020 and −10 @proximal −0.016, within tolerance; outliers are single records where one eCB step (dpre 0.29) differs at W ≈ θ_eCB.
- Model vs data (same in live and offline, so a rule limit on split1, not a simulator issue): low-frequency post-before-pre LTD too shallow (0.1 Hz −10/−25, 10 Hz −10, 20 Hz −10: model 0.85-0.94 vs data 0.57-0.69); 40/50 Hz −10 LTP too small (1.13/1.10 vs 1.51/1.70); S07 pair mglu_block 1.71 vs 2.13; Letzkus 3AP −10 distal (validation) 0.82 vs 1.42.
- Figure: results/S1_V5c_full.png/.pdf; tables results/S1_V5c_full_{targets,records}.csv. Interim (596 records): results/S1_V5c_partial.*.

## Runs
- pilot 22144694 (8:54, 5.24 GB at 12 workers, 80% CPU), compare 22144695.
- serial full 22145079 cancelled at 596 records (54 min) for speed; 7 shards x 32 CPU 22147501-07: 4:20-9:33, 13.1-15.0 GB of 28G (next: 19G, 0:15); merge + compare 22147508: 2:00, 20.2 GB of 26G.
