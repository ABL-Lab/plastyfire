# Does the doublet keep-list filter bias the burst/train targets? (agent F, 2026-10-02)

Status: list-level counts below are final. The per-bin Letzkus table is produced by Slurm job 22291773
(`/scratch/dhuruva/doublet_bias/count.py` -> `count.out`, 1 CPU 1G 0:15, was pending at write time). Fill section 2 from it.

## 1. What the filter is (run_stage_fit.sh l. 25-52)
- Lists: `/scratch/dhuruva/s2g0321/doublets/{l5,l23l5,l23l23}_keep_pairs.txt` (non-strict). They drop a pair only when the
  post cell has no exactly-one-AP c_post stimulus (`*_cpost_spikes.csv`, blank n_spikes), i.e. it doublets at threshold.
- The strict lists additionally drop pairs failing the prefire guardrail in any protocol. They are NOT what the fit uses.
- Pair level, not synapse level: a pair is dropped for all protocols, including bursts where doublets are tolerated.

## 2. Counts (pairs; synapses from ebner/pair_geometry_L23PC_L5TTPC*.csv)
| pathway | template | kept | excluded | strict-excluded | note |
|---|---|---|---|---|---|
| L5-L5 (Markram, Sjostrom, Letzkus-L5 none) | 24 | 22 | 2 (192879-186028, 207453-189325) | 3 | |
| L2/3->L5 (Letzkus 3AP, Ebner) | 120 | 102 | 18 (15%) | 78 (kept 42) | |
| L2/3-L2/3 (Zilberter train10, 5AP) | 120 | 120 | 0 | 5 | default filter is a no-op |

Letzkus pairs/synapses before the filter (pair_geometry): distal (rise>2.7 ms) 80 pairs / 501 syn; proximal 40 / 298 syn.
Per-bin before/after (rise split, path_mean bins 0-200/200-300/300-430/430-600/>600 um, geo split, all_protocols only):
see `/scratch/dhuruva/doublet_bias/count.out` (job 22291773). FILL: excluded pairs/syn per bin, mean path_mean,
n_apical_gt450/600, post_apical_max, epsp of excluded vs kept.

## 3. Verdict
- Zilberter train10 50 Hz and 5AP rows (L2/3-L2/3): NO bias from the default filter, 120/120 pairs kept.
  Caveat: the strict list would remove 5 pairs (e.g. 155875-145161, 161901-15591 deliver 392-393 of 400 spikes in
  train10, i.e. missed spikes, not doublets), so those runs are guardrail-marginal, not filtered.
- L5-L5: 22/24 kept, no meaningful bias.
- Letzkus 3AP 200 Hz: 18 of 120 pairs (15%) removed. 102 kept is ample in total (Letzkus n = 7-10 per cell), so the risk is
  not sample size overall but composition. Bias is plausible for the strongest-Ca sites: doublet-at-threshold cells are the
  large-soma, strongly coupled (bursting) TTPCs (Chagnac-Amitai 1990; Schaefer 2003), the ones most likely to have the biggest
  distal Ca spikes, which is the quantity the "3AP -10 distal LTP" target needs (DISTAL_LTP.md; the target is already dropped).
  Whether it is yes/no per distance bin is decided by section 2: bias = YES if the excluded pairs' share is clearly higher
  in the >430 um / apical-far bin than in proximal, or if the >600 um bin falls below ~8 pairs (only 11 pairs
  had path_mean>600 um before filtering, so each loss there is large).
- Filtering is inconsistent with the protocol: the exclusion criterion is a single-AP 2 ms c_post, but Letzkus bursts
  already tolerate late spikes (--allow-late-spikes). Doublet cells are valid for 3AP protocols (DOUBLETS_FITTING.md 1.4).

## 4. Recommendation
1. Do not change the L5-L5 and L2/3-L2/3 filter (no effect / negligible).
2. For Letzkus 3AP 200 Hz (proximal and distal) use a doublet-tolerant variant: keep all 120 pairs for those rows
   (new keep list `l23l5_keep_pairs_letzkus.txt` = all 120, passed only to the Letzkus targets), keep the strict filter for
   single-AP and Ebner rows, per-protocol `late_spikes_ok` (option (d) in DOUBLETS_FITTING.md), and report the doublet-cell
   stratum separately. c_post for those cells stays single-bAP Ca where available; bap_ref falls back to the pathway median
   for them (nmiss), which needs reporting.
3. Required pre-check: the 18 excluded pairs need valid prefire/npz for the 3AP protocols; if they lack them, a prefire
   run (heavy-run sizing, log first) is needed before this variant can be used.
4. Interim: if the section-2 table shows no excess loss in the distal bins, keep the current filter and only log the 18.
