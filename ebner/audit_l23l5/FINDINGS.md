# T28: why delta-cooker over-potentiates L2/3 -> L5 TTPC (S&H 2006, 5+5 @ 50 Hz, pre +10 ms)

Read-only audit, login node only, no Slurm. Scripts 01..07 and all csv here. Model = fit_results/delta-cooker.json, full rho ODE.
Reproduced: L23->L5 **1.65 +- 0.07** (119 pairs; 1.56 without the variance factor). Same model, same readout, L5->L5 Ebner set
(111 pairs, 07_readout_l5l5.py): **1.46 +- 0.05** (paper L5->L5 1.40). Paper L23->L5 1.06 (distal 0.86, proximal ~1.21).
So L5->L5 is fine and L23->L5 is +0.19 above it, +0.59 above data.

## Ranking (by evidence)

1. **All-or-none LTP at synapses the bAP reaches (structural, rule + readout).** 03_readout.py, 05_transitions_tm.py.
   Synapses with c_post above the resting floor (<2e-4; 313/794 = 39%): effcai peak / theta_p median 2.62, 100% have peak > theta_p
   and **P(rho0=0 -> 1) = 1.00** (L5->L5: 0.87). Each flipped weak synapse is x3.5 (gmax_p/gmax_d = 2, Use_p/Use_d = 0.87/0.5). Pair
   ratio by fraction of c_post-floor synapses: none floor (n=15) 2.23; <50% (24) 1.94; 50-99% (48) 1.79; **all floor (n=32) 0.95
   (median 0.69)**. The all-floor pairs are at the S&H distal value (0.86); every pair that has any synapse in bAP reach is 1.8-2.2.
   By mean path: <200 um 1.86, 200-300 1.82, 300-450 1.69, >450 1.21 (data: proximal ~1.21, distal 0.86), i.e. +0.5 at every distance
   and LTP still at 300-450 um. Floor fraction vs path (spearman) only 0.68: 47% of synapses at 200-300 um and 29% at 300-450 um still
   have bAP Ca. Letzkus 1 AP: same split, c_post-floor synapses P(1->0) 0.82 (LTD, right sign) vs non-floor P(0->1) 0.64 (pooled 1.04
   vs 0.72 data). The thresholds are in units of c_pre/c_post and were fit on L5->L5; there is no L2/3->L5-specific lever, and
   bAP-reached synapses saturate far above theta_p.
2. **rho0 assignment differs between the two pathways (edge artefact, L23->L5 only).** 01_dump_edges.py, 06_rho0_split_cf.py.
   L5->L5 edges: rho0 = 1 for conductance >= per-connection median (create_dhuruva_edges.py:192, only for L5:A/B combos at :75-80);
   mean conductance rho0=1 2.64 vs rho0=0 1.21 nS, 58% of pairs exactly 50/50. L23->L5: rho0_GB/Use_d/Use_p/gmax_d/gmax_p are the
   original Ecker-circuit values, untouched (identical to the original edges.h5 for all 798 synapses; our file only replaces theta_d/theta_p
   by -1); rho0 is independent of conductance (0.48 vs 0.48 nS), only 15% of pairs are 50/50 (8% all 0 or all 1). Random rho0 puts large synapses
   in the depressed state as often as small ones, so the x3.5 flips of the weak ones move the EPSP more. Linear EPSP proxy
   (sum gmax*Use, final states held fixed, NOT a re-simulation): actual rho0 1.48 +- 0.07, re-split by the L5->L5 rule 1.22 +- 0.05
   (-0.26). The proxy underestimates the full readout by ~0.1 (L5->L5 proxy 1.36 vs 1.46), so take -0.2 to -0.3 as the size.
3. **Variance factor in the readout**: (1 + min(cv2, 0.25)) = 1.05 on average (1.03 median), +0.09 (batch_v2.py:144). Small pairs
   (median 5 synapses) only.
4. **Release / NMDA / spine volume: not a cause (first order).** 02_synapse_table.py, 04_pre_vs_post_drive.py, regression in transcript.
   Per-synapse effcai/c_pre (peak): 3.34 vs 3.05 mean, median 2.60 vs 2.74. After controlling c_post-floor and apical, the L23 offset is
   +9% (log 0.09); floor alone is -34%. Share of effcai peak from the pre->post windows is 0.36-0.43 in both sets. Mean-field
   Tsodyks-Markram over the 5 pre spikes: U 0.50/D 671/F 17 (L23) vs U 0.38/D 352/F 8 (L5): release per spike [0.50 .30 .14 .07 .05] vs
   [0.37 .25 .17 .12 .09], decay-weighted sum 0.85 vs 0.83 c_pre-equivalents: higher Use is offset by slower recovery, no net
   depression difference. Nonfloor synapses do get 26% more peak/c_pre than L5->L5 (4.5 vs 3.4; pre-window-only peak > theta_p for 40% vs
   14%), a weak pointer at higher Use (0.50 vs 0.38) that the regression does not isolate. Needs a sim to test (see open points).

## Q1 Recipe and edges
- The circuit edges are NOT built from recipe_dhuruva.csv. All three recipes have identical L2/3 -> L5:A/B rows (line 28 etc.):
  gsyn 0.47, nrrp 1.5, dtc 1.74, u 0.50, d 671, f 17, gsynSRSF (NMDA) 0.70, synapse_type E2_L23PC_L5TTPC, spinevol 0.0433.
  They differ on L5:A -> L5:A (line 154): dhuruva gsyn 1.60 nrrp 2.9 u 0.50 d 672.7 f 15.9 NMDA 0.5625 spinevol 0.100 (E2);
  mod gsyn 1.94 nrrp 2.8 u 0.38 d 671 f 25.7 NMDA 0.55 spinevol 0.163 (E2_L5TTPC); andras = mod with d 365.6.
- Edges file (data/dhuruva_modified_edges.h5, used by data/dhuruva_delta_circuit_config.json): L5->L5 (892 syn): u 0.38, d median 352
  (mean 362), f median 8.1, gsyn 1.96, nrrp 2.77, volume_CR 0.162: matches recipe_andras (Ecker/Chindemi circuit), not recipe_dhuruva.
  L23->L5 (798 syn, 120 pairs, all syn_type_id 120, 81% apical): gsyn 0.48, nrrp 1.54 (median 1), u 0.50, d 671, f 17, volume_CR 0.044,
  conductance_scale_factor 0.70, u_hill 2.79, k_gsyn 2 and k_u 0.2 relations hold (Use_d = u^5 at rho0 1, Use_p = u^0.2 at rho0 0).
  theta_d/theta_p = -1 for every L23->L5 synapse (set for L5->L5). Per-synapse: edges_l23l5.csv, edges_l5l5.csv.
- Stage-by-stage consequence: EPSP 0.33 mV unitary (ebner/pair_geometry csv), 4x smaller gsyn, 3.7x smaller spine volume, but median
  c_pre 0.053 vs 0.057 (c_pre is a per-synapse calibration and absorbs gsyn/nrrp/volume).

## Q2 NMDA ratio
- Set by the edge property conductance_scale_factor, exactly 0.70 on all 798 L23->L5 synapses, 0.55 on all 892 L5->L5. bluecellulab
  sets gmax_NMDA = gmax0_AMPA * conductance_scale_factor (.venv/.../bluecellulab/synapse/synapse_types.py:284); the mod default
  (glusynapse_v2/mod/GluSynapse_v2.mod:151, 0.55 nS) is overridden. Recipe source column gsynSRSF (0.70 vs 0.55), carried into the
  circuit by the original build; our edge script does not touch it (MODIFIED_FIELDS, create_dhuruva_edges.py:69).
- 0.70 is the low end of the ~0.7-1 range for L2/3 -> L5, and already 27% above L5->L5. A larger ratio would raise NMDA Ca,
  but c_pre (measured per synapse at Use 1, rho0 1, simulator_edges.py:271-277) scales with gmax_NMDA too, so theta moves with it; net
  effect is second order (only the bAP-gated supralinearity). Not tested with a sim.

## Q3 Which part of the rule
- Thetas: theta_d = 1.03 c_pre + 1.11 c_post, theta_p = 1.47 c_pre + 1.29 c_post (batch_v2.py:195-199); theta_p/theta_d median 1.43, theta_d
  0.058-0.064, theta_p 0.083-0.088 (effcai units). c_pre median 0.053 (0.0268-0.225), c_post median 1e-4, mean 0.0056; c_post at floor:
  481/794 = 60.6% (S&H protocol table), 38.7% is the cache number of the geometry run (different selection). L5->L5: c_pre median 0.057, c_post floor 35%.
- Distal (floor) synapses, S&H: peak > theta_d 96%, peak > theta_p 79%, peak/theta_p 1.41 (peak/c_pre 2.1 vs theta_p/c_pre 1.47), time
  above theta_p 0.6 s of 152 s; outcome P(0->1) 0.34, P(1->0) 0.63: net LTD in state, as wanted.
- Is distal LTP from c_pre alone? Partly: peak/c_pre = 2.1 >= 1.47 for floor synapses. Lower bound of a pre-only train (effcai
  filter of spine Ca inside the 10 ms pre->post windows): exceeds theta_d at 24.5% (L5->L5 14%), theta_p at 3.5% (both). So
  c_pre-driven Ca alone reaches the LTD band at L23 distal synapses but almost never LTP; the extra comes from bAP-gated NMDA after the AP.
- rho0 split: 51.6% of synapses rho0=1 (L5->L5 52.8%); final state (S&H): rho_f has almost no memory of rho0 at floor synapses (P(0->1) 0.34 vs
  P(1->0) 0.63), full reset at non-floor ones (100% end potentiated). Mean rho_f 0.53.
- 02_synapse_table.py/syn_table_*.csv: per-synapse td, tp, peak, time above thresholds, rho_f.

## Q4 Other structure
- 81% of L23->L5 synapses are apical (path median 309 um) vs 23% for L5->L5: the floor fraction (61% vs 35%) is what makes L23 *less*
  potentiating at the synapse level; it is the non-floor third that does the damage (item 1).
- No evidence for a Use/D/F or NMDA-ratio driven difference (item 4). Spine volume, gsyn, nrrp are absorbed by c_pre.
- Full-protocol Ca model differences (bAP attenuation in the delta emodel, where c_post leaves the floor) is the lever most worth a
  direct look: compare the delta emodel bAP amplitude vs path distance with Stuart & Sakmann / Larkum data, since 47% of synapses at 200-300 um
  and 29% at 300-450 um still see bAP Ca here.

## Open points / suggested next steps (need the user)
1. Re-run the offline S&H + Letzkus eval with L23->L5 rho0 re-assigned by the per-connection median rule (new edges file or an
   offline rho0 override plus basis recomputation; the proxy says -0.2..-0.3).
2. Check c_post vs synapse path distance (bAP reach) in the delta emodel; the floor cut at ~300 um decides the proximal/distal split.
3. Only a simulation can test NMDA ratio / Use directly (one L23 pair, nscale 0.7 vs 0.55, Use 0.5 vs 0.38): not run.
