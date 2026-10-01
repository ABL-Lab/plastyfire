# v5: Chindemi 2022 plus one shared spine VDCC-Ca pool

Goal (user): "simpler plasticity that can explain all". Rule: only a00, a01, a10, a11, γd, γp are fitted. Every other
constant is a measurement (E) or an independent calibration (C). Inputs are synapse-local and Ca-based, and one rule
applies to every pathway.

Code (v4 files untouched):
- `gpu_v5_rho.py`: `GPUModelV5(GPUModelV4)` with a new shared-pool kernel. With `v5_mode` 0 it runs the v4 path unchanged.
- `fit_v5.py`: copy of `fit_v4n.py` (`--extra`, `--joint`), with the v4 A_mglu / A_NO slots removed from the DE vector
  (held at exactly 0), and `--check-v4`.
- `run_fit_v5.sh`: copy of `run_fit_v4.sh`. The v5 constants are in FILTERS, and the outputs are `results/v5_<TAG>`.

## 1. Equations

Per synapse, with c* = effcai (τ* 278.3 ms), c_pre and c_post as in Chindemi:

- θd = a00 c_pre + a01 c_post, θp = a10 c_pre + a11 c_post (ρ-γ = 1, Chindemi's form).
- **Shared pool** (the v4 C1 gate pool, unchanged): dV/dt = −V/τE1 + (−I_VDCC)/i_scale. It is fed only by the
  synapse's own spine-head VDCC Ca current.
- **Potentiation gate:** pot = Θ(c* − θp) Θ(V − θ_V), dep = Θ(c* − θd). As in C1, a crossing of θp that fails the
  gate counts as depression.
- **ρ:** τ_ind dρ/dt = −ρ(1−ρ)(ρ*−ρ) + γp(1−ρ) pot − γd ρ dep (1−pot).
- **eCB-LTD, on the same pool with the same threshold:** at each own pre-spike arrival (count c) with V > θ_V,
  d ← d_min + (d − d_min)(1 − A_eCB)^c. With A_eCB = 1, one gated pre spike sets d = d_min.
- **Readout:** U_SE = min(1, U0(1+d)), unchanged.
  - AM251 / ifenprodil (mglu_block): no eCB step.
  - NMDAR block: ρ frozen and d = 0 (eCB-LTD needs preNMDAR, Sjöström 2003).
  - NO block: same as control, because there is no NO arm.

Reading: V is the occupancy of a sensor tethered at the spine VDCCs (nanodomain, CA_DECODE §1a), with τE1 as its
memory. When it is active it does two things:
1. It licenses LTP when c* is also high.
2. It drives Ca-triggered eCB (2-AG) synthesis. Retrograde CB1 action then needs the synapse's own pre spike
   (preNMDAR coincidence).

The order "post Ca before pre" comes from the pool's state at the moment the pre spike arrives. There is no K
crossing, no (1−b) weight, no T trace and no NO trace.

What v4 C1 had and v5 removes:
- eCB chain: K, the (1−b) event weighting, S, T, τ_T, θTe, θTg (replaced by the pool and θ_V), and A_mglu (→ A_eCB = 1).
- NO chain: N, τ_NO, θNOi, Z, τ_Z, θZ, A_NO.
- ρ-γ.

## 2. Parameters

| param | status | value | anchor (preparation) |
|---|---|---|---|
| a00, a01, a10, a11 | fitted (Chindemi) | DE 0–5 | Chindemi 2022, 10.1038/s41467-022-30214-w |
| γd, γp | fitted (Chindemi) | 20–250, 100–600 | GAMMA_LIT.md |
| θ_VDCC (theta_V) | **fitted, flagged** | DE 0–50 | Calibration C-V against Kampa 2006 was null (22135045: og-delta has no 200 Hz supralinearity; geometry sets the gate). PARAM_ANCHORS option (ii). In v5 the one threshold serves both the gate and eCB, so v5 has one flagged constant where v4 had θ_VDCC plus θTe/θTg/τ_T. |
| d_min | E | −0.29 | CB1-agonist (ACEA) LTD ceiling 71 ± 3 %, occluding tLTD. Rat V1 L5→L5 pairs, P12–21, 32–34 °C (Sjöström 2003, 10.1016/S0896-6273(03)00476-8; PARAM_ANCHORS) |
| τE1 | E (bracket) | 100 ms | eCB drive memory: synthesis/release < 50 ms at 37 °C and 75–190 ms at 22 °C (rat CA1, Heinbockel 2005, 10.1523/JNEUROSCI.2078-05.2005); L4→L2/3 eCB-LTD coincidence scale ≈ 125 ms (rat S1, Bender 2006, 10.1523/JNEUROSCI.0176-06.2006). Kept at 100 ms so that the M1 limit is v4 C1 exactly. |
| A_eCB | structural limit | 1 | tLTD is frequency-independent (0.1–20 Hz ≈ 0.7) and occluded by AEA, so the protocols saturate (Sjöström 2003 Fig 8). PARAM_ANCHORS §3 option A. |
| i_scale | units | 1e-5 nA | only θ_V·i_scale enters |
| ρ-γ | structural | 1 | Chindemi / GB form (C_REF drops out) |
| τ_ind, ρ*, τ* | Chindemi | 70 s, 0.5, 278.3 ms | Chindemi |

Count: **7 fitted** (6 Chindemi + θ_V flagged), **2 anchored** (d_min, τE1) and 1 structural limit (A_eCB = 1). v4 C1
has 18 fitted.

Expected costs, judged from the C1Ajn_s5 residuals:
- **No NO arm.** Step pair: control 1.62 ± 0.07 and NO block 1.36 ± 0.07 get the same prediction, which gives
  χ² ≥ ≈ 6.9 on those two (v4: 3.5).
- **eCB is also active during 20–50 Hz +10 trains**, because the pool is above θ_V there. ρ-LTP must reach about 1.9
  so that 50 Hz +10 control ≈ 1.35 against ifenprodil 1.59 ± 0.19.
- **Expected gains.** 20 Hz −10 (v4 0.92 vs 0.65, AM251 1.20 vs 1.02) and 10 Hz −10 (0.78 vs 0.57) should improve,
  because the (1−b) weight that suppressed eCB in fast trains is gone.

## 3. Validation and fits

Joint L5 + L2/3→L5 fits, 39 targets (30 + 9). Setup:
- DROPT: Letzkus 3AP −10 distal (control and nmdar_block), as in MINIMAL_LADDER.md.
- Default L5 dirs and groups, no sj04.
- FITGAMMA=1.

Resources: 2g MIG, 1 CPU, 50G, def-emuller.

| job | tag | what | time | result |
|---|---|---|---|---|
| 22135813 | V5chk | `--check-v4`, MAXITER 0, seed C1Ajn_s5 (a's, γ's, θ_V; ρ-γ 1, A_mglu = A_NO = 0) | 0:15 | pending |
| 22135814 | V5_s5 | v5 fit, seeded C1Ajn_s5, MAXITER 300, afterok 22135813 | 0:45 | pending |
| 22135815 | V5_s6 | v5 fit, unseeded (pop[0] a's from delta-cooker), MAXITER 300, afterok 22135813 | 0:45 | pending |

The check passes only if the log line `CHECK-V4 OK` is present. Otherwise the job exits non-zero and the fits never
start. For each model, l5 and l23, it asserts to 1e-9 in χ²:
- M1: the v5 kernel with A_eCB 0 equals the v4 vamp kernel (C1 gate, eCB/NO off).
- M0: the v5 kernel with A_eCB 0 and θ_V 0 equals the v4 rates kernel (vamp 0) and the v4 vamp kernel at θ_V 0.

The log also prints the v4 M0 total next to `results/v4_LM0chk.json` (ladder 22135018, same a's and γ's), and the
seed's χ² under the full v5 rule (`results/v5_V5chk.json`).

How to compare (n 39, AIC = χ² + 2k, BIC = χ² + k ln 39):
- v5 has k = 7.
- It beats C1Ajn_s5 (χ² 51.98, k 18, AIC 87.98, BIC 117.92) if χ² < 73.98 (AIC) or < 92.27 (BIC).
- The same bars apply against LM1 (also k 7). It must beat LM1 on χ² directly, because the eCB arm adds no free
  parameter.

## 4. Variants v5b / v5c (added after MECH_NECESSITY.md §E)

The pool test in MECH_NECESSITY §E (job 22135544) found that the (1−b) event weight is what keeps eCB off in fast
+10 trains: the weighted 50 Hz +10 pool is 18× lower (3.9 → 0.22 θ_VDCC), while the −10 protocols keep most of
theirs. The plain v5 pool has no such weight (the §2 risk). Both variants are switched by `v5_mode` 2, and the
default v5 (`v5_mode` 1) is unchanged, bit for bit.

- **v5b (`v5_mode` 2):** the eCB trigger reads W, a copy of V with the same τE1 and i_scale, whose input is weighted:
  dW/dt = −W/τE1 + (1 − min(b, 1))(−I_VDCC)/i_scale.
  - An own pre spike with W > θ_V sets d ← d_min (A_eCB 1).
  - The gate still reads the unweighted V with the same θ_V.
  - k = 7, the same as v5.
- **v5c:** v5b with a separate θ_eCB on W (`--free-filters theta_V,theta_eCB`, log box 0.05–50, seeded at θ_V
  5.50). θ_eCB is **fitted and flagged** (no anchor). k = 8.

**Origin and anchor of b.** b is the synapse's own glutamate-bound NMDAR state, exactly as v4 t_drive 4 defines it
(`model_v2.glu_weight`):
- +1 at each own pre arrival (prespike + NetCon delay);
- decays with τ = tau_d_NMDA = 70 ms;
- capped at 1;
- arrivals at or before the time count.

70 ms is the NMDAR decay constant of the synapse mod file. It is not a fitted v4 value; it is anchored through
Chindemi's synapse calibration (PARAM_ANCHORS, hidden constants). There is no K level: the weight multiplies the
continuous current instead of K-crossing events.

One numerical difference from v4: the kernel applies each arrival at its grid sample (`arrival_index`, first sample at
or after the arrival, 0.25 ms grid in the spike windows), whereas `glu_weight` uses exact arrival times.

| job | tag | variant | free (k) | seed | after |
|---|---|---|---|---|---|
| 22135953 | V5b_s5 | v5b | 7 | C1Ajn_s5 | afterok 22135813 |
| 22135954 | V5c_s5 | v5c | 8 | C1Ajn_s5 (θ_eCB start 5.50) | afterok 22135813 |

Both jobs also run `--check-v4` on the weighted kernel: with A_eCB 0, the M1 and M0 limits must equal the v4 kernels,
asserted before the DE. Sizing is 2g, 50G, 0:45, as for the v5 fits.

## 5. v5-E2 (v5_mode 3) and the 3-pathway fits

Input: L23L23_DIAG.md E2 (run 2, 22135934–36, C1Ajn_s5 params, no refit). In v4 terms, own VDCC events feed the eCB
drive S1 only while c_VDCC ≤ θ_V. That gives L2/3→L2/3 / L5 / L2/3→L5 = 98.06 / 50.13 / 5.08 = 153.3 (base 162.2).
It is the only uniform variant that helps all three pathways.

**v5 form (no new parameter).** v5-E2 is v5b, except that the current of an own VDCC excursion feeds W only if the
excursion *began* while the pool was unprimed:

dW/dt = −W/τE1 + g_E (1 − min(b, 1))(−I_VDCC)/i_scale

- **The latch g_E.** g_E = [V ≤ θ_V], evaluated at the onset of each own VDCC excursion and held until the next onset.
  An onset is a K_ca crossing (the extracted cev events), exactly as v4 t_drive 4 defines its events.
- **Trigger.** At an own pre arrival with W > θ_V, d ← d_min.
- **Why this is the faithful equivalent.** In E2 an event's whole S1 impulse counts or not, depending on c_VDCC at the
  event. Here an excursion's whole charge counts or not, depending on V at its onset.
- **Effect.** A single bAP into a quiet pool feeds W in full, so single −10 eCB is kept. In a train, later bAPs find
  V > θ_V and add nothing; W then decays with τE1, so train −10 eCB is cut.
- **Rejected alternatives:**
  - Gating every step by [V ≤ θ_V] caps W at about θ_V, so the trigger W > θ_V would never fire.
  - Using the upward θ_V crossing of V would be a different rule (eCB on priming, not from unprimed events).

K_ca is v4's structural detection floor (PARAM_ANCHORS C-K). It is not a new parameter.

As PARAM_ANCHORS and L23L23_DIAG note, the biological case is weak: 2-AG synthesis rises with Ca. E2 is a
phenomenological "eCB from unprimed bAPs" gate. k = 7 (6 Chindemi + θ_V flagged).

**Fits: 3 pathways, 54 targets (30 L5 + 9 L2/3→L5 + 15 Zilberter L2/3→L2/3).**
- `EXTRA l23l23:paired_l23l23:glusynapse_v2/extracted/zilberter_l23l23_delta-prefire-vseg-rs:basis_results_edges_zilberter_l23l23_delta_rs`, plus the same DROPT as before.
- 3g MIG, 1 CPU, 79G (63.0 GB measured in 22135207 / 22134473, +25%).
- Fits get 0:45. Basis: 3-pathway M0 fit 22135207 ran 11:41 at 3.0 s per generation; v5 has 56 members × 300
  generations ≈ 18 min, plus about 2.5 min load, +50% ≈ 31 min.
- v5b was the arm chosen: none of v5 / v5b / v5c had finished at submission, so the default rule applied.

| job | tag | variant | seed | time | after |
|---|---|---|---|---|---|
| 22136256 | V5e2zchk | E2, MAXITER 0, `--check-v4` (M1/M0 limits on 3 pathways) + E2 score at the seed | C1Ajn_s5 | 0:15 | — |
| 22136257 | V5bz_s5 | v5b, 3 pathways | C1Ajn_s5 | 0:45 | afterok 22136256 |
| 22136258 | V5e2z_s5 | E2, 3 pathways | C1Ajn_s5 | 0:45 | afterok 22136256 |
| 22136259 | V5e2z_s6 | E2, 3 pathways | unseeded | 0:45 | afterok 22136256 |

How to compare (n 54, BIC = χ² + k ln 54, ln 54 = 3.989), against C1Ajn_s5 on 3 pathways: 162.24, k 18, AIC 198.24,
BIC 234.04. With k 7, a v5 variant wins if χ² < 184.24 (AIC) or < 206.12 (BIC).

V5chk 22135813 PASSED: M1 and M0 diff 0.0 on both pathways, and v4 M0 equals v4_LM0chk.json exactly
(1865.645030116). It ran 2:44, 31.6 GB MaxRSS, 98% CPU. Under the plain v5 rule the seed scores χ² 339.9
(a starting point only, not a fit).

## 6. Results (39 targets) and follow-ups

| fit | variant | χ² total | L5 / L2/3→L5 | k | AIC | BIC | job (seff) |
|---|---|---|---|---|---|---|---|
| V5_s5 | v5 | 169.86 | | 7 | 183.9 | 195.5 | 22135814 (12:07, 32.5 GB) |
| V5_s6 | v5, unseeded | 165.98 | | 7 | 180.0 | 191.6 | 22135815 (7:14, 27.2 GB) |
| V5b_s5 | v5b | 123.90 | 112.7 / 11.2 | 7 | 137.9 | 149.5 | 22135953 (13:03, 27.2 GB) |
| **V5c_s5** | **v5c** | **85.13** | **83.87 / 1.26** | **8** | **101.1** | **114.4** | 22135954 (13:57, 39.1 GB) |
| ref C1Ajn_s5 | v4 C1 | 51.98 | 46.89 / 5.08 | 18 | 87.98 | 117.92 | |
| ref M2a | v4 anchored eCB | 102.9 | | 8 | | | |

V5c_s5 fitted values: θ_V 0.244, θ_eCB 6.86, γd 167, γp 258. It beats C1Ajn_s5 on BIC with 8 parameters instead of 18.

θ_V = 0.244 lies 0.5 % above the box floor. The floor is 0, which turns the gate off (`thV > 0` in the kernel), so
the gate is close to always open. V5c_nogate_s5 tests this: θ_V = 0 and k = 7, with θ_eCB still separate. If its χ²
is within about 2 of V5c_s5, drop the gate.

Note: E2 latches "unprimed" against θ_V. At θ_V ≈ 0.24 almost every pool is primed, so in v5c-E2 the fit has to move
θ_V up for E2 to act.

Follow-ups (2g 50G 0:30 for 39 targets; 3g 79G 0:45 for 54 targets):

| job | tag | what | seed |
|---|---|---|---|
| 22136662 | V5c_s6 | v5c, basin check | unseeded |
| 22136663 | V5c_s7 | v5c, MAXITER 600 (0:45), convergence check | V5c_s5 |
| 22136664 | V5cz_s5 | v5c, 3 pathways (54 targets) | V5c_s5 + C1Ajn_s5 |
| 22136665 | V5ce2_s5 | v5c-E2 (`v5_mode` 3 + free θ_eCB), k 8 | V5c_s5 |
| 22136666 | V5ce2z_s5 | v5c-E2, 3 pathways | V5c_s5 |
| 22136667 | V5c_nogate_s5 | v5c with θ_V = 0 (gate off), k 7 | V5c_s5 |

V5e2zchk 22136256 PASSED: `CHECK-V4 OK` on 3 pathways (M1 371.570167348 and M0 2501.145914843 equal), 5:17, 54.0 GB.

## 7. Preparation-matched Letzkus split (GEOM_L23 = ebner/pair_geometry_L23PC_L5TTPC_geo.csv)

In this split a pair is distal when at least half of its synapses are apical beyond 430 µm. Rescoring without refit:
L2/3→L5 χ² rises from 5.08 to 26.1 for C1Ajn_s5 and from 16.4 to 28.3 for V5_s6.

All three jobs are v5c on 39 targets (2g, 50G).

| job | tag | what | seed | time |
|---|---|---|---|---|
| 22136762 | V5c_georescore | MAXITER 0 rescore of V5c_s5 on the geo split | V5c_s5 | 0:15 |
| 22136763 | V5c_geo_s5 | v5c fit, geo split | V5c_s5 + C1Ajn_s5 | 0:30 |
| 22136764 | V5c_geo_s6 | v5c fit, geo split | unseeded | 0:30 |
