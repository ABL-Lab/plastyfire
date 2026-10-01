# Synapse-local eCB T drive for the Sjöström 2003 window (2026-09-30, test only, no model code changed)

**Setup.** The extracted npz hold only effcai, shaft_cai and vdcc. They have no local V, no ica_NMDA and no cacr (the pkl traces do have ica_NMDA), so I ran new single-cell sims. `local_t/sims.py` uses the diag_burst setup and keeps full per-synapse traces of v, shaft cai, cai_CR, effcai, ica_VDCC and ica_NMDA, for og and for ljp_VDCC=5 (s5).
- Coverage: all 24 subset24 pairs and **191 synapses** (143 basal, 48 apical; 32 at <60 µm, 95 at 60–150 µm, 64 at >150 µm).
- Cost: jobs 22072098 + 22072313, 2.4 CPU·h total, 1 CPU / 1.1G / 0:15 per pair, MaxRSS ≤ 0.89 GB.

**How I read the window.** The gate reads T at the pre arrival.
- The +10/+25 single pairings are therefore LTD-free for any post-driven drive.
- The real tests: −10/−25/−50 single and burst −120/−200 must give LTD; −100/−120/−200 single must not.
- Score (`local_t/analyze.py`, `check_hp.py`): the fraction of synapses where **one common parameter set and θ_Tg** puts every LTD row ≥ 2θ_Tg and every no-LTD row ≤ θ_Tg. It is also "strict": the own EPSP must not gate a second pre spike 50 ms later.

| Drive (local input) | Extra params | Window works, all / >150 µm |
|---|---|---|
| reference: `t_drive 2`, somatic APs (not allowed) | – | 1.00 / 1.00 |
| **D1: local-V event + own-pre veto** (θ_V 3 mV, w 15 ms) | θ_V (+ w fixed) | **0.99 / 0.97** |
| D1, θ_V 4 / 6 mV | | 0.96 / 0.89; 0.95 / 0.86 |
| D1 rest-free, v − v̄ with v̄ τ_b 300 ms (10 ms) | τ_b fixed | 0.99 / 0.97 (0.96 / 0.89) |
| A1: local-V event, no veto | θ_V | 0.30–0.63 strict; 0.99 if the EPSP is ignored |
| A2: V time above threshold / A3: linear low-pass of pos(V−θ) | θ_V | 0.13 / 0.05 |
| B1 / B2: shaft Ca event / two-stage | θ | 0.14 / 0.00 |
| spine Ca (cacr) event | θ | 0.15 |
| C1: ica_VDCC event (og or ljp5); with veto | θ_i | 0.63 (0.64); 0.88 (0.70 at >150 µm) |
| C2: ica_VDCC two-stage, og or ljp5 | θ_i | 0.03 |

**Why the others fail.** Every amplitude-based drive (V area, V linear, shaft or spine Ca, VDCC charge, with or without ljp5) scales with the bAP, which drops about 10× from <60 µm to >250 µm. So no common threshold works across synapses.
- ljp5 lifts VDCC Ca but not the spread, and leaves every event-based score unchanged.
- An amplitude-free count of local depolarising events works. It inherits the `t_drive 2` burst supralinearity (S += 1, τ_E1 100 ms, pos(S − θ_Te), then τ_T).
- The veto is needed because the unitary EPSP exceeds the local bAP peak at 20% of synapses (distal, clustered contacts; up to +6 mV).

**D1 matches `t_drive 2` row by row.** At θ_Te 0.85 and τ_T 40, median T is:
- LTD rows: 0.86 / 0.75 / 0.40 / 9.6 / 1.3, against 0.73 / 0.79 / 0.42 / 10.2 / 1.4 for `t_drive 2`;
- no-LTD rows: 0.115 / 0.070 / 0.009, against 0.12 / 0.07 / 0.01.

So the fitted θ_Te, τ_T and θ_Tg should carry over.

**D1 improves the train rows.** Median gate over the 5 arrivals, +10 / −10:

| Train | D1 | `t_drive 2` | Data |
|---|---|---|---|
| 10 Hz | 0.00 / 1.00 | 0.65 / 1.0 | 1.16 / 0.57 |
| 20 Hz | 0.00 / 1.00 | 0.8 / 1.0 | 1.31 / 0.65 |
| 50 Hz | 0.00 / 0.74 | 0.8 / 1.0 | 1.57 / 1.70 |

At 50 Hz −10, rho LTP still has to win, as it does now.

**Where D1 fails:** 2 of 191 synapses (basal 332 µm, apical 769 µm), where the bAP is < 2 mV and no local-V drive can see it. At θ_V ≥ 4 mV, 5 tuft synapses (685–724 µm, bAP 6–9 mV) also lose one burst AP to the afterdepolarisation.

**Biology (one line each)**
- Event count on local V: 2-AG synthesis (PLCβ/DAGLα) needs a fast local depolarisation plus VDCC nanodomain Ca. Each transient saturates the nanodomain, so the drive counts events rather than bAP amplitude.
  **Corrected 2026-09-30 (ECB_TRIGGER_LIT.md): this rationale is wrong.** eCB production is Ca-triggered (EGTA ≈ BAPTA block, so not a nanodomain; Ca uncaging alone suffices; subthreshold 250 ms steps replace spikes at L5-L5). D1 is a phenomenological proxy for a VDCC Ca event, not a mechanism.
- Supralinearity, pos(S − θ_Te): cooperative PLCβ activation, the same as in `t_drive 2`.
- Veto for 15 ms after an own pre arrival: while the own NMDARs are glutamate-bound, the depolarisation's Ca goes to the NMDAR/CaMKII (LTP) side, so pre-before-post makes no eCB. This one is phenomenological; w is fixed, not fitted.

**Recommendation: `t_drive 3` = D1.** It reads only the local segment voltage v at the synapse, a slow local baseline v̄ of it, and the synapse's own arrival times. It adds 0 free parameters (θ_V 3 mV, w 15 ms and τ_b 300 ms fixed), or 1 if θ_V is freed.

**Changes needed (not made, not submitted)**
1. `plastyfire/simulator_edges.py`: add a trace var `v_seg` that records `seg._ref_v` at the synapse location, as diag_burst does. Pass it via `--trace-vars cai_CR,shaft_cai,ica_VDCC,v_seg`.
2. `extract_v2.py`: add `vev`, the per-synapse upward crossing times of v − v̄ > θ_V, padded (n_syn, n_max). Store the times rather than the trace, because the windowed grid is too coarse for a 1 ms bAP.
3. `batch_v2._sig`: for `t_drive 3`, return per-synapse impulses at `vev`, dropping those within w of that synapse's `arr`. `model_v2.feature_T` then takes the existing `t_drive 2` branch (make it `in (2, 3)`).
4. Mod, for validation later: STATE v̄ (v̄' = (v − v̄)/τ_b). In BEFORE STEP, if v − v̄ crosses θ_V upward and t − t_lastpre ≥ w, then S += 1. Set t_lastpre = t in NET_RECEIVE.

**Re-prefire cost for paired_l5 on subset24.** The cost is the sims; the extra recording is negligible.
- Sjöström 2003, 11 protocols: measured 22057545, 27 CPU·h used (3 × 32 CPUs, 19–26 min, 143–154 GB).
- r50 bursts: 22068340, 6.6 CPU·h (3 × 16 CPUs, 12–15 min, 64 GB).
- Sjöström 2001, 10 protocols, and Markram, 7 protocols: estimated ≤ 24 and ≤ 17 CPU·h, scaled per protocol·pair from 22057545. These are shorter protocols, so this is an upper bound.
- **Total ≈ 60–75 CPU·h. This is over 30 CPU·h, so it must be logged in DECISIONS before submitting.** Extraction adds < 0.2 CPU·h (22068341).

**Caveats:** these are diag sims, one run per protocol, with no minis or background noise and a 0.1 ms delay. A 3 mV threshold assumes a quiet in-vitro membrane, so check it on the stochastic prefire traces once `v_seg` exists.

## Implementation (2026-09-30, `t_drive 3` in the offline pipeline; nothing committed, GPU path untested)
- **Files:** `plastyfire/simulator_edges.py` (`v_seg` trace, recorded only if named in `--trace-vars`; default output unchanged), `model_v2.py` (`theta_V` 3, `w_V` 15, `tau_b` 300 in DEFAULTS; `v_events`, `veto`; `feature_T` takes `t_drive in (2, 3)`), `extract_v2.py` (`vev` n_syn x n_max, NaN pad, from the full 0.025 ms trace; raw v dropped), `batch_v2.py` (`v_impulses` in `_sig`, `vev` loaded), `jax_v2.py` (per-lane impulse array for `t_drive 3`; not run, no jax on the login node).
- **Detector:** v_bar is an EMA from v[0]; one event per excursion, re-armed 2 mV below theta_V (`V_HYST`, as `local_t/analyze.HYST`).
- **Tests:** `test_t_drive3_veto`, `test_t_drive3_window`, `test_t_drive3_equals_t_drive2` added; `pytest tests/test_units.py -k "not jax"` gives 27 passed. GPU equivalence tests not run.
- **Pilot:** job 22073044 (2 pairs x 11 Sjostrom 2003 + 2 r50, 25 synapses, hash `delta-prefire-vseg`): 17:25, MaxRSS 102 GB (3.9 GB/sim) of 200G, CPU eff 49% (single wave, tail-bound), 526 s/sim (1.3x the same protocols without v_seg). Extract 22073365: 1:48, 4 CPUs, MaxRSS 94 GB (23.6 GB/worker, request was 90G, so use 30 GB/worker). Dirs `extracted/sj03_delta-prefire-vseg`, `extracted/sj03r50_delta-prefire-vseg`.
- **Events:** bAP events start ~3 ms before the soma spike time (window -8..+30 ms). Every synapse gets exactly 1 event per post AP (50/50 at 0.1 Hz, 75/75 and 250/250 for bursts). Post-before-pre +50: 50/50. Pre-before-post +25: only 17 of 50 survive (EPSP excursion not yet re-armed, or vetoed), which is the desired LTD-free row. Own EPSPs give 47 raw events per synapse in pre-only, all removed by the veto.
- **Noise:** 0 false events in 144,000 synapse-seconds of baseline outside [-8, +30] ms of a post AP and inside the veto (0 of 325 synapse-protocols). No background noise or minis exist in these sims, so the 3 mV threshold is untested against real membrane noise. On real npz, t_drive 3 gates match t_drive 2 (theta_Te 0.85, tau_T 40, theta_Tg 0.1): -25 0.52 vs 0.55, -50 0.26 vs 0.28, burst -120 1.00 vs 1.00, -200 0.78 vs 0.82, +25/+50/-120/-200 single 0.
- **Full run (24 pairs, NOT submitted, over 30 CPU h):** sim work 63 CPU h (sj03+r50 42.8, sj01 14.1, Markram 6.4), charged about 72 CPU h. Part A `sbatch glusynapse_v2/run_prefire_vseg_sj03_full.sh` (3 x 32 CPU, 195G, 0:45, ~48 CPU h; 9 sj03 + 2 r50 protocols, the nreps-15 bursts dropped) and part B `sbatch glusynapse_v2/run_prefire_vseg_sj01_markram_full.sh` (3 x 16 CPU, 110G, 0:45, ~24 CPU h; part B has no v_seg measurement, sized from prior logs x1.32). Both write `simulations/<pair>/<proto>/{bluecellulab_results,simulation_edges}_delta-prefire-vseg` (about 0.2 GB/sim, ~135 GB total). Extract afterwards to `extracted/{sj03,sj03r50,sj01,markram}_delta-prefire-vseg` with `extract_v2.py --param-hash delta-prefire-vseg --window`: about 2.4 CPU h, e.g. 8 workers, 240G, 0:30.

## Mod implementation (2026-09-30, `t_drive_GB = 3` in `mod/GluSynapse_v2.mod`)
- **Build:** `mod_build_td3/` (copy of `build/*.mod` plus the new mod, `nrnivmodl .`). `build/` is untouched; tests pick a build with `MECH_DIR=glusynapse_v2/mod_build_td3`.
- **GLOBALs:** `t_drive_GB` (0), `theta_Te_GB` 0.5, `Te_scale_GB` 1, `tau_E1_GB` 100, `theta_Tg_GB` 0, `theta_V_GB` 3, `w_V_GB` 15, `tau_b_GB` 300, `hyst_V_GB` 2. `dpre_min_GB` already existed.
- **States:** `S_GB' = -S/tau_E1`, `vbar_GB' = (v - vbar)/tau_b` (starts at v). `T' = -T/tau_T + T_drive()`, where T_drive is `pos(S - theta_Te)/Te_scale` for t_drive 3, else the old `pre_drive`.
- **Detector:** `WATCH (v - vbar > theta_V)` flag 6 and `WATCH (v - vbar < theta_V - hyst)` flag 7, set in flag 1 only if t_drive is 3. WATCH is root-found under CVODE and checked every step at fixed dt, unlike a BREAKPOINT check plus net_send.
- **Event:** flag 6 with `armed`: armed = 0, then `S += 1` if `t - tlast >= w_V`. Flag 7 re-arms. A vetoed crossing still uses up the arm, as in `model_v2.v_events` then `veto`.
- **Veto clock:** `tlast = t` in NET_RECEIVE, which is prespike + NetCon delay, as `batch_v2.v_impulses` (prespike + delay). Events in [arrival, arrival + w_V) are dropped.
- **Gate:** `dpre -= A_mglu*tanh(pos(T - theta_Tg))*(dpre - dpre_min)`. T >= 0, so with theta_Tg = 0 this is the old tanh(T); the old suite and a t_drive 0 dpre are bit-identical to `build/`.
- **Gate test** (`tests/test_t_drive3_gate.py`, 2 ms steps to -10 mV, 0.025 ms): post-before-pre dpre -0.0716, pre-before-post +0 (vetoed), burst of 5 at 100 Hz -0.768, veto edge (event at +15.6 ms) -0.210. Offline and NEURON agree to 1e-13 relative.
- **Residual:** only event-time discretisation. Both sides put the event on the first 0.025 ms sample after the crossing, and the Euler vbar differs from the offline EMA by about 1e-4 of the excursion, which cannot move a 60 mV event. Near-threshold excursions (within ~0.01 mV of 3 mV) could differ by one sample. CVODE gives the same dpre within 0.4% (event time is root-found).
- **Not covered:** cells with many segments and noise. The 2 mV hysteresis and the 3 mV threshold are as in the offline model, and GPU/jax stay unchanged.
