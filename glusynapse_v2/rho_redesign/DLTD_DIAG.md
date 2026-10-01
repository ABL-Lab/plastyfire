# Why v4_C1Ajs_s5 cannot make Sjostrom 2004 dLTD

Diagnosis only; nothing in the model or fitters was changed.

- Script: `rho_redesign/diag_dltd.py`. It replays the gpu_v4_rho kernel (vamp 1, t_drive 4) on the CPU, one record at a time.
- Job 22132795: 0:12, 910 MB, 100% CPU. Log: `logs/diag_dltd_22132795.out`.
- Data: 24 pairs and 191 synapses per protocol.
  - sj04 = `sjostrom04_dltd_step250ms` (sj04_delta-prefire-vca).
  - dt-10 = `sjostrom_0.1hz_dt-10ms` (ebner_delta-prefire-vca). This is the protocol that works.

## 1. The replay reproduces the fit

| protocol | condition | fit csv | replay | data |
|---|---|---|---|---|
| sj04 | control | 0.8872 | 0.8872 | 0.69 ± 0.04 |
| sj04 | mglu_block | 0.8872 | 0.8872 | 1.06 ± 0.05 |
| dt-10 | control | 0.8343 | 0.8343 | 0.69 ± 0.07 |
| dt-10 | mglu_block | 1.0207 | 1.0207 | 1.01 ± 0.04 |

The replay matches the fit to 4 decimals, so the per-synapse numbers below can be trusted.

## 2. Pre side: the eCB gate

The table gives medians over the 191 synapses, with fractions of synapses where marked. K = 2.469e-8 nA.

| quantity | sj04 | dt-10 |
|---|---|---|
| VDCC events per pairing before the pre arrival | **0** (0% of syn) | 1 (94% of syn) |
| VDCC events per pairing after the arrival | 1 (99.5% of syn), at +1.5 ms | 0 (26% have any) |
| event weight 1 − b | **0.06** | 1.00 |
| peak −I_VDCC / K | **4.8** | 1190 |
| −I_VDCC / K, 50–250 ms after arrival | 3.4 (never below 0.5 K, so no re-arm) | 0.19 |
| S max | 0.17 (θTe 0.125) | 1.00 |
| T max | 0.33 | 13.6 |
| max over arrivals of T − θTg (θTg 7.21) | **−7.21** (0% > 0) | +0.42 (70% > 0) |
| max T − θTg in the 400 ms after an arrival | −6.9 (4% > 0) | +6.4 |
| d0 (control dpre) | **0** | −0.241 (= dpre_min) |
| d1 (mglu_block dpre) | 0 | 0 |

**The "94% before / 70% xg > 0" figures come from the dt-10 summary line (log line 545), not from sj04.** The sj04 line (log line 273) has 0% of synapses with an event before the arrival and 0% with xg > 0. So the no-events-before-arrival idea holds for sj04, and dt-10 is the control that shows the gate opening.

**d1 = 0 in both protocols by construction.** d1 is the mglu_block accumulator: the A_mglu step is never applied to it, and the NO branch is inactive here. The eCB effect lives only in d0:
- dt-10: d0 is driven to dpre_min, so control 0.83 vs AM251 1.02;
- sj04: d0 = 0, so control = mglu_block.

Why sj04 gets no eCB:
1. **Ordering.** The yaml puts the pre spike at step onset. The step's single K-crossing comes 1.5 ms after the arrival, but the A_mglu step reads T at the arrival, when T = 0. The previous pairing was 10 s earlier, so nothing carries over.
2. **Weight.** That crossing falls inside the synapse's own glutamate window, so 1 − b ≈ 0.06. This is the Nevian order built into t_drive 4.
3. **One event per step.** −I_VDCC stays at about 3.4 K for the whole 250 ms, so it never drops below 0.5 K to re-arm. A 250 ms depolarisation counts as one weak event, while each bAP counts as one full event.

K is not the problem: the step already crosses it at 99.5% of synapses.

## 3. Post side: the ρ rule (sj04)

The table gives c* (effcai) time totals over 50 pairings, in ms, as medians with means in brackets.

| quantity | sj04 | dt-10 |
|---|---|---|
| synapses with peak c* > θd | **136 / 191** | 55 / 191 |
| synapses with peak c* > θp | 102 / 191 | — |
| time c* > θd | 553 [1143] | 0 [192] |
| time c* > θp (raw) | 56 [626] | 0 [47] |
| depression actually applied (dep and not pot) | 553 [1140] | 0 [160] |
| potentiation actually applied (after the V gate) | **0 [3]** | 0 [32] |
| ρ flips 1→0 / 0→1 (of 99 / 92) | **14 / 0** | 1 / 2 |

- The step plus NMDA (glutamate bound while the cell sits depolarised for 250 ms) lifts c* above θd at most synapses, higher than a bAP 10 ms before the pre spike does.
- The VDCC-amplitude gate (θV = 6.7) blocks potentiation, because the step's VDCC current is about 250× below a bAP's. Under C1 (vamp_mode 1) a blocked potentiation step counts as depression, so nearly all of the time above θp also depresses.
- The 14 flipped synapses spent 1702 ms in depression on average (minimum 1206 ms). The 85 that stayed at ρ = 1 spent 402 ms (maximum 1224 ms). So the flip threshold is about 1.2 s.
- These 14 flips alone give 0.887. AM251 does not touch them, which contradicts the data (1.06).

## 4. What Sjostrom 2004 says dLTD requires

| requirement | evidence (page, figure) |
|---|---|
| CB1 receptors | AM251 "also blocked dLTD (P = 0.81)" (p3340 text, Fig 3B) |
| presynaptic NR2B NMDARs | ifenprodil "completely abolished dLTD" (p3340, Fig 3B). It does not act on postsynaptic receptors at this age. |
| no postsynaptic LTD component | both blockers give about 1.0–1.1 (Fig 3B); expression is presynaptic by CV and STD (Fig 2, p3340) |
| depolarisation threshold | a single connection's EPSP is not enough (p3339, Fig 1B circles); the step goes to −52 mV (Methods p3339) |
| timing | pre 25–125 ms *before* the depolarisation: no LTD. Pre 50–400 ms *after* it, or coincident: LTD (p3341 text, Fig 4 legend) |
| what the target pools | Fig 1B triangles = "all timings resulting in LTD in Fig. 4 pooled" (p3339). So 0.69 pools coincident and post-before-pre pairings. |
| Ca source | not tested: no BAPTA, VDCC blocker, mGluR antagonist or APV in this paper. Low-threshold dendritic Ca channels are only suggested (p3341). |

**Yaml timing recheck.** Sjostrom2004_L5TTPC_L5TTPC.yaml puts the pre spike at step onset (dt 0, "ASSUMED"). That is the least favourable of the pooled timings. The paper's own model (p3341) has eCB released by the depolarisation and still present when the pre spike arrives up to 400 ms after it ends. So the current target compares a pooled post-before-pre-weighted mean against one onset-coincident simulation.

## 5. Candidate fixes

All candidates are synapse-local and uniform across pathways. None is implemented. They are ranked by how well the pharmacology supports them and by risk.

| # | fix | sj04 effect | effect on fitted protocols | support |
|---|---|---|---|---|
| 1 | **Protocol / target, not the rule.** Add sj04 ids with the pre spike inside or after the step (e.g. +100 ms after onset, and 50 / 200 ms after the step ends) and score the pooled 0.69 against them, not against dt 0 alone. | Estimated from the fitted τE1 = 100 and τT = 26, with the onset event at weight 1 because b = 0 when the pre spike comes later. T reaches θTg from about 5 to about 110 ms after the onset crossing, so pre at +100 ms gives full eCB LTD, like dt-10. Pre more than 110 ms after onset (any post-step timing) still gives T ≈ 0. | none (new ids only) | Fig 1B pooling, Fig 4 |
| 2 | **Duration-coded VDCC events.** While −I_VDCC stays above K, re-arm on a refractory time τr (about 20–30 ms) as well as below 0.5 K. | A bAP excursion (a few ms) is still exactly one event, so the per-bAP unit is kept (DECISIONS 2026-09-30). The 250 ms step gives about 10 events and S of about 4 at step end. S then stays above θTe for about 350 ms, giving the 50–400 ms post-step window of Fig 4. | Sj07 step200ms_pair gains eCB LTD (risk to 1.62; the mglu gap 2.11 vs 1.68 already shows eCB acting there). Bursts change only if −I_VDCC stays above 0.5 K between spikes (dt-10 falls to 0.19 K, so single APs are unchanged). The 0.1 Hz −120/−200 runs (single AP, one excursion) are unchanged. This needs a refit, plus a check of excursion widths in 20–50 Hz bursts. | p3341: the depolarisation releases eCB via dendritic Ca channels; CB1 blocked |
| 3 | **Post side: C2 instead of C1** (vamp_mode 2: a blocked potentiation step is neutral). This already exists in the kernel; it is a refit, not code. | It removes the "blocked pot → dep" time; for the flipped synapses that is about 775 of 1702 ms, which puts them below the ~1.2 s flip threshold, so most of the 14 flips should go. The prediction then moves towards 1.0 for both arms (AM251 target 1.06). | Every L5 LTD target with drug data is mglu-sensitive (dt-10, 10 Hz −10, 20 Hz −10 all about 1.0 under mglu_block), so post ρ-depression is barely needed there. But C2 changed the earlier ROUND2 fits; refit and compare χ² with C1Ajs_s5. | AM251 and ifenprodil abolish all dLTD, so there is no postsynaptic LTD (Fig 3B, Fig 2) |
| 4 | **Gate integrated over the own glutamate window** (A_mglu rate ∝ b · tanh(pos(T − θTg)) rather than a step at the arrival). | It helps the dt 0 coincident case only together with #2: the onset event still has 1 − b ≈ 0.06, so later step events are needed. | Risk for pre-before-post (+10 ms) LTP and no-LTD targets. The 1 − b weighting limits it, but it shifts the timing window. | the presynaptic NMDAR + CB1 coincidence model (p3341) |
| 5 | Lower K, or K relative to the own bAP peak | none: the step already crosses K at 99.5% of synapses. A relative K (the bAP is about 1200 K, the step 4.8 K) would *remove* the step event. | changes every protocol's event set | **reject** |
| 6 | Graded VDCC charge or Ca drive replacing events | could work in principle | already ruled out (DECISIONS 2026-09-30: the per-bAP unit must hold to ±4%) | low |

**Supported by the pharmacology: #3 together with #1 or #2.** Fig 3B demands two things:
1. the eCB / A_mglu pathway must produce all of the depression (#1 or #2);
2. the AM251-insensitive ρ depression must disappear (#3).

Fixing only the pre side leaves the 14 flips, so AM251 would stay at about 0.89. Fixing only the post side gives about 1.0 for both arms.

**Suggested order:**
1. Re-simulate sj04 at the Fig 4 timings (#1); this needs no rule change.
2. Refit with vamp_mode 2 (#3).
3. Only if pre timings more than 110 ms after onset must give LTD, prototype #2 against excursion widths in the burst protocols.

## 6. Timing re-simulation

Proposal only. targets.py is unchanged.

**Protocols.** Four new ids in Sjostrom2004_L5TTPC_L5TTPC.yaml (same 250 ms step, threshold minus 20 pA, 50 pairings at 0.1 Hz). The pre time is `dt` with `dt_ref: stim`, where pre = step onset minus dt, so dt < 0 means the pre spike comes after the onset.

| id | dt (ms) | pre spike |
|---|---|---|
| sjostrom04_dltd_step250ms (existing) | 0 | at onset (coincident) |
| ..._pre+100 | -100 | 100 ms after onset, inside the step |
| ..._post+50 | -300 | 50 ms after step end |
| ..._post+200 | -450 | 200 ms after step end |
| ..._pre-100 | +100 | 100 ms before onset: no-LTD control, data ~1.0 |

**Jobs** (24 subset24 pairs, og-delta, hash delta-prefire-vseg): simwriter 22133660 (2G, 0:15, 8 CPU), prefire 22133663 (70G, 1:15, 12 workers), extract 22133664 (91G, 0:15, 4 workers), chained with afterok. Output: extracted/sj04_delta-prefire-vca/<pair>__<protocol>.npz, 96 new files.

**Scoring proposal.**
- LTD arm: the pooled 0.69 is scored against the unweighted mean of the 4 LTD timings (existing id at dt 0, pre+100, post+50, post+200), averaged over pairs. Fig 1B pools all LTD timings of Fig 4, and the paper does not give per-timing n, so equal weights are the neutral choice. Also report each timing separately to see which ones fail.
- Control: pre-100 is scored against ~1.0 (no LTD). It is not part of the pooled mean.
- AM251 arm: the same four LTD timings with the eCB / CB1 pathway blocked, mean scored against the AM251 target 1.06 (Fig 3B). The control needs no separate AM251 arm (already ~1.0). Which switch implements AM251 in the replay is the one used for the existing sj04 AM251 prediction; reuse it unchanged for each new id.
