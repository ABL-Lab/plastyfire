# A synapse-local Ca drive for the eCB-LTD trace T (2026-09-30, offline numpy on local_t/out, no model code changed)

## Recommendation
- **Drive T from spine VDCC Ca-influx events, and replace the 15 ms hand veto with the synapse's own glutamate-bound state.**
  - S += (1 − b(t)) at each upward crossing of −ica_VDCC through K. That is C1, i.e. the n → ∞ limit of candidate 1's saturating sensor.
  - b = the own B_NMDA state (it jumps at the own release and decays with the existing tau_d_NMDA of 70 ms), normalised so that one release gives 1 and capped at 1.
  - The rest of the rule is unchanged: τ_E1 100, pos(S − θ_Te), τ_T, and the gate at the own pre spike.
- **Scores: 0.90 of synapses at ljp 0 and 0.93 at ljp 25**, with every synapse under 250 µm at ≥ 0.90.
  - Connections: majority 0.96 / 1.00, mean T 0.92 / 0.88.
  - These are the same numbers as the hand-veto C1 on every scored row.
- **New free parameters: 1 (K, the influx threshold).** The veto window w is gone, so the net change against C1 + hand veto is −1.
  - `ljp_VDCC = 25` is optional for this drive: it adds only +0.03. Keeping ljp 0 avoids the effcai refit.
- **Main negative result: no linear or integrated Ca drive can pass the current rule, even when it is perfectly normalised per synapse.** A passing drive must deliver the same unit per bAP to within ±4%, and only per-event saturation does that.
  - With the own bAP charge as exact normaliser (NRM), the peak of S still spreads 0.82–0.97 (p10–p90, 60–250 µm). This is because Ca entry lasts several ms and its time course differs across synapses (only 51–85% of the charge has entered 6 ms after the stimulus at 60–250 µm). The score is 0.42–0.58.
- **Biology, stated honestly.**
  - The drive is Ca-triggered and uses VDCC Ca only, which matches the ECB_TRIGGER_LIT Verdict. It reads influx across the whole spine, not a nanodomain.
  - Its per-event saturation contradicts Lenz & Alger's linear scaling (Verdict c, Q5). So it is still a phenomenological "a VDCC Ca event occurred" proxy (Verdict: D1 is acceptable only as such a proxy), now expressed in Ca rather than V.
  - The veto is the mechanism that Q5 asks for: Ca must come *before* glutamate/mGluR→PLC (Nevian & Sakmann 2006), and Ca entering while the own glutamate is bound does not prime.
  - Caveats: mGluR involvement at L5-L5 is untested. B_NMDA serves only as a "glutamate is bound here" clock; the NMDAR conductance and its Ca are not used, since postsynaptic NMDARs are not needed.
- **Losing synapses beyond 250 µm is acceptable.**
  - 22–43% of those synapses have no VDCC event at all: silent, no LTD.
  - The first crossing is late (median 5.4–6.4 ms, p90 10–13 ms), so the −10 ms row fails at 43–65% of them.
  - That is consistent with no bAP-tLTD at distal inputs (Sjöström & Häusser 2006; Verdict last bullet).
- **Optional, +2 fixed parameters: own-history K** (candidates 1+2, K_i = 0.32 × the synapse's own recent bAP influx peak).
  - Scores 0.94 at ljp 0 and at ljp 25, [1.00 / 1.00 / 0.98 / 0.57], connections 1.00 / 1.00. It is flat for k 0.2–0.5 and n 4–8 (0.91–0.95), and it is independent of ljp.
  - With a realistic history that also contains the own (vetoed) EPSP peaks, it falls to 0.86 / 0.92. So it is not worth the extra state unless 150–250 µm matters.

## Candidates (per synapse: all [<60 / 60–150 / 150–250 / >250 µm]; per connection: majority / mean T; ljp 0 | ljp 25)
| candidate | new params | synapses, ljp 0 | synapses, ljp 25 | conn ljp 0 / ljp 25 | biology (ECB_TRIGGER_LIT) |
|---|---|---|---|---|---|
| ref D1: local-V event + hand veto (VDCC_DECODE) | θ_V, w | 0.99 [1.00/1.00/1.00/0.91] | – | 1.00/1.00 | V is upstream only; proxy (Verdict) |
| ref C1: influx event + hand 15 ms veto | K, w | 0.90 [1.00/0.96/0.90/0.52] | 0.93 [0.97/0.99/0.95/0.57] | 0.96/0.92, 1.00/0.88 | Ca-triggered, VDCC; per-event saturation, not linear |
| C1 without veto | K | 0.47 [0.78/0.49/0.22/0.35] | 0.35 [0.44/0.39/0.17/0.35] | 0.46/0.21, 0.33/0.33 | own EPSP Ca then primes T |
| **4: C1 + own-glutamate veto (1 − B_NMDA)** | **K** | **0.90 [1.00/0.96/0.90/0.52]** | **0.93 [0.97/0.99/0.95/0.57]** | **0.96/0.92, 1.00/0.88** | **Ca before glutamate (Nevian & Sakmann 2006; Q5)** |
| 4b: veto while B_NMDA > 0.8 | K, θ_b | 0.90 (same bins) | 0.93 (same) | same | ≡ a 15.6 ms window; θ_b is w in disguise |
| 4c: veto + S reset at own pre | K | 0.90 (same) | 0.93 (same) | same | primed state consumed at the read (Q5); lowers train gates (below) |
| 1: Hill onset on influx, S += pos(ΔHill(q;K,n)), n 8 | K (n fixed) | 0.90 [0.97/0.96/0.90/0.52] | 0.89 [0.97/0.95/0.88/0.52] | 0.96/0.83, 0.92/0.79 | cooperative Ca sensor, saturates per event; continuous C1 |
| 1, n 16 / n 4 | K | 0.92 / 0.88 | 0.90 / 0.86 | 0.96/0.92, 0.96/0.92 (n 16) | n 16 is not biological; n 4 costs 2–4% |
| 1: Hill onset on volume VDCC Ca c | K, n | 0.02 | 0.64 [1.00/0.83/0.27/0.00] | 0.00/0.00, 0.71/0.58 | volume-averaged (EGTA≈BAPTA), but the 12 ms Ca tail saturates bursts |
| 1: Hill rate (time above K), on c / on q | K, n | 0.00 / 0.08 | 0.01 / 0.09 | ≤ 0.17 | time above K ∝ log(amplitude) |
| linear influx integral, population scale | θ | 0.02 | 0.03 | ≤ 0.04 | the literal Q5 sensor (Lenz & Alger slope ~1) |
| 2: linear influx / own bAP charge (exact, NRM) | τ_M, M0 | 0.50 [0.91/0.56/0.32/0.00] | 0.42 [0.72/0.49/0.27/0.00] | 0.54/0.71, 0.50/0.67 | metaplastic gain; not in the literature file |
| 2: linear Ca integral / own (NRMc) | τ_M, M0 | 0.00 | 0.00 | 0.00 | ∫Ca builds over ~30 ms, so T ≈ 0 at −10 |
| **1+2: Hill onset, K_i = 0.32 × own bAP peak, n 8** | k (τ_K, K0 fixed) | 0.94 [1.00/1.00/0.98/0.57] | 0.94 [1.00/1.00/0.98/0.57] | 1.00/1.00, 1.00/1.00 | threshold relative to own Ca (Nevian: a threshold exists, no absolute value) |
| 1+2, history incl. own vetoed EPSP peaks | same | 0.86 [1.00/0.95/0.71/0.57] | 0.92 [1.00/0.97/0.93/0.57] | 0.92/0.96, 0.96/0.96 | 8–10% of synapses have vetoed EPSP peak > bAP peak |
| 3: other T nonlinearity (log, √, step, Hill, S/S_slow) | 0–1 | λ-tolerance ≤ 1.23-fold | (need ≥ 250-fold) | – | no-go, see §3 |

Train gates (median gate over 5 arrivals; not scored). At 20 Hz −10, where the data want LTD:
- hand veto 1.00;
- 1 − B_NMDA: 0.65 / 0.69 (a bAP 40 ms after the own pre is weighted 0.44);
- with reset: 0.40 / 0.32.

The 50 Hz −10 train is 0.75–0.96 in every case, and rho LTP must still win there, as it does now. The +10 trains are 0 for every vetoed drive.

## Findings per question
1. **Saturating Ca sensor.**
   - **On influx q: one K works across 2.4–3.8 decades.** It is best at K = 0.1–0.3% of the median bAP peak with n ≥ 8, and every synapse with q_pk ≳ 1.5K (n 8) gets a unit to within 4%.
   - It is **identical to C1** (the step limit), which scores as well or better. So the Hill adds nothing but smoothness; the drive is extracted as impulses anyway.
   - **The sub-K distal range goes silent (no LTD)**, not wrong. Beyond 250 µm, 13–43% of synapses are silent, and another 22–30% fail in other ways, mostly late crossings that miss the −10 row.
   - **On volume Ca c: no single K works.** The Ca tail (τ 12 ms) keeps proximal spines above K between burst APs (50 ms apart), so the later burst APs add nothing, while distal peaks are below K. The pass band is about 0.5 decade (estimated from the tail ratio), and the best score is 0.64 at ljp 25.
2. **Own slow-history normalisation.**
   - **Unstable to history.** After a pre-only baseline (the standard in-vitro recording), the own vetoed EPSP Ca becomes the unit, so pre-only trains are LTD-drive-equivalent at 100% of synapses.
   - **A floor M0 fixes this only by re-introducing amplitude coding.** At p25 of the bAP charge, 21–38% still give EPSP ≥ 0.5 unit, and 25% of bAP drives fall below the 4% band.
   - **With no bAPs the gain 1/M diverges up to 1/M0.** A rate-weighted mean also makes λ depend on the cell's firing rate.
   - **Applied to K of the Hill sensor instead (1+2), it is stable.** Each event adds at most 1 whatever K is. In silence K decays to K0, the veto keeps the EPSP out (E row passes at 100% for K0 = 1e-4 to 1e-2 × median), and the next bAP gives about 1.
3. **Different T nonlinearity: a no-go by construction** (`tol3.py`, ideal impulses at the somatic AP times).
   - For any static g(S) of one trace S = λΣe^{−t/τ_E1}, S stays above threshold for τ_E1·ln(λ/θ). So the −50 vs −100 row split tolerates at most e^{50/τ_E1}.
   - Measured widest λ range: 1.08 at τ_E1 100 and 1.23 at 200 (lin, √, step, log all alike); Hill-in-S 1.06.
   - Shorter τ_E1 widens it but kills burst summation (burst APs are 50 ms apart; the read comes 117/197 ms after the last AP).
   - λ-free ratios (S_100/S_2000) also cancel the burst count, so they pass nothing.
   - Needed: 250-fold (ljp 25) to 7000-fold (ljp 0). A count-preserving, amplitude-free step *is* per-event saturation, so the fix belongs in the drive, not in T.
4. **Mechanistic veto: works, with 0 new parameters.**
   - (1 − B_NMDA) with the existing tau_d 70 ms reproduces every scored row of the hand 15 ms veto for all event drives. B_NMDA jumps at release, so a crossing on the EPSP foot is already blocked.
   - Veto on the conductance time course instead (dual exponential, 1.6 ms rise): unchanged for the multiplicative form, but it drops the threshold form to 0.67–0.74, because the EPSP crossing in the first 0.5 ms escapes. Use B, not g.
   - The reset-at-pre variant changes no scored row but halves the 20 Hz −10 train gate, so do not use it.
   - Release failures give b = 0 in the mod. That is correct: no glutamate means no veto. The offline pipeline should use release times if it has them, otherwise arrivals.

## Implementation, if adopted (not done)
- **Offline** (`batch_v2._sig`, t_drive 3 path): impulses at the upward crossings of −ica_VDCC above K, each weighted by 1 − Σ exp(−(t − t_arr)/70) (capped at 0..1), instead of dropping those within w.
  - This needs the `ica_VDCC` trace (already a trace var), not `v_seg`.
- **Mod:**
  - In BEFORE STEP, on an armed upward crossing of −ica_VDCC through K: S += 1 − min(1, B_NMDA·Nrrp/factor), so that one released vesicle gives b = 1.
  - Re-arm below K (hysteresis optional).
  - The prefire must be re-run (cost ≈ 60–75 CPU·h, as in LOCAL_T_DRIVE.md; this must be logged in DECISIONS first).

**Files:** `cell_audit/ca_drive/`:
- scripts: `common.py`, `cands.py`, `extra.py`, `robust.py`, `silence.py`, `tol3.py`;
- results: `out/*.json`.

The replay cache (140 MB) was written to the session scratch dir via `$CA_CACHE`; rebuild it in ~20 s.

**Compute:** login node only, about 25 CPU·min in total (numpy, 1 thread per process). No Slurm job.
