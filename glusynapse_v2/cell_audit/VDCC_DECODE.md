# Why synapse-local Ca cannot drive the eCB-LTD window (2026-09-30, analysis only, no model code changed)

## Verdict
- **The main cause is the rule, not the channels.** The fitted window (t_drive 2 form: S += drive per bAP, τ_E1 100, pos(S − θ_Te), τ_T, one common θ_Tg, LTD rows ≥ 2θ_Tg) passes only if **every synapse gets the same per-bAP drive to within ±4%**.
  - The widest pass band is a 1.08-fold range of λ (per-bAP drive / fitted unit), for any θ_Te in 0.3–1.5 and τ_T 20–40 ms; τ_T 80 ms gives none (`tolerance.json`, `figs/tolerance.png`). This is because T ∝ ~(λ − θ_Te)², so a 10% change in λ moves the −50 and −100 rows 5–25-fold.
  - An event count gives exactly 1 per bAP, which is why D1 and `t_drive 2` pass. An amplitude-coded signal cannot: even the bAP **voltage** amplitude spans 1.0 decade (p5–p95), and 2.4–3.1-fold IQR inside each distance bin.
- **The channels then make it far worse.** A 10× fall in V becomes a 10^3.8 fall in spine VDCC Ca (Ca ∝ amp^3.5; e-fold per 8.7 mV of peak V). That is the m² foot of the **HVA-only** GluSynapse R-type (m² half-activation at +2.4 mV).
- **Total spine Ca (cai_CR) has an NMDA floor.** The own-EPSP spine Ca exceeds the bAP spine Ca at 49% of synapses; at 150–250 µm the bAP gives 0.75% of the EPSP value.
- **Volume-averaged integrals also fail on timing.** ∫Ca builds over ~30 ms, so at the −10 ms read-out (7 ms after the bAP) T ≈ 0. Score 0.00 in every distance bin, even with oracle per-synapse normalisation.
- **User's hypothesis (shifted curves): half right.**
  - The spine VDCC *is* right-shifted: its gating comes from 110 mM Ba²⁺ patches with no surface-charge correction.
  - The dendritic Ca_HVA2 and Ca_LVAst are *not* right-shifted against their sources. Their changes (see §1) all *reduce* availability at rest.
  - The documented correction (−25 mV) compresses the spine Ca spread from 3.8 to 2.4 decades. It rescues only **event-coded** Ca drives, never integrated ones.
- **Best Ca-based drive: a spine-VDCC Ca-influx event (threshold crossing of ica_VDCC) plus the own-pre veto, with ljp_VDCC = 25.**
  - Synapses: 0.92 (<60 / 60–150 / 150–250 / >250 µm: 0.97 / 0.99 / 0.93 / 0.57).
  - Connections: majority of contacts pass 1.00, mean-T 0.88.
  - For comparison, D1 scores 0.99 (1.00 / 1.00 / 1.00 / 0.91), 1.00 and 1.00.
  - The loss is confined to >250 µm, where Sjöström & Häusser 2006 suggest no bAP tLTD. This drive is still an event code, not the volume-averaged linear sensor of ECB_TRIGGER_LIT.md.

## 1. Channels as used (celsius 34; `vdcc_decode/channels.py`, `figs/channels.png`, `channels.json`)
| channel | m² V½ | h V½ | h∞(−78) | τm at −30 / −44 | m²·h∞(rest) at −44 / −27 / 0 mV | vs source |
|---|---|---|---|---|---|---|
| spine VDCC (GluSynapse R-type) | **+2.4** | −39 | 0.99 | 1 ms fixed (a guess), no q10 | 3e-4 / 0.010 / 0.42 | M&J 1995 values, measured in 110 mM BaCl₂, cell-attached |
| spine VDCC, ljp 25 | −22.6 | −64 | 0.82 | 1 ms | 0.033 / 0.30 / 0.76 | 110 mM Ba → 2 mM Ca surface-charge correction |
| dend Ca_HVA2 | −29.3 | −70.8 | 0.59 | 5.8 / ~7 ms | 0.0026 / 0.37 / 0.58 | Reuveni 1993 rates (= Hay Ca_HVA); **h shifted −12.3 mV** (unverified); τ ×1.18 (q10 to 36 °C) |
| dend Ca_LVAst | −34.8 | **−90.1** | **0.13** | 3.5 / ~8 ms | 0.016 / 0.10 / 0.13 | Hay 2011 file: A&J 1996 curves shifted −10 mV ("LJP"). A&J give inactivation V½ −80 in 2 mM Ca, which leaves 0.42 available at rest |

- **Source check.** Magee & Johnston 1995 recorded dendritic Ca channels in 110 mM BaCl₂ ([doi](https://doi.org/10.1113/jphysiol.1995.sp020862)).
  - High Ba²⁺ shifts activation positive by ~32 mV (100 mM Ba vs physiological; Smith et al. 1993, [doi](https://doi.org/10.1085/jgp.101.5.767)).
  - Ca²⁺ shifts it back +10 mV relative to equimolar Ba²⁺ (Zamponi & Snutch 1996, [doi](https://doi.org/10.1007/BF02207290)).
  - So the correction to 2 mM Ca is **−20 to −32 mV**. Chindemi's values carry only the m² correction, and s5 (−5 mV) covers a fraction of it. The spine VDCC has no low-threshold (T-type) component.
- **Densities** (S/cm², `SSCx-AAD-delta-emodels/cADpyr_L5TPC.hoc`; antic-delta differs only in basal NaTg/Ka):

| | soma | basal | apical | axon |
|---|---|---|---|---|
| Ca_HVA2 | 8.4e-4 | 2.39e-3 | 2.39e-3 | 5.1e-4 |
| Ca_LVAst | 2.2e-4 | 2.39e-3 | 2.39e-3 | 9.4e-4 |

  - Basal Ca_LVAst is **not zero**: the basal dendrites carry 24 pS/µm² of each type, uniform, which is more than Hay 2011 (no basal Ca channels).
  - Dendrites use `cad` (kE 62, γ₀ 0.24, τ ≈ 66·diam ms). Soma and axon use CaDynamics_DC0. GluSynapse only READs cai, so shaft cai carries no synaptic Ca.
  - **Missing dendritic Ca: no.**
- **Kinetic filter.** Both dendritic channels activate with τm 3.5–8 ms at −30 to −45 mV, much longer than a 1–2 ms bAP, so a bAP opens them far below m∞. This adds to the steepness. These τm values are unverified against 34 °C data.

## 2. Decoding on local_t/out (og, 191 synapses; `decode.py`, `decode.json`, `figs/decode.png`)
**Replay validation.** Gating was replayed open-loop on the recorded local V, with g and K fitted per synapse.
- Recorded VDCC charge is reproduced to 1.00 [1.00–1.01], and spine Ca peak to 0.996 [0.98–1.00], for og and for s5 (ljp 5, not fitted).
- g·K^(2/3) is constant to 0.9% CV, as the volume scaling predicts.
- Closed-loop ljp 25 (job 22088814, 2 pairs, 31 synapses) against the replay: VDCC charge median 0.99, worst 0.82 (EPSP); spine Ca peak 0.97–0.99. V changes ≤ 6 mV at clustered synapses.

Medians per bin, bAP-only (ap1) and own EPSP (epsp):

| bin (µm) | n | bAP amp (mV) / peak V | EPSP amp | spine VDCC charge, bAP/EPSP | spine cai_CR pk bAP (µM) / ratio bAP:EPSP | shaft cai pk (µM) / ratio |
|---|---|---|---|---|---|---|
| <60 | 32 | 77 / −0.8 | 5.0 | 5100× | 0.83 / 3.1 | 0.42 / 21000 |
| 60–150 | 95 | 48 / −30 | 9.2 | 245× | 0.043 / 2.1 | 0.077 / 1560 |
| 150–250 | 41 | 30 / −48 | 22 | **1.7×** | 0.0037 / **0.0075** | 0.0035 / 1.2 |
| >250 | 23 | 10 / −67 | 3.1 | 9× | 0.0006 / 5.8 | 0.00015 / 12 |

- **Attenuation.** V falls 10× (1.0 decade, p5–p95). VDCC charge falls 3.9 decades, spine Ca 4.0 and shaft Ca 4.4. The log-log slope against bAP amplitude is 3.5 / 3.6 / 4.1, i.e. e-fold per 8.7 / 8.4 / 7.9 mV of peak V.
- **Does the ~9 mV e-fold slope match the HVA foot? Yes.** R-type m∞² with k 9.5 goes e-fold per 5 mV at −30 and per 14 mV at 0 mV; the chord from −44 to 0 mV is 6 mV. bAP broadening with distance flattens it towards the 8.7 mV observed. This matches the parallel bap_ca_map (e-fold per 8–10 mV).
- **EPSP floor.** The own EPSP beats the bAP at 19% of synapses for VDCC current and for V (clustered distal contacts), and at 49% for total spine Ca (NMDA).
- **Shaft cai is driven locally.** Replayed HVA2+LVA charge vs shaft cai: log r = 0.98; the `cad` fit with diam free gives R² 0.997 (diam median 0.64 µm). LVA carries 79% of it (24% at <60 µm, 98% beyond 150 µm), yet 87% of LVA is inactivated at rest.
- **Which mechanism kills the signal, ranked:**
  1. amplitude coding versus the knife-edge rule;
  2. the HVA-only, Ba-shifted spine VDCC steepening the 10× V attenuation to ~10⁴×;
  3. the NMDA floor (total spine Ca);
  4. dendritic kinetics and inactivation (shaft).
- bAP attenuation itself (Na/K) is within data (BASAL_NAKV.md) and is not the limiter.

## 3. Counterfactuals (`counterfactual.py`, `ttype.py`, `tolerance.py`, `conn_check.py`; score as check_hp.py, strict, MARGIN 1)
**Drive shapes tested.** Linear or two-stage integrals of VDCC-only spine Ca (NMDAR share excluded, as ECB_TRIGGER_LIT Q5 asks), of total spine Ca, and of shaft cai. Grids span threshold 0–0.3× the median bAP peak, θ_Te 0.1–1.3 (0.3–1.0 in 0.01 steps for the oracle), τ_T 20–80, with or without the 15 ms own-pre veto. Event counts are included for reference.

| signal (variant) | best integrated drive | event + veto: all [<60 / 60–150 / >150] | conn: T summed / majority |
|---|---|---|---|
| spine VDCC current, og | 0.02 | 0.88 [1.00 / 0.96 / 0.70] | 0.46 / 0.92 |
| spine VDCC current, ljp 25 | 0.03 | **0.92 [0.97 / 0.99 / 0.80]** | 0.50 / **1.00** |
| VDCC-only spine Ca, ljp 0 / 15 / 25 / 30 | 0.00 | 0.43 / 0.62 / 0.76 / 0.80 | 0.33–0.42 |
| total spine Ca (cai_CR) | 0.00 | 0.43 | 0.33 |
| shaft cai: as used / HVA h unshifted / LVA at A&J / both | 0.00 | 0.08 / 0.08 / 0.09 / 0.10 | ≤ 0.08 |
| spine T-type (Ca_LVAst kinetics), calibrated to R0 at <60 µm | 0.00 | 0.87 (Ca event 0.49) | 0.42 / – |
| reference D1 (V event + veto) | – | 0.99 [1.00 / 1.00 / 0.97] | 0.58 / 1.00 |

- **Fits restricted to <60 µm or <150 µm, with bin-normalised thresholds.** Integrated drives still score 0.00. With exact per-synapse normalisation (oracle) they score 0.00 with all rows and 0.005 without the −10 row. Timing plus the ±4% tolerance, not distal synapses alone, is what kills them.
- **The summed-T connection score fails even for D1 (0.58)**, because contact counts range from 4 to 18 and the rule cannot absorb a 4.5× scale. Use per-contact majority or mean T instead (D1 1.00 / 1.00).
- **T-type in the spine (coordinator's hypothesis).** Density was calibrated so the <60 µm bAP Ca equals R0's 0.82 µM.
  - At 40–60 mV bAPs, spine Ca is 0.22 µM with Ca_LVAst as used, 0.035 at A&J values, 0.38 with R ljp 25, and 0.035 with R0.
  - Beyond 150 µm every variant gives ≤ 0.08 µM.
  - Spread: 3.6 decades for the T-type (e-fold 9.9 mV) and 4.7 at A&J values, against 3.8 for R0 and 2.4 for R ljp 25.
  - A T-type moves the foot left, but it is just as steep (k 6) and 87% inactivated at rest, so it does not flatten the distance profile.
- **Proposed fix, if a Ca-based drive is required.**
  1. Set `ljp_VDCC_GluSynapse = 25` (−20 to −32 mV is documented), with `gca_bar_VDCC` × 0.31 so the proximal bAP Ca stays at Chindemi's 1.4 µM (×0.18 keeps og's 0.82 µM). This changes effcai everywhere and needs a refit.
  2. Drive T from **spine VDCC influx events** (upward crossings of −ica_VDCC above ~0.3% of the median bAP peak), with the 15 ms own-pre veto. This is D1 with ica_VDCC in place of v: 0.92 vs 0.99 of synapses, all lost synapses beyond 150 µm.
  3. A volume-averaged, linear Ca drive (the literature sensor) would need the rule changed, not the channels. Candidates are a faster S (no τ_E1 build-up at −10 ms) and a per-synapse normalisation, e.g. a slow local Ca baseline or sliding threshold, so that θ_Te is expressed in units of the synapse's own bAP. This is not tested here.

**Files:** `cell_audit/vdcc_decode/` (channels, replay, decode, counterfactual, ttype, tolerance, oracle_check, conn_check, closed_loop(+_compare).py; run_closed_loop.sh; *.json; figs/channels, decode, ttype, tolerance.png).
**Compute:** everything above ran on the login node (numpy, ~10 min CPU in total), plus one Slurm array, 22088814. That array had 2 tasks, each 1 CPU / 1100M / 0:15, and measured 1:57 and 1:14 elapsed, MaxRSS 0.64 and 0.93 GB, CPU efficiency 98% and 95%: 0.05 CPU·h in total.
