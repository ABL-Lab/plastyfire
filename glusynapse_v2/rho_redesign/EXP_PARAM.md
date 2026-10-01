# v5c in experimentalist quantities (Cpre, Cpre_APV, Cpre_Mg0, Cpost), 2026-10-01

User request: write the best rule (v5c) in the language of experimentalists rather than c_VDCC pool units, using
per-synapse Ca measurements:
- **Cpre**: one EPSP at 1 mM Mg (Chindemi's c_pre).
- **Cpre_APV**: one EPSP with NMDARs blocked (gNMDA 0). This is the VDCC (+AMPA-depolarisation) Ca of an EPSP.
- **Cpre_Mg0**: one EPSP at Mg 0. This is the full-NMDA EPSP.
- **Cpost**: one bAP (Chindemi's c_post).

Each comes as a peak of effcai_GB (cpre, cpre_apv, cpre_mg0, cpost) and as a peak of free spine Ca in µM (`*_cai`).
The measurement (A4) goes to `/scratch/dhuruva/split1/cexp/<pathway>.csv` (l5l5, l23l5, l23l23), keyed by
pre_gid, post_gid, syn_id.

Derived per synapse:
- NMDA share of one EPSP: (Cpre − Cpre_APV)/Cpre.
- Mg-block depth (NMDA Ca held back by Mg at rest): (Cpre_Mg0 − Cpre)/Cpre.

Reference fit: `results/v5_S1_V5c_s6_39.json` (split1, 39 targets, χ² 87.84, k 8, BIC 117.1):
- θ_V 1.031: LTP gate on the own-VDCC pool V.
- θ_eCB 20.68: eCB trigger on the (1−b)-weighted pool W.
- γd 34.8, γp 507.

**Unit bridge (C, ANCHOR_V5C §2).** One event with peak free spine VDCC Ca C_cai leaves a pool level of
P(C) = 48.3 · C_cai, read 10 ms after the event (tight: L5 basal q10–q90 44.7–60.0; L2/3→L5 46.6). So θ/P(C) is a
dimensionless "how many such events".
- The bridge was calibrated on bAPs.
- For Cpre_APV (a slower EPSP VDCC transient) it must be checked against A4's vdcc_q: the pool from charge alone is
  vdcc_q / i_scale.

## 1a. Pure relabelling (no refit)

Absolute values via the bridge:
- θ_V = 1.031 → **21 nM** peak free spine VDCC Ca. This is a floor ("the spine sees its own VDCC Ca at all");
  C-Vf shows it is identified (×0.5 / ×2 give Δχ² +11 / +15 on og-delta).
- θ_eCB = 20.68 → **0.43 µM** peak free spine Ca, with no glutamate bound. It sits inside the "submicromolar"
  Ca-assisted PLCβ range (Hashimotodani 2005, via review). That is a bracket only; no lock without the figure values.

Per synapse, once cexp exists, the computation is:
- r_V,X = θ_V / P(X) and r_E,X = θ_eCB / P(X), for X ∈ {Cpre_APV, Cpre, Cpre_Mg0, Cpost}.
- Report q10/25/50/75/90 per pathway and P(r < 1), i.e. the share of synapses where a single such event already
  crosses the threshold.

Readings:
- r_V,Cpre_APV: how many APV-EPSPs' worth of VDCC Ca the gate needs.
- r_E,Cpost: how many bAPs the eCB trigger needs, at +10 ms.
- P(r_E,Cpost < 1) is the share of synapses where a single bAP 10 ms before the pre spike gives tLTD. On og-delta
  this was 0.44.

Implemented in `fit_v6.py` (function `relabel`). Run it as MAXITER 0 with `CEXP_DIR=/scratch/dhuruva/split1/cexp` and
no scales. Results:
- The log gets `relabel <pathway> <event>: ...` lines, plus the NMDA-share and Mg-depth quantiles.
- The per-synapse table goes to `/scratch/dhuruva/split1/cexp_relabel/<tag>_<model>.csv`.
- The same run also re-asserts the v5c repro.
- It also prints the consistency of cexp cpre/cpost against the npz c_pre/c_post, which should be ≈ 0.

**Result (22156526, real cexp, unique synapses: L5 n 178, L2/3→L5 n 700; NaN 0; REPRO bit-identical).**

The cai bridge is unusable. Every `*_cai` peak includes the 0.070 µM resting floor (48.3 × 0.070 = 3.39), so APV
EPSPs and most bAPs sit on that floor. The numbers below therefore use the **charge bridge**: P = 1e5 · vdcc_q_X,
the pool peak of one event.

Pool of one event, q10 / q50 / q90:

| event | L5 | L2/3→L5 |
|---|---|---|
| EPSP (1 mM Mg) | 0.25 / 0.72 / 298 | 0.085 / 0.145 / 0.34 |
| EPSP + APV | 0.25 / 0.77 / 280 | 0.090 / 0.152 / 0.36 |
| EPSP, Mg 0 | 0.44 / 16.5 / 505 | 0.088 / 0.19 / 0.90 |
| bAP | 0.34 / 3.77 / 134 | 0.15 / 0.30 / 33 |

The v5c thresholds as multiples of one event (q10 / q50 / q90, then the fraction of synapses where a single event
crosses the threshold):

| ratio | L5 | L2/3→L5 |
|---|---|---|
| θ_V / P(bAP) | 0.008 / **0.27** / 3.1 (0.66) | 0.03 / **3.45** / 7.0 (0.26) |
| θ_V / P(APV-EPSP) | 0.004 / 1.34 / 4.1 (0.47) | 2.9 / 6.8 / 11.5 (0.03) |
| θ_eCB / P(bAP) | 0.15 / **5.5** / 62 (0.39) | 0.63 / **69** / 141 (0.14) |

What this shows:
- **The VDCC Ca of an EPSP does not depend on NMDA.** The APV and 1 mM charges are the same. In effcai, by contrast,
  an EPSP's Ca is 99.8 % NMDA (L5 q10 0.76).
- **Mg-block depth** (Cpre_Mg0 − Cpre)/Cpre is 15 / 60 / 86 at L5 and 64 / 89 / 104 at L2/3→L5.
- **Cpost/Cpre** is 0.0004 / 0.008 / 0.29 at L5 and 0.001 / 0.002 / 0.38 at L2/3→L5.

Reading:
- θ_V = 1.03 units = 1.03e-5 nA·ms ≈ **32 Ca ions through the spine's own VDCCs** (PARAM_ANCHORS units).
  That is 0.27 of a median L5 bAP but 3.5 of a median L2/3→L5 bAP (mostly apical).
- θ_eCB ≈ 0.43 µM-equivalent ≈ 640 Ca ions. That is 5.5 bAPs at the median L5 synapse and 69 at L2/3→L5, so a
  single bAP triggers eCB only at the 39 % (L5) / 14 % (L2/3→L5) of synapses with large VDCC Ca. These are proximal
  and basal.
- Per-synapse VDCC charge spans about 3 decades with location. The **absolute** v5c thresholds use that spread as a
  location filter: proximal synapses gate and trigger, distal ones do not. §3 shows that normalising per synapse
  removes exactly this.

## 1b. Rewrite candidates (every threshold a measured per-synapse quantity × one coefficient)

In every candidate, Chindemi is kept: θd = a00·Cpre + a01·Cpost and θp = a10·Cpre + a11·Cpost, with a00..a11, γd, γp
free. What the four a's mean:
- **a00:** LTD starts once c* exceeds a00 single-EPSP Ca plus a01 single-bAP Ca.
- **a10, a11:** the same reading for LTP.

The new coefficients are k_V and k_E. Where the user rule needs an anchor, it is listed; otherwise the coefficient is
**flagged**. d_min, τE1, A_eCB and τ_d_NMDA stay as in v5c.

| cand. | gate (LTP) | eCB trigger | θd / θp | k | anchor status |
|---|---|---|---|---|---|
| E0 relabel | V > θ_V (21 nM equiv.) | W > θ_eCB (0.43 µM equiv.) | Chindemi | 8 | both flagged; θ_eCB bracket only (submicromolar PLCβ) |
| **E1** | V > k_V · P(Cpre_APV) | W > k_E · P(Cpost) | Chindemi | 8 | k_V flagged; k_E anchor candidate (tLTD window) |
| **E2** | V > k_V · P(Cpost) | W > k_E · P(Cpost) | Chindemi | 8 (6 if both anchors hold) | k_V anchor candidate (LTP frequency threshold); k_E as E1 |
| **E3** | as E2 | as E2 | θp = a10 · **Cpre_Mg0** + a11 · Cpost, θd Chindemi | 8 | as E2; a10 stays Chindemi-free |

What each new coefficient means:
- **k_V (E1):** LTP needs the pairing to give own-VDCC Ca at least k_V× that of a single EPSP with NMDARs blocked
  (APV).
  - There is no paper value. The test would be how many EPSPs (APV, no bAP) must summate in 100 ms before LTP is
    licensed, which has not been measured as such.
  - **Flagged.**
- **k_V (E2/E3):** LTP needs at least k_V bAPs' worth of own-VDCC Ca within the pool memory (τE1 100 ms).
  - With +10 ms pairings at f Hz, the pool peaks at 1/(1 − e^(−1/(f·τE1))) single-event levels: 1.58 at 10 Hz,
    2.54 at 20 Hz, 4.5 at 40 Hz.
  - LTP frequency threshold between 10 and 20 Hz (Sjöström 2001, 10.1016/S0896-6273(01)00542-6, rat V1 L5→L5) gives
    k_V ≈ 1.6–2.5.
  - Those protocols are fit targets, so this is a consistency bracket only. A lock needs an independent preparation's
    frequency threshold.
  - **Flagged, with that candidate.**
- **k_E (E1–E3):** eCB-LTD fires at an own pre spike if the glutamate-free VDCC Ca in the pool is at least k_E× that
  of one bAP read 10 ms after it.
  - Equivalently, a single bAP gives tLTD for pre spikes up to T_LTD = 10 ms + τE1·ln(1/k_E) after it.
  - Anchor candidate: the independent eCB coincidence window, T ≈ 125 ms (rat S1 L4→L2/3, Bender 2006,
    10.1523/JNEUROSCI.0176-06.2006), gives k_E = e^(−1.15) = **0.32** (C/E, confidence L, other pathway).
  - The L5 window (−10/−25 LTD, −100/−120 none) gives 0.41–0.92, but it is in the targets.
- **a10 with Cpre_Mg0 (E3):** the LTP threshold counts the EPSP at its unblocked NMDA Ca.
  - With a10 = a11 = 1, LTP needs more Ca than the linear sum of a Mg-free EPSP and a bAP, which is a strict
    supralinearity criterion.
  - NMDA dependence then enters through the Mg-block depth Cpre_Mg0 − Cpre of each synapse.
  - Free (Chindemi slot).

Structural notes:
- **Synapse-local only.** E1–E3 use only the synapse's own Ca measurements as constants, which is the same status as
  Chindemi's c_pre/c_post. The drives are unchanged: own VDCC Ca, own pre spikes, effcai.
- **Uniform across pathways.** The rule is uniform. Pathway differences enter only through each synapse's own
  measured Cs.
- **Main risk.** The tLTD magnitude changes character:
  - In v5c, a global θ_eCB lets only part of the synapses trigger on a single bAP (P ≈ 0.44), which grades the
    population LTD.
  - With θ_E ∝ Cpost, the single-bAP trigger becomes all-or-none per synapse across the population (k_E < 1: all;
    > 1: none), so the population LTD is about d_min (−0.29) or 0.
  - The −0.29 ceiling is close to the L5 tLTD data (≈ 0.7), so this may be acceptable. The refit decides.
- **E1 degeneracy.** E1 degenerates if Cpre_APV ≈ 0 at most synapses (AMPA-only EPSPs open few VDCCs): the gate
  threshold → 0 loses the floor role. In that case E2 replaces E1 in the refit pair.

Original plan (superseded, see §2 bridge note and §3): **E1 and E3**. Together they cover all three new quantities, and E2 is the nested
control if E3 wins. Replace E1 by E2 if the median P(Cpre_APV) < 0.1·θ_V (degenerate). For each: seeded from
v5_S1_V5c_s6_39 (thresholds converted at the median synapse), plus one unseeded basin check. Compare on BIC
(n 39) against 117.1; with k 8 the bar is χ² < 87.84.

## 2. Implementation (new files only)

- `gpu_v6_rho.py`:
  - `GPUModelV6(GPUModelV5)`. Its kernel is a copy of the v5 kernel with θ_V,i = θ_V·uV[i] and θ_eCB,i = θ_eCB·uE[i]
    (the E2 latch also reads θ_V,i). The gate on/off switch is still the global θ_V > 0.
  - `set_scales(uV, uE, xd, xp)`; xd/xp replace c_pre in θd/θp (rho_gamma 1 only).
  - uV = uE = 1 is v5c exactly.
- `fit_v6.py`: copy of fit_v5.py (including A7's WEIGHTS).
  - Options: `CEXP_DIR` (or STANDIN), `SCALE_V` / `SCALE_E` / `THETA_PRE` (json {col: coef}, column `one` = 1),
    `BOX`, `CEXP_MAP`.
  - Seed conversion: by the median u, applied to θ's and to a00/a10.
  - Relabel report (§1a).
  - MAXITER 0 repro assert at 1e-6 (exit non-zero on a diff).
- `run_fit_v6.sh`: copy of run_fit_v5.sh, output `results/v6_<TAG>`. `STANDIN_SMOKE=1` adds a stand-in run
  (E3 spec, MAXITER 2, output on /scratch).
- **Bridge used in the refits: the charge (q), not cai.** A4's csv has the own VDCC charge of each event (vdcc_q,
  vdcc_q_apv, vdcc_q_mg0, vdcc_q_post, nA·ms). The charge is the pool's own input, so 1e5 · vdcc_q_X is that event's
  pool peak with no bridge constant.
  - k = 1 means "as much own-VDCC Ca influx as one such event".
  - Pilot rows (L5 basal) show why cai is not used:
    - 48.3 · cpost_cai differs from 1e5 · vdcc_q_post by ×1.4–2.2, so cai is not proportional to the pool input.
    - cpre_apv_cai ≈ 0.07 µM, which is probably the resting level (to be checked).
    - In effcai units Cpre_APV/Cpre ≈ 6e-4, which is degenerate (A4).
  - In charge units Cpre_APV is not degenerate: vdcc_q_apv ≈ vdcc_q ≈ 3e-6 nA·ms, i.e. 0.3 pool units. The EPSP's
    VDCC charge barely depends on NMDA. So E1 is kept in its charge form.
- **Undefined measurements (uniform rule).**
  - Affected synapses: cpost and vdcc_q_post are NaN at the 8 post cells that cannot fire one AP. That is L5 pair
    207453-189325, 13 of 191 synapses, and 7 L2/3→L5 post cells, 56 of 798 synapses. The npz c_post is invalid at
    the same cells.
  - Rule: any NaN scale gets the pooled median of the finite scales over all pathways, i.e. the v5c absolute
    threshold at those synapses. For Xd/Xp it is c_pre × the pooled median X/c_pre.
  - The count is logged as `NaN -> median rows`. The Chindemi θd/θp at those synapses keep the npz c_post, as in
    every earlier fit.
- Candidate env settings (split1 39-target env as in run_split1_refits.sh, `SET={"v5_mode": 2}`,
  `FREE=theta_V,theta_eCB`, so θ_V/θ_eCB are k_V/k_E; `CEXP_DIR=/scratch/dhuruva/split1/cexp`; files L5L5/L23L5/L23L23
  are the CEXP_MAP default):
  - E1: `SCALE_V='{"vdcc_q_apv": 1e5}' SCALE_E='{"vdcc_q_post": 1e5}'`
  - E2: `SCALE_V='{"vdcc_q_post": 1e5}' SCALE_E='{"vdcc_q_post": 1e5}'`
  - E3: E2 + `THETA_PRE='{"p": {"cpre_mg0/cpre": 1.0}}'`. Xp = npz c_pre × the cexp ratio Cpre_Mg0/Cpre, because cexp cpre differs from the cache pair by pair (L2/3→L2/3 48 of 120 pairs > 2 % off; A4 CEXP.md). The Chindemi c_pre/c_post stay the cache values. The ratio is 62–104 (median), so the seed's a10 is divided by about 80.
  - Box: k_V 0–50 (linear) and k_E 0.05–50 (log) are the v5 boxes; override with BOX if a fit hits an edge.
- With the charge bridge the readings in §1b become:
  - "k_V × the VDCC influx of one APV-EPSP (E1) / one bAP (E2, E3)";
  - "k_E × the VDCC influx of one bAP". The pool peaks at the event, so the tLTD window edge is
    T = τE1 · ln(1/k_E), and Bender's 125 ms gives k_E ≈ 0.29.

## 3. Jobs

| job | tag | what | size | result |
|---|---|---|---|---|
| 22155608 | S1_V6r | MAXITER 0 rescore of v5_S1_V5c_s6_39, v6 kernel at uV = uE = 1, then the STAND-IN smoke | 2g, 40G, 0:15 | **REPRO OK, bit-identical**; 5:30, 37.40 GB |
| 22156526 | S1_V6relabel | MAXITER 0, real CEXP, no scales: §1a numbers (log `relabel ...`, /scratch/dhuruva/split1/cexp_relabel/v6_S1_V6relabel_*.csv) + repro | 2g, 47G, 0:15 | pending |
| 22156527 / 22156528 | S1_V6E2_s5 / _s6 | E2 refit, seeded v5_S1_V5c_s6_39 / unseeded | 2g, 50G, 0:30, afterok 22156526 | pending |
| 22156559 / 22156560 | S1_V6E3_s5 / _s6 | E3 refit, seeded / unseeded; Xp = npz c_pre × cexp (Cpre_Mg0/Cpre). 22156529/30 were cancelled while pending: they used the absolute cexp Cpre_Mg0 | same | pending |
| 22156531 / 22156532 | S1_V6E1_s5 / _s6 | E1 (charge form) refit, seeded / unseeded | same | pending |

Compare on n 39, k 8 against v5_S1_V5c_s6_39 (χ² 87.84, BIC 117.1). If a fit lands near the
anchor candidates (k_E ≈ 0.29, k_V ≈ 1.6–2.5), the next step is a k-6 validation refit with those values fixed.

A4 facts (CEXP.md) to use in the analysis:
- The NMDA share of single-EPSP Ca is 0.998–0.9995, and Cpre_APV sits at the resting floor (0.070 µM). So E1 is
  degenerate in effcai and cai. It runs only in its charge form: vdcc_q_apv is non-zero, about 0.3 pool units at the
  L5 pilot rows. Treat E1 as an exploratory gate.
- Cpost/Cpre depends strongly on location: L5 basal 0.033 / apical 0.002; L2/3→L5 basal 0.14 / apical 0.002;
  L2/3→L2/3 basal 0.009 / apical 0.23.

### 3a. Results (39 targets; v5c reference 87.84 = L5 73.17 + L2/3→L5 14.66, BIC 117.1; k 8 for all)

| fit | χ² total | L5 / L2/3→L5 | BIC | k_V | k_E | γd / γp | job (seff) |
|---|---|---|---|---|---|---|---|
| E2_s5 | 110.34 | 66.72 / 43.62 | 139.6 | 0.59 | 0.168 | 37 / 487 | 22156527 (12:27, 25.9 GB) |
| **E2_s6** | **97.00** | **62.80** / 34.20 | 126.3 | 0.88 | 0.168 | 44 / 558 | 22156528 (12:24, 25.9 GB) |
| E1_s5 | 119.65 | 84.12 / 35.52 | 149.0 | 2.54 | 0.75 | 25 / 588 | 22156531 (13:57, 32.8 GB) |
| E1_s6 | 98.78 | 93.73 / **5.05** | 128.1 | 1.29 | 0.64 | 119 / 264 | 22156532 (13:03, 26.7 GB) |
| E3_s5 / s6 | 1035.7 / 1097.9 | — | — | — | — | — | 22156559/60 (3:06/2:24): failed, see below |

All three rewrites lose to v5c. E2 is better than v5c on L5 alone (62.8 vs 73.2) but much worse on L2/3→L5.
E1_s6 is the reverse.

**Why per-synapse gate scaling hurts.** Per-target Δz² of E2_s6 against v5c, read from the results csvs:

| target | condition | data | v5c | E2_s6 | Δz² |
|---|---|---|---|---|---|
| Letzkus 3AP 200 Hz +10, distal (L2/3→L5) | control | 0.79 | 0.84 | 1.04 | **+16.1** |
| Sjöström 0.1 Hz −120 | control | 1.05 | 0.93 | 0.72 | **+9.6** |
| Sjöström07 step pair | control | 1.62 | 1.51 | 1.40 | +7.2 |
| Letzkus 3AP −10, proximal | control | 0.89 | 0.81 | 0.78 | +5.7 |
| Sjöström 20 Hz −10 | mglu_block | 1.02 | 1.11 | 1.18 | +3.4 |
| Sjöström 0.1 Hz −10 / −25 | control | 0.69 / 0.65 | 0.88 / 0.89 | 0.72 / 0.72 | −7.3 / −6.3 |
| Sjöström 20 Hz −10 / 10 Hz −10 | control | 0.65 / 0.57 | 0.94 / 0.86 | 0.84 / 0.76 | −6.0 / −3.6 |

1. **The gate loses its floor.** With θ_V,i = k_V · P(bAP)_i, the threshold falls to about 0.1 units at distal
   apical synapses, whose bAP charge has q10 0.15. Those synapses then pass the gate and potentiate (Letzkus distal
   +10: 1.04 instead of 0.79). v5c's absolute floor of about 32 Ca ions kept them out. The gate's job in v5c is
   "does this spine have VDCC Ca at all?", which is an absolute quantity. A ratio to the synapse's own bAP cannot
   express it.
2. **eCB becomes all-or-none across the population.** With θ_eCB,i = k_E · P(bAP)_i (k_E 0.17), every synapse has
   the same single-bAP window, T = τE1 · ln(1/k_E) ≈ 180 ms:
   - Gain: −10 and −25 now give full tLTD (0.72), which is the main improvement.
   - Loss: −120 depresses too (0.72 vs 1.05).
   - In v5c, the spread of P(bAP) relative to one global θ_eCB gave each synapse its own window, so the population
     window was graded (−10 > −25 > −120). The k_E fit (0.17) is in the range of the Bender-window value (0.29).
3. **E1** (gate ∝ APV-EPSP VDCC charge) keeps a floor-like gate on L2/3→L5 (θ_V/P(APV) q10 2.9, 5.05 there) but
   fails on L5, where the APV-EPSP charge spans 3 decades as well.

**E3 diagnosis.** Xp = c_pre · Cpre_Mg0/Cpre has a ratio median of 78 (L5 q10–q90 16–88; L2/3→L5 64–105).
- Only the seed was converted (a10/78). The DE box for a10 stayed [0, 5], so a10 ≈ 0.02 sits at the edge, and random
  members have θp up to 400 × c_pre.
- The converted seed also fails fit_v2's admissibility rule (≥ 15 % of synapses with θp below their peak) by 1.1 %,
  because of the ratio spread.
- So **0 of 64 members were admissible in every generation** (log: `0 admissible`). DE minimised the 1e3 penalty and
  converged on it after 78 generations. The χ² of about 1000 is a rule with no LTP, not a fit.
- **Fix (E3n):** normalise the ratio in the spec: `THETA_PRE={"p": {"cpre_mg0/cpre": 0.01282289}}` = 1/77.9855. Then
  Xp/c_pre has median 1, and a10 keeps Chindemi's scale and box. Experimental meaning: a10 · Cpre_Mg0/78.
- Jobs: rescore 22159388 (S1_V6E3n_r, MAXITER 0, 47G 0:15) and refit 22159389 (S1_V6E3n_s5, seeded, 41G 0:30).

**E3n result: dropped.** Rescore 22159388 gave 321 (2:46, 28.8 GB). The refit 22159389 gave 3652 (4:04, 37.2 GB);
its DE stayed on the admissibility penalty, with ≤ 13 admissible members. With Cpre_Mg0 in θp, the potentiation
rule cannot be satisfied on this data. Cpre_Mg0 stays a reported quantity only (Mg-block depth, §1a).

### 3b. Joint winner V7vn2 (A1): `results/v7_W2_N3.json`, 65 targets, k 8

The joint set is L5 40, L2/3→L5 9 and L2/3→L2/3 16 targets.

| fit | χ² total | L5 | L2/3→L5 | L2/3→L2/3 | BIC |
|---|---|---|---|---|---|
| **V7vn2 W2_N3** | **276.39** | 128.39 | 19.12 | 128.88 | **309.8** |
| W2_N2 (seeded) | 276.94 | | | | |
| V7v joint J2 | 294.34 | | | | |
| v5c params, rescored | ≈ 437 | | | | |

On 39 targets, V7vn2 gives 68.15 / 68.12 (k_E 0.173 / 0.182) and V7v 67.77; v5c gives 87.84.

Worst residuals of W2_N3 (z):
- Zilberter 1AP +10: +3.9
- Zilberter 5AP 20 Hz −10: −4.1
- Zilberter train10 +4 (last) / +5: −4.0 / +5.2
- Sjöström 2001 40 Hz dt 0: +4.7
- Sjöström 2001 50 Hz −10: −3.3
- Sjöström 2001 0.1 Hz +10: +3.4

So the L2/3→L2/3 timing at +4/+5/+10 and the Sjöström 2001 high-frequency LTP are where the rule still misses.

## 4. The best rule in experimentalist language (V7vn2, `results/v7_W2_N3.json`)

**Per-synapse measurements (calibration, Chindemi status).** Each is one event, measured at the synapse:
- **Cpre:** one presynaptic spike at 1 mM Mg, the EPSP Ca.
- **Cpost:** one bAP alone, the bAP Ca.
- **P1:** the own-VDCC Ca influx of that bAP (cexp vdcc_q_post; 1 pool unit = 1e-5 nA·ms ≈ 31 Ca ions).

**Rule.**
- **ρ (LTP / LTD, Chindemi).** c* is the spine Ca integrated over 0.28 s.
  - If c* > θd = 1.22·Cpre + 2.41·Cpost, the synapse depresses; if c* > θp = 1.79·Cpre + 2.73·Cpost, it potentiates.
  - Read: "LTP needs integrated Ca above 1.8 EPSPs plus 2.7 bAPs".
  - Rates γd 47.2 and γp 344 (s⁻¹ scale, Chindemi τ_ind 70 s, ρ* 0.5).
  - **Free:** a00, a01, a10, a11, γd, γp.
- **LTP licence.** Potentiation counts only while the spine's own VDCCs have let in at least θ_V = 3.58 units
  ≈ **110 Ca ions** within the last ~100 ms (τE1). A potentiating crossing without the licence counts as depression.
  - In bAP units (§1a bridge, rescaled to θ_V 3.58), θ_V is about 1 bAP at the median L5 synapse (q10–q90 0.03–11)
    and 12 bAPs at the median L2/3→L5 synapse (0.1–24); L2/3→L2/3 was not in the relabel.
  - So it is a floor that proximal spines pass with one bAP and distal spines only with bursts.
  - **Flagged** (fitted, no anchor; Kampa C-V null).
- **eCB-LTD (presynaptic).**
  - **Trigger:** at each own presynaptic spike, the eCB step fires if the synapse's glutamate-free VDCC Ca pool W
    exceeds θ_eCB,i = k_E·P1_i. W is the VDCC influx weighted by (1 − NMDAR glutamate-bound state b), with
    τ_d_NMDA 70 ms.
  - k_E = **0.487**: a single bAP up to τE1·ln(1/k_E) ≈ **72 ms** before the pre spike gives tLTD.
  - **Veto:** the step is cancelled if the synapse's own VDCC Ca influx over [t_pre, t_pre + 25 ms] exceeds θ_eCB,i.
    In other words, a postsynaptic spike shortly after the pre spike (pre-before-post) prevents eCB-LTD.
  - **Effect:** U_SE × (1 + d_min), with d_min = −0.29.
  - Status of the constants:
    - k_E is **anchored**: it lies inside A1's independent Sjöström 2003 bracket, 0.37–0.78.
    - d_min −0.29 is **anchored**: the CB1-agonist ceiling, Sjöström 2003.
    - τ_d_NMDA 70 ms is **anchored**: the synapse mod-file NMDAR decay.
    - τE1 100 ms is an E bracket (Heinbockel 2005 / Bender 2006).
    - A_eCB 1 is a structural limit.
    - The veto window of 25 ms is A1's (see A1's notes for its status).
- **Count:** k = 8 fitted. These are the 6 Chindemi parameters plus θ_V (flagged) and k_E (fitted, inside an
  independent bracket).

**What an experimentalist would measure to test this rule.**
- **Calibration.** At identified spines, image the Ca transient of:
  - one EPSP at 1 mM Mg (Cpre);
  - the same EPSP in APV, to give its VDCC part;
  - the same EPSP at Mg 0, for the Mg-block depth (predicted ≈ 60–90× Cpre in effcai terms, but it is not used by the
    rule);
  - one bAP alone (Cpost, and with Ni²⁺/Cd²⁺ the VDCC part P1).
- **Location predictions.** The model predicts that P1 falls by about 3 decades from proximal basal to distal apical
  spines. The predicted single-bAP tLTD window edge then scales with the synapse's own P1 only through k_E, so it is
  ≈ 70 ms everywhere.
- **Key new prediction: the 25 ms veto.** A post spike within 25 ms after a pre spike should abolish eCB-LTD for that
  pairing even when a bAP preceded the pre spike. Test with post-pre-post triplets, with AM251 as the control.
- **LTP licence.** It predicts that distal spines with little VDCC Ca need ≥ 10 bAPs per 100 ms (bursts) before
  pairing can potentiate, while proximal spines need about one. Block VDCCs locally (Ni²⁺) and compare.

## 4b. Shaft-Ca licence instead of spine VDCC ions (G8X_G1a, `results/v7_G8X_G1a_u6.json`, v7x timing, 65 targets)

This version replaces the ion-count LTP licence with a measurable dendritic Ca concentration. No new spine channel is
added. Everything else in §4 is unchanged (ρ thresholds, eCB trigger, the 25 ms veto, d_min); here the trigger is read at
the exact pre arrival (v7x).
- **LTP licence (shaft).** Potentiation counts only if the free Ca in the dendritic shaft under the synapse has
  risen ≥ **0.15 µM above rest** at some point in the last **100 ms**. A potentiating crossing without the licence
  counts as depression.
  - θ_G 0.15 µM is **anchored**: bAP-evoked dendritic Ca transients, Helmchen 1996. It is not fitted.
  - When θ_G is left free it fits to 0.168 µM (288.58 / BIC 322.0), close to the anchor.
  - Median single-bAP shaft ΔCa within 100 µm: L5 0.53 µM and L2/3→L2/3 0.23 µM, so proximal synapses pass with one
    bAP. Distal synapses need bursts.
  - Its Spearman correlation with the spine VDCC P1 is 0.86, so the shaft gate ranks synapses almost like the ion gate.
- **Fitted parameters:** γd 35.2 and γp 370.7. k_E is 0.496, which stays inside the Sjöström 2003 bracket.
- **Count:** k = **7**, the 6 Chindemi parameters plus k_E. With the anchored shaft licence there is no flagged parameter.
- **Score:** χ² 292.21 / BIC 321.4. Per pathway: L5 150.7, L2/3→L5 17.5, L2/3→L2/3 124.1.
  - The spine-VDCC version on the same timing (W2_X_s7) gives 278.0 / BIC 311.4 with k 8.
  - So the anchored shaft gate costs about 14 χ² (about 10 BIC) for one parameter fewer and a rule stated in µM.
- **Experimental test:** image the shaft next to the spine with a low-affinity dye during the pairing. The model predicts that
  pairings whose shaft ΔCa never reaches 0.15 µM within the last 100 ms give no LTP, whatever the spine Ca.
- **Live BCL (GluSynapseV8, gate_src 2): PASS 68/68**, live χ² 298.1 vs offline 292.2 (BV/BCL_G8X_G1a_u6.md).
