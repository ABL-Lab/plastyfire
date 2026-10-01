# Nevian & Sakmann 2006: the extracellular preparation and the mGluR pathway

Nevian T, Sakmann B (2006) J Neurosci 26:11001, doi:10.1523/JNEUROSCI.1749-06.2006 (text: `nevian/…md`).
Status: design note (2026-09-30). No code changed. The model, jax and batch files are untouched.

## 1. Data (checked against the paper text)
Prep: rat barrel-cortex L2/3 PCs, P13–15, 32–36 °C, 2 Ca / 1 Mg, **bicuculline 10 µM** (GABA-A blocked, GABA-B not).
Input: a pipette near basal dendrites, 50–150 µm from the soma (~50 µm in Fig 1A), 3–10 µA / 100 µs, set for a
**single-component** compound EPSP of 1–3 mV (Fig 1B example: 2 → 3.5 mV). Protocol: 60 pairings at 0.1 Hz; readout at 20–40 min.
Δt' = first AP − EPSP onset, and Δt = nearest AP − EPSP onset.

| Protocol (post APs) | Control | Drug results (ratio ± SEM, n) |
|---|---|---|
| EPSP only / 3AP50 only | 0.97±0.08 (5) / 1.05±0.2 (6) | |
| 3AP50 Δt' −90 / ≤−90, ≥+150 | ~1.0, n.s. (4–6) [yaml 1.00±0.09, 6] | |
| **3AP50 Δt' −50 (Δt −10)** | **0.67±0.06 (9)** [yaml 0.68±0.05, 10] | MCPG 1.06±0.16 (7); AM251 1.09±0.14 (5); U73122 1.11±0.07 (4); heparin 0.68±0.12 (9, no block); post MK-801 **0.60±0.11 (7, no block)**; D-APV 0.98±0.11 (6, block); nimodipine 0.70±0.10 (8); Ni²⁺ 0.70±0.03 (4); nimo+Ni 1.00±0.05 (3); EGTA 1 mM / BAPTA 1 mM block (IC50 0.39 / 0.36 mM) |
| 3AP50 Δt' −30 | n.s. (8) [yaml 0.98±0.12] | |
| 3AP50 Δt' −10 | 1.5±0.3 (6) [yaml 1.42±0.19, 12] | |
| **3AP50 Δt +10** | **2.08±0.25 (11)** [yaml 2.01±0.22] | MCPG 1.79±0.33 (9); D-APV 1.08±0.11 (3); MK-801 0.95±0.19 (6); nimodipine 1.92±0.31 (7); 2 mM EGTA/BAPTA block (IC50 ~0.5 mM) |
| 3AP50 Δt' +50 | ~1.0 (Fig 2B) [yaml 0.92±0.11, 4] | |
| 1AP −10 / +10 | 0.80±0.07 (5) / 1.04±0.08 (10) | Ni²⁺ at −10: 1.05±0.16 (3, block) |
| 2AP50 −10 / +10 | 0.72±0.12 (5) / 1.95±0.31 (9) | |
| 3AP20 −10 / +10 | 0.72±0.14 (7) / 1.09±0.27 (5) | |
| 3AP100 −10 / +10 | 0.52±0.12 (4) / 2.29±0.48 (7) | **MCPG at −10: 1.53±0.24 (11), so LTD flips to LTP** |

Imaging: MCPG and U73122 leave spine Ca unchanged (MCPG at 100 Hz −10: 88±6 %; U73122: 99±11 %).
Same peak Ca can give LTP or LTD (Fig 8). The [yaml] values in `nevian_protocols.py` differ slightly from the text
for 3 rows (−50, +10, Δt'−10). All differences are within 1 SEM, so they are probably figure-derived. Fix only if Nevian is scored.
**Pharmacological signature:** LTD needs mGluR, PLC, CB1 and non-postsynaptic NMDARs (D-APV blocks it but post MK-801 doesn't). It does not need IP3R.
That is Sjöström's L5 downstream path (CB1 + presynaptic NMDAR) with an **mGluR/PLC trigger** added.

## 2. What the pipette delivers, compared with a unitary pair (Q1)
- **Recruited volume.** From the current–distance relation r ≈ √(I/K), with K ~1–4 mA/mm² (Stoney 1968, doi:10.1152/jn.1968.31.5.659;
  Histed 2009, doi:10.1016/j.neuron.2009.07.016), 3–10 µA activates axons within **r_act ≈ 30–90 µm**. Axons are excited before somata.
  Every fibre fires within <1 ms of the pulse, so the input is **synchronous** (onset jitter under 1 ms).
- **Recorded cell.** In our circuit a contact gives ~0.25 mV (L2/3→L2/3 unitary mean 1.13 mV from ~4–5 synapses, `ebner/things_to_note.md`).
  Literature: 1.0±0.7 mV from ~2.8 mostly basal contacts (Feldmeyer 2006, doi:10.1113/jphysiol.2006.105106).
  So 1–3 mV takes **~5–20 active contacts from ~3–20 fibres**, mostly within r_act of the tip. That is 1–5 unitary pairs' worth of input to the recorded cell.
- **Neuropil (the real difference).** The pulse fires every axon in that volume, not just the ones contacting the recorded cell.
  Take ~0.7–1×10⁹ synapses/mm³, so ~4×10⁵ within 50 µm. The recorded cell has roughly 300–1000 spines in that volume, and 5–20 of them are driven, so about 1–5 % of the local synapses are active.
  The same fraction of the neuropil gives **~10³–10⁴ co-active release sites, with the nearest one ~1.5–2.5 µm away** (0.55·n^−1/3).
  A unitary pair has 4–5 release sites spread over the whole dendritic tree, tens to hundreds of µm apart, in an otherwise silent slice.
  So the two preparations deliver similar EPSPs and similar numbers of contacts on the recorded cell. The ambient glutamate differs by orders of magnitude.

## 3. Why mGluR shows up here but not in L5 pairs (Q2)
1. **Spillover and cooperativity (main factor).** Group I mGluR (mGluR5) sits perisynaptically, in an annulus just outside the PSD
   (Luján 1996, doi:10.1111/j.1460-9568.1996.tb01611.x). It is activated when glutamate pools from many nearby sites, and uptake limits
   that pooling (Arnth-Jensen 2002, doi:10.1038/nn825; Scanziani 1997, doi:10.1038/385630a0). mGluR-driven eCB release needs
   *spatially dense* activation: dispersed activation of the same number of synapses fails (Marcaggi & Attwell 2005, doi:10.1038/nn1458; Brown 2003, doi:10.1038/nn1126).
   Sjöström 2007 (local PDF, doi:10.1016/j.neuropharm.2006.07.021) makes the same argument: extracellular stimulation "should generate more glutamate spill-over than the activation of unitary connections".
2. **The sequence detector is PLCβ.** PLCβ needs both Gq (mGluR) and Ca (Hashimotodani 2005, doi:10.1016/j.neuron.2005.01.004).
   With post before pre, bAP/VDCC Ca is already high when the glutamate arrives. That explains sequence sensitivity, the Ni/nimodipine and EGTA block, and why post MK-801 doesn't block.
3. **The L5 evidence is weaker than we've assumed.** Sjöström 2003 tested LY341495 only against *ACEA-induced* depression, which asks whether presynaptic autoreceptors carry the activity dependence.
   The mGluR dependence of eCB *production* in L5 tLTD was never tested in the local papers (2003, 2004, 2007). "L5 tLTD is mGluR-independent" is therefore not established.
4. **Cell type and pathway also matter; age doesn't.** Ages match (P13–15 vs P12–21). mGluR-dependent LTD does occur in **unitary** L4–L4 pairs (Egger 1999, doi:10.1038/16026),
   and extracellular L4→L2/3 tLTD needs mGluR, PLC and CB1 (Bender 2006, doi:10.1523/JNEUROSCI.0176-06.2006; the IP3R dependence there differs).
   So spillover amplifies the mGluR route, but it is not the only way to engage it.
5. **Our model forces a cell-type difference anyway** (`nevian_gate_check.py`, analytic, seconds). The L5-fitted `t_drive 2` eCB gate
   (paired_l5 td2 fits s1–s3) is **saturated at 1.00** for every Nevian −Δt burst, *including Δt' −90*, where the data show no change. It is also mGluR-independent by construction, so MCPG couldn't block it.
   Sjöström's window is broad: a 5-AP burst ending 200 ms before the pre spike still depresses. Nevian's is narrow: a 3-AP burst ending 50 ms before does nothing.
   One AP-driven eCB source can't give both. **In the Nevian prep the AP-only eCB route has to be ≈0, and LTD has to come from a fast Ca×mGluR coincidence.**

## 4. Other effects of the preparation (Q3)
- **GABA-A** is blocked, so the pilot's exc-only input is right. **GABA-B** is intact, but a single pulse at 0.1 Hz does little; ignore it.
  Bicuculline raises network excitability, but the "single-component EPSP" criterion excludes polysynaptic EPSPs.
- **Neuromodulatory fibres** (ACh, NA, 5-HT) are also recruited. Muscarinic M1 is Gq→PLCβ and can drive 2-AG as well. MCPG's complete block says mGluR is *necessary*; ACh is at most a modulator. Don't model it.
- **Astrocytes.** Astrocytes carry the eCB→glutamate→pre-NMDAR step in L4→L2/3 tLTD (Min & Nevian 2012, doi:10.1038/nn.3075). That fits D-APV blocking while post MK-801 doesn't, and it sits downstream of T, so it needs no new term.
- Citations are from memory: the PubMed tool was down this session. Only the Sjöström 2003/2007 quotes were read from local PDFs. Check the DOIs before they go into the paper.

## 5. The pipette mode (other session, read only) (Q4)
What it does (`plastyfire/pipette.py`, pilots `pip0-{143065,18559,22785}`):
- It picks a seeded EXC basal synapse at 50–150 µm path distance as the tip.
- It recruits **whole pre cells**, ranked by their nearest basal synapse, until the EPSP is ~2 mV (accepted 1–3 mV). All of each cell's synapses are active.
- Inhibition is excluded, release is stochastic, and GluSynapse uses circuit rho0.
- Result: **2–3 pre cells, 4–10 synapses, EPSP 1.87–2.41 mV**. 3–6 synapses fall in the 50–150 µm range.
- Measured here (edges h5): the synapses lie **0–300 µm (Euclidean) from the tip**, and only 2–4 per cell are within 20 µm. Their **delays are the circuit soma→synapse delays, 0.8–7.0 ms** (L4 pre 89346 in 18559: 6.7–7.0 ms).

Assessment:
- The recorded-cell side is roughly right: synapse count and EPSP are at the low end of the 5–20 contact estimate. Whole-cell recruitment is correct physics.
- Too few fibres are recruited, because our L2/3 connections are multi-synaptic. That is a circuit property; accept it.
- Needed changes (for the owning session):
  - (a) **Delays.** The spike starts at the tip, so use one common latency (~1 ms) for every group synapse, not the pre-soma delay. The 6 ms spread smears onset at Δt = ±10 ms.
  - (b) **Log per-synapse distance to the tip** in `pipette_<post>.json`; the mGluR covariate g_i needs it.
  - (c) The neuropil glutamate can't come from simulated neighbours: GluSynapse has no spillover. It enters only through g_i (§6). No extra cells are needed.
  - (d) **Baseline drift.** EPSP-only drifts rho to ~0.8 (data 0.97±0.08). This has to be fixed or normalised before any −Δt LTD is read.
- With (a)–(d), the pipette mode is the minimal faithful simulation: an L2/3 PC, synchronous clustered basal input of 1–3 mV, no inhibition, and Nevian timings.

## 6. Recommended design (Q5): T stays the common eCB readout; mGluR is a second, spillover-gated source
The existing presynaptic step is kept: at each arrival, `dpre −= A_mglu·tanh(pos(T − θ_Tg))·(dpre − dmin)`.
It stands for CB1 + pre-NMDAR, and `nmdar_block` (D-APV) and AM251 act on it. What's added is a second source of T at the arrival:

    S_M' = −S_M/τ_M,   S_M += 1 at each post AP          (fast post-AP Ca proxy at the spine, τ_M = 20 ms, fixed)
    at arrival j of synapse i, before the dpre step:
        T ← T + A_M · g_i · S_M(t_j⁻)                      (mGluR–PLCβ coincidence: glutamate × Ca already there)
    u_T (t_drive 2) → γ_ct · pos(S − θ_Te)/Te_scale       (AP-driven eCB, γ_ct = 1 for L5 post, 0 for L2/3 post)
    g_i = exp(−d_i² / 2 r_act²) for pipette input (d_i = synapse–tip distance, r_act = 50 µm); g_i = 0 for unitary pairs

- The jump comes before the gate because 2-AG production (tens of ms) is shorter than the presynaptic coincidence window (pre-NMDAR, ~100 ms). Folding both into the arrival event keeps the exact event-by-event solution.
- **Parameters.** A_M is the only free one. τ_M, r_act and γ_ct are fixed: τ_M from spine Ca decay, r_act from the stimulus current, γ_ct as a hypothesis switch forced by §3.5. g_i is a covariate computed from geometry.
- **Conditions.**
  - `mglur_block` (MCPG, U73122): A_M = 0.
  - `cb1_block` (AM251): A_mglu = 0, which is today's `mglu_block`. The Sjöström AM251 rows map here unchanged.
  - `post_nmdar` (MK-801): rho frozen, T intact, so it predicts LTD 0.60 kept and LTP lost, as the data show.
  - `nmdar_block`: unchanged.
- **Paired data alone.** Every paired target has g = 0, so A_M drops out exactly. The paired_l5 fit and its 14 parameters are unchanged (backward compatible; add a unit test).
- **What Nevian can do.** It can only *validate*. A_M is scanned in 1-D and never fitted, and the output is a curve: does any A_M reproduce the pattern?
- **Timing check with S_M** (arrival values; gate = tanh(A_M g S_M − θ_Tg) with A_M g = 1, θ_Tg = 0.28):

  | Protocol | S_M | Gate |
  |---|---|---|
  | 3AP100 −10 | 1.20 | 0.73 |
  | 3AP50 −10 | 0.91 | 0.56 |
  | 2AP50 −10 | 0.83 | 0.50 |
  | 3AP20 −10 | 0.66 | 0.36 |
  | 1AP −10 | 0.61 | 0.32 |
  | Δt' −90 | 0.12 | 0 |
  | all +Δt | 0 | 0 |

  The order matches the measured LTD depth (0.52 < 0.67 < 0.72 ≈ 0.72 < 0.80, and ~1.0 at −90), and all +Δt protocols are untouched, as MCPG-insensitive LTP requires.
- **Known risks.**
  - Δt' −30 (0.98) gets gate 0.50, so it needs rho LTP from the trailing AP to cancel.
  - MCPG 100 Hz → 1.53 needs rho LTP with all APs *before* the EPSP.
  - The L5 A_mglu (0.006–0.01 per pairing) caps Use at ~0.66–0.79 after 60 pairings, while the 100 Hz −10 target is 0.52.
  - Nevian +10 burst LTP remains the known cell-model limit (DECISIONS 09-29).
- **Rejected.**
  - Validation only with no mechanism: it leaves MCPG and Δt' −90 unexplained.
  - A pure cooperativity gate on the existing T: the L5 path already saturates in Nevian, so MCPG couldn't block it.
- **Testable prediction:** Nevian protocols on unitary L2/3→L2/3 pairs (g = 0) give no −Δt LTD in the model.

## 7. Test plan (small; no fitting)
| # | What | Cost |
|---|---|---|
| T0 | `nevian_gate_check.py` (done): the L5 t_drive 2 gate saturates on all Nevian −Δt protocols including −90; S_M ordering as above | 0 |
| T1 | Unit test: g = 0 ⇒ dpre identical to current model_v2 on the paired_l5 records (numpy, login node) | login node, minutes |
| T2 | Pipette prefire, 16 control Nevian protocols × the 3 pilot posts, **with fix (a)**. Run by the ebner session or with its OK. Sized from the measured pip0-18559 prefire (207 s per protocol, 2 workers, MaxRSS 15 GB): one job per post, `--cpus-per-task=2 --mem=19G --time=1:25:00` (16 × 207 s = 55 min, +50 %) | 3 jobs × 2 CPU × ~1 h ≈ **5.5 CPU·h** |
| T3 | Offline (numpy model_v2 + analytical_method rho, login node): A_M ∈ [0, 5] (21 values) × r_act ∈ {30, 50, 90} µm × γ_ct ∈ {0, 1} × conditions {control, mglur_block, cb1_block, post_nmdar, nmdar_block}. Report the χ² curve over the 21 Nevian rows and the pass/fail checks in §6 | login node, <30 min |
| T4 | Only if T3 passes: 30 pipette groups for the variance. Same per-post sizing as T2 (≈ 55 CPU·h); flag in DECISIONS first, per the >30 CPU·h rule | later |
