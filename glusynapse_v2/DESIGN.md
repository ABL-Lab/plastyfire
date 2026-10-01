# GluSynapse_v2 — presynaptic plasticity extension (design spec, T8)

Status: draft 2026-09-27. Owner session: plasty_model_research. Board: tsk project `glusynapse_v2`, thread `glusynapse-v2` (T8–T19).

## 1. Why a v2

* The current rule (Chindemi 2022, `GluSynapse.mod` as loaded from
  `DEES_cell_packages/other_mods/modified_mechanisms/GluSynapse.mod`) has **one** bistable efficacy ρ
  that drives both gmax_AMPA (post) and Use (pre) through Use_d/Use_p, gmax_d/gmax_p.
* Chindemi 2022 states the gap itself (Discussion): with a single ρ the model "could not reproduce
  presynaptic-only LTD (Sjöström 2003)" and suggests "modeling separate pre- and post-synaptic efficacies".
* Ebner 2019 shows that four pathways fit Sjöström 2001, Nevian & Sakmann 2006 and Letzkus 2006 in one
  L5 PC model: pre-LTD (mGluR→PLC→eCB→CB1), pre-LTP (L-VGCC→NO→GC, gated by preceding presynaptic activity),
  post-LTD and post-LTP (NMDAR Ca²⁺, phosphatase/kinase). Our ρ rule already covers the two post pathways.
* Tong 2020 (Padamsey 2017 data): ΔPr = η(H − Glu). Hebbian/NO signal raises Pr **only if presynaptic
  activity precedes NO release by ~7–10 ms**; glutamate release (pre-NMDAR) lowers Pr. NO release requires
  strong dendritic depolarisation / L-VGCC (dendritic spikes, bursts).
* Tong 2021: NO from strongly depolarised, clustered-active dendrite **lowers** Pr at nearby
  (≤4 µm) *inactive* terminals (hetLTD; needs NO + CaMKII, abolished by L-NAME; independent of post changes).

Together: NO sign depends on whether the terminal was active just before the NO signal; mGluR/eCB gives
post-before-pre LTD; both act on release probability, independently of ρ.

## 2. Minimal-change principle

1. v2 = v1 + one presynaptic modifier state `dpre_GB` and three traces. Nothing in the existing ρ, calcium,
   AMPA/NMDA, VDCC or TM-release code changes.
2. With all v2 amplitudes = 0 (the defaults), v2 is **bit-identical** in behaviour to v1 (regression test, T14).
3. New POINT_PROCESS name `GluSynapse_v2`, own directory `glusynapse_v2/mod/`, own compiled library.
   The shared `DEES_cell_packages/x86_64` is never overwritten.
4. edges.h5 fields (Use_d/Use_p, gmax_d/p, theta_d/p, rho0, …) keep their meaning.

## 3. Equations

Existing (unchanged): ρ rule with `dep*(1-pot)` gating, `effcai` from spine `cai_CR`, gmax_AMPA and Use relax
to their ρ-dependent targets with τ_exp.

New signals (all dimensionless except where noted):

* `ca_sh = max(shaft_cai − ca_sh_rest, 0)` (mM): the section Ca from CaDynamics (`USEION ca READ cai`, already
  read by v1 as `shaft_cai`). It is the dendritic Ca-spike / bAP-Ca signal; spine `cai_CR` is not used here,
  so the pre pathways are driven by dendritic (VGCC) Ca, not by NMDAR Ca, as in Ebner and Padamsey.

**Pre-LTD (mGluR → eCB, Ebner E = D·T, Nevian 2006 burst LTD)**

    T' = −T/τ_T + [ca_sh − θ_T]+ / ca_scale              (postsynaptic VGCC-Ca trace, "Ca bound to PLC")
    at each presynaptic spike:  dpre ← dpre − A_mglu · tanh(T)   (mGluR activation coincides with T)

Only post-before-pre (T already elevated when the pre spike arrives) gives LTD; pre-before-post is invisible
to this pathway (consistent with Nevian 2006 and Sjöström 2003).

**Pre NO pathway (Padamsey 2017 / Tong 2020 pre-LTP; Tong 2021 hetLTD)**

    N' = −N/τ_NO + [ca_sh − θ_NO]+ / ca_scale            (NO synthesis needs strong dendritic Ca, θ_NO > θ_T)
    Z' = −Z/τ_Z,     Z ← Z + 1 at each presynaptic spike (presynaptic activity/coincidence trace)
    dpre' += (A_NO · tanh(N) · (tanh(Z) − z_het)) / (1e3 · τ_ind)

`z_het = 0` gives pure pre-LTP (Ebner X = Z·N). `z_het > 0` lets NO depress terminals that were silent
(Tong 2021). Every synapse of a paired connection is active in our protocols, so `z_het` is fixed at 0 for
the first fits and kept only as a hook for later heterosynaptic work.

**Expression**

    dpre ∈ [dpre_min, dpre_max]    (hard clip; default [−0.8, 1.0])
    Use_target = min(1, (Use_d + ρ (Use_p − Use_d)) · (1 + dpre))
    Use_GB' = (Use_target − Use_GB) / (1e3 · τ_exp)       (same filter as v1)

Continuous (not bistable) presynaptic state: Chindemi 2022 notes that presynaptic changes in hippocampus
"support multiple stable states, or evolve on a continuum".

## 4. New parameters (all GLOBAL unless noted)

| name | default | meaning | fit? |
|---|---|---|---|
| `ca_sh_rest_GB` | 6.5e-5 mM | CaDynamics minCai | fixed |
| `ca_scale_GB` | 1e-3 mM | normalisation of the Ca drive | fixed |
| `theta_T_GB` | 1e-4 mM | VGCC-Ca threshold, pre-LTD trace | yes |
| `tau_T_GB` | 50 ms | PLC/eCB trace | yes |
| `A_mglu_GB` | 0 | pre-LTD per pre spike | yes |
| `theta_NO_GB` | 5e-4 mM | dendritic Ca-spike threshold for NO | yes |
| `tau_NO_GB` | 100 ms | NO lifetime | yes |
| `tau_Z_GB` | 10 ms | presynaptic trace (Padamsey 7–10 ms window) | yes (5–30) |
| `A_NO_GB` | 0 | pre NO rate | yes |
| `z_het_GB` | 0 | heterosynaptic NO depression offset | fixed 0 |
| `dpre_min_GB`, `dpre_max_GB` | −0.8, 1.0 | bounds on dpre | fixed |
| `dpre0_GB` (RANGE) | 0 | initial dpre | fixed |

Defaults for the thresholds are placeholders and are set from the recorded `shaft_cai` distribution
(bAP-only vs Ca-spike peaks) before fitting.

## 5. Analytical-method consequences

* **Inputs**: each record needs `effcai` (have it), a decimated `shaft_cai`, and the presynaptic spike times.
  `extract.py` gets a `--with-shaft` flag (T9). The Ebner T5 run records cai_CR, shaft_cai and ica_VDCC.
  For the Markram 10Hz_* delta-prefire set, check whether shaft_cai exists; if it doesn't, a re-run is needed.
* **Exactness**: as for v1, the traces come from prefire runs, where Use/gmax don't feed back into the calcium
  during induction. The offline dpre integrator mirrors §3 exactly and is gated against BCL-v2 prefire (T15).
* **EPSP readout** (T16): the current basis stores ρ∈{0,1} per synapse, moving gmax and Use together. For v2
  the recommended first step is Option B: scale synapse i's basis contribution by its Use_target ratio
  (the C01 test pulses are sparse, so the first-pulse EPSP is ~linear in Use). Check it against a small
  2-D basis (Option A) on a few pairs before trusting it.
* **Objective**: sum over protocols of ((pred − target)/SEM)², for Markram 1997 plus the Ebner set
  (Sjöström 2001 frequency, Nevian 2006, Letzkus 2006; targets from ebner T7).

## 6. Expected qualitative behaviour (test cases for T14/T17)

| protocol | post ρ (v1) | pre dpre (v2) | net |
|---|---|---|---|
| Markram 10 Hz −10 ms | LTD | mGluR LTD (post before pre) | LTD |
| Markram 10 Hz +10 ms | LTP | small (single bAP-Ca, below θ_NO) | LTP |
| Nevian 3 AP 50 Hz, burst before EPSP (Δt′ = −50) | weak | mGluR LTD dominates | LTD |
| Sjöström ≥40 Hz, ±10 ms | LTP | NO pre-LTP (bursts give dendritic Ca) | LTP at both signs |
| Letzkus 3 AP 200 Hz, distal | depends on Ca spike | NO pre-LTP after a Ca spike, gated by Z | timing reversal |

## 7. Open decisions

* Use_target multiplicative (chosen: backward compatible, keeps Chindemi's coupled component) vs. moving Use
  entirely to an independent ρ_pre (cleaner separation, breaks v1 equivalence). Revisit if the fit rails dpre.
* mGluR drive per spike vs. per *released* vesicle. Per spike is chosen: it's deterministic, so the offline
  model stays exact. Per-release is the more biological choice and is kept as a later option.
* Apical/basal (Letzkus proximal vs distal) needs section labels in the extraction (T10).

## 8. Implementation notes (2026-09-27)

* **Offline model** `model_v2.py`:
  - dpre is solved event by event. The filter pass costs about 3 s per 45 Markram records; an amplitude
    change costs about 0.3 s.
  - Gate: `tests/test_offline_gate.py` matches the NEURON mod to ≤0.03 % on the 0.025 ms, 0.25 ms and
    windowed grids.
  - Spike arrival = prespike + edge `delay`.
* **Batch** `batch_v2.py`:
  - Takes records of any grid. Its v1 path reproduces analytical_method Batch: rho_fast and rho_full
    agree to 6e-6, and the EPSP ratios are identical.
  - Conditions: `mglu_block` sets A_mglu = 0; `post_nmdar` freezes rho at rho0.
* **Readout check** (`check_readout.py`):
  - Setting Use directly with the v1 mechanism is invalid, because Use_GB relaxes back to its ρ target
    (τ_exp = 100 s) during the 2 min C01 run. The measured ×0.70 at s = 0.5 is exactly that relaxation.
  - The check therefore uses the drop-in v2 library with dpre0 = s − 1.
  - The C01 ISI is 4 s against Dep ≈ 400 ms, so release should be linear in Use.
  - All-potentiated superposition overestimates the EPSP by about 10 % (6.30 predicted vs 5.73 measured,
    pair 180351-198084). This comes from the v1 basis and affects v1 equally.
* **Letzkus proximal/distal:**
  - Letzkus split inputs by somatic uEPSP rise time (median 2.7 ms), not by distance.
  - In silico: a median split of pairs by the basis EPSP rise time, or by mean synapse path distance.
    Pass it as `loc={pair: "proximal"|"distal"}`.
  - L5TTPC→L5TTPC contacts are mostly basal/oblique, so "distal" pairs may be few.
* **Calibration evidence (Markram 10 Hz delta-cooker traces):**
  - Shaft Ca peaks at the synapse span 0.07–8 µM (median ≈ 0.5 µM) with a single bAP, so θ_NO must sit
    well above the Markram range.
  - With τ_T = 50 ms, T carries over the 90 ms from post to the next pre spike in the +10 ms protocol, and
    mGluR LTD hits LTP protocols too. The Markram +5/+10 targets force τ_T short and/or θ_T high.
* **Fitting strategy:**
  - The post a-params (4) and the pre amplitudes (A_mglu, A_NO) are cheap: an inner DE.
  - The filter parameters (θ_T, τ_T, θ_NO, τ_NO, τ_Z) cost one feature pass per set: an outer coarse
    grid or DE with caching.
