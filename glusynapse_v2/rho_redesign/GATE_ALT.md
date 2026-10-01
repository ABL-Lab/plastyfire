# Replacing the LTP licence θ_V (V7vn2 W2_N3) with a quantity experimentalists image (A13, 2026-10-01)
Today: potentiation counts only while Vg (own spine VDCC influx, leaky τE1 100 ms) ≥ θ_V 3.58 ≈ 110 Ca ions; fitted, no anchor.
Its job (EXP_PARAM §3a, §4): an **absolute location floor**. Proximal spines pass with ~1 bAP, distal apical spines need bursts
(Letzkus 200 Hz +10 distal must stay 0.79). A per-synapse relative gate removes the floor and lost (E2: χ² 97 vs 88).
So the replacement must be an absolute, synapse-local Ca quantity that falls steeply with distance for single bAPs and
recovers with bursts/Ca spikes. Literature via PubMed, DOI links below.
## Literature anchors (via PubMed)
- Single bAP, L5, fura-2: peak Δ[Ca] 267 ± 109 nM basal, 128 ± 25 nM proximal apical, 43 nM soma; 5 APs at 10-12 Hz: 710 nM
  proximal apical; amplitude falls toward the bifurcation (Schiller 1995, [doi](https://doi.org/10.1113/jphysiol.1995.sp020902)).
- Dye-free extrapolation, proximal apical L5/CA1: single-AP transients **150-300 nM**, decay < 100 ms at 37 °C; during trains
  [Ca] reaches a plateau (τ < 200 ms) linear in AP frequency (Helmchen 1996, [doi](https://doi.org/10.1016/S0006-3495(96)79653-4)).
- Distal: single APs give little or no tuft Ca in L2/3 (Waters 2003, [doi](https://doi.org/10.1523/JNEUROSCI.23-24-08558.2003)); L5 bursts above
  a critical frequency (60-200 Hz) evoke Ca spikes with the largest Ca 400-700 µm from soma (Larkum 1999, [doi](https://doi.org/10.1073/pnas.96.25.14600));
  L2/3 critical frequency 130 Hz (Larkum 2007, [doi](https://doi.org/10.1523/JNEUROSCI.1717-07.2007)); distal L5 basal: 200 Hz bursts supralinear,
  ~100 Hz critical (Kampa & Stuart 2006, [doi](https://doi.org/10.1523/JNEUROSCI.3062-05.2006)).
- LTP requirement: L5→L5 LTP needs 200 Hz bursts (single APs, 50 Hz bursts fail); Ni²⁺ blocks the supralinear basal Ca and LTP
  (Kampa 2006, [doi](https://doi.org/10.1113/jphysiol.2006.111062)). Distal L2/3→L5: single-AP +10 shifts to LTD with distance, bursts/BAC
  depolarisation restore LTP (Letzkus 2006, [doi](https://doi.org/10.1523/JNEUROSCI.2650-06.2006)). Distal sign set by bAP spread, switched by
  cooperativity (Sjöström & Häusser 2006, [doi](https://doi.org/10.1016/j.neuron.2006.06.017)). Distal CA1 LTP needs dendritic spikes (Golding
  2002, [doi](https://doi.org/10.1038/nature00854)); single-burst LTP needs dendritic spikes + L-type (Remy & Spruston 2007, [doi](https://doi.org/10.1073/pnas.0707919104)).
  Distal basal synapses fail to potentiate with global activity (Gordon 2006, [doi](https://doi.org/10.1523/JNEUROSCI.3502-06.2006)).
- Against spine total Ca: equal volume-averaged spine Ca gives LTP or LTD (Nevian & Sakmann 2006, [doi](https://doi.org/10.1523/JNEUROSCI.1749-06.2006));
  CA_DECODE: spine ca_pk / ca_int do not separate (AUC 0.56 / 0.35 vs vd_int 0.20). Photolysis: brief high Ca → LTP, long modest → LTD
  (Yang 1999, [doi](https://doi.org/10.1152/jn.1999.81.2.781)), no µM gate value usable here.
## Candidates
| candidate | local? | experimentalist units | location dependence (single bAP / burst) | anchor (value, source) | risk |
|---|---|---|---|---|---|
| 1 shaft peak Δ[Ca] at synapse, hold 100 ms | yes (shaft Ca at site) | µM, "dendritic Ca transient per AP" | basal/prox apical 0.13-0.3 µM per bAP; tuft ≈ 0; bursts/Ca spike → µM | **0.15 µM** = one proximal bAP (Helmchen 1996 lower bound; Schiller 1995) | model shaft Ca scale/decay may differ from imaging; shaft cai may include synaptic Ca if the mod writes ica |
| 1b shaft mean Δ[Ca] (100 ms low-pass) | yes | µM plateau, "Ca vs AP frequency" | same, plus frequency code (plateau ∝ f) | single-bAP G_peak = 0.15 µM × f(τc) → needs model τc; calibrate | anchor depends on model shaft decay |
| 2 spine total [Ca], peak or ∫ | yes | µM in spine | NMDA-dominated (99.8 % of EPSP effcai) → tracks c*, weak location filter | none; fails CA_DECODE and Nevian 2006 | redundant with θp; 0.070 µM floor (EXP_PARAM §1a) |
| 3a c* alone (higher θ) | yes | Cpre/Cpost units | none beyond c* (NMDA distal Ca is large) | none | cannot keep Letzkus distal +10 LTD (NMDA-dependent, APV 1.03) |
| 3b c* > θp AND dendritic Ca event | yes, via 1 | "θp AND ≥ X µM local Ca" | = candidate 1 | = candidate 1 | this is the proposed combination, not a separate quantity |
| 4 bAP / burst count | **no** (cell-level) | spikes, Hz | none per se | Kampa 2006 (200 Hz) | forbidden as a drive; allowed only as read through 1/5 |
| 5 Ca spike / plateau via shaft Ca | yes | µM, "supralinear burst Ca" | single bAP fails everywhere; 200 Hz burst passes basal-distal and apical | ≈ 0.5 µM: above ~2 summed proximal bAPs, below 5-AP train 0.71 µM (Schiller 1995) | flagged (bracket, no lock); may kill L5 20-50 Hz LTP (Sjöström 2001) |
## Ranked top 3 (all keep the V7vn2 combination: licence L(t) multiplies the potentiation lane; a θp crossing without L counts as depression)
Notation: c_sh(t) = cai of the dendritic segment that carries the synapse (shaft, not spine), Δc_sh = c_sh − c_sh,rest, with c_sh,rest
the per-synapse value just before the first stimulus. pot(t) = H(c* − θp) · L(t). θp, θd, γ, eCB unchanged.
1. **G1 peak shaft Ca (recommended).** L(t) = H( max_{s ∈ [t − T_h, t]} Δc_sh(s) − θ_sh ), T_h = 100 ms (τE1 bracket, Heinbockel 2005 / Bender 2006),
   θ_sh = **0.15 µM** above rest (anchored: physiological single-AP dendritic transient, Helmchen 1996 150-300 nM; Schiller 1995 basal 267 nM).
   Experimentalist reading: "LTP needs integrated spine Ca above θp *and*, within the last 100 ms, a dendritic Ca transient at the synapse at least as
   large as one bAP gives near the soma (≥ 0.15 µM)." Trains integrate through the shaft's own Ca decay. If anchored, k drops 8 → 7.
   Fallback if model proximal single-bAP Δc_sh is outside 0.15-0.30 µM: θ_sh = model median single-bAP Δc_sh peak at synapses ≤ 100 µm path
   distance (independent calibration, like c_post; uniform). Validation run: θ_sh free once (flagged), check it lands in 0.1-0.3 µM.
2. **G2 mean shaft Ca.** dG/dt = (Δc_sh − G)/τG, τG = 100 ms (= τE1); L(t) = H(G − θ_G), θ_G in µM. G at steady state = mean Δ[Ca] (Helmchen
   1996 plateau ∝ AP frequency), so it is the closest analogue of Vg. θ_G = G_peak of one bAP at the proximal calibration set (independent
   calibration; for exponential decay τc, G_peak = 0.15 µM · r^{1/(1−r)}, r = τc/τG, e.g. 0.25 × 0.15 µM at τc 50 ms). Use if G1 fails on 20-50 Hz trains.
3. **G5 dendritic Ca event.** Same as G1 with θ_sh = 0.5 µM (bracket 0.3-0.7: above 1-2 summed proximal bAPs, below a 5-AP train, Schiller 1995;
   Kampa 2006 single APs and 50 Hz bursts fail, 200 Hz supralinear passes). **Flagged** (no lock). Prediction: Sjöström 0.1 Hz +10 (z +3.4) improves;
   risk: Sjöström 2001 20-50 Hz LTP and Markram 10 Hz +10 lose LTP. Score as a fixed-param rescore against G1 only.
## What A14 must find in the records (or add)
- **Shaft cai at the synapse segment** (sec(x).cai of the parent dendrite at the synapse location), same grid as effcai (0.25 ms spike windows);
  per-synapse resting value before the first stimulus.
- Whether GluSynapse spine Ca (NMDA/VDCC) is written into segment cai (mod WRITE ica / USEION ca). If yes, shaft cai carries synaptic Ca and must
  be flagged (it then mixes NMDA into the gate).
- Emodel CaDynamics per region (decay, gamma, depth) for og-delta / antic-delta, and dendritic Ca channel presence in basal, apical trunk, tuft;
  a flat shaft cai anywhere kills G1/G2 there.
- Calibration run (cexp-like): single isolated bAP peak Δc_sh and decay τc per synapse, with path distance and basal/apical/tuft label; report
  q10/50/90 per pathway and the proximal (≤ 100 µm) median against 0.15-0.30 µM. Also peak Δc_sh for a 3AP 200 Hz burst (Letzkus distal sites).
- Correlation of single-bAP Δc_sh with vdcc_q_post (cexp) across synapses: if ρ_Spearman ≳ 0.9, G1 should reproduce V7vn2's location filter.
