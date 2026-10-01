# Anchoring every non-Chindemi parameter (plan, 2026-10-01)

User rule (2026-10-01): only Chindemi's own rule parameters are fitted to the STDP targets: **a00, a01, a10, a11, γd, γp**.
Every other parameter is either
- **E**: an experimental number from a paper, with preparation, species, age and temperature; or
- **C**: locked by a calibration simulation against data that are *not* fit targets (the model is Chindemi's check of
  the NMDA Ca current against GHK, and of spine Ca transients against Sabatini 2002).

"Fixed at our fitted value" is not allowed. PARAM_REDUCTION.md R1–R3 is rejected (γp "at Chindemi", sloppy
parameters at fitted values); its sources and its identifiability analysis (job 22133878) are reused here.
Reference fit: `results/v4_C1Ajn_s5.json` (χ² 51.98 over 39 targets, 18 free).

Unit conversions used below. They are the model's own spine constants: volume 0.087 µm³, free fraction η = 0.04,
τCa = 12 ms (Sabatini 2002, via Chindemi), i_scale = 1e-5 nA.
- **c_VDCC (θ_VDCC units).** One unit of c_VDCC (= θ_VDCC units) = 1e-5 nA·ms of own-VDCC charge = 31 Ca ions =
  0.60 µM total Ca entering the spine head (about 0.024 µM free).
- **VDCC current to steady free Ca.** A sustained own-VDCC current I gives a steady free spine Ca rise of
  I·η·τCa/(2F·vol) = 0.029 µM per fA.

## 1. Table

Values are best fit (C1Ajn_s5). Confidence (conf.) is H, M or L.

| param | fitted | meaning in measurable terms | anchor | value or calibration | source (DOI) | conf. |
|---|---|---|---|---|---|---|
| a00, a01, a10, a11, γd, γp | 1.13, 1.17, 1.17, 3.00, 51.7, 158 | Chindemi threshold factors and ρ rates | **free** | fitted | 10.1038/s41467-022-30214-w | — |
| **dpre_min** | −0.245 | Maximal fractional drop of U_SE (release probability) that eCB/CB1 can produce, i.e. the CB1-agonist LTD ceiling | **E** | **−0.29 ± 0.03**. Basis: ACEA 125 nM + 15–30 Hz pre firing gives 71 ± 3 % (n 14); AEA gives 79 ± 2 % (n 5, → −0.21). Occlusion: tLTD then AEA 77 ± 8 % vs tLTD 73 ± 3 %, so tLTD sits at the same ceiling. Expression is presynaptic (CV, STD), and the agonist LTD is BAPTA-insensitive (no ρ part), so ratio = 1 + d. Prep: rat V1 L5→L5 pairs, P12–21, 32–34 °C. Readout check C-D. | 10.1016/S0896-6273(03)00476-8 (Figs 4, 8A; text p644) | H (value), M (ACEA vs AEA spread −0.21 to −0.29) |
| **A_mglu** | 0.31 | Fraction of the remaining distance to d_min covered per fully gated own pre spike. 1/A ≈ 3 gated spikes reach 63 % of the maximal eCB-LTD. | **none** | No measurement of tLTD vs number of pairings at L5 (none found). Qualitative only: tLTD is frequency-independent at 0.1–20 Hz (≈0.7, Sjöström 2003 Fig 8B) and occluded by AEA, so the protocols saturate. Options in §3. | 10.1016/S0896-6273(03)00476-8 | — |
| **θTe** | 0.22 | Offset on the event trace S: S must exceed it before eCB production starts. One full-weight event exceeds it for τE1·ln(1/θTe) = 150 ms. | **none** (qualitative E: θTe < 1) | A single bAP is enough for tLTD (L5 0.1 Hz single-AP pairings; L2/3 single-AP tLTD needs Ni²⁺-sensitive VDCC Ca), so θTe < 1. No number. | 10.1016/S0896-6273(03)00476-8; 10.1523/JNEUROSCI.1749-06.2006 | — |
| **τ_T** | 30.4 ms | eCB synthesis/release lag: the onset delay from Ca event to eCB at the terminal | **E (upper bound only)** | eCB synthesis + release complete within 75–190 ms of the stimulus at 22 °C and < 50 ms at 37 °C (rat CA1 → GABAergic, caged Ca/glutamate/AEA). At 32–34 °C this gives τ_T ≲ 50–75 ms; the fitted 30 ms is inside. No lower bound. | 10.1523/JNEUROSCI.2078-05.2005 | M (bound) |
| **θTg** | 6.18 | Gate threshold on T at the own pre spike. Together with θTe and τ_T it sets S* = θTe + θTg/τ_T ≈ 0.43: the gate opens while S > S*. For one full-weight Ca event that is a pre spike ≈ 10–85 ms after the bAP. | **none** | The window itself (−10/−25 LTD, −100/−120 none) is in the fit targets. Independent window: L4→L2/3 eCB-LTD coincidence time scale ≈ 125 ms (rat S1, Bender 2006), which gives S* = e^(−125/100) = 0.29. Different pathway. Options in §3. | 10.1523/JNEUROSCI.0176-06.2006 | — |
| **τE1** (fixed 100 ms) | 100 ms | Memory of the eCB-production drive after a Ca event. It also serves as the memory of the C1 gate sensor. | **E (bracket)** | Synthesis/release lag as for τ_T (75–190 ms at 22 °C, < 50 ms at 37 °C). DSI outlasts the Ca rise (Wang & Zucker 2001). FAAH block widens the L5 tLTD window to −400 ms, so degradation, not Ca decay, ends the eCB signal (Sjöström 2003 Fig 9). DSI decay τ ≈ 20 s (Lenz & Alger 1999) is the CB1 effect lifetime, not this drive. | 10.1523/JNEUROSCI.2078-05.2005; 10.1111/j.1469-7793.2001.t01-1-00757.x; 10.1111/j.1469-7793.1999.00147.x | L–M |
| **K** (fixed, 2.469e-8 nA = 0.2 % of the median bAP VDCC peak) | — | Smallest own-VDCC current excursion that counts as one Ca event (detection floor) | **C (insensitivity check)** | Structural. The event must be VDCC Ca (Ni²⁺, Nevian 2006), and a 250 ms subthreshold step must count (Sjöström 2004; it crosses K at 99.5 % of synapses, DLTD_DIAG). **C-K:** re-score at K/2 and 2K from the extracted cev_lo / cev_hi (K_mult 0.5 / 2, MAXITER 0). If Δχ² < 1, K is declared a non-parameter. | 10.1523/JNEUROSCI.1749-06.2006; 10.1152/jn.00376.2004 | — |
| **A_NO** | 4042 | Rate of the NO-driven rise of U_SE toward dpre_max while NO and pre activity coincide. Per ms at full gate: A_NO/τ_ind = 0.058 ms⁻¹. | **C (conditional)** | **C-NO2:** replay the NO-photolysis protocol: 30–60 single pre spikes at 5 Hz, NO uncaged 7–10 ms after each, ΔPr +0.29 ± 0.07 (n 10). This needs θZ < 1 (below); with the fitted θZ 1.75 a single spike gives no NO-LTP. If ΔPr/Pr0 saturates at dpre_max, it gives only a lower bound. | 10.7554/eLife.29688 (Fig 5F–H) | L (hippocampus) |
| **θNOi** | 1.88e-4 nA (188 fA) | Own-VDCC current above which nNOS makes NO. 188 fA alone would hold free spine Ca at ≈ 5.4 µM. | **C (arithmetic)** | **C-NO1:** nNOS half-activation is 0.2–0.3 µM free Ca (rat brain synaptosomal NOS: 0.3 µM for NO, 0.2 µM for citrulline; NOS1 in cell lysates 200–300 nM). Above the 70 nM rest that is a rise of 0.13–0.23 µM, so **θNOi ≈ 4.5–8 fA**, 25–40× below the fit. Caveat: nNOS is tethered to NMDAR/PSD-95, and its Ca may be a nanodomain value, not this volume average. In CA1 the NO for pre-LTP needs L-VGCC Ca (Padamsey 2017), which supports the VDCC drive. | 10.1111/j.1471-4159.1992.tb10090.x; 10.1074/jbc.271.37.22679; 10.7554/eLife.29688 | M (EC50), L (mapping) |
| **τ_NO** | 6.99 ms | NO lifetime at the synapse (consumption) | **E** | **τ_NO = 6.7 ms (5–10)**. NO inactivation by tissue is Michaelis–Menten with Vmax 1–2 µM/s and Km ≈ 10 nM. At ≤ 10 nM it is first order, k = Vmax/Km = 100–200 s⁻¹, τ = 5–10 ms, t½ < 10 ms. Prep: rat cerebellar slices (diffusion–inactivation fit to cGMP); temperature not checked. | 10.1113/jphysiol.2006.118380 | M–H |
| **τ_Z** | 48.7 ms | Memory of the terminal for its own spikes when NO arrives | **E (lower bound only)** | Pre spikes 7–10 ms *before* NO photolysis potentiate (+0.29). Pre spikes 7–10 ms *after* NO do nothing (−0.01 ± 0.04, n 8). So τ_Z ≳ 10 ms and NO is short-lived. Prep: rat CA3–CA1 organotypic (P7–8, 7–14 DIV) and acute slices, 31–33 °C. | 10.7554/eLife.29688 | L (prep) |
| **θZ** | 1.75 | Number of recent own pre spikes needed for NO to act. The fitted 1.75 means ≥ 2 spikes within ≈ τ_Z. | **E (hippocampus) or none** | One pre spike per pairing suffices for NO-LTP (60 pairings at 5 Hz, ΔPr +0.19 to +0.38), so θZ < e^(−10/τ_Z). Simplest: θZ = 0, as in the Tong 2020 form tanh(Z). At L5, NO-LTP was tested only with high-frequency pairing (Sjöström 2007), so a ≥ 2-spike rule is neither tested nor contradicted there. | 10.7554/eLife.29688; 10.1016/j.neuropharm.2006.07.021 | L |
| **θ_VDCC** (JSON theta_V) | 5.50 (≈ 170 Ca ions ≈ 3.3 µM total entry per 100 ms) | Own-spine VDCC Ca (100 ms sensor memory) needed for potentiation | **C failed → none for this emodel** | **C-V (job 22135045) RESULT: no lock.** og-delta has no 200 Hz supralinearity at basal L2/3->L5 synapses (dV3/(dV1·linear sum) median 0.86 at 200 Hz vs 0.87 at 50 Hz; Kampa: supralinear only at 200 Hz). Youden J ≈ 0.02 (S&H 50 Hz) and 0.00 (L5 Sjöström 50 Hz) at every level; at 5.50, P(open) is 0.44 for 200 Hz vs 0.52–0.63 for 50 Hz. Basal dV3 spans q10 0.05 to q90 100 across synapses, so geometry, not frequency, sets the gate. Intended design: lock θ_VDCC between the c_VDCC of 3 bAPs at 50 Hz and 3 bAPs at 200 Hz at basal L5 synapses. Basis: in L5 basal dendrites, 200 Hz AP bursts give supralinear Ca that low Ni²⁺ blocks, and LTP, while 50 Hz bursts give neither (Kampa 2006). In the apical dendrite the critical frequency is 60–200 Hz (Larkum 1999). Kampa's LTP outcome is not a fit target, and the lock uses the Ca observation. | 10.1113/jphysiol.2006.111062; 10.1073/pnas.96.25.14600 | none: the model cannot reproduce the observable |
| **rho_gamma** | 0.52 | Exponent of the c_post scaling of θd/θp: how the thresholds follow the single-bAP Ca at each synapse | **none** | No experimental counterpart (Chindemi and GB use 1). Options in §3. | — | — |
| **C_REF** (fixed 0.01 mM) | — | Reference Ca at which c_q = c_post | **none** | It only exists when γ ≠ 1; it drops out at γ = 1. | — | — |
| **i_scale** (fixed 1e-5 nA) | — | Units constant | **not a parameter** | C1 gate: only θ_VDCC·i_scale enters, so θ_VDCC is reported as charge. eCB (t_drive 4): unused. NO: it sets the N scale inside tanh(N), so it merges with θNOi. Lock both by C-NO1, or replace tanh(N) by a Hill of the VDCC-driven Ca at the nNOS EC50 so that i_scale disappears. | — | — |

Hidden constants that also need a line in the paper:
- dpre_max = 1.0 (largest NO-driven fractional rise of U_SE). E-candidate: Padamsey 2017 ΔPr +0.29 to +0.38 at
  low-Pr synapses. Baseline Pr must be read from their Fig 5G before it is converted to a fraction.
- τ_glu = 70 ms (own glutamate-bound state). This is the synapse model's NMDAR decay, so it is already anchored through
  Chindemi's synapse calibration.
- ca_hyst 0.5, Te_scale 1, θ_N 0: structural (one event per excursion; absorbed into θTg / tanh scale).
- τ_ind, ρ*, τ*: Chindemi.

## 2. Calibrations, ranked (value / cost)

Each cost below is sized from a measured job named in brackets.

1. **C-V θ_VDCC vs Kampa 2006. DONE, job 22135045 (0:53, 1.16 GB of 1.17, CPU 60 %): NULL. No 200 Hz supralinearity (0.86 vs 0.87 at 50 Hz), J ≈ 0, so θ_VDCC cannot be locked with og-delta. Options: (i) an emodel with a basal Ca spike or supralinear burst Ca (new name, new mod + emodel dir), then rerun C-V; (ii) flag θ_VDCC as fitted (tight, 4.8–6.3); (iii) drop the C1 gate. C-V2 would not help (the bAP-only 200 Hz trace is already in the data).** Original design: Script: `run_calib_vdcc_kampa.sh`, with the numpy code
   inline. Output: `/scratch/dhuruva/param_anchors/calib_vdcc_kampa_syn.csv` and the log
   `logs/calib_vdcc_kampa_22135045.out`.
   - Request: 1 CPU, 1.2G, 0:15 (from diag_dltd 22132795: 0:12, 910 MB), ≤ 0.25 CPU·h.
   - Inputs: existing extracted −ica_VDCC on L5 TTPC dendrites:
     - Letzkus 1 AP;
     - Letzkus 3 AP @ 200 Hz at −500 ms, which is bAP-only;
     - S&H 50 Hz (same L2/3→L5 synapses);
     - Sjöström 0.1 Hz and 50 Hz on L5→L5.
   - Observable: the increment of c_VDCC after the 1st, 3rd and 5th AP of each burst, per synapse.
   - Lock: θ* maximises P(200 Hz open) − P(50 Hz open) at basal synapses (Youden). The job also reports the
     supralinearity dV3 / (dV1 × linear sum).
   - Check afterwards:
     - (a) Is the model supralinear at 200 Hz and linear at 50 Hz, as Kampa found? If not, og-delta lacks the basal
       Ca spike, and θ* is set only by the τE1 summation difference (weak; record it as a model limit).
     - (b) θ* and its J ≥ max − 0.05 interval vs the fitted 5.50 (the job also prints P(open) at 5.50).
     - (c) The 50 Hz traces carry EPSP Ca (S&H +10, Sjöström −10). If that confound matters, run C-V2.
   - Value: high. θ_VDCC is the tightest new parameter (×1.22 across fits) and carries the L2/3 transfer.
2. **C-τNO and C-NO1. Arithmetic, done in the table, 0 cost.**
   - τ_NO = 6.7 ms (H&G 2006).
   - θNOi = 4.5–8 fA (nNOS EC50).
   - Their effect is only known after the validation refit (item 7).
3. **C-D dpre_min readout check.**
   - Request: 1 CPU, 1.2G, 0:15 (diag_dltd basis).
   - Method: with ρ frozen at ρ0 and d = −0.29 at every synapse, compute the EPSP ratio from the existing per-pair EPSP
     basis (model_v2.epsp_ratio_v2, the fit readout).
   - Pass: the ratio is 0.71 ± 0.01. The U_SE cap and MVR binomial could bend it; if it misses, lock dpre_min to the
     d that gives 0.71.
4. **C-K insensitivity of K.**
   - Request: two MAXITER-0 evaluations (SET K_mult 0.5 and 2: cev_lo / cev_hi already extracted), each on a 2g MIG
     slice, 1 CPU, 49G, 0:15 (from 22133806: 2:41, 39.05 GB).
   - Outcome: if Δχ² < 1 both ways, K is reported as a detection floor with no fitted content.
   - Value: medium.
5. **C-NO2 A_NO vs NO photolysis (Padamsey 2017).**
   - Method: offline replay of the dpre ODE with an N pulse (saturating, decay τ_NO) and Z from single spikes 7–10 ms
     before it, 30–60 pairings at 5 Hz. Match ΔPr/Pr0.
   - Request: 1 CPU, 1G, 0:15.
   - Preconditions: θZ < 1 must be adopted, and baseline Pr must be read from Fig 5G.
   - Value: low–medium (hippocampal prep).
6. **C-V2 bAP-only burst frequency series. NOT NEEDED: C-V is null on the bAP-only (−500 ms) 200 Hz trace itself.**
   - Protocol: 3 APs at 50, 100, 150 and 200 Hz, 10 bursts at 0.1 Hz, no pre. Run on the 24 subset L5→L5 pairs with
     og-delta, with new protocol ids in a new yaml.
   - Simwriter: 2G, 0:15 (22127910: 1:39, 1.43 GB).
   - Prefire: 96 runs of ≈ 1/5 of sj04's 500 s bio, so ≈ 60–80 s each.
     - 8 workers × 5–6 GB → 50G, 0:30 (22127911: 5–6 GB/worker, 261–390 s per 500 s-bio run), ≈ 2–3 CPU·h.
   - Extraction: size from the seff of 22128325.
   - Value: it removes the EPSP confound and finds the model's own critical frequency.
7. **Validation refit after the locks** (not run: no plasticity fits in this task).
   - Seeded from C1Ajn_s5, with only the 6 Chindemi parameters plus whatever §3 leaves flagged.
   - Request: 2g MIG, 1 CPU, 50G, 0:45 (22132568: 39.5 GB, 25 min), plus one unseeded basin check.
   - Nested variants for the §3 options: 3–5 refits, ≈ 2–4 MIG-h.

## 3. No anchor: options for the user

| param | option A (remove / simplify) | option B | recommendation |
|---|---|---|---|
| A_mglu | A_mglu = 1. Each fully gated pre spike moves d all the way to d_min, graded only by tanh(T − θTg). Justification: occlusion plus frequency independence say tLTD saturates. | Keep as a declared fitted rate: the eCB analogue of γd, a presynaptic rate with no Chindemi counterpart. | B for now, flagged. A is outside the DE range (0.11–0.76), so test it in the item-7 nested refit. |
| θTe, τ_T, θTg | Collapse to one threshold on S: gate = tanh-free step H(S − S*), with S* = exp(−W/τE1). τ_T takes the Heinbockel bound (≤ 50 ms), or the T stage is dropped. W comes from the independent L4→L2/3 window (≈ 125 ms, Bender 2006), giving S* = 0.29. This needs a kernel variant in a new file (gpu_v4_rho.py is not edited). | Keep S* (one combination) as a flagged fitted parameter. θTe and τ_T are then pure reparameterisation, so fix θTe = 0 and τ_T at the bound. This removes two parameters without using fitted values. | B: 3 → 1 flagged. Option A depends on accepting another pathway's window under the uniform rule. Ask the user. |
| rho_gamma, C_REF | γ = 1 (Chindemi's form), so C_REF drops out. The unseeded C1Ajn_s6 found γ 0.98 at χ² 63.8, against 52.0 seeded at 0.52. | Keep γ flagged as fitted. | Test γ = 1 in the item-7 refit. If Δχ² ≤ ~4 + 2 (AIC for one parameter), adopt it. |
| θZ (if the hippocampal anchor is rejected) | θZ = 0 (Tong 2020 form) | Flag as fitted ("≥ 2 spikes", untested at L5) | Ask the user: take the hippocampal E value or flag it. |

**Resulting fitted set.**
- *Strict* (every lock above accepted, options B): a00, a01, a10, a11, γd, γp plus flagged A_mglu, S*, rho_gamma, θ_VDCC (C-V null),
  and θZ if the hippocampal value is rejected. That is **10–11 free**, down from 18.
- *If the §3 simplifications A are accepted:* **6 free**, the Chindemi set.

## 4. Risks of each lock

Judged from the DE population (Δχ² ≤ 4 ranges, PARAM_REDUCTION §2a) and the C1Ajn_s5 residuals.

| lock | inside the Δχ² ≤ 4 range? | targets likely to change |
|---|---|---|
| dpre_min −0.29 | No (−0.267 to −0.215) | The eCB-LTD misses should *improve*: 0.1 Hz −10 0.80 vs 0.69, −25 0.77 vs 0.65, 10 Hz −10 0.78 vs 0.57, 20 Hz −10 0.92 vs 0.65. Risk: burst −120 / −200 (0.77 / 0.82 vs 0.79) and sj07 pre-only (0.90 vs 1.00) over-depress. Net small. Use −0.21 (AEA) as the sensitivity bound. |
| τ_NO 6.7 ms | Yes (3.8–11.8) | Negligible; it trades with γp, which stays free |
| θNOi 4.5–8 fA (×1/25–1/40) | No | **High.** NO-LTP in every protocol with large VDCC Ca and pre activity. 20 Hz −10 control (0.92 vs 0.65) and its AM251 arm (1.20 vs 1.02, already z 2.6) get worse; 50 Hz −10 (1.22 vs 1.70) improves. A_NO would have to fall. Likely needs the Hill reformulation, or a decision that nNOS reads nanodomain Ca (EC50 not transferable). |
| θZ → 0 | No (1.37–2.14 in every fit) | **High.** Single-spike protocols gain NO-LTP: 0.1 Hz +10 (1.06 vs 0.97, already z 2.3), Letzkus 1AP +10. τ_Z ↔ θZ correlate (+0.62), so τ_Z cannot compensate if it is also locked. |
| τ_Z ≥ 10 ms | Bound only | none |
| θ_VDCC | C-V gave no lock (above), so it stays flagged | L2/3 +10 targets depend on it (they are the reason for C1). Dropping the gate returns the pre-C1 transfer failure (χ² 502.9 over 9). |
| A_mglu = 1 | No (0.15–0.54) | Faster LTD per pairing. Pairs with θTg; if S* stays free it can partly compensate. |
| S* from W = 125 ms | Effective S* now ≈ 0.43 (W ≈ 85 ms) | A wider window: 0.1 Hz −100/−120 (no-LTD targets) may depress; burst −120 / −200 gain LTD |
| γ = 1 | No (0.49–0.64) | L2/3 +10 targets (the pre-gate td4_s2 transfer gave χ² 502.9 over 9); a11 compensates (r −0.43) |

The eCB and NO parameters have no measurement in L5 or L2/3→L5 pairs other than the fit targets. Most E values come
from other preparations:
- hippocampus: Heinbockel, Padamsey, Wang & Zucker;
- cerebellum: Hall & Garthwaite;
- enzymology: the nNOS EC50.
Each such value carries the preparation listed in the table.

## 5. Calibration results (coordinator decisions 2026-10-01; lines appended by the jobs)

Accepted locks: dpre_min −0.29 (subject to C-D below), τ_NO 6.7 ms, θNOi from the nNOS EC50, i_scale = units.
Test variants (ladder): θZ = 0 and A_mglu = 1. S* stays fitted, flagged as the one unanchored eCB parameter; S* = 0.29
from Bender 2006 needs the user's OK. Jobs: C-D and C-NO2 run `run_calib_anchors.sh` (MODE=CD / NO2). C-K runs
`run_fit_v4.sh` with MAXITER=0, all 18 C1Ajn_s5 values fixed, K_mult 0.5 / 2 (cev_lo / cev_hi), compared with
χ² 51.9757 at K_mult 1. Submitted: C-D 22135560 (1 CPU, 1500M, 0:15), C-NO2 22135563 (1 CPU, 1G, 0:15), C-K 22135619 (K/2) and 22135620 (2K) (2g MIG, 49G, 0:15, def-emuller). Fill seff here.
- C-K K_mult 2.0 (job 22135620, MAXITER 0, all 18 fixed at C1Ajn_s5; baseline chi2 total 51.9757 at K_mult 1): chi2 L5 47.49 over 30 targets (gamma_d 51.704, gamma_p 158.394, vamp (1, 5.500694980600818));chi2 L2/3->L5 5.18 over 9 targets;chi2 total 52.6724 over 39 targets, 0 free;
- C-K K_mult 0.5 (job 22135619, MAXITER 0, all 18 fixed at C1Ajn_s5; baseline chi2 total 51.9757 at K_mult 1): chi2 L5 54.38 over 30 targets (gamma_d 51.704, gamma_p 158.394, vamp (1, 5.500694980600818));chi2 L2/3->L5 5.38 over 9 targets;chi2 total 59.7662 over 39 targets, 0 free;
- C-D (22135560; 0:05, 248 MB; the compute finished, only the flock append failed, so this line was added by hand): mean ratio over the 24 L5 pairs is 0.719 / 0.765 / 0.800 at d −0.29 / −0.245 / −0.21. d giving 0.71 (ACEA) = −0.299, d giving 0.79 (AEA) = −0.220. **dpre_min locked at −0.30** (agrees with the paper anchor −0.29).
- C-NO2 (22135563; 0:02, 21 MB; append failed as above): with the fitted θZ 1.75, a single pre spike gives no NO-LTP (gate integral 0), which contradicts Padamsey 2017. With θZ 0, the A_NO giving ΔPr +0.29 after 60 pairings is 32–1123 depending on τ_Z (10 / 48.7 ms), N0 and Pr0 0.3–0.6. The fitted 4042 saturates d = 1. So A_NO is not locked by this calibration (range spans 35×), which supports dropping NO (MECH_NECESSITY §NO).
- C-K: K/2 → χ² 59.77 (+7.79), 2K → 52.67 (+0.70). K is not a non-parameter at the low side; the detection floor matters when lowered. Keep K at its v4 value as the detection-floor definition and flag it; v5 has no K.
- Fix for the job scripts: `flock` with `-c` takes one string argument; quote the whole command.
