# Reducing the free parameters of the joint fit (from C1Ajn_s5)

Start: `results/v4_C1Ajn_s5.json`, χ² 51.98 over 39 targets (L5 46.89/30, L2/3→L5 5.08/9), 18 free.
Fixed already: τ_ind 70 s, ρ* 0.5, τ* 278.3 ms, τ_E1 100 ms, i_scale 1e-5, C_REF 0.01 mM.

## 1. Literature anchors (PubMed / Europe PMC survey; every DOI checked)

| param | fitted | literature value / range | preparation | source |
|---|---|---|---|---|
| d_min (max fractional U_SE drop, eCB-LTD) | −0.245 | tLTD 0.73 ± 0.03, expressed presynaptically (CV, STD index) → ≈ −0.27; Sjöström 2007 split (ctrl 1.62, NO block 1.36, AM251 2.13) → eCB factor 0.76 → −0.24; dLTD 0.65–0.70 | rat V1 L5→L5 pairs, P12–21, 32–34 °C | 10.1016/S0896-6273(03)00476-8; 10.1016/j.neuropharm.2006.07.021; 10.1152/jn.00376.2004 |
| τ_NO (NO trace) | 6.99 ms | NO t½ < 10 ms in tissue at ≤ 10 nM (consumption Vmax 1–2 µM/s, Km ≈ 10 nM) → τ ≤ 14 ms; older 0.5–5 s values are superseded | rat cerebellar slices (diffusion-model fit); reviews | 10.1113/jphysiol.2006.118380; 10.1016/j.niox.2009.07.002; 10.1111/j.1460-9568.2008.06285.x |
| γp (ρ potentiation rate) | 158.4 | Chindemi 2022 216.2 (same c* family, τ* 278 ms, τ_ind 70 s); γd/γp 0.13–0.54 across GB2012/2016, Deperrois 2020 (GAMMA_LIT.md) | in silico, fitted to rat L5 slice STDP | 10.1038/s41467-022-30214-w; GAMMA_LIT.md |
| γd | 51.7 | Chindemi 101.5; Deperrois V1 31.98 (γ/τ 0.16–3.95 s⁻¹) | model fits | GAMMA_LIT.md |
| τ_E1 (eCB drive, already fixed 100 ms) | — | L5→L5 tLTD at −10/−25, not at −120/−200 ms, rescued by FAAH block; L4→L2/3 LTD timescale ≈ 125 ms | rat V1 L5; rat S1 barrel | 10.1016/S0896-6273(03)00476-8; 10.1523/JNEUROSCI.0176-06.2006 |
| τ_Z, θZ (pre spike count for NO-LTP) | 48.7 ms, 1.75 | Hippocampus: one pre spike leading NO by 7–10 ms suffices (ΔPr +0.29); L5 NO-LTP only tested with 40–60 Hz bursts, so a ≥ 2-spike requirement is untested, not contradicted | rat CA3–CA1 organotypic; rat V1 L5 | 10.7554/eLife.29688; 10.1016/j.neuropharm.2006.07.021 |
| θ_VDCC (spine VDCC Ca gate) | 5.50 | No absolute value. Requirement measured: L5→L5 basal LTP needs 200 Hz bursts with supralinear Ni²⁺-sensitive Ca (50 Hz fails); distal L2/3→L5 LTP abolished by 100 µM Ni²⁺ (0.99), LTD not (0.81) | rat L5 / L2/3→L5 slices, 32–36 °C | 10.1113/jphysiol.2006.111062; 10.1523/JNEUROSCI.2650-06.2006 |
| a00, a01, a10, a11 | 1.13, 1.17, 1.17, 3.00 | Chindemi basal b00 1.002, b01 1.954, b10 1.159, b11 2.483 (apical 1.127, 2.456, 5.236, 1.782): fitted model quantities, units of the c* peaks of isolated pre / post spikes | in silico | 10.1038/s41467-022-30214-w |
| ρ-γ (c_post exponent) | 0.52 | none; Chindemi and GB use 1 (linear). Indirect: distance-dependent LTP→LTD shift rescued by bursts / depolarisation | rat L5, L2/3→L5 | 10.1016/j.neuron.2006.06.017; 10.1523/JNEUROSCI.2650-06.2006 |
| A_mglu, θTe, θTg, τ_T (eCB gate) | 0.31, 0.22, 6.2, 30 ms | none. Constraint only: tLTD is frequency-independent (0.1 Hz singles 0.65 vs 20–30 Hz bursts 0.70), so one post→pre pairing must cross θTg without summing pairings | rat L5 pairs | 10.1152/jn.00376.2004; 10.1146/annurev.physiol.010908.163149 |
| A_NO, θNOi | 4042, 1.9e-4 | none (model units). The NO-LTP size is anchored instead: ≈ ×1.19 (1.62/1.36) | rat V1 L5 | 10.1016/j.neuropharm.2006.07.021 |

Measurable re-expressions (not used yet): θ_VDCC could be set between the filtered own-VDCC drive of 3 APs at 50 Hz and
3 APs at 200 Hz at a proximal basal L5 spine (Kampa 2006); the eCB gate could be set so a single 0.1 Hz post→pre
pairing just crosses θTg. Both would need a scan of the extracted traces, not a literature number.

## 2. Identifiability from the fits

Analysis job 22133878 (1 CPU, numpy; `/scratch/dhuruva/param_reduction/ident.py`; seff: 39 MB, <1 s).
- **DE population (C1Ajn_s5 ckpt, 144 members, nit 150).** 12 / 33 / 126 members lie within Δχ² 2 / 4 / 10 of the best.
  The late population has partly collapsed, so these widths are lower bounds on the true sloppiness.
- **Tight (width ≤ 0.07 of the box at Δχ² ≤ 4, and stable across fits):** a00 (1.11–1.18), a10 (1.15–1.23),
  θ_VDCC (4.8–6.3), γd (46–62), d_min (−0.267 to −0.215).
  The stiffest direction is a00 − a10, i.e. the c_pre parts of θd and θp are pinned against each other.
- **Sloppy (Δχ² ≤ 4 spans ×2–5):** A_NO 1675–8985, τ_NO 3.8–11.8 ms, A_mglu 0.15–0.54, θTe 0.11–0.34, τ_T 19–50 ms,
  γp 126–208, a01 0.67–1.41, τ_Z 40–74 ms, θZ 1.56–2.10, ρ-γ 0.49–0.64.
- **Largest-variance directions (sloppy combinations):** τ_NO with γp; θZ with τ_Z (r = +0.62: a longer count trace
  with a higher count threshold gives the same gate); ρ-γ with γp. eCB: θTe–θTg (−0.49) and θTe–τ_T (+0.48), so the
  eCB trace offset, gate and time constant trade against each other. a11–ρ-γ (−0.43) is the expected cq scaling.
- **Across fits** (best values; same 39-target data: Ajn_s5 seeded, Ajn_s6 unseeded; Ajd 36 targets; Aj4/Ajs + sj04):

| param | Ajn_s5 | Ajn_s6 | Ajd_s5 | Ajd_s3 | Ajd_s4 | Aj4_s7 | Ajs_s5 | spread |
|---|---|---|---|---|---|---|---|---|
| a00 | 1.127 | 1.248 | 1.125 | 1.146 | 1.165 | 1.178 | 1.146 | ×1.11 |
| a01 | 1.168 | 0.579 | 1.268 | 1.035 | 0.877 | 0.371 | 0.886 | ×3.4 |
| a10 | 1.171 | 1.298 | 1.259 | 1.209 | 1.240 | 1.214 | 1.219 | ×1.11 |
| a11 | 3.005 | 1.844 | 2.902 | 2.909 | 2.398 | 3.001 | 3.043 | ×1.65 |
| A_mglu | 0.308 | 0.158 | 0.401 | 0.345 | 0.229 | 0.131 | 0.449 | ×3.4 |
| A_NO | 4042 | 5311 | 5197 | 2563 | 1254 | 1867 | 6364 | ×5.1 |
| θTe | 0.225 | 0.370 | 0.165 | 0.165 | 0.528 | 0.040 | 0.125 | ×13 |
| τ_T (ms) | 30.4 | 45.5 | 28.8 | 27.3 | 17.9 | 21.0 | 26.2 | ×2.5 |
| θTg | 6.18 | 3.99 | 6.66 | 6.36 | 2.54 | 7.81 | 7.21 | ×3.1 |
| d_min | −0.245 | −0.244 | −0.230 | −0.225 | −0.267 | −0.260 | −0.241 | 0.041 |
| θNOi | 1.88e-4 | 4.22e-4 | 2.35e-4 | 1.85e-4 | 2.31e-4 | 7.7e-5 | 2.67e-4 | ×5.5 |
| τ_NO (ms) | 6.99 | 12.8 | 9.02 | 6.70 | 14.6 | 4.67 | 8.58 | ×3.1 |
| τ_Z (ms) | 48.7 | 49.7 | 42.0 | 51.8 | 41.7 | 68.6 | 57.3 | ×1.65 |
| θZ | 1.75 | 1.73 | 1.45 | 1.52 | 1.37 | 2.14 | 1.66 | ×1.56 |
| θ_VDCC | 5.50 | 6.05 | 6.00 | 5.93 | 6.12 | 6.15 | 6.69 | ×1.22 |
| ρ-γ | 0.523 | 0.977 | 0.454 | 0.599 | 0.754 | 0.575 | 0.527 | ×2.2 |
| γd | 51.7 | 89.4 | 52.9 | 54.6 | 57.6 | 55.6 | 49.6 | ×1.8 |
| γp | 158 | 298 | 190 | 208 | 217 | 220 | 243 | ×1.9 |
| χ² total | 51.98 | 63.81 | 52.99 | 54.42 | 66.98 | 105.7 | 85.14 | |

(Aj4_s5 is identical to Ajn_s5: its DE never left the seed.) The `--tie {"theta_NOi": "theta_Ti"}` in run_fit_v4.sh is
a no-op at t_drive 4 (θTi is not in the parameter dict, fit_v2.tie copies only when it is), so θNOi is an
independent free parameter.

## 3. Proposal

| param | fitted | literature | identifiable? | decision (variant that first fixes it) |
|---|---|---|---|---|
| a00 | 1.127 | Chindemi b00 1.00 (fit) | yes (stiffest) | free |
| a01 | 1.168 | Chindemi b01 1.95 (fit) | sloppy (×3.4), trades with a11, ρ-γ | free (fixing it and ρ-γ together is risky) |
| a10 | 1.171 | Chindemi b10 1.16 (fit) | yes | free |
| a11 | 3.005 | Chindemi b11 2.48 (fit) | moderate | free |
| ρ-γ | 0.523 | none (default 1 breaks L2/3) | sloppy 0.49–0.64 | fix 0.5 (sloppy; R3) |
| θ_VDCC | 5.50 | requirement only | yes (×1.2) | free (core gate) |
| γd | 51.7 | 32–101 in fits | yes | free |
| γp | 158.4 | Chindemi 216.2 | sloppy 126–208 | fix 216.2 (lit; R1) |
| A_mglu | 0.308 | none | sloppy | free (sets the eCB-LTD speed) |
| θTe | 0.225 | none | sloppy (×13 across fits) | fix 0.2245 (sloppy; R2) |
| τ_T | 30.4 ms | none (tens of ms) | sloppy 19–50 | fix 30.44 (sloppy; R2) |
| θTg | 6.18 | none | moderate | free |
| d_min | −0.245 | −0.24 to −0.27 | yes | fix −0.25 (lit; R1) |
| A_NO | 4042 | none | sloppy (×5) | fix 4042 (sloppy; R2) |
| θNOi | 1.9e-4 | none | moderate | free |
| τ_NO | 6.99 ms | ≤ 14 ms (t½ < 10 ms) | sloppy 3.8–11.8 | fix 7.0 (lit; R1) |
| τ_Z | 48.7 ms | ≥ 10 ms | sloppy, tied to θZ (r 0.62) | fix 48.74 (sloppy; R2) |
| θZ | 1.75 | 1 spike (hippocampus); L5 untested | moderate | free |

Groups:
- **Fixed at the literature value (R1, 15 free):** d_min −0.25, τ_NO 7 ms, γp 216.2.
- **Tied:** θNOi stays independent; the existing θNOi→θTi tie is a no-op (section 2). τ_Z is effectively tied to θZ by fixing τ_Z.
- **Fixed at the fitted value because sloppy (R2, 11 free):** R1 + θTe 0.2245, τ_T 30.44 ms, τ_Z 48.74 ms, A_NO 4042.
- **R3 (10 free):** R2 + ρ-γ 0.5.

Free in R3 (10): a00, a01, a10, a11, θ_VDCC, γd, A_mglu, θTg, θNOi, θZ.

## 4. Refits

`fit_v4.py --fix-params` (env `FIX` in run_fit_v4.sh) removes the named parameters from the DE vector; the default `{}`
is unchanged. Repro (MAXITER=0), both χ² 51.9757 over 39 = C1Ajn_s5 exactly:
- 22133806: all 18 fixed, no seed (0 free). 2:41, 39.05 GB MaxRSS, 93% CPU.
- 22133808: 6 fixed, the other 12 from the C1Ajn_s5 seed. 2:42, 39.11 GB.

Refits: JOINT=1 FITGAMMA=1, SET vamp_mode 1, FREE=theta_V,rho_gamma, DROPT Letzkus 3AP −10 distal (39 targets),
seeds C1Ajn_s5 + C1Ajd_s5; 2g MIG, 1 CPU, 50G, 0:45 (MaxRSS 39.5 GB from 22132568 / 39.1 GB from the repros + 25%; 25 min + 50%).

| job | tag | free | fixed | seed |
|---|---|---|---|---|
| 22134084 | PR1_s5 | 15 | R1 | seeded |
| 22134085 | PR2_s5 | 11 | R2 | seeded |
| 22134086 | PR3_s5 | 10 | R3 | seeded |
| 22134087 | PR3_s6 | 10 | R3 | unseeded basin check |

Comparison (n = 39; χ² = −2 log L + const, so AIC = χ² + 2k, BIC = χ² + k ln 39). Full: AIC 87.98, BIC 117.92.
A reduced fit is preferred if its χ² is below: R1 (k 15) 57.98 (AIC) / 62.96 (BIC); R2 (k 11) 65.98 / 77.62;
R3 (k 10) 67.98 / 81.28.

| tag | χ² L5 / L2/3 / total | k | AIC | BIC |
|---|---|---|---|---|
| C1Ajn_s5 | 46.89 / 5.08 / 51.98 | 18 | 87.98 | 117.92 |
| PR1_s5 | pending | 15 | | |
| PR2_s5 | pending | 11 | | |
| PR3_s5 | pending | 10 | | |
| PR3_s6 | pending | 10 | | |
