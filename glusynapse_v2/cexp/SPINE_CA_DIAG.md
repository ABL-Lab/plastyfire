# Spine Ca diagnosis after E1 (A17, 2026-10-01)

Inputs: CHINDEMI_SPINE_CA.md, CEXP_LIT.md, CEXP.md, measure_cexp.py, `plastyfire/plastyfire/GluSynapse.mod` (the one compiled in DEES x86_64; mg is RANGE there), `/scratch/dhuruva/split1/cexp{,_ljp25}/*.csv`. Paper values via PubMed (DOIs given).

Measured diagnostics: `CEXP_DIAG=1` runs in `/scratch/dhuruva/split1/cexp_diag/{ljp25,ljp25_js,ljp25_l23}/parts/`, 3 pairs each (L5: 33 synapses; L2/3→L2/3: 13). Jobs 22184271/2/3. "Counted" means read off the sorted csv at 2 decimals by hand; no tool computation.

**What each condition measures** (measure_cexp.py): every synapse of the pair gets `rho0=1, Use_p=1`, so released = Nrrp: a full-RRP EPSP on all contacts at once, not a single vesicle or the mean trial. The cell sits at emodel rest (no holding current). APV: gmax_NMDA=0. Mg0: mg=0 on the pair's synapses. Post: one somatic bAP. Recorded: effcai peak after 1000 ms, cai_CR peak after 990 ms (no baseline subtracted), ica charges integrated over **510 ms** (990–1500 ms).

## Measured (diag jobs) vs A17 predictions
| quantity | prediction | measured | verdict |
|---|---|---|---|
| L5 v_rest at synapse | −75 to −78 mV | **−77.8 to −79.1 mV** (33 syn) | confirmed |
| L5 pre_eff0 (effcai pedestal) ≈ cpre_apv, proximal | ≈ | eff0 0.76–1.6e-3; **80–87 % of cpre_apv** at 15–64 µm (e.g. 1.28e-3 of 1.55e-3); cai0 72.8–75.9 nM | confirmed |
| L5 J&S: Cpre ×, Mg0/pre ÷ | ×2.7–3, ÷2.7–3 | proximal effcai ×2.9–2.95 (free cai 1.06→2.88 µM); distal shared-branch ×1.9; apical 357 µm ×1.2; Mg0 unchanged, so Mg0/pre ÷2.9 proximal | confirmed (proximal) |
| L2/3 v_rest | about −77 mV | **−81.9 to −82.2 mV** (13 syn); cai0 71.8–73.5 nM | refuted (5 mV more negative) |
| local EPSP at synapse (pre_vpk − v_rest) | "a few mV" | L5 proximal **+9 to +12 mV** (to −66…−69.5); shared-branch distal to −3…−26 mV; L2/3 +7 to +16 mV | corrected |

## (a) Root causes, ranked
**1. Mg0/pre ×60–100 is the Mg block at the resting operating point: a model property.**
- B = 1/(1+e^(−0.072V)·Mg/2.552) gives 1/B = 109 at −78 mV, 145 at −82, and 50 at −67 (the L5 proximal EPSP peak).
- The measured proximal effcai Mg0/pre is 78–82 in these 3 pairs (L5), between 1/B at rest and at peak. The pooled L2/3 value of ~100 is likewise 1/B(−82) partly relieved by a +7–16 mV EPSP (so my earlier "L2/3 ⇒ −77 mV" inference was wrong).
- Where the local EPSP reaches −3…−26 mV (2–4 contacts of a pair sharing a branch), the fold is 3–5×: L5 82.6 µm 9.4→35.7 µM, 95.1 µm 11.0→54.8 µM.
- Ranked: (1a) **operating point**: −78/−82 mV rest, no spine-head voltage beyond the dendritic +9–12 mV. (1b) **Block steepness**: Jahr & Stevens 1990 (0.062/3.57, [doi](https://doi.org/10.1523/JNEUROSCI.10-09-03178.1990); standard form, PMC text empty) measured ÷2.9 → about 27×, still not few-fold. (1c) **Measurement**: effcai vs peak free Ca (58 vs 33–89 pooled); full RRP lowers the fold; Mg0 free Ca of 35–136 µM is beyond any indicator's range. (1d) **g_NMDA / Pf** cancel in the ratio; absolute check 1 nS·0.05·60 mV = 3 pA → 4 µM/ms·12 ms = 49 µM vs 35–80 observed; Pf at −70 mV = 8 % of NMDA current (Schneggenburger GHK, as Chindemi).
- The "few-fold" target has **no source** (CEXP_LIT row 5).

**2. The 2–7 % APV-resistant fraction is mostly the E1 window-current pedestal plus distance pooling.**
- Measured: the proximal cpre_apv (1.1–1.7e-3) is 80–87 % resting pedestal (eff0). Pedestal-corrected proximal APV share = (1.55−1.28)/(53.6−1.28) = **0.3–0.6 %**.
- Distal shared-branch synapses reach 40–56 % ((9.18−0.07)/(16.39−0.07) at 142 µm). The medians mix the two.
- E1 window current: in the exponential regime it scales by e^(50/9.5)·0.31·(0.80/0.98) ≈ ×49, vs vdcc_q_post ×37–43 at apical sites, where there is no bAP. **That ×37–43 is window current, not bAP Ca.** The vdcc_q-based V_rest estimate (−78 mV) matched the measured −77.8…−79.1.
- At the measured proximal EPSP peak (−67 mV), VDCC m∞² is ×10 its rest value, but only for the few-ms AMPA phase. A 10–30 % share needs ×20–60, since m∞² ∝ e^(2ΔV/9.5). That means **+14 to +19 mV beyond the present +9–12 mV**: a spine-head EPSP (Harnett 2012, R_neck ~500 MΩ, [doi](https://doi.org/10.1038/nature11554)).
- gca is not the knob: bAP Ca pins it (1.5 µM at < 60 µm), and ×20–60 would give 30–90 µM. Volume is not either (∝V^(−1/3)). NMDA Ca is not too high: proximal cpre is about 1 µM at full RRP, vs Sabatini 2002's 0.7 ([doi](https://doi.org/10.1016/s0896-6273(02)00573-1)).
- Cortical data: "almost entirely" NMDAR (Kovalchuk 2000 [doi](https://doi.org/10.1523/JNEUROSCI.20-05-01791.2000); Nevian & Sakmann 2004 [doi](https://doi.org/10.1523/JNEUROSCI.3332-03.2004)). The 0.1–0.3 range is from CA1.

**3. L5 basal bAP/EPSP (free Ca) passes Koester at < 80 µm.** `dist_um` is the path distance; L5L5 basal has **58/143 synapses < 80 µm**. Counted medians:
- **1.30 at < 60 µm** (n 28), **0.84 at 60–80 µm** (n 29), **≈1.0 at < 80 µm** (n 57). Koester & Sakmann 1998 report "comparable up to 80 µm" ([doi](https://doi.org/10.1073/pnas.95.16.9596)).
- Beyond 100 µm, 0.02–0.5: the bAP peaks at −24…−54 mV there (measured post_vpk), while shared-branch cpre reaches 5–17 µM. **The pooled 0.46 mixes distances.**
- Other pathways: L2/3→L5 basal cpre 0.5–1, cpost 2–3 µM to 130 µm (≈ 2.6). L2/3→L2/3 basal cpost 1–2 µM to 40 µm, 0.2–0.6 beyond 60.
- Caveat: the trial-mean EPSP is about U× our full-RRP cpre, so the like-for-like ratio is up to 1/U higher.

**4. Artefacts vs properties.**
- **Artefacts:** the effcai pedestal (2–3 % of proximal cpre, 6 % of proximal cpost, but **≥50–75 % of cpost at apical/distal sites**, e.g. 1.22e-3 of 1.60e-3 at 357 µm); the 510 ms window-current integral (vdcc_q, ×37–43); the 0.070 floor instead of the measured baseline (72–76 nM); distance pooling; full-RRP all-contact cpre; effcai-vs-peak weighting.
- **Model properties:** the Mg-block curve and the −78/−82 mV rest; no spine-head voltage; the E1 window current (+2–6 nM resting cai, inside Sabatini's 70 ± 29 nM).

## (b) Fixes within the constraints (no new channel; paper globals or measurement changes)
- **M1 (measurement, recommended):** subtract eff0/cai0 (now recorded with `CEXP_DIAG=1`; default output unchanged); charges as q − I(989)·window over 1000–1300 ms; report Δpeak free Ca next to effcai; bin at < 60 µm (Chindemi/Sabatini) and < 80 µm (Koester).
- **M2 (measurement, optional):** compare at the paper's holding potential via somatic DC. L5–L5 pairs: −60 ± 2 mV (Markram 1997, [doi](https://doi.org/10.1113/jphysiol.1997.sp022031)). Predicted at −60 vs −78: NMDA Ca ×3.6 (1/B 30.5 vs 109); VDCC ×44·0.48 ≈ ×21 (APV share about 2.6 %); resting cai about 165 nM.
- **P1 (global, anchored, not recommended, now measured):** Jahr & Stevens block (slope_NMDA 0.062, scale_NMDA 3.57). It gives Mg0/pre about 27× but proximal free cpre of 2.3–3.1 µM vs Sabatini's 0.7. That breaks Chindemi's validation; V-C&R 2003 ([doi](https://doi.org/10.1152/jn.01038.2002)) is the cortical fit.
- **Outside the constraints:** only a spine-neck voltage (A15 E2, R_neck 500 MΩ, new mod, no channel) gives the +14–19 mV a 10–30 % APV share needs. That is the user's call.

## (c) Downstream
- **M1:** cexp only; no prefire or refit. The pedestal also sits in prefire effcai and in c_pre/c_post, so θ = a·c_pre + b·c_post absorbs it. But **at apical/distal synapses, c_post is mostly pedestal**: their θ terms in c_post carry no bAP information. In the ljp25 refit, check θ_d > eff0 per synapse, and consider c_post − eff0 in the fit-threshold definition (a measurement change; needs a re-cache, not a re-prefire).
- **M2:** cexp only. If applied to the prefire protocols: c_pre about ×3.6, all 8 thresholds move, so a new prefire and refit.
- **P1:** new prefire, cexp and refit; c_pre ×2.9 proximal and ×1.2–1.9 distal, so θ_d/θ_p shift non-uniformly and pairing supralinearity changes.
- **E2:** new mod, prefire, cexp and refit.

## (d) Recommendation (unchanged by the data, one addition)
Keep E1 and Chindemi's Mg block; adopt M1 now. Accept Mg0 ×60–100 and the 0.3–0.6 % proximal APV share as model properties, consistent with the measured −78/−82 mV operating point and with cortical (not CA1) data. Drop the uncited few-fold Mg0 target; bAP/EPSP passes Koester at < 80 µm. **New:** before the ljp25 refit, report c_post − eff0 per synapse (apical/distal c_post is ≥50 % pedestal). Ask the user about E2 only if the CA1 APV share is required.
