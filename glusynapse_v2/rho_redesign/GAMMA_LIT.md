# γd / γp in published calcium-threshold (Graupner–Brunel-type) rules

Rule form: `τ dρ/dt = −ρ(1−ρ)(ρ*−ρ) + γp Θ(c−θp)(1−ρ) − γd Θ(c−θd) ρ`. γ is dimensionless and τ is in seconds,
so while c is above threshold the rate is **γ/τ (s⁻¹)**. Averaged over a protocol it is Γ = γ·α, where α is the fraction of time
above threshold, and τ_eff = τ/(Γp+Γd) (GB2012 main text; Graupner 2016 Methods). Θ is a Heaviside on the
*instantaneous* Ca trace (GB-type: c in units of θd = 1) or on the *leaky-integrated* c* (Chindemi: τ* ≈ 278 ms, per-synapse
θ = a·Cpre + b·Cpost). With a slow c*, more time is spent above threshold per event, so absolute γ values are only roughly
comparable between the two families. γd/γp and γ/τ are the comparable quantities.
"No cubic" means the paper drops the bistable term (continuum of states). Ratios and γ/τ are my arithmetic from the listed numbers.

| Paper / set | γd | γp | τ (s) | ρ* | γd/γp | γd/τ, γp/τ (s⁻¹) | Ca / thresholds, notes | Source |
|---|---|---|---|---|---|---|---|---|
| GB2012 DP (generic) | 200 | 321.808 | 150 | 0.5 | 0.62 | 1.33, 2.15 | τCa 20 ms, Cpre 1, Cpost 2, θd 1, θp 1.3; ratio set by "no change at large Δt" (main text) | SI Table S1 ‡ |
| GB2012 DPD / P / D (generic) | 250 / 160 / 500 | 550 / 257.447 / 550 | 150 | 0.5 | 0.45 / 0.62 / 0.91 | — | illustrative curve-shape sets, not data fits | SI Table S1 ‡ |
| GB2012 DPD′ / D′ (generic) | 50 / 60 | 600 / 600 | 150 | 0.5 | 0.08 / 0.10 | — | θp 2.5 / 3.5; not data fits | SI Table S1 ‡ |
| GB2012 cortical slices (Sjöström 2001, V1 L5) | 331.909 | 725.085 | 346.3615 | 0.5 | 0.46 | 0.96, 2.09 | τCa 22.69 ms, Cpre 0.5618, Cpost 1.2396, θd 1, θp 1.3, D 4.61 ms, σ 3.35 | SI Table S2 ‡; same numbers in Higgins 2014 Table 1 |
| GB2012 hippocampal slices (Wittenberg & Wang) | 313.0965 | 1645.59 | 688.355 | 0.5 | 0.19 | 0.45, 2.39 | τCa 48.8 ms, Cpre 1, Cpost 0.276 (Cpost also in Fig. 3 legend), θp 1.3 | SI Table S2 ‡ |
| GB2012 hippocampal cultures (Wang 2005) | 61.141 | 113.6545 | 33.7596 | 0.5 | 0.54 | 1.81, 3.37 | τCa 11.95 ms, Cpre 0.58, Cpost 1.76, θp 1.3 | SI Table S2 ‡ |
| Higgins, Graupner & Brunel 2014, in vitro | 331.909 | 725.085 | 346.3615 | 0.5 | 0.46 | 0.96, 2.09 | GB2012 cortical set reused; in vivo set rescales only Cpre/Cpost (0.337/0.744) | Table 1, 10.1371/journal.pcbi.1003834 |
| Graupner, Wallisch & Ostojic 2016, linear Ca | 137.7586 | 597.08922 | 520.76 (τcb) | no cubic | 0.23 | 0.26, 1.15 | V1 L5 regular + jittered pairs (Sjöström 2001); τCa 22.27 ms, Cpre 0.844, Cpost 1.621, θp 2.009 | Table 2, 10.1523/JNEUROSCI.0104-16.2016 |
| Graupner 2016, nonlinear Ca (n = 2) | 111.82515 | 894.23695 | 707.02 (τcb) | no cubic | 0.13 | 0.16, 1.26 | τCa 18.93 ms, θp 4.998 | Table 2 (same DOI) |
| Deperrois & Graupner 2020, visual cx: +STD / no STD / nonlinear | 111.32 / **31.98** / 183.51 | 564.39 / 161.99 / 1000 (at bound) | 299.9 / 79.98 / 525.9 | no cubic | 0.20 / 0.20 / 0.18 | 0.37, 1.88 / 0.40, 2.03 / 0.35, 1.90 | Sjöström 2001 L5–L5 V1; θp 1.63 / 1.63 fixed / 2.31; fit bounds γd 20–1000, γp 100–1000 (authors' code `params.py`) | Table 1, 10.1371/journal.pcbi.1008265 |
| Deperrois & Graupner 2020, somatosensory cx: +STD / no STD / nonlinear | 176.54 / 105.05 / 157.34 | 579.58 / 406.98 / 518.17 | 143.1 / 26.6 / 196.8 | no cubic | 0.30 / 0.26 / 0.30 | 1.23, 4.05 / 3.95, 15.3 / 0.80, 2.63 | Markram 1997 L5–L5 | Table 1 (same DOI) |
| Chindemi et al. 2022 (final) | 101.5 | 216.2 | 70 † | 0.5 | 0.47 | 1.45, 3.09 | c* integrator τ* 278.318 (bounds printed "(150, 350) s"; Ecker 2025 gives ms); 8 a/b threshold factors; fit bounds γ 1–300 | Table 2, 10.1038/s41467-022-30214-w |
| Chindemi bioRxiv v1 (2020) | 71.1 | 225.8 | 70 † (assumed) | 0.5 | 0.31 | 1.02, 3.23 | τ_effcai 314.4; same 1–300 bounds | Table 3, 10.1101/2020.04.19.043117 |
| Ecker et al. 2025 (SSCx network, Chindemi rule) | 101.5 (alt. 100) | 216.2 (alt. 450) | 70 (stated) | 0.5 | 0.47 (0.22) | 1.45, 3.09 | alternative set 100/450, τ* 200 ms, in a supplementary-figure legend | Eq. 1 + fig. suppl., 10.7554/eLife.101850 |
| BBP GluSynapse.mod defaults (code, not a paper) | 100 | 450 | 70 (tau_ind_GB) | 0.5 | 0.22 | 1.43, 6.43 | tau_effca 200 ms; overridden by circuit config | neurodamus-models `common/mod/GluSynapse.mod` |
| **Our joint fit** | **53** | **190** | 70 | 0.5 | **0.28** | **0.76, 2.71** | c* as in Chindemi | — |

‡ The GB2012 SI PDF (10.1073/pnas.1109359109) could not be retrieved (captcha). These values come from the authors' code
(github.com/mgraupe/CalciumBasedPlasticityModel, `synapticChange.py`, sets named after the paper's cases). The cortical set is
independently confirmed by Higgins 2014 Table 1, which cites GB2012. Check the other sets against Table S1/S2 before quoting them.
† Chindemi 2022 gives no number for τ (grep of the main text and SI is negative). 70 s is from GluSynapse.mod and from Ecker 2025,
which states τ = 70 s with γd 101.5 / γp 216.2.
**Not located:** γp 199.77 does not appear in Chindemi 2022 (main text, SI or peer-review file), in the bioRxiv v1 or in GluSynapse.mod. Its origin is unverified.
Shouval 2002 and Rachmuth 2011 use the Ω(Ca)/η(Ca) rule, which has no γd/γp analogue, so they are not tabulated (not checked in detail).

## Recommendation
1. γd < γp in every data-fitted set. Published γd/γp ranges from 0.13 to 0.54 (Graupner 2016 Table 2; GB2012 Table S2), and our 0.28 is mid-range. It sits next to Chindemi v1 (0.31), Deperrois somatosensory (0.26–0.30) and Graupner 2016 linear (0.23).
2. γd = 53 at τ = 70 s gives 0.76 s⁻¹. That is inside the published γd/τ range of 0.16–3.95 s⁻¹, and Deperrois 2020 V1 (no STD) fits γd = 32 (Table 1). So 50 is not a lower edge in the literature, and the bound is truncating the fit.
3. Proposed bounds: **γd 20–250, γp 100–600**. The γd floor of 20 matches Deperrois 2020's fit bounds (20–1000; code `params.py`) and leaves room below the 32 fitted there. The ceilings cover GB2012 cortical at τ = 70 (γ/τ ≈ 0.96 and 2.1 s⁻¹ → 67 and 147) up to Deperrois somatosensory (γp/τ ≈ 4 s⁻¹ → 284) with margin, and stay near Chindemi's search box (1–300).
4. Optionally constrain γd/γp to 0.1–0.65, the range of all GB2012/2016 and Deperrois sets, or fit the ratio directly, since GB2012 fixes it by the large-Δt balance condition.
5. Caveat: our Θ acts on the integrated c* (as in Chindemi), not on raw Ca, so compare via γ/τ and γd/γp rather than raw γ.
