# Experimental-style calcium measures per synapse (delta-split1)

Per synapse, no .mod edits (mg and gmax_NMDA are RANGE variables set per synapse):
- Cpre: one pre spike, 1 mM Mg (Chindemi c_pre).
- Cpre_APV: same, gmax_NMDA = 0.
- Cpre_Mg0: same, mg = 0 (full NMDA, no Mg block).
- Cpost: one bAP (soma TStim, amplitude from spike_threshold_finder), 1 mM Mg.

`effcai` = effcai_GB peak (Chindemi units; what the thresholds use). `_cai` = cai_CR peak in uM (free spine Ca, absolute; resting floor 0.070 uM).
Charge `vdcc_q`, `nmda_q` = inward integral of ica_VDCC / ica_NMDA from 990 ms (pC), for the pre spike at 1 mM; `_apv`, `_mg0`, `_post` variants in the csv.
Files: /scratch/dhuruva/split1/cexp/{L5L5,L23L5,L23L23}.csv (syn_id = global edge id, as in the cache and extracted npz).
Code: measure_cexp.py, run_measure_cexp.sh, summarize_cexp.py. Summary log: logs/cexp_summ_22156396.out.

## Medians (IQR in the log), effcai units

| pathway | loc | n syn | cpre | cpre_apv | cpre_mg0 | cpost | mg0/cpre | cpost/cpre | NMDA share 1-apv/cpre |
|---|---|---|---|---|---|---|---|---|---|
| L5->L5 | basal | 143 | 0.060 | 1.0e-4 | 3.5 | 1.9e-3 | 62 | 0.033 | 0.998 |
| L5->L5 | apical | 48 | 0.058 | 5.9e-5 | 3.4 | 1.1e-4 | 64 | 0.002 | 0.999 |
| L2/3->L5 | basal | 154 | 0.047 | 3.8e-5 | 4.6 | 7.3e-3 | 102 | 0.14 | 0.999 |
| L2/3->L5 | apical | 644 | 0.054 | 5.2e-5 | 4.6 | 9.7e-5 | 88 | 0.002 | 0.999 |
| L2/3->L2/3 | basal | 387 | 0.044 | 1.9e-5 | 4.3 | 4.1e-4 | 101 | 0.009 | 0.9995 |
| L2/3->L2/3 | apical | 102 | 0.046 | 2.0e-5 | 4.3 | 8.8e-3 | 104 | 0.23 | 0.9995 |

Free spine Ca (cai_CR, uM): cpre 0.87-1.6 (L5 basal 1.6, L2/3->L5 basal 0.87); APV 0.070-0.076 (= resting, ratio to cpre 0.07); Mg0 51-72 (mg0/cpre 33-89); cpost 0.07 at distal synapses up to 0.64-0.99 (L2/3->L5 basal, L2/3->L2/3 apical medians) and 0.21 for L5 basal.
Charges (pC): nmda_q 0.0057 (L5), 0.0015 (L2/3->L5), 0.0024 (L2/3->L2/3); vdcc_q 7e-6, 1.4e-6, 1.1e-6, i.e. VDCC is ~1e-3 of NMDA for a single pre spike at 1 mM. Under APV the VDCC charge is unchanged (ratio 0.85-1.06).

## Reading
- A single pre spike in this model is >99.8 % NMDAR-dependent (effcai); APV leaves only the VDCC-driven residue at the resting floor in cai_CR.
- Removing the Mg block raises Cpre 60-100x (effcai), 30-90x (free Ca).
- Cpost is strongly distance dependent: negligible (<1e-3 x cpre median, resting floor in uM) at apical/distal synapses, 0.03-0.2 x cpre at basal ones; the wide IQRs reflect this. Cpost is NaN (no single-AP stimulus up to 5 nA at 1.5/3/5 ms) for 6.8 % of L5 synapses (13, one pair 207453-189325), 7.0 % of L2/3->L5 (56, 7 pairs), 0 in L2/3->L2/3. The same cells failed in the original cache build.

## Check against the cache (c_pre / c_post of the fit, defit2 globals)
- cpost: ratio 0.989-1.031 (all pathways). L5 cpre: ratio 1.000 on all 24 pairs.
- cpre differs pair-specifically in L2/3 pathways: L2/3->L5 median 0.996, 13/120 pairs off by >2 % (worst pair 153238-206803, 0.72, min 0.61, max 1.25 over synapses); L2/3->L2/3 median 0.995, 48/120 pairs off by >2 % (range 0.19-1.74; e.g. 164624-145443 0.70, 2285-8885 1.33). The offset is not tied to location (basal/apical medians 0.99-1.00) and weakly to distance (corr -0.13 / -0.29). Not diagnosed; probable cause (unverified): stochastic release / per-synapse seed differing between the cache build and these runs, since pairs with low gmax_NMDA and L2/3 pairs are the ones that move. Use these tables for relative quantities (ratios within one run); for thresholds the fitter keeps the cache c_pre.

## Literature (all numbers UNVERIFIED, from memory, not checked)
- Nevian & Sakmann 2006 (J Neurosci): spine Ca transients evoked by single/few EPSPs in L2/3 and L5 pyramidal spines were largely NMDAR (APV-sensitive); the model's >99 % APV-sensitive fraction is at the upper end of that picture, and the model has no APV-resistant (VDCC/other) component for a single EPSP. Exact fractions not checked.
- Kampa et al. 2004 (J Physiol; pdf in $R): NMDAR-dependent spine Ca with Mg-block relief by back-propagating APs; fraction not checked.
- The Mg0 ratio (60-100x) has no experimental counterpart checked.

## Sizing (measured)
Arrays (CHUNK=12 pairs, 1 CPU): 1.5-8.7 min per task, ~40 s per pair, CPU eff ~97 %, MaxRSS up to 2.09 GB (limit 2G; use 2.7G next time), ~1.5 CPU-h total.
