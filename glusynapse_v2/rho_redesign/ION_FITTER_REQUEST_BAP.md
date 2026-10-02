# Request to ion_fitter: L5 TTPC emodel with a data-like distal apical bAP (draft, 2026-10-02)

From the glusynapse_v2 orchestrator. The user approved the request and asked that it be very specific. Background: EMODEL_LIMITS.md (D1, F1, F4), DISTAL_BAP_LIT.md.

## Why
The plasticity rule reads only synapse-local Ca. At the real L2/3->L5 synapse sites, split2 gives the following bAP amplitudes:

| Distance | Model bAP | Letzkus 2006 Fig 6D |
|---|---|---|
| < 200 um | 74 mV | |
| 200-450 um | 13 mV | |
| > 450 um | 2 mV | 45-65 mV |

The distal/proximal single-AP spine Ca ratio is 3e-4 to 3e-3, against ~0.1-0.3 expected. split2 fails bap_620 in 30/30 cells, and its best feasible value is ~8 mV.

The Letzkus 1AP +10 (0.72 +- 0.03) and 3AP +10 distal (0.79) LTD rows cannot be fitted without the distal bAP. A data-like distal bAP would also relax the Markram vs L2/3->L5 trade-off.

## What to build
- **Name:** `delta-split3-ap`. It gets a new emodel dir and circuit config. Any new or changed channel gets a new .mod and SUFFIX, e.g. `NaTg_ap`.
- **Cell:** L5TPC (cADpyr_L5TPC) only.
- **Start from:** split2-bk25 (the locked-in L5 emodel).
- **Free parameters:** apical only.
  - apical NaTg density and its vshiftm/vshifth;
  - apical Ka_kampa density, uniform or with at most a weak increasing gradient (L5 trunk: I_A +2.3 pA per 100 um, Bekkers 2000; total K falls distally, Korngreen & Sakmann 2000). Do not use the CA1 Hoffman gradient: it makes the distal bAP worse;
  - apical cm on the trunk (the cm 1 lever named in BAP_REPORT).
- **Fixed:** everything else at split2-bk25 values: soma, axon, basal, apical Ca densities, apical kBK and Ih.

## Targets (all must be met)
1. **Distal bAP** (single AP, measured on the main trunk; Letzkus 2006 J Neurosci, doi:10.1523/JNEUROSCI.2650-06.2006, Fig 6D):

| Distance | bAP amplitude |
|---|---|
| 100 um | ~100 mV |
| 400 um | 60-65 mV |
| 620 um | 45 +- 10 mV (the bap_620 feature) |

   The ratio of the 620 um amplitude to the 100 um amplitude must be 0.45-0.65.
2. **Ni-sensitive burst Ca along the trunk** (Fig 6C; 3 APs at 200 Hz; the Ni bound is LVA = 0):

| Distance | control/Ni burst integral |
|---|---|
| 200-250 um | ~1.15-1.2 |
| 300-400 um | ~1.3-1.45 |
| 450-660 um | ~1.5 |

3. **200 Hz burst fidelity:** 3 x 2 ms pulses at 3-5 nA, 200 Hz, must give 3 APs in every cell (Letzkus 2006 Methods). split2 fails this in 201 of 300 tasks. This is a soma/AIS objective, and doublets at threshold stay allowed.
4. **Keep:** all current split2-bk25 somatic features within their present scores, BAC firing 8/8, and the long apical plateau (no plateau-width cap; only a runaway guard).

## Constraints (user rules)
- No apical Ca hotspot zone: apical Ca stays uniform. An apical != basal split is allowed.
- Doublets are allowed.
- If uniform channels cannot reach target 1, stop. Report the best achievable profile as a known limit, and don't add a zone or a distance-dependent Ca.

## Step 1: feasibility probe (before any full fit)
- Run a quick probe on 8-10 L5 posts: sweep apical NaTg (x1-x4), the shift and trunk cm 1.
- Report the bAP at 100/400/620 um and the 3AP200 fidelity.
- Size it from 22064266 (1.0 GB), so pilot 1.3G 0:15. Take the time limit from the pilot's seff.

## Deliverables
- The emodel dir and circuit config `delta-split3-ap`.
- A table of bAP amplitude at 100/200/400/620/800 um and Fig 6C control/Ni ratios for 10 cells.
- The 200 Hz fidelity count over all 144 L5 post cells.
- The somatic feature scores against split2-bk25.

## Downstream cost (ours, after delivery)
Prefire, cexp and basis for L5->L5 + L2/3->L5 cost about 308 CPU-h, logged in DECISIONS before submission.

## Part 2 (same L5 package, user 2026-10-02): basal Ca, name `delta-split3-bca`
Deliver parts 1 and 2 as one package so that the downstream pipeline runs once. A suitable combined name is `delta-split3-apbca`; choose one name and tell us.

**Why.**
- Basal Ca_HVA2 = Ca_LVAst = 0.0023897659621395684. These are og-delta values, never refitted (L5TPC.hoc l. 320-321; json `og_delta_basal`).
- The basal 1-AP Ca ratio of >130 to <130 um is 0.32 in the model, against 0.95 +- 0.20 in Kampa & Stuart 2006.
- The shaft-Ca LTP licence is blocked in 509/1725 L5 potentiation rows. This caps 40/50 Hz +10 LTP at ~1.4 against 1.53/1.57 in the data.

**Free parameters.** Basal Ca_HVA2 and Ca_LVAst, untied. Keep the basal bAP voltage: 48 mV at 100-150 um already matches Nevian 2007. If an R-type channel is needed, it is a new mod `CaR_dend` with its own SUFFIX.

**Targets.**

| Target | Value | Source |
|---|---|---|
| 1-AP Ca, >130 / <130 um | 0.95 +- 0.20 | Kampa & Stuart 2006, doi:10.1523/JNEUROSCI.3062-05.2006 |
| Critical frequency | ~100 Hz | Kampa & Stuart 2006 |
| VSD 3rd/1st AP | 2.11 +- 0.22 | Kampa & Stuart 2006 |
| 200 Hz burst Ca | 4.6 +- 0.5x | Kampa, Letzkus & Stuart 2006, doi:10.1113/jphysiol.2006.111062 |
| Same burst in Ni | 3.5 +- 0.6x | Kampa, Letzkus & Stuart 2006 |
| Basal Ca spikes | none | Nevian 2007, doi:10.1038/nn1826 |

Basal kBK stays at the bk25 value.

## Part 3 (separate L2/3 package, user 2026-10-02): `l23delta-d1`
New dir `SSCx-AAD-l23delta-d1-emodels`, for L2TPC and L3TPC.

**Why.**
- The L2/3 hocs use the L5 og-delta dendrite values byte-for-byte, identical in og-delta, split2 and bk25:
  - apical NaTg 0.0229, Ka 0.0228, kBK 0.00356;
  - basal NaTg 0.003, Ka 0.002, HVA = LVA = 0.00239, kBK 0.0427.
- Shaft Ca in 50 Hz trains reaches ~16 uM in thin L2/3 dendrites, against ~1 uM expected.

**Free parameters.** Own apical and basal Na, Ka and Ca densities (split allowed). Basal kBK set like bk25 (0.0025) unless a check shows it matters.

**Targets.**

| Target | Value | Source |
|---|---|---|
| bAP at 150 um | ~75 mV | Waters 2003, doi:10.1523/JNEUROSCI.23-24-08558.2003 |
| bAP attenuation length | 317 um | Waters 2003 |
| Apical Ca | 10 % of max at 65 % of the soma-pia distance | Waters 2003 |
| Critical frequency | 130 Hz | Larkum 2007, doi:10.1523/JNEUROSCI.1717-07.2007 |
| Dendritic spike half-width | 2-4 ms | Larkum 2007 |
| Sag | little (low Ih) | Larkum 2007 |
| Basal spine Ca vs shaft | ~ equal, slight decline with distance | Koester & Sakmann 2000, doi:10.1111/j.1469-7793.2000.00625.x |

Guard, not a fit value: the shaft Ca plateau in a 50 Hz train should be ~1 uM and grow linearly with frequency. This is inferred from Helmchen 1996 (150-300 nM per AP, decay < 100 ms).

Keep the current L2/3 somatic features within their present scores.

**Downstream cost.** The L2/3->L2/3 chain is ~71 CPU-h on our side.

## Still open (not in this request)
- An age-matched (3-6 week) cell for Letzkus. Letzkus' cells have Rin 21 MOhm, against 142-191 MOhm for the juvenile cells. Tell us whether your pipeline can produce an older-age variant, but don't build it yet.
