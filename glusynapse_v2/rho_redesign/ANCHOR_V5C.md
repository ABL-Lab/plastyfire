# Anchoring the two v5c thresholds (θ_V, θ_eCB), 2026-10-01

Target: v5c (V5_DESIGN §4, `results/v5_V5c_s7.json`) with only a00, a01, a10, a11, γd, γp fitted. The two remaining
flagged constants are θ_V = 0.245 (potentiation gate on the spine VDCC-Ca pool V) and θ_eCB = 6.91 (eCB trigger on the
(1−b)-weighted pool W). Both are in the pool units of PARAM_ANCHORS line 15: 1 unit = 1e-5 nA·ms of own-VDCC charge,
which is 31 Ca ions or 0.60 µM total Ca (≈ 0.024 µM free) in 0.087 µm³. The pool memory is τE1 = 100 ms.

Script: `run_calib_vdcc_kampa_split1.sh`, a copy of `run_calib_vdcc_kampa.sh` with two modes: MODE=V for C-V and
MODE=E for C-E. Both run on the delta-split1 extracted traces (`extracted/ebner_l23l5_delta-split1-prefire-vseg-rs`,
`extracted/ebner_delta-split1-prefire-vca`). Section types are read from the pathway's own edges file. Outputs:
`/scratch/dhuruva/param_anchors/calib_kampa_split1_{V,E}_syn.csv` and `logs/calib_kampa_split1_<job>.out`.
Request: 1 CPU, 1500M, 0:15 (basis 22135045: 0:53, 1.16 GB).

## 1. θ_V by C-V on delta-split1 (job 22140071)

Observable (Kampa 2006, 10.1113/jphysiol.2006.111062, rat L5 basal dendrites): with 3 APs, a 200 Hz burst gives
supralinear Ca that low Ni²⁺ blocks, plus LTP; a 50 Hz burst gives linear Ca and no LTP. The LTP outcome is not a
fit target, and the lock uses only the Ca observation.

Lock: θ_V* is the pool level V after the 3rd AP that separates 200 Hz from 50 Hz at basal L2/3→L5 synapses. Two
criteria are computed:
- Youden J (unpaired), as in the og-delta run;
- the paired bracket F(θ) = P(V3,50 Hz < θ < V3,200 Hz) on the same synapses. This is new: og-delta failed because
  geometry spreads V over 3 decades, and a paired criterion removes that spread.

The job also reports:
- the pool supralinearity V3 / (V1 × linear sum);
- the shaft-Ca ratio at 200 vs 50 Hz (shaft_cai, Kampa's own observable);
- P(open) at θ_V 0.245 (v5c) and at 5.50 (C1).

Data caveat: split1 extraction 22136807 (another session's pipeline) timed out after 309 of its 339 files, and its
step 2 never ran. As a result, the bAP-only arm `letzkus_3ap_200hz_dt-500ms` is missing and
`letzkus_3ap_200hz_dt+10ms` has 51 files (og-delta: 56). This run therefore uses the 200 Hz +10 arm, which is
EPSP-matched to the 50 Hz comparator (S&H 50 Hz +10, one EPSP 10 ms before the 1st AP). The script picks up the dt−500
arm automatically once it has been extracted; rerun with MODE=V then.

Decision rule:
- **(a) Lock.** Split1 is supralinear at 200 Hz only (SI clearly > 1 at 200 Hz, ≈ 1 at 50 Hz), and F* ≥ 0.5 or J ≥ 0.3.
  Then θ_V is set to θ* and locked in v5c units (C). This needs a validation refit with θ_V fixed, which is expected to
  cost χ², because the fit sits at the box floor (0.245, gate ≈ always open).
- **(b) Not separable.** Split1 is supralinear but the gate is not separable (F*, J small). Then θ_V cannot be a global
  constant, and the gate is dropped (θ_V = 0, test V5c_nogate_s5 22136667: drop if Δχ² ≲ 2).
- **(c) No supralinearity.** Same outcome as og-delta (null).

**Result: (b/c), no lock (22140071: 0:46, 1.46 GB of 1.5 G, CPU 61 %).**

Data used: basal L2/3→L5 n 68 (200 Hz +10) and 128 (S&H). Results:
- **The spine VDCC pool has no 200 Hz supralinearity.** The pool SI is 0.69 at 200 Hz [q10 0.49, q90 1.29] and 0.78
  at 50 Hz. The median V3 is 12.6 in both arms.
- **Paired, per synapse.** P(V3,200 > V3,50) = 0.56, and F* = 0.074. Youden J = 0.02 against S&H and 0.00 against L5
  Sjöström 50 Hz −10, the same as og-delta.
- **The v5c gate is open in both arms.** At θ_V 0.245 it is open at 84 % of basal synapses for 200 Hz and 84 % for
  50 Hz (95 % for L5 50 Hz −10). In practice the gate is a no-op.
- **Shaft Ca does change with frequency, but the test cannot tell whether it is supralinear.** The shaft_cai rise
  after AP3, 200 vs 50 Hz, has median 1.56 [q10 1.12]. Linear summation with a shaft τ of about 20 ms would also give
  about 1.59, and a 1-AP-normalised shaft SI needs the bAP-only arm.

Reading: split1's dendritic Ca spike, which is what passes Kampa at cell level, does not reach the synapse's own spine
VDCC current, and V is fed only by that current. So θ_V cannot be anchored by Kampa on V.

The gate cannot be dropped. V5c_nogate_s5 (θ_V = 0) has χ² 161.47 against 84.69 for V5c_s7 (BIC 187.1 vs 114.0).

At 0.245, θ_V is not a frequency gate but a floor. It closes potentiation at synapses whose own VDCC Ca is almost nil:
- In the pool distribution, basal synapses have V3 q10 ≈ 0.1 and q50 ≈ 12, while apical synapses have q50 ≈ 0.1.
- In one-bAP free-Ca terms (§2 mapping), 0.245 is about 5 nM peak, or about 8 Ca ions.

Proposal: anchor θ_V the way K was handled (C-K).
- **Floor test (C-Vf).** Rescore V5c_s7 with MAXITER 0, all parameters fixed, at θ_V × 0.5 and × 2.
- **If Δχ² < 1 both ways:** θ_V is a structural detection floor with no fitted content. It means "the synapse sees own
  VDCC Ca", and its value is the gap in the pool distribution.
- **If not:** θ_V stays flagged.
- **Sizing,** from the C-K MAXITER-0 runs 22135619/20 (2g MIG): 2g, 1 CPU, 49G, 0:15.
- **Kampa route:** a Kampa-anchored gate is only possible on shaft Ca, which is allowed as a synapse-local Ca drive.
  That needs a new kernel file, and a C-V rerun on shaft Ca with the dt−500 arm.

## 2. θ_eCB: literature anchor and C-E (job 22140072)

Ca thresholds for postsynaptic eCB release (via PubMed):

| source | prep | trigger | [Ca]i threshold |
|---|---|---|---|
| Brenowitz & Regehr 2003, 10.1523/JNEUROSCI.23-15-06373.2003 | rat cerebellar Purkinje cells | Ca alone (depolarisation) | half-max ≈ 15 µM |
| Wang & Zucker 2001, 10.1111/j.1469-7793.2001.t01-1-00757.x | rat CA1 pyramidal cells, NPE photolysis | Ca alone (DSI and PSI) | half-max 3.6–3.9 µM, linear in [Ca]i |
| Hashimotodani 2005, 10.1016/j.neuron.2005.01.004 | cultured hippocampal neurons, PLCβ1 | Gq receptor plus Ca coincidence | physiological [Ca]i (sub-µM); the exact values must be read from the figures (no PMC full text) |
| Maejima 2001, 10.1016/s0896-6273(01)00375-0 | rat cerebellar Purkinje cells, mGluR1 | receptor-driven, no Ca rise needed | none (no Ca threshold) |

Which threshold applies:
- Single bAPs drive L5 tLTD (Sjöström 2003), and they never reach a 4–15 µM bulk Ca level. The Ca-alone
  (DSI/DSE) thresholds are therefore upper bounds.
- The candidate anchor is the coincidence (PLCβ) threshold of Hashimotodani 2005. tLTD needs a coincident
  glutamate signal, and the model's mglu_block arm assumes the same.
- Confidence is L: the preparation is cultured hippocampal neurons, and the value has not been read from the figure yet.

C-E maps a [Ca]i threshold onto W with the model's own spine Ca (Chindemi constants, b = 0 as in Ca-clamp
experiments). There are two mappings, and they differ by τE1/τCa ≈ 8:
- **Steady (Ca clamp, photolysis for seconds):** W_ss/ΔCa_ss = τE1/(i_scale · K_Ca · τCa) ≈ 345 units per µM. On this
  scale θ_eCB 6.91 ≈ 0.020 µM.
- **Transient (one bAP):** V(t_bAP + 10 ms) divided by the peak free VDCC ΔCa of the same bAP, per synapse. The primary
  arm is L5 Sjöström 0.1 Hz −10, where the pre spike arrives at +10 ms. The job prints θ_eCB* for [Ca]i thresholds from
  0.1 to 15 µM, and the fraction of single-bAP synapses that would trigger.

The steady mapping is the faithful one for clamp-type data. With it, 0.1–0.3 µM gives θ_eCB ≈ 35–100, about 5–15×
the fit, so single-bAP tLTD would fail at many synapses. The transient mapping lands near the fit if ΔCa_thr ≈ 0.2 µM.

Decision rule:
- If the job shows that θ_eCB*(steady, Hashimotodani range) leaves P(single-bAP trigger) ≳ 0.5 at L5, lock it (C) and
  run the validation refit.
- Otherwise θ_eCB has no transferable anchor. In that case either (i) report it flagged as the one eCB constant, or
  (ii) reformulate the trigger as a Hill function of free spine Ca (PLCβ EC50) in a new kernel file. That needs the
  user's OK.

**Result (22140072: 0:19, 0.97 GB, CPU 74 %).**
- **Model spine Ca for one bAP.** At L5 basal synapses (n 135) the peak free VDCC ΔCa has median 0.21 µM
  [q10 0.002, q90 2.7].
- **Transient mapping.** It is tight and the same in both pathways: 48.3 units per µM (L5 basal q10–q90 44.7–60.0;
  L2/3→L5 basal 46.6). The steady mapping is 349.8 units per µM.
- **What the fitted θ_eCB 6.91 corresponds to.** Peak ΔCa 0.143 µM (transient) or 0.020 µM (steady).
- **Lit thresholds converted to θ_eCB** (in brackets: P that one bAP triggers at L5, at +10 ms):

| [Ca]i threshold | θ_eCB, transient mapping | θ_eCB, steady mapping |
|---|---|---|
| 0.1 µM | 4.8 (0.47) | 35 (0.30) |
| 0.2 µM | 9.7 (0.43) | 70 (0.21) |
| 0.3 µM | 14.5 (0.41) | 105 (0.14) |
| Wang & Zucker 3.6–3.9 µM | 174–189 (0.04) | 1259–1364 (0.006) |
| Brenowitz & Regehr 15 µM | 725 (0.006) | 5246 (0) |

At the fit, P = 0.44.

Proposal: lock θ_eCB = 48.3 × Ca_thr (C, confidence L), using the transient mapping. The reason for the transient
mapping is that the tLTD trigger is a single-bAP Ca event, and PLCβ reads instantaneous Ca. The steady mapping
describes sustained clamp data, not this protocol class.
- **Ca_thr.** It comes from Hashimotodani 2005 (PLCβ1 coincidence, physiological [Ca]i). The value must be read from
  the figures. For 0.1–0.3 µM, θ_eCB = 4.8–14.5, which brackets the fit.
- **Ca-alone thresholds are excluded.** They would leave single-bAP tLTD at ≤ 4 % of synapses.
- **Validation refit:** θ_eCB fixed at 48.3 × Ca_thr, and θ_V set as the C-Vf floor test decides. If both locks hold,
  only the 6 Chindemi parameters are free.

## 3. Next steps

1. Done: both logs have been read and seff is in the script header. The next MODE=V run needs `--mem=1900M`.
2. Ask the split1 owner to rerun `/scratch/dhuruva/split1/jobs/ext_l23l5.sh` (`--skip-existing`; 30 step-1 files plus
   the step-2 nulls with dt−500 remain). It needs more than 0:15; the 22136807 rate suggests ≈ 0:30. Then rerun MODE=V
   for the bAP-only arm.
3. For each lock that is accepted: a v5c validation refit with that threshold fixed, seeded from V5c_s7 plus one
   unseeded basin check (2g, 50G, 0:30, as for the V5c follow-ups).
4. Read the Hashimotodani 2005 [Ca]i values from Figs 2–4 (paper PDF) before any θ_eCB lock.

## θ_eCB literature status (orchestrator, 2026-10-01)
- Hashimotodani 2005 (Neuron 45:257) full text not reachable (publisher 403, no PMC copy). Secondary sources (Ohno-Shosaku & Kano 2014 review, PMC4237895) give only "submicromolar" for Ca2+-assisted RER and "micromolar" for CaER.
- Consistency bracket only (not a lock): submicromolar 0.1-1 µM × 48.3 units/µM = θ_eCB 4.8-48; fit 6.91 (0.14 µM) lies inside. Lock needs the paper's clamped [Ca]i (figure values); ask the user for the PDF.

## C-Vf result
- θ_V profile (22140352-55, 13:38-14:12, 38.3-39.1 GB of 50G): χ² at θ_V ×0.25/0.5/1/2/4 = 152.3/95.9/84.7/100.0/109.6 (θ_eCB refit 49.2/4.65/6.91/14.4/17.3). Not flat: θ_V is identified (Δχ² +11 / +15 at ×0.5 / ×2, approx 95% range 0.2-0.3), so it is a real fitted threshold and stays flagged (C-V null on og-delta and split1).
