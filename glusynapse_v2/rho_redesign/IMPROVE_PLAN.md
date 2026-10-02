# Plan to improve fit accuracy from S1C (2026-10-02)

The user picked S1C. This plan pulls together three diagnoses:
- S1C_L5_FAILURE.md (L5 validation misses);
- S1C_L23_FAILURE.md (both L2/3 pathways);
- STAGE3_DESIGN.md (the new core list, weights, acceptance tests and decision tree).

## Where S1C stands

| | fitted chi2 | validation chi2 |
|---|---|---|
| S1C | 10.3 over 7 L5 shape targets | L5 258 / 33, L2/3->L5 355 / 9, L2/3->L2/3 336 / 16 |

Markram +/-10 is hit. Both L2/3 pathways fail on sign: they predict LTP where the data show LTD. They do not fail by flattening. On the 20 core targets below, a flat model (1.0 everywhere) scores 498 and S1C scores 571.

## Why it fails (short)

1. **No postsynaptic rho-LTD.**
   - gamma_d sits at its lower bound (20.8), and the depression band closes wherever c_post is small.
   - All LTD is therefore eCB, which works only post-before-pre and bottoms out at 0.71.
   - So every pre-before-post LTD target fails: Letzkus 1AP +10, Zilberter 1AP +10 and the 5AP rows.
   - This is **parameter-limited**. The unseeded S1C run has the same L5 shape but an open band, and scores 108 against 355 on L2/3->L5.
2. **The eCB veto is too wide** (25 ms integral, plus (1-b) weighting).
   - At 20 Hz and above, every pre spike has a bAP within 25 ms of it, so eCB never fires.
   - So 20 Hz dt 0/+25 and 40 Hz dt 0 give LTP against LTD in the data.
   - The depolarising step in S07 also vetoes eCB, so its mGluR-block arm is too weak.
   - This is **structural**: refitting cannot fix it.
3. **LTP is capped at about 1.3-1.4** by the licence (50 Hz -10: 1.70 measured; S07 pair).
   - It is partly structural.
   - NO (S2C) raises S07, but the S2C validation was scored at the wrong A_NO/d_NO_max, a bug that is now fixed and being rescored.
4. **L2/3->L2/3 train LTD under AM251 and with mGluR** needs CB1-independent postsynaptic LTD.
   - No current mechanism gives it. This is **structural**.
5. **Input limits (emodel, not the rule).**
   - The distal bAP Ca at L2/3->L5 synapses is about 6x smaller than at L5 synapses.
   - The licence never opens at distal Letzkus synapses.

## Steps

### Step 0: running now (no decision needed)

Jobs:
- 22289873 and 22289874: the S1C rule refitted on 7 L5 + 23 L2/3 targets (seeded and unseeded). Do L2/3 targets in the fit keep Markram?
- 22289864: eCB variants (veto_T 15, unweighted W, A_eCB 0.2).
- 22289865: LTP-cap variants (gamma_p x2, licence off).
- 22290236-39: corrected S2C/S2N validation rescores.

### Step 1: stage 3 with L2/3 in the fit (20 core targets, STAGE3_DESIGN.md section 2)

**Core targets:**
- **L5 (10):** the 7 shape targets, plus 20 Hz -10 control, 20 Hz -10 mGluR block and S07 step pair.
- **L2/3->L5 (4):** Letzkus 1AP +10, 3AP +10 proximal, 3AP +10 distal and 3AP -10 proximal.
- **L2/3->L2/3 (6):** Zilberter 1AP +10 and 1AP -10, train10 +4, train10 -10 (control and mGluR block), and 5AP 10 Hz +10.

These are the worst-fit targets that carry distinct information. 18 are non-null and only 2 are nulls, so they can't pull the fit towards 1.

**Weights** (through `--weights`):
- a 0.05 SEM floor, so very tight SEMs don't dominate;
- the null drug arm at 0.5;
- Markram +10 at 2.

**Runs:**
- 3C: the S1C rule, 7 free, seeded from s1C_s and s2E_s, plus an unseeded run.
- 3N: with NO, 9 free, seeded from the corrected S2C_s.

**Sizing:** each fit 91G 0:45:00 on a GPU slice; validation rescore 84G 0:15:00.

**Hard pass:**
- Markram +/-10 within 1 SEM;
- at least 16 of the 18 non-null cores on the correct side of 1;
- slope of pred on data at least 0.6;
- prediction range at least 0.7.

The per-pathway limits are in STAGE3_DESIGN.md section 4.

### Step 2: eCB veto fix (new kernel gpu_v11_rho.py; v10 is not edited)

The change:
- Use a veto window (t_a+3, t_a+17] ms, read as a transient peak, in place of the 25 ms integral.
- Drop the (1-b) weighting.
- No new free parameter: both edges are bracketed by Sjostrom 2001 (+10 LTP, +25 LTD).

What it should fix: 20/40 Hz dt 0/+25 and the S07 eCB depth. A graded A_eCB is added only if the 22289864 diagnostic shows end-of-train leak.

The BCL then needs a new mod: GluSynapseV10.mod, with a new name.

### Step 3: LTP ceiling

- If 22289865 shows saturation, test the low-pass shaft licence G2. It is already in v10, and its threshold comes from the single-bAP calibration, with no new parameter.
- Otherwise keep NO (3N) if it buys more than 6 chi2 (BIC) on the shared cores.

### Step 4: L2/3->L2/3 postsynaptic LTD (only if steps 1-3 leave it failing)

- Add the new arm MVD (mGluR-gated VDCC depression; working name vdep, approved by the user 2026-10-02): depression from the synapse's own VDCC Ca within its own glutamate window. Per the review it is a band (theta_MVD_lo, theta_MVD_hi), not one-sided; kernel gpu_v11_rho.py, mod GluSynapseV10.mod (new SUFFIX).
- One parameter, anchored to Zilberter's D890 Ca ratio (0.37).
- Synapse-local and uniform.
- Zilberter 1AP +/-10 stays a uniformity conflict, because AM251 shows the equivalent L5 single pairings must not change. Accept about 20 chi2 or drop them; that is the user's decision.

### Step 5: inputs

- If distal Letzkus LTD still fails, the cause is the weak distal bAP Ca of the emodel. That goes to ion_fitter (new emodel name), not to a new rule parameter.
- Before that, test one re-parameterisation with no new parameter: seed a00 about 1.2, a01 about 2.

### Step 6: validation and BCL

- Rescore everything with plot_stage_fit.py (STAGE 3).
- Run the V9/V10 BCL last on the winner.

## Decision points for the user

- Zilberter 1AP +/-10: accept the conflict, or fit them at reduced weight?
- MVD approved (2026-10-02). Zilberter 1AP +/-10 kept in the core at weight 0.5 (orchestrator decision).
