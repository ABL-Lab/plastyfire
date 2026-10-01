# Which targets need each mechanism beyond Chindemi 2022 (MECH_NECESSITY)

Goal (user): "simpler plasticity that can explain all". Rule: only a00, a01, a10, a11, gamma_d and gamma_p are free. Every
added parameter needs a paper measurement or an independent calibration. Drives are synapse-local and Ca-based, and the
rule is uniform across pathways.

Script: `mech_necessity.py`, read-only.
- `--floors`: groups, chi2 per group, pharmacology floors. Uses CSVs only.
- `--pool`: the gate's own spine VDCC-Ca pool, read at own pre arrivals.

Jobs:
- 22135544: both parts. 0:52, 2.09 GB, 85% CPU. Log: `logs/mech_necessity_22135544.out`.
- 22135545: fill-in. Runs `--floors` after the ladder fits 22135019-22 and 22135205-8 (afterany).
- 22135725: sj04 timing pool. Runs after the extraction 22133664.

Reference fits:
- C1Ajn_s5: the current rule, 18 free, 39 targets.
- C1Ajz_pilot (22134473): the same parameters scored on all 54 targets.

## 0. Bottom line

| mechanism | needed by (forced by the data) | chi2 floor without it | new params | verdict |
|---|---|---|---|---|
| **M0: Chindemi core** | all targets | – | 6 free | keep |
| **G: spine VDCC-Ca potentiation gate** (C1) | 6 distal/proximal L2/3->L5 targets | ~60 (joint A0 L2/3 chi2 66-75 vs 5.1 with the gate) | theta_VDCC (free; only the requirement is anchored) + tau_E1 100 ms (needs an anchor) | **keep** |
| **E: eCB / CB1 presynaptic LTD** | 4 L5 control/AM251 pairs (8 rows) + sj04 | 45.0 (+33.4 with sj04) | theta_eCB + A_eCB. dpre_min -0.25 is anchored (Sjostrom 2003/2007). It reuses the gate's pool, tau_E1, i_scale and the (1-b) weight | **keep, simplified from 5 to 2 params** |
| NO-dependent presynaptic LTP | 1 row (Sj07 step-pair, L-NAME) | 6.4 | 4 free (A_NO, theta_NOi, tau_Z, theta_Z) + tau_NO anchored | **drop** (A_NO = 0). Costs AIC 8 / BIC 16 for at most 6.4 chi2. The code path stays for heterosynaptic / cluster work |
| rho_gamma | none | – | 1 with no anchor | keep at 1 (Chindemi) |

So the proposal is 6 Chindemi + theta_VDCC + theta_eCB + A_eCB = **9 parameters**, against 18 in C1Ajn_s5. tau_E1 and
dpre_min are anchored constants. One variant has 8: it ties theta_eCB to theta_VDCC, one VDCC-tethered sensor read by both
arms. It is worth one fit.

Two targets no monotone local rule can meet. They are listed in section 6:
- Zilberter train-LTD under AM251. It is a uniformity conflict and needs a user decision.
- sj04. It is probably an emodel limit.

## 1. The 54 targets grouped by phenomenon

Each target is counted once, in its primary group (`G` in mech_necessity.py; secondary tags are given there). The chi2 is
C1Ajz_pilot, i.e. C1Ajn parameters on all 54 targets.

| group | targets (data) | L5 / L2/3->L5 / L2/3->L2/3 | chi2 now |
|---|---|---|---|
| **LTP**: pre-before-post and its frequency dependence | Markram 10 Hz +5 (1.20), +10 (1.20); Sj +10 at 0.1 Hz (0.97, no LTP), 10 (1.16), 20 (1.31), 40 (1.53), 50 Hz (1.57); Sj -10 at 40 (1.51) and 50 Hz (1.70) (LTP at high f for either order); Sj07 step-pair control (1.62); Z 1AP +10 (**0.64, LTD**), Z 5AP 10 Hz +10 (**0.76, LTD**), Z 5AP 20 Hz +10 (1.07) | 10 / 0 / 3 | 23.9 |
| **LTD**: post-before-pre | Markram 10 Hz -10 (0.79); Sj 0.1 Hz -10 (0.69), -25 (0.65); Sj 10 Hz -10 (0.57), 20 Hz -10 (0.65); Z 1AP -10 (0.56), Z 5AP 20 Hz -10 (0.93) | 5 / 0 / 2 | 32.9 |
| **NMDAR block** | APV: Sj 0.1 Hz -10 (1.07), 20 Hz -10 (1.02), 50 Hz +10 (1.04); Letzkus 3AP +10 distal (1.03) | 3 / 1 / 0 | 2.3 (constant) |
| **AM251 (eCB)** | Markram 10 Hz -10 (1.07, Zilberter L5); Sj 0.1 Hz -10 (1.01), 20 Hz -10 (1.02): LTD abolished. Sj07 pair (2.13): LTP enhanced. Sj 50 Hz +10 (1.59): no effect. Z train +4 last (1.73 vs 1.49): n.s. Z train -10 last (0.73 vs 0.72): **no effect** | 5 / 0 / 2 | 34.5 |
| **NO** | Sj07 pair with L-NAME (1.36 vs 1.62) | 1 / 0 / 0 | 2.9 |
| **distal vs proximal** | S&H 50 Hz +10 distal (0.86), all (1.06); Letzkus 1AP +10 (0.72); 3AP +10 prox (1.30) / distal (0.79); 3AP -10 prox (0.89). Dropped: 3AP -10 distal (1.42) + APV | 0 / 6 / 0 | 3.2 |
| **burst / train** | Sj burst 5x20 Hz -120 (0.79), -200 (0.79); Z train10 +4 last (1.49), -4 last (0.99), -10 last (0.72), +5 before the 1st AP (0.97); train4 +4 (0.76), train8 +4 (1.15) | 2 / 0 / 6 | 57.1 |
| **no-change controls** | Sj 0.1 Hz -120 (1.05), -200 (1.05); Sj07 pre-only, post-only (1.00); Letzkus nopost (0.98), -500 (0.99); Z pre-only (0.98), train post-only (1.03) | 4 / 2 / 2 | 5.4 |
| **subthreshold-depolarisation LTD** (planned) | sj04: 4 LTD timings pooled (0.69 +- 0.04), AM251 + ifenprodil (1.06 +- 0.05), pre -100 control (~1.0) | 4-target set | not scored |

Totals: 13 + 7 + 4 + 7 + 1 + 6 + 8 + 8 = 54, and the chi2 sums to 162.2 (= 46.9 + 5.1 + 110.3). Most of the current misses
are on L2/3->L2/3: bursts 56.9, AM251 27.4, LTD 16.1.

## 2. Can pure Chindemi (M0) do each group, and what is logically necessary?

M0 is Chindemi 2022: rho with theta_d/theta_p = a * c_pre + a * c_post on effcai, with Ca from NMDA + VDCC, and no pre
side. Under M0 the drug switches act as follows:
- mglu_block and no_block give exactly the control prediction.
- nmdar_block freezes rho, so it gives 1.013 (L5) or 1.032 (L2/3).

The M0 fits (LM0_s5/s6, LM0z_s5/s6) are still running. So far only the MAXITER=0 checks exist, at C1Ajn's a's:
- LM0chk: 1866 over 39 targets.
- LM0zchk: 2501 over 54 targets.
These checks only show that eCB/NO are switched off. They are not the M0 optimum. Section 7 says what to fill in.

| group | M0 in principle? | necessary extra | evidence |
|---|---|---|---|
| LTP / frequency, L5 | **yes**: Chindemi fitted these sets with this structure | none | Sj 40/50 Hz -10 LTP (1.51/1.70) is missed by every variant tried (td3, td4, A0, C1: 1.10-1.25). No mechanism in this list addresses it; see section 6 |
| pre-before-post LTD, L2/3->L2/3 (Z 1AP +10 0.64, 5AP 10 Hz 0.76) vs L5 0.1 Hz +10 (0.97) | probably yes, through the per-synapse c_pre/c_post scaling of the thresholds | none proven; open until LM0z/LM1z | own pool at arrival ~0.003 theta_VDCC (22135544), so no eCB can act here; it must be rho LTD |
| post-before-pre LTD: sign | yes (rho between theta_d and theta_p) | – | |
| post-before-pre LTD: L5 pharmacology | **no**: control = AM251 in M0 | **CB1-dependent presynaptic factor (E)** | floors: Markram 10 Hz -10 10.9, Sj 0.1 Hz -10 15.8, Sj 20 Hz -10 10.5 (section 3) |
| NMDAR dependence | yes: block removes NMDA Ca. In the code it is a rho/dpre freeze | none | constant 2.31 chi2 in every model. The L5 APV result (pre NR2B, Sjostrom 2003) agrees with E acting only at own pre arrivals. The Zilberter APV-resistant LTD is out by user decision (non-NMDA Ca source) |
| AM251-sensitive LTD / LTP | **no** | **E** | 4 L5 pairs: floor 45.0 summed. sj04: floor 33.4 |
| AM251-insensitive arms | yes, trivially | – (and E must stay **silent** here) | Sj 50 Hz +10, Z train +4, Z train -10: floors 0.004 / 0.80 / 0.014. C1Ajn's eCB fires at Z train -10 (AM251 pred 1.06 vs 0.73, z 4.7, 22 chi2), so the current eCB form costs more than M0 here |
| NO-dependent LTP | no, for the L-NAME arm only | an NO-like presynaptic LTP | 1 row, floor 6.44 once E is in. Below the AIC/BIC cost of 4 params, so **not necessary** for the objective |
| distal vs proximal | **no** | **source-specific Ca signal on LTP (G)** | Same protocol, different sign: 50 Hz +10 gives 1.57 at L5->L5 vs 1.06 / 0.86 at L2/3->L5. Among synapses that cross theta_p, effcai/theta_p is identical between L2/3 should-depress and L5 should-potentiate (AUC 0.52; ROUND2). No sensor on total spine Ca separates them (CA_DECODE: window confound). Own VDCC charge does (AUC 0.20; 3AP distal vs prox 0.21). Requirement anchored by Letzkus Ni2+ (distal LTP), Zilberter D890 (train-LTP 0.57) and S&H bAP boosting |
| bursts, L5 (-120/-200 LTD, single -120/-200 none) | yes in principle (effcai summation), or through the dose of E | none beyond E | burst pool at arrival is only 1.6-1.7x the single-AP pool (burst -200 median 0.107 < single -120 0.134 theta_VDCC). So the burst/single split must come from 5 vs 1 arrivals per pairing (dose), not from the threshold |
| trains, L2/3->L2/3: AP-count threshold (4: 0.76, 8: 1.15, 10: 1.49) | yes in principle (effcai and the gate pool both sum over the train) | none new (G gives the VDCC requirement, D890) | not fitted yet; C1Ajn params give 0.69 / 0.85 / 0.89 |
| trains: order inside the train (Z12, pre before the 1st AP, 0.97 vs Z11 before the 10th, 1.49) | **doubtful** for any rule on effcai | a rule-form issue, not a new mechanism | effcai is a slow low-pass, so it is high late in the train in both cases. Z11 needs NMDA Ca to coincide with a high VDCC pool. C1 reads effcai, not the NMDA transient (C1Ajn: 1.16, z 3.1). Revisit only if LM1z leaves it as a top miss |
| sj04 subthreshold LTD | no (AM251 and ifenprodil abolish it) | **E**, with no AM251-insensitive post LTD (C2-like), at the Fig 4 timings | DLTD_DIAG. Strength is limited by the emodel (section 6) |

## 3. Pharmacology floors (data only, no model)

When a model cannot tell drug arms apart, its best chi2 for a set of arms is the inverse-variance pooled value:
floor = min over p of sum ((m_i - p) / s_i)^2.

| target | arms M0 cannot tell apart | data | floor M0 | floor M0 + E |
|---|---|---|---|---|
| Markram 10 Hz -10 | control = AM251 | 0.79 +- 0.03 / 1.07 +- 0.08 | **10.91** | 0 |
| Sj 0.1 Hz -10 | control = AM251 | 0.69 +- 0.07 / 1.01 +- 0.04 | **15.75** | 0 |
| Sj 20 Hz -10 | control = AM251 | 0.65 +- 0.09 / 1.02 +- 0.07 | **10.53** | 0 |
| Sj 50 Hz +10 | control = AM251 | 1.57 / 1.59 | 0.00 | 0 |
| Sj07 step-pair | control = AM251 = L-NAME | 1.62 / 2.13 / 1.36 | **14.29** | 6.44 (control = L-NAME) |
| Z train +4 last | control = AM251 | 1.49 / 1.73 | 0.80 | 0 |
| Z train -10 last | control = AM251 | 0.72 / 0.73 | 0.01 | 0 (only if E is silent there) |
| sj04 dLTD (planned) | control = AM251 | 0.69 +- 0.04 / 1.06 +- 0.05 | **33.39** | 0 |
| nmdar_block x4 | frozen in every model | | 2.31 | 2.31 |

What the floors imply:
- **M0 on 54 targets is at least 54.6** (52.30 + 2.31) before any control target is missed; 88.0 with sj04.
- **M0 on the 39-target ladder set is at least 53.8**, above C1Ajn's 51.97. So no M0 fit can match the current chi2.
  By AIC, the 6-parameter M0 is not excluded (AIC >= 65.8 vs 88.0 for C1Ajn), so the ladder fits decide that.
- **E is worth up to 45.0 chi2 (78.4 with sj04)** for 2-3 params, so it is clearly justified.
- **NO is worth at most 6.44.** It costs 4 unanchored params: AIC +8, BIC +16 (ln 54 = 3.99). Net negative even if it
  fit perfectly.

## 4. Simplest form of each kept mechanism

### G: potentiation gate (unchanged from C1)

The rule is c_VDCC' = -c_VDCC/tau_E1 + (-ica_VDCC)/i_scale, and pot = H(c* - theta_p) H(c_VDCC - theta_VDCC).
- theta_VDCC converges to 5.5-6.1 in every joint fit, so it is well constrained by the data.
- Its anchor is a requirement only (Ni2+ / D890), not a value. The candidate calibration is Zilberter Fig 5: the D890
  Ca reduction (train Ca x 0.37 +- 0.04) together with the 4-vs-8 AP threshold. This needs checking before it counts as
  independent of the fitted ratios.
- tau_E1 = 100 ms still needs an anchor (MINIMAL_LADDER).

### E: eCB-LTD reusing the gate's pool

At each own pre arrival, the step is dpre -= A_eCB when c_eCB(arrival) > theta_eCB. The floor dpre_min = -0.25 is
anchored. The pool is

c_eCB' = -c_eCB/tau_E1 + (1 - b)(-ica_VDCC)/i_scale

- It is the same current, tau_E1 and i_scale as G.
- b is the synapse's own glutamate-bound NMDA state (t_drive 4), so it adds no parameter. It encodes "VDCC Ca before own
  glutamate" (the Nevian order).
- The own-arrival condition carries the pre NR2B coincidence (Sjostrom 2003). It also keeps E dead under APV through the
  existing freeze.

Against the current eCB (theta_Te, tau_T, theta_Tg, A_mglu, dpre_min), the free count goes from 5 to 2. A kernel change
is needed (gpu_v4_rho eCB branch); not done here, as briefed.

**Pool test (22135544, C1Ajn_s5 params, theta_VDCC 5.50; 191 L5 / 115 L2/3->L2/3 synapses).** "w" is the (1-b)-weighted
pool. The column "w frac > 1" is the fraction of arrivals above 1 theta_VDCC.

| protocol | class | w pool at arrival, median (q25-q75), in theta_VDCC | w frac > 1 |
|---|---|---|---|
| Sj 0.1 Hz -10 | eCB+ | 0.36 (0.05-4.9) | 0.39 |
| Sj 20 Hz -10 | eCB+ | 1.09 (0.13-5.7) | 0.36 |
| Markram 10 Hz -10 | eCB+ | 0.82 (0.13-5.7) | 0.39 |
| Sj07 step-pair | eCB+ (AM251 enhances LTP) | 1.36 (0.19-4.8) | 0.30 |
| sj04 step, pre at onset | eCB+ | 0.009 (= resting level) | 0.00 |
| Sj 50 Hz +10 | eCB- | 0.22 (0.04-1.8); plain pool 3.9 | 0.19 |
| **Z train10 -10 last** | **eCB-** | **1.24 (0.24-6.5)** | **0.52** |
| Sj 0.1 Hz -120 / -200 | no change | 0.13 / 0.07 | 0.30 / 0.24 |
| Sj burst -120 / -200 | LTD | 0.23 / 0.11 | 0.36 / 0.28 |
| Sj 0.1 Hz +10, Z 1AP +10 | no eCB expected | 0.009 / 0.003 | 0.00 |

Findings:
1. **(1 - b) is what makes the order work.** The 50 Hz +10 pool drops 18x (3.9 to 0.22) and 10 Hz +10 drops 10x. The -10
   protocols keep most of their pool. So eCB-insensitivity at 50 Hz +10 comes from the weight, not from a threshold.
2. **Timing selectivity from an absolute threshold is weak.** Synapse-to-synapse spread (q25-q75 ~100x) swamps the
   e^-1.1 decay from -10 to -120 (frac > 1: 0.39 vs 0.30).
   - A threshold relative to the synapse's own single-AP pool would make the timing exact (-10 gives 0.90, -120 0.30,
     -200 0.14), in the same way Chindemi scales theta by c_post.
   - It needs one more per-synapse precomputed quantity (pool of one post AP, e.g. from `sjostrom_post_only_0.1hz`) but
     no extra parameter.
   - The L5 single -120/-200 (1.05 +- 0.10) constraint is soft, so try absolute first.
3. **The L5 burst -200 LTD vs single -120 no-change has to come from dose.** The burst pool is lower than the single-AP
   pool. The difference is 5 arrivals per pairing vs 1. So A_eCB is a real (non-saturating) parameter. An anchor
   candidate is the per-pairing LTD build-up in Sjostrom 2003, if printed.
4. **Tying theta_eCB = theta_VDCC is plausible on these numbers.** The eCB+ medians are 0.4-1.4 theta_VDCC. So variant
   E-tied (k = 8) is worth one fit.
5. **Z train -10 last has the largest pool of all (frac 0.52), yet AM251 has no effect there.** No monotone threshold on
   any VDCC pool can turn L5 -10 on and this off; see section 6.

## 5. Count of targets per mechanism vs parameter cost

| mechanism | rows that require it | rows it shapes | params (free) | chi2 at stake |
|---|---|---|---|---|
| M0 | – | 54 | 6 | – |
| G | 6 (DIST) | Z trains x6, Sj07, L5 +10 bursts | 1 (+ tau_E1 anchor) | ~60 (A0 vs C1 joint fits) |
| E | 8 + sj04 (2 scored rows) | all L5 -10 LTD, bursts, 50 Hz +10 | 2 (1 if tied) | >= 45.0 (+ 33.4) |
| NO | 1 | – | 4 | <= 6.4 |

## 6. Constraints no local rule can meet with this emodel (or this rule set)

1. **Letzkus 3AP -10 distal LTP (1.42) and its APV arm.**
   - It needs a tuft Ca spike or post-burst plateau. og-delta basal/oblique distal dendrites have none (DISTAL_LTP.md,
     22115348).
   - The one cell with a plateau gives 2.08 under the unchanged rule.
   - Emodel limit; already dropped (DROPT).
2. **sj04 dLTD (probably emodel).**
   - The 250 ms step at -52 mV gives -I_VDCC of about 3.4 K (DLTD_DIAG), i.e. a steady pool of about 3.4 x 2.469e-8 x
     100 / 1e-5 = 0.84, which is **0.15 theta_VDCC**.
   - That is the same as the median single-bAP pool 120 ms after the AP (0.13), which must give no LTD.
   - So in og-delta, any VDCC-Ca eCB drive that keeps single -120 silent leaves sj04 marginal, at most at the "pre inside
     the step" timings.
   - The paper points to low-threshold dendritic Ca channels (p3341), which this emodel barely opens below spike
     threshold.
   - Per-synapse check: 22135725, on the Fig 4 timings once 22133664 has extracted them.
3. **Z train-LTD under AM251 (0.73 = control) next to L5 10 Hz -10 under AM251 (1.07), both in Zilberter 2009.**
   - This is not an emodel limit. It is a uniformity conflict.
   - L5 tLTD is presynaptic and CB1-dependent. L2/3->L2/3 train-LTD is postsynaptic and CB1-independent, and its pool is
     the largest of all.
   - No synapse-local Ca quantity here separates them: both have 1 pre spike per pairing, and the order and b are the same.
   - Options, which need a user decision:
     - (a) keep both and accept about 20 chi2 at Z19;
     - (b) move Z18/Z19 to validation as "no CB1 at L2/3 PC -> L2/3 PC terminals" (a pathway property, against the
       uniformity rule);
     - (c) look for a pre-terminal-local discriminator, such as pre NR2B. None is known in our inputs.
4. **Sj 40/50 Hz -10 LTP (1.51/1.70): suspected emodel limit, not proven.**
   - All rule variants since td3 under-predict it (1.10-1.25).
   - High-frequency LTP for either order needs supralinear dendritic Ca during 50 Hz bursts (Sjostrom 2001
     cooperativity).
   - Diagnose it as in DISTAL_LTP before calling it an emodel limit.
5. **Not representable now** (simulation or condition limits, not rule limits): Z5 step + AP, Z6 extracellular, Z7/Z26
   Mg-free, Z24 D890, Z25 BAPTA, Z22/Z23 post-mGluR block, and Hardingham 2007. Zilberter APV arms Z20/Z21 are out by
   user decision.

Letzkus 1AP +10 (0.72) is reachable. C1's blocked-potentiation-to-depression conversion gives 0.73. Under C2 it would be
capped at about 0.9 (ROUND2).

## 7. To fill in when the ladder fits finish

Fill-in job 22135545 writes `logs/mech_necessity_fill_22135545.out`: chi2 per group per tag, and per-target pred / z
for LM0_s5/s6, LM1_s5/s6 (39 targets) and LM0z_s5/s6, LM1z_s5/s6 (54 targets). Read:
- **M0 eCB + LTD groups vs the floor.** If M0's chi2 there is near the floor (L5 37.2 + Sj07 14.3), post rho alone fits
  the control LTD and E is needed only for the drug arms.
- **M1 - M0 in DIST.** This quantifies G (expected: tens of chi2). Also check whether G moves the Z train / BURST group.
- **L2/3->L2/3 under M0z/M1z.** Do the train-count targets (Z11, Z16, Z17) and the anti-Hebbian single pairings (Z1, Z3)
  fit with Chindemi thresholds alone? This settles the open rows of section 2.

Rerun commands, if the log is missing:
- Floors: `sbatch --account=rrg-emuller -c1 --mem=1G -t 0:15:00 --wrap "source .venv/bin/activate && python glusynapse_v2/rho_redesign/mech_necessity.py --floors"`
- Pool: same, with `--pool --mem=2600M`.

Next rung after that is M2 = M1 + E, as in section 4: k 9, and k 8 with theta_eCB tied to theta_VDCC; NO off. It needs
the eCB branch in gpu_v4_rho to read the (1 - b) c_VDCC pool, and that kernel change is the fitter owner's to make.
