# Protocol audit: bath Ca/Mg, spike-evoking pulse, pairing count and rate per source paper (2026-10-02)

Question. Do our STDP protocols use each source paper's own extracellular Ca/Mg, post-spike pulse, pairing count and pairing
rate, and what must change? (User approval 2026-10-02: "set Ca, pulse shape and pairing counts per paper".) The emodel and its
age stay fixed, and the rule reads only Ca. Every change below is a bath or stimulus condition: it changes the Ca the rule sees
and adds no free parameter.

Short answer. Temperature, Mg and the pairing counts and rates match every paper, except one doubtful case (the Sjostrom 2003
burst count). Two mismatches matter:
1. **Ca_o.** Everything runs at 2.0 mM, but Sjostrom 2001 used **2.5 mM**. This affects its 17 L5->L5 control ids.
2. **Post pulse.** All L5 Sjostrom-lineage and Markram ids use **3 ms** pulses at the lowest firing amplitude on a coarse grid.
   The papers used **5 ms** pulses at fixed suprathreshold currents (0.8-1.8 nA). Letzkus uses 2 ms, as we do, but at 3-5 nA,
   while we search 1-8 nA.

The machinery to change Ca_o exists only in part. At 2.5 mM, release probability goes up x1.88 and NMDA Ca +22 %, while VDCC Ca
rises only about 3 %. The net sign on the targets therefore needs a pilot run.

## 1. How the conditions enter the pipeline (read-only code check)

- **Ca_o is hard-coded.** `extracellular_calcium: 2.0` and `GluSynapse.cao_CR: 2.0` are fixed in simwriter `_sim_config`
  (plastyfire/plastyfire/simwriter.py:757, 785). There is no per-yaml key.
  - simulator_edges.py:640-649 can override both.
  - The flag `--extracellular-calcium` exists in pairrunner_edges_fit.py:79 and precompute_cpre_cpost.py:283.
  - It is **absent from run_de_fit2_pool.py**, which is the prefire runner of the split1/split2 chains (spine/run_s2g0321.sh `pf()`).
- **NMDA Ca in bluecellulab.** Ca_o sets cao_CR_GluSynapse (bluecellulab/simulation/neuron_globals.py:29-33). That sets:
  - the NMDA Ca fraction `Pf = 4cao/(4cao+120/1.38)*0.6` (glusynapse_v2/build/GluSynapse.mod:324): 0.0505 at 2.0 mM and
    0.0619 at 2.5 mM (**+22 %**);
  - the spine-VDCC Nernst E_Ca (mod:329).
- **Release in bluecellulab.** Ca_o scales Use_d/Use_p by a constrained Hill function (bluecellulab/synapse/synapse_types.py:190-206,
  263-266). With u_hill 2.79 (ebner/audit_l23l5/FINDINGS.md:46) the factor is 1.000 at 2.0 mM and **1.876 at 2.5 mM**. Release
  saturates at Use = 1.
- **`cao0_ca_ion` is never set.** Dendritic Ca_HVA2/Ca_LVAst use a linear driving force with a Nernst eca (build/Ca_HVA2.mod:65,
  Ca_LVAst.mod:42), so shaft Ca stays at 2 mM. Setting cao0 to 2.5 would add 2.95 mV to E_Ca at 34 C, which is about **+3 %**
  VDCC current near the AP peak (inferred). Chindemi 2022 likewise changed VDCC only through E_Ca (Methods "L5-TTPC to L5-TTPC
  in low calcium", step b).
- **The c_pre/c_post cache is keyed by (pre, post) only** (simulator_edges.py:739). A 2.5 mM run with the existing cache keeps
  the 2 mM thresholds. The c_pre mini-sim forces Use = 1 (simulator_edges.py:283-285), so c_pre tracks Pf but not release.
- **The offline readout uses unscaled Use_d/Use_p** (glusynapse_v2/model_v2.py:71, 268-269; gpu_v3.py:331). The EPSP bases were
  simulated at 2.0 mM.
- **Mg** is 1 mM, the GluSynapse default (build/GluSynapse.mod:134). **Temperature:** no `celsius` is written (simwriter.py:784-790),
  so bluecellulab uses its default of 34.0 C (sonata_simulation_config.py:263-265).
- **Pulse.** Each yaml sets `width` and the grid `amp_min/amp_max/amp_lev`. `amp: find` takes the lowest grid amplitude that gives
  exactly one AP per pulse (simwriter.py:291-305). Pre spikes are replayed spike times, so the pre-cell pulse never matters.
- **Readout.** The prefire covers the pairing plus 1 s (simwriter.py:829). The model reads the induced state, not a time window.

## 2. Per-paper tables

Yaml shorthands (all in plastyfire/configs):

| shorthand | file |
|---|---|
| SAB | Sabrina_L5TTPC_L5TTPC_STDP.yaml |
| E5 | Ebner2019_L5TTPC_L5TTPC.yaml |
| E235 | Ebner2019_L23PC_L5TTPC.yaml |
| X5 | L5extra_split1.yaml |

The live copies are in /scratch/dhuruva/split1/configs (e.g. l5_sj01.yaml:42 has width 3.0).

**Markram et al. 1997, Science (doi 10.1126/science.275.5297.213); companion J Physiol 500:409 (doi 10.1113/jphysiol.1997.sp022031)**

| item | paper | ours | verdict |
|---|---|---|---|
| Ca/Mg | not verified (Science paywalled; PMC1159394 is a scan) | 2.0/1 | unverified (Chindemi 2022 used 2.0) |
| temp | 32-34 C (J Physiol abstract) | 34 C | match |
| pulse | 5 ms, ~1 nA per EBNER_HOW.md:123; not re-verified | 3 ms; find 0.5-3.0 nA, 10 levels (SAB:22-27) | mismatch if 5 ms is confirmed |
| pairings | 5+5 APs, 10 bursts every 4 s | nreps 10, T 4000 (SAB:28,30; X5:49-53) | match |
| readout | maximum deviation, 10-50 min (EBNER_HOW.md:124) | induced state | not modelled |

**Sjostrom, Turrigiano & Nelson 2001, Neuron 32:1149 (doi 10.1016/s0896-6273(01)00542-6), pp. 1161-1162 (lab PDF)**

| item | paper | ours | verdict |
|---|---|---|---|
| Ca/Mg | ACSF "MgCl2, 1; ... CaCl2, 2.5" | 2.0/1 (simwriter.py:757, 785) | **mismatch** |
| temp | 32-34 C | 34 C | match |
| pulse | "5-ms-long current injections (0.8-1.5 nA)" | 3 ms; find 0.5-3.0, 10 levels (E5:38-42; X5) | **mismatch** |
| pairings | 0.1 Hz: 50; >= 10 Hz: 5+5 APs x 15 at 0.1 Hz | nreps 50 / 15, T 10000 (E5:51-60; X5:55-61) | match |
| readout | mean from 10 min to end (>= 40 min) | induced state | not modelled |

**Sjostrom et al. 2003, Neuron 39:641 (doi 10.1016/s0896-6273(03)00476-8), p. 652 (local PDF)**

| item | paper | ours | verdict |
|---|---|---|---|
| Ca/Mg | "2 CaCl2 (unless otherwise specified)", 1 MgCl2 | 2.0/1 | match |
| temp | 32-34 C | 34 C | match |
| pulse | "as previously described (2001)", so 5 ms (inferred) | 3 ms (Sjostrom2003 yaml:48) | mismatch (inferred) |
| pairings | 0.1 Hz single pairings; Fig 9C burst count not stated | 50 at 0.1 Hz; bursts r50 (Sjostrom2003b yaml:56-57, the fitted ids) | 0.1 Hz match; burst: see note |
| readout | from 15 min after induction | induced state | not modelled |

Two notes on Sjostrom 2003:
- The drug targets #4-9 reuse the Sjostrom 2001 ids sjostrom_{0.1hz_dt-10, 20hz_dt-10, 50hz_dt+10} (targets.py:73-82). These
  must stay on 2.0 mM traces.
- Sjostrom 2004 (Methods p. 3339) says 20 Hz 5-AP bursts "were paired 15 times every 10 s, as previously described (Sjostrom et
  al. 2003, 2001)". So the lab convention is 15 bursts, and r50 is an assumption (inferred mismatch).

**Sjostrom et al. 2004, J Neurophysiol 92:3338 (doi 10.1152/jn.00376.2004), pp. 3338-3339.** Paper:
- 2 mM Ca / 1 mM Mg, 32-34 C;
- pre spikes by 5 ms pulses at 0.9-1.5 nA (replayed in our runs);
- 250 ms subthreshold post step, 50-60 pairings at 0.1-0.2 Hz;
- readout from 15 min.

Ours: 50 pairings at 0.1 Hz with a step search (Sjostrom2004 yaml:43-45). **Match.**

**Sjostrom et al. 2007, Neuropharmacology 52:176 (doi 10.1016/j.neuropharm.2006.07.021), p. 177.** Paper:
- 2 mM Ca / 1 mM Mg, 32-34 C;
- 200 ms steps at 0.5-1.8 nA (mean 1.2), 30 times every 10 s;
- readout from 15 min.

Ours: 30 steps every 10 s, step amplitude for 10-12 APs on a 0.1-3.0 nA grid (Sjostrom2007 yaml:54-56, 62-63). **Match.**

**Sjostrom & Hausser 2006, Neuron 51:227 (doi 10.1016/j.neuron.2006.06.017), p. 236 (local PDF)**

| item | paper | ours | verdict |
|---|---|---|---|
| Ca/Mg, temp | 2/1 mM; 32-35 C | 2.0/1; 34 C | match |
| pulse | "5 ms long current injections (1.0-1.8 nA)" | width 5 (E235:52); find 0.5-3.0, 10 levels | width match; amplitude window mismatch |
| pairings | 5 APs at 50 Hz, +10 ms, 15 times every 10 s | 15 x 10 s | match |

**Letzkus, Kampa & Stuart 2006, J Neurosci 26:10420 (doi 10.1523/JNEUROSCI.2650-06.2006), p. 10421 (local PDF)**

| item | paper | ours | verdict |
|---|---|---|---|
| Ca/Mg, temp | 2 CaCl2 / 1 MgCl2; 34-35 C | 2.0/1; 34 C | match |
| pulse | "brief somatic current injections (2 ms; 3-5 nA)" | width 2; find 1-8 nA, 150 levels (E235:59-67; nulls yaml:55-58) | width match; **amplitude mismatch** |
| pairings | "100-200 times at 1 Hz" | 100 at T 1000 | match (lower bound) |
| readout | mean at 20-30 min | induced state | not modelled |

**Zilberter et al. 2009, Cereb Cortex 19:2308 (doi 10.1093/cercor/bhn247), p. 2309 (local PDF)**

| item | paper | ours | verdict |
|---|---|---|---|
| Ca/Mg, temp | 2 CaCl2 / 1 MgCl2; 32-34 C | 2.0/1; 34 C | match |
| pulse | not stated | 3 ms, marked ASSUMED (Zilberter2009 yaml:3, 43) | not verifiable |
| pairings | "stimulated 40 times, every 5 s" | nreps 40, T 5000 (yaml:44-46) | match |
| readout | "5 min after the conditioning until the end" | induced state | not modelled (includes early potentiation) |

**Egger, Feldmeyer & Sakmann 1999, Nat Neurosci 2:1098 (doi 10.1038/16026)** (validation only). According to the PubMed abstract,
these are **L4 spiny stellate** pairs at P14, with 5 pre APs at 10 or 20 Hz. configs/L23L23extra_split1.yaml:50-52 labels them
L2/3 PC -> L2/3 PC, which is a pathway mismatch. Bath and pulse were not verified (paywalled).

## 3. Effect of each mismatch on the Ca the rule sees

**A. Ca_o 2.5 mM for the Sjostrom 2001 controls.**

Affected targets: core #15, #16, #20, #23, #24; validation #17-19, #21-22; the L5extra sjs01 rows.

Three routes change synapse-local Ca:
- NMDA Ca per release: +22 % (Pf).
- Release probability: x1.88, capped at Use = 1. That gives more NMDA Ca events at 0.1 Hz and more depletion within 10-50 Hz trains.
- VDCC and shaft Ca: about +3 %, and only if cao0_ca_ion is set.

The direction depends on how thresholds are set (all inferred):
- **Thresholds recomputed at 2.5 mM** (Chindemi 2022 step d). The Pf gain cancels and the release gain remains. Plasticity at
  0.1 Hz gets stronger: the -10 target (0.69) is helped and the +10 target (0.97) is at risk. The effect at 10-50 Hz is small and uncertain.
- **Thresholds fixed at 2 mM.** NMDA Ca rises 22 % against fixed thresholds, so everything moves towards LTP. This helps
  40/50 Hz -10 (targets 1.51, 1.70) and hurts the LTD targets #16 and #20 (0.69, 0.65).
- **Readout at 2.5 mM.** Use_p saturates at 1, so the range of presynaptic potentiation shrinks. Example for u = 0.5: Use_p goes
  from 0.87 to 1 and Use_d from 0.031 to 0.059. The predicted LTP drops, which works against the 1.70 target.

Data check (stated): at 50 Hz +10 the 2 mM control of Sjostrom 2003 Fig 7D is 1.64 +- 0.16 (n 16; targets.py:79-80), and the
2.5 mM Sjostrom 2001 value is 1.57 +- 0.26. So for this target the data show no large Ca_o effect.

What must be rerun:
- a new cache at 2.5 mM (if thresholds are recomputed);
- the prefire of the 17 ids (new hash, same workdirs);
- extraction;
- optionally the L5 EPSP basis at 2.5 mM.
- cexp does not need a rerun: its anchors are calibrations at the reference bath.

Code changes:
- `--extracellular-calcium` (plus optional `cao0_ca_ion`) in run_de_fit2_pool.py;
- a per-protocol Use scale in the readout;
- a separate `--extra` model in the fit (own extracted dir, cache and basis).

Kernel cost: keeping the 3 drug ids at 2.0 mM adds about 12 % to the L5 kernel work; keeping all 10 adds 34 % (TARGETS.md cost column).

**B. 5 ms post pulses at the paper currents.**

Which ids: all L5 Sjostrom-lineage ids, Markram if 5 ms is confirmed, and the S&H and Letzkus amplitude windows.

Effects (inferred unless marked):
- At threshold, a 5 ms pulse fires later and keeps injecting for about 2-3 ms after the AP. This leaves a somatic and proximal
  after-depolarisation.
- Markram 1997 (J Physiol) places 63 % of L5 contacts on basal dendrites at about 82 um. Pre-before-post pairings there would get
  more NMDA unblock, and the residual depolarisation summates at 10-50 Hz. Sjostrom 2001 Fig 5D/E (stated) shows LTP growing
  sigmoidally with that depolarisation.
- Net direction: more LTP at +10 for high frequencies; little change at -10 and 0.1 Hz.
- Doublet risk rises in bursting split2 cells (DOUBLETS_FITTING.md 1.1, 2). The user rule "doublet = one firing" applies.
- Letzkus: per the yaml note E5:82-83, low currents give late Ca-burst spikes. Moving to 3-5 nA probably removes them. That
  changes distal Ca after the burst and the number of guardrail failures.

What must be rerun:
- simwriter recalibration (new stimulus keys and new ids, e.g. `_w5`);
- prefire;
- extraction.
- The basis (pre-only test pulses) and the c_post cache (it has its own single-AP search) stay valid.

**C. Sjostrom 2003 bursts at r15 instead of r50.** Our integrators are slow. With 3.3x fewer pairings, the Fig 9C burst LTD
(target 0.79, core #14) gets weaker. Rerun the prefire and extraction of the 2 ids.

## 4. Cost (scaled from MEASURED lines; split1 values, split2 still pending at spine/run_s2g0321.sh:11)

| change | stage | measured basis | estimate |
|---|---|---|---|
| A | cache | cache5 22223209: 8 CPU, 2:17, 5.1 G (run_ljp25g0321.sh:3) | 0.3 CPU-h |
| A | prefire, 17 ids x 24 pairs | pfA 22223212: 3x16 CPU, max 31:48 (Sj01 ~89 % of bio time); pfx 22223215: 3x16 CPU, 19:53 | ~36 CPU-h |
| A | extraction | ext_sj01 ~13.5 min on 4 CPU (submit_s2g0321.sh:66), plus ext_l5x | ~2 CPU-h |
| A (opt.) | L5 EPSP basis | basis5 22136796: 24 x 12 CPU, max 18:11, 9.9 GiB | <= 87 CPU-h |
| B, L5 | simwriter | sims on 8 CPU: mk 0:41, sj01 2:40, sj03 0:03, r50 0:18, extra 2:46 (run_s2g0321.sh:9) | ~1 CPU-h |
| B, L5 | prefire, all pulse ids | pfA + pfB 22223213 (3x32 CPU, 29:46, 189 G) + pfx | ~90 CPU-h (A rides along) |
| B, L5 | extraction | ext_l5 22223221: 1:04:06 on 4 CPU, 57 G | ~5 CPU-h |
| B, L2/3->L5 | sim + prefire + ext | sim l23l5 5:06; pf23_5 22223216 (32 CPU, 36:27, 169.5 G); pf23_5n 22223217 (36:33); ext23_5 20:54 | ~35 CPU-h |
| C | prefire + ext, 2 ids | share of pfB (r15 = 300 s of bio time) | ~3 CPU-h |

Totals: A+B+C come to about 135 CPU-h without the basis and about 220 CPU-h with it. Each prefire step exceeds 30 CPU-h, so it
must go in the heavy-run log before submission. Outputs go to /scratch (inode quota). Use the +25 % memory values next to the
measured lines; several of those runs were at their limit.

## 5. Ranked changes (XS < 5, S < 40, M < 100, L > 100 CPU-h)

1. **[S + code] Ca_o 2.5 mM pilot on the 17 Sjostrom 2001 ids.** Keep the 2.0 mM traces for the 3 drug ids. Steps:
   - plumb the flag into run_de_fit2_pool.py;
   - build a 2.5 mM cache;
   - apply the Use x1.876 scale in the readout;
   - evaluate the current best parameters without a refit.
2. **[M, same run as 1] 5 ms post pulses for every L5 Sjostrom-lineage id**, and for Markram once 5 ms is confirmed.
   Amplitude: the lowest single-AP amplitude inside 0.8-1.5 nA, searched in 0.05 nA steps. Log the cells that have no valid
   amplitude in that window.
3. **[S] L2/3->L5 amplitude windows.** Letzkus: 2 ms at 3-5 nA (1AP, 3AP and nulls). S&H: 5 ms at 1.0-1.8 nA.
4. **[XS] Sjostrom 2003 Fig 9C bursts at 15 repetitions** (lab convention). Keep r50 as a sensitivity run.
5. **[M] L5 EPSP basis at 2.5 mM**, only if step 1 shows that Use saturation moves a 2.5 mM target by more than ~0.5 SEM.
6. **[0] Readout windows.** Record Markram's maximum-deviation readout and Zilberter's 5-min start in the targets; the model has
   no post-induction time course to match.
7. **[0] Egger yaml pathway label:** it is L4 spiny stellate (validation only).

No change needed for temperature (34 C), Mg (1 mM), or the pairing counts and rates of Markram, Sjostrom 2001/2004/2007,
S&H 2006, Letzkus and Zilberter.

## 6. Decisions for the user

- **(a) Thresholds at 2.5 mM.** Two options:
  - recompute c_pre/c_post at the bath Ca (Chindemi 2022, step d), with a00..a11 defined under the experiment's conditions
    (recommended);
  - keep the 2 mM thresholds (run as a sensitivity check).
- **(b) Bath Ca on the dendritic VDCCs** via `cao0_ca_ion`. This is a bath condition, not a channel change, and adds about +3 %.
  Recommended yes, to match Chindemi step b.
- **(c) Amplitude policy.** Either the lowest firing current inside the paper window (recommended) or a fixed current in the
  middle of the window.
- **(d) Burst count.** r15 (lab convention) or r50.

## Open points

- Markram 1997 Ca/Mg and pulse are not verified (Science is paywalled; the PMC copy is a scan).
- Sjostrom 2003 re-analysed 37 pairs from Sjostrom 2001 by CV (LTD_LIT.md:9), so some of its numbers may come from 2.5 mM cells.
  The drug and timing targets come from new 2003 experiments ("2 CaCl2 unless otherwise specified").
- The per-cell calibrated amplitudes (single_cells pkls) were not inspected; that needs a small Slurm job. The share of Letzkus
  cells already inside 3-5 nA is unknown.
- Zilberter's pulse width is not in the paper, so 3 ms remains an assumption.

References (metadata and DOIs via PubMed):
- Markram 1997 Science, 10.1126/science.275.5297.213
- Markram 1997 J Physiol, 10.1113/jphysiol.1997.sp022031
- Sjostrom 2001 Neuron, 10.1016/s0896-6273(01)00542-6
- Sjostrom 2003 Neuron, 10.1016/s0896-6273(03)00476-8
- Sjostrom 2004 J Neurophysiol, 10.1152/jn.00376.2004
- Sjostrom 2007 Neuropharmacology, 10.1016/j.neuropharm.2006.07.021
- Sjostrom & Hausser 2006 Neuron, 10.1016/j.neuron.2006.06.017
- Letzkus 2006 J Neurosci, 10.1523/JNEUROSCI.2650-06.2006
- Zilberter 2009 Cereb Cortex, 10.1093/cercor/bhn247
- Egger 1999 Nat Neurosci, 10.1038/16026
- Chindemi 2022 Nat Commun (papers/ local md, Methods "L5-TTPC to L5-TTPC in low calcium")
