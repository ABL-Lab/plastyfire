# What triggers postsynaptic eCB production for tLTD: Ca or voltage? (lit check, 2026-09-30)

Sources: local PDFs in `plastyfire/papers/` (read in full: Sjöström 2003, 2004, 2007; Nevian & Sakmann 2006; Letzkus 2006; Sjöström & Häusser 2006). Everything else was checked against PubMed metadata and abstracts, and DOIs were checked there too. "Abstract only" means I could not read the full text. **UNVERIFIED** marks a claim I could not check.

## Verdict
- **eCB production is Ca-triggered.** Every tLTD/DSI/DSE preparation that was tested is blocked by postsynaptic Ca chelation. Uncaging Ca with no depolarisation is enough to evoke DSI. In L2/3 the Ca must come from VDCCs, and it must arrive *before* mGluR→PLC activation.
- **No paper I found shows Ca-independent, depolarisation-driven eCB release.** The one Ca-independent route is strong mGluR1→PLCβ4 activation (Maejima 2001). It is driven by agonist, not by voltage.
- **Voltage acts upstream only.** It opens VDCCs, and possibly biases the agonist affinity of mGluR1a/mAChR. That GPCR voltage-sensitivity has been shown only in *Xenopus* oocytes and has never been tested for eCB release or tLTD.
- **So the user's criterion ("accept voltage only if eCB release works via voltage") is not met.** A local-V event counter (D1) can be kept only as a phenomenological proxy for "a VDCC Ca event occurred", not as the mechanism.
- **Line 49 of LOCAL_T_DRIVE.md is wrong on three points:**
  - (a) A *fast* depolarisation is not needed. At L5-L5, 250 ms subthreshold steps substitute fully for spikes (Sjöström 2004), and Ca uncaging alone works (Wang & Zucker 2001).
  - (b) The Ca sensor is **not in a VDCC nanodomain**. EGTA blocks as well as BAPTA does (Nevian 2006: IC50 0.39 vs 0.36 mM). DSI is "activated by micromolar [Ca2+]i acting far from sites of Ca2+ entry" (Wang & Zucker 2001).
  - (c) DSI scales ~linearly with integrated Ca influx (log-log slope ~1, Lenz & Alger 1999). It does not saturate per event.
  - Also, at L5-L5 the pharmacology points to an FAAH-degraded eCB (AEA-like). PLCβ/DAGLα were never tested at L5-L5.
- **The distal failure of Ca drives may be correct biology, not a bug.** Most L5-L5 contacts are proximal/basal, and the recorded tLTD connections were proximal. Distal L5 inputs follow different rules (Sjöström & Häusser 2006; Letzkus 2006).

## Q1. Ca, voltage, or both?
**Evidence for Ca (strong, direct, many preparations)**
- **L5-L5 tLTD:** 10 mM postsynaptic BAPTA blocks tLTD but not CB1-agonist LTD (Sjöström 2003, Fig. 4F; conc. from Methods).
- **L2/3 tLTD:** 1 mM EGTA or 1 mM BAPTA blocks LTD. Half-block is at 0.39 mM EGTA and 0.36 mM BAPTA (Nevian & Sakmann 2006, Fig. 4C-D).
- **L4→L2/3 tLTD:** requires "calcium from voltage-sensitive channels and IP3 receptor-gated stores" (Bender 2006, abstract).
- **DSI (CA1):** blocked by 10 mM EGTA or BAPTA. DSI tracks total Ca influx with a log-log slope of ~1 (Lenz & Alger 1999). Ca uncaging with nitrophenyl-EGTA (PSI) mimics DSI, and both are linear in [Ca]i (Wang & Zucker 2001).
- **DSE (cerebellum):** abolished by postsynaptic BAPTA (Kreitzer & Regehr 2001). Depolarised hippocampal neurons release AEA and 2-AG "in a Ca2+-dependent manner" (Wilson & Nicoll 2001, abstract).
- **PLCβ1 as coincidence detector:** receptor-driven eCB release depends on physiological [Ca]i and is strongly boosted by depolarisation-induced Ca. It is absent in PLCβ1 KO (Hashimotodani 2005, abstract). The quantitative Ca dependence is UNVERIFIED (no full text).
- **Enzyme:** DGLα KO abolishes retrograde eCB suppression in cerebellum, hippocampus and striatum, so 2-AG made by DGLα is the retrograde messenger there (Tanimura 2010, abstract). Whether DGLα itself is Ca-dependent is UNVERIFIED.

**Evidence for voltage (weak, indirect)**
- **Sjöström 2004 (L5-L5 dLTD):** subthreshold 250 ms steps to about −52 mV, paired with pre spikes, give LTD that is CB1- and NR2B-dependent and indistinguishable from tLTD. BAPTA was *not* tested. The authors attribute it to Ca through low-threshold Ca channels ("which may provide sufficient Ca2+ influx for endocannabinoid release"). This fits Ca and voltage equally well.
- **Ohana 2006:** mGluR1a apparent glutamate affinity *increases* with depolarisation. Ben-Chaim 2006: m1/m2 muscarinic receptors carry gating charge. Both were shown in oocytes. Neither tested neurons, eCB or LTD, and both need agonist.
- **Maejima 2001:** mGluR1-driven eCB release in Purkinje cells without a postsynaptic Ca rise. This is Ca-independent but *agonist*-driven, not voltage-driven.
- **Maejima 2005:** weak mGluR1 plus mild depolarisation gives 2-AG. The paper frames it as "Ca2+-assisted" (submicromolar Ca), i.e. voltage acting through Ca.

**Which is stronger:** Ca, by a wide margin. The voltage evidence is either Ca-confounded (Sjöström 2004, Maejima 2005) or heterologous and agonist-dependent (Ohana, Ben-Chaim).

## Q2. Ca source, compartment, numbers
- **VDCC, subtype depends on the protocol (L2/3, Nevian 2006):**
  - Single-AP tLTD (−10 ms) is abolished by 50 µM Ni2+ (T-type).
  - Burst (3 AP, 50 Hz) tLTD survives Ni2+ alone and nimodipine alone, but is blocked by Ni2+ plus 10 µM nimodipine.
  - Nimodipine cuts the spine Ca transient for post→pre to 69% but leaves LTD intact.
  - Heparin (IP3R) has no effect. MCPG, U73122 and AM251 each block LTD.
  - Order matters: VDCC Ca must precede mGluR activation.
- **Stores:** needed at L4→L2/3 (Bender 2006, abstract: IP3R-gated stores), not needed at L2/3 (Nevian 2006, heparin). The specific VSCC subtype in Bender 2006 is UNVERIFIED (full text not accessible; Nevian cites it for T-type).
- **NMDAR:** postsynaptic NMDAR not required for tLTD. Evidence: BAPTA spares CB1-agonist LTD and NMDAR dependence is presynaptic (Sjöström 2003); MK-801 loaded pre vs post (Rodríguez-Moreno & Paulsen 2008); Bender 2006.
- **DSI (CA1):** N-type (ω-conotoxin GVIA) almost abolishes it. L-type contributes only with unclamped spikes. Ni2+ 100 µM and ω-agatoxin TK have no effect (Lenz, Wagner & Alger 1998).
- **Compartment:** Nevian measured spine Ca, and LTD tracked volume-averaged spine Ca. EGTA ≈ BAPTA means the sensor is far from the channel mouth (not a nanodomain). Rancz & Häusser 2006: local dendritic Ca spikes trigger eCB release with spatial specificity (abstract).
- **Numbers (all verified from abstracts or full text):**

| System | Ca for half-max eCB effect | Source |
|---|---|---|
| CA1 DSI and uncaging PSI | ~3.6–3.9 µM [Ca]i | Wang & Zucker 2001 |
| Purkinje DSE/DSI | ~15 µM | Brenowitz & Regehr 2003 |
| Purkinje, with weak mGluR1 | submicromolar | Maejima 2005 |
| PF + CF, local mGluR | mGluR "locally reduced the dendritic calcium levels required" | Brenowitz & Regehr 2005 |
| L2/3 tLTD | LTD threshold ≈ ½ the LTP threshold (relative ΔG/R, no absolute value) | Nevian 2006 |

- **Time course:** DSI decays with τ ≈ 20 s, independent of buffer (Lenz & Alger 1999), and outlasts the Ca rise (Wang & Zucker 2001). The L5 tLTD window widens with FAAH block (AA-5-HT), with AEA-uptake block (AM404) and with postsynaptic bursting (Sjöström 2003, Fig. 9). So eCB degradation, not Ca decay, sets the window.

## Q3. Pathway-specific postsynaptic requirements
- **L5→L5, Sjöström 2003:**
  - Postsynaptic BAPTA (10 mM) blocks tLTD.
  - CB1 required (AM251). NR2B/presynaptic NMDAR required.
  - AEA and ACEA mimic tLTD and occlude it.
  - The mGluR antagonist LY341495 was tested **only on ACEA-LTD**, not on tLTD, so a postsynaptic mGluR/PLC role at L5-L5 is **untested**.
  - No VDCC blocker, PLC or DAGL test.
- **L5→L5, Sjöström 2004:** subthreshold depolarisation substitutes for spikes. CB1- and NR2B-dependent. The single-connection EPSP alone is not enough. The synapses were proximal (20–80% rise time 1.3 ms), and the authors note that "synapses more distally ... may not be able to undergo dLTD".
- **L5→L5, Sjöström 2007:** AM251 during high-frequency pairing enhances LTP, so a presynaptic eCB-LTD runs alongside LTP. No postsynaptic Ca/VDCC/PLC pharmacology.
- **L4→L2/3, Bender 2006 (abstract):** group I mGluR + VSCC Ca + IP3R stores + CB1 + presynaptic NMDAR. Postsynaptic NMDAR not needed. Coincidence window ~125 ms for LTD vs ~25 ms for LTP.
- **L2/3 (extracellular stim 50–150 µm on basal dendrites), Nevian 2006:**
  - EGTA/BAPTA-sensitive.
  - T-type VDCC for single APs; T+L for bursts.
  - MCPG, U73122 and AM251 block. Heparin has no effect.
  - Ca must precede mGluR activation.
  - This is validation only (not a paired recording).

## Q4. Distance
- **Markram 1997 (L5-L5 anatomy):** 63% of contacts are basal at 82 ± 35 µm, 27% are apical oblique at 145 ± 59 µm, and the mean per connection is 147 µm. Sjöström & Häusser 2006 say "the vast majority of L5-to-L5 synapses are on the proximal, basal dendritic tree", with L5-L5 rise time 2.2 ± 0.2 ms.
- **Sjöström & Häusser 2006:**
  - AP-EPSP pairing (pre→post) gives LTP at proximal inputs and **LTD at distal ones** (rise time > 3 ms; 80%, n = 28). Distal LTD switches to LTP with cooperativity or dendritic depolarisation that boosts the bAP.
  - Unpaired 50 Hz distal EPSPs with no bAP give **AM251-sensitive LTD** (71%, n = 10). So distally, synaptic depolarisation or Ca can drive eCB-LTD *without* a bAP.
  - Post→pre tLTD was not mapped against distance.
- **Letzkus 2006 (L2/3→L5):**
  - Single APs at +10 give distance-dependent LTD.
  - With bursts, **distal inputs get LTP at −10** (blocked by APV and Ni2+) and LTD at +10 (APV-sensitive, Ni2+-insensitive).
  - eCB was not tested.
  - So distally the post→pre LTD window disappears or reverses.
- **Froemke 2005 (L2/3 visual, apical):** the LTD window is *broader* distally, and its mechanism is AP-induced NMDAR suppression, not eCB (CB1 not reported in the abstract).
- **Implication:**
  - The data do not require every synapse at >150 µm to show the Sjöström 2003 post→pre window. The measured window comes from mostly proximal contacts.
  - At distal sites the literature shows *different* rules (LTP at −10 with bursts; eCB-LTD driven by local EPSP trains).
  - A Ca drive that loses the bAP signal distally is consistent with this. Scoring "fraction of synapses" over all 191 is stricter than the data.
  - Score at the **connection level** (somatic EPSP change summed over a pair's contacts) instead.

## Q5. Model implication (synapse-local only)
- **Drive the eCB trace with synapse-local VDCC Ca**, spine or shaft volume-averaged, not a nanodomain. Use `ica_VDCC`-derived Ca, or cai_CR with the NMDAR share excluded, since postsynaptic NMDAR is not needed.
  - Integrate linearly above a threshold (Lenz & Alger slope ~1). Do not count events.
  - A threshold is justified (Nevian: an LTD threshold exists, at about ½ the LTP one). No absolute value transfers from the data, so it stays a fitted parameter.
- **Coincidence with the synapse's own glutamate (mGluR→PLCβ), order-sensitive:**
  - VDCC Ca must *precede* the pre spike (Nevian 2006; Hashimotodani 2005, PLCβ1 detector; Maejima 2005, Ca-assisted).
  - This replaces the phenomenological 15 ms own-pre veto with a mechanism: a PLC priming trace P(t) set by recent local Ca, read out when the pre spike arrives. That is what the gate already does (T read at pre arrival).
  - Caveat: mGluR involvement at L5-L5 is untested (above), so keep it pathway-uniform and generic.
- **Time course:** the eCB/T trace should decay on a time scale set by degradation (the tLTD window widens with FAAH block). Induction window of tens of ms at L5-L5 (LTD at −25/−50, none at −100 in the fitted data), ~125 ms at L4→L2/3 (Bender 2006). DSI's τ ≈ 20 s is the retrograde-effect decay, not the induction window.
- **What to test next (no new data needed):**
  - (i) Re-score B1/C1 Ca drives at the connection level on proximal-dominated pairs.
  - (ii) Check that the chosen drive gives LTD for Sjöström 2004 dLTD (250 ms subthreshold step + pre) and for distal unpaired 50 Hz (Sjöström & Häusser 2006). These are the two protocols that discriminate the drives. D1 with its own-pre veto probably fails the second (UNVERIFIED, not simulated).
- **If D1 is kept anyway,** document it as "a VDCC-opening proxy, phenomenological". Do not cite 2-AG/DAGLα/nanodomain as its basis.

## Papers
| Paper | Key finding used | DOI |
|---|---|---|
| Sjöström, Turrigiano, Nelson 2003 Neuron 39:641 | L5-L5 tLTD: post BAPTA 10 mM blocks; CB1 + pre NMDAR; FAAH block widens window | 10.1016/s0896-6273(03)00476-8 |
| Sjöström, Turrigiano, Nelson 2004 J Neurophysiol 92:3338 | Subthreshold depolarisation substitutes for spikes; CB1/NR2B; proximal synapses | 10.1152/jn.00376.2004 |
| Sjöström, Turrigiano, Nelson 2007 Neuropharmacology 52:176 | AM251 unmasks LTP during HF pairing; no post Ca pharmacology | 10.1016/j.neuropharm.2006.07.021 |
| Nevian & Sakmann 2006 J Neurosci 26:11001 | L2/3 tLTD: EGTA≈BAPTA; T-VDCC (1 AP), T+L (burst); mGluR/PLC/CB1; no IP3R | 10.1523/JNEUROSCI.1749-06.2006 |
| Bender, Bender, Brasier, Feldman 2006 J Neurosci 26:4166 | L4→L2/3 tLTD: mGluR + VSCC Ca + IP3R stores, CB1 (abstract only) | 10.1523/JNEUROSCI.0176-06.2006 |
| Rodríguez-Moreno & Paulsen 2008 Nat Neurosci 11:744 | tLTD needs presynaptic, not postsynaptic, NMDAR (L4→L2/3) | 10.1038/nn.2125 |
| Sjöström & Häusser 2006 Neuron 51:227 | Distal inputs: LTD at pre→post; AM251-sensitive LTD from 50 Hz EPSPs without bAP | 10.1016/j.neuron.2006.06.017 |
| Letzkus, Kampa, Stuart 2006 J Neurosci 26:10420 | Distal L2/3→L5: LTP at −10 (bursts), LTD at +10; location-dependent rules | 10.1523/JNEUROSCI.2650-06.2006 |
| Froemke, Poo, Dan 2005 Nature 434:221 | Distal LTD window broader; NMDAR-suppression mechanism | 10.1038/nature03366 |
| Markram et al. 1997 J Physiol 500:409 | L5-L5 contacts: 63% basal (82 µm), 27% oblique (145 µm) | 10.1113/jphysiol.1997.sp022031 |
| Hashimotodani et al. 2005 Neuron 45:257 | PLCβ1 = Ca-dependent coincidence detector for eCB release (abstract) | 10.1016/j.neuron.2005.01.004 |
| Maejima et al. 2001 Neuron 31:463 | mGluR1→eCB release without postsynaptic Ca rise (agonist-driven) | 10.1016/s0896-6273(01)00375-0 |
| Maejima et al. 2005 J Neurosci 25:6826 | Three modes; Ca-assisted mGluR1-PLCβ4 (submicromolar Ca) makes 2-AG | 10.1523/JNEUROSCI.0945-05.2005 |
| Tanimura et al. 2010 Neuron 65:320 | DGLα-made 2-AG mediates retrograde suppression | 10.1016/j.neuron.2010.01.021 |
| Wang & Zucker 2001 J Physiol 533:757 | Ca uncaging (no depolarisation) gives DSI-like PSI; half-max 3.6–3.9 µM; far from channels | 10.1111/j.1469-7793.2001.t01-1-00757.x |
| Lenz & Alger 1999 J Physiol 521:147 | DSI ∝ Ca influx (slope ~1); 10 mM EGTA/BAPTA block; τ ≈ 20 s | 10.1111/j.1469-7793.1999.00147.x |
| Lenz, Wagner, Alger 1998 J Physiol 512:61 | DSI via N-type (and L-type with spikes); Ni2+ no effect | 10.1111/j.1469-7793.1998.061bf.x |
| Brenowitz & Regehr 2003 J Neurosci 23:6373 | Purkinje: ~15 µM Ca for half-max eCB inhibition | 10.1523/JNEUROSCI.23-15-06373.2003 |
| Brenowitz & Regehr 2005 Neuron 45:419 | Local mGluR lowers the Ca needed; supralinear local Ca | 10.1016/j.neuron.2004.12.045 |
| Kreitzer & Regehr 2001 Neuron 29:717 | DSE abolished by postsynaptic BAPTA | 10.1016/s0896-6273(01)00246-x |
| Wilson & Nicoll 2001 Nature 410:588 | Ca-dependent AEA/2-AG release mediates DSI | 10.1038/35069076 |
| Rancz & Häusser 2006 J Neurosci 26:5428 | Local dendritic Ca spikes trigger eCB release (abstract) | 10.1523/JNEUROSCI.5284-05.2006 |
| Ohana et al. 2006 J Biol Chem 281:24204 | mGluR1a/mGluR3 agonist affinity is voltage-sensitive (oocytes) | 10.1074/jbc.M513447200 |
| Ben-Chaim et al. 2006 Nature 444:106 | m1/m2 mAChR gating charge coupled to agonist affinity (oocytes) | 10.1038/nature05259 |
