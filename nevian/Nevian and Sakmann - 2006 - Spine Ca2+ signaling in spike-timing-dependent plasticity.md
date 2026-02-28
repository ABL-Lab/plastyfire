The Journal of Neuroscience, October 25, 2006 • 26(43):11001–11013 • 11001

Cellular/Molecular

# Spine Ca<sup>2+</sup> Signaling in Spike-Timing-Dependent Plasticity

**Thomas Nevian and Bert Sakmann**
Department of Cell Physiology, Max-Planck Institute for Medical Research, D-69120 Heidelberg, Germany

Calcium is a second messenger, which can trigger the modification of synaptic efficacy. We investigated the question of whether a differential rise in postsynaptic Ca<sup>2+</sup> ([Ca<sup>2+</sup>]<sub>i</sub>) alone is sufficient to account for the induction of long-term potentiation (LTP) and long-term depression (LTD) of EPSPs in the basal dendrites of layer 2/3 pyramidal neurons of the somatosensory cortex. Volume-averaged [Ca<sup>2+</sup>]<sub>i</sub> transients were measured in spines of the basal dendritic arbor for spike-timing-dependent plasticity induction protocols. The rise in [Ca<sup>2+</sup>]<sub>i</sub> was uncorrelated to the direction of the change in synaptic efficacy, because several pairing protocols evoked similar spine [Ca<sup>2+</sup>]<sub>i</sub> transients but resulted in either LTP or LTD. The sequence dependence of near-coincident presynaptic and postsynaptic activity on the direction of changes in synaptic strength suggested that LTP and LTD were induced by two processes, which were controlled separately by postsynaptic [Ca<sup>2+</sup>]<sub>i</sub> levels. Activation of voltage-dependent Ca<sup>2+</sup> channels before metabotropic glutamate receptors (mGluRs) resulted in the phospholipase C-dependent (PLC-dependent) synthesis of endocannabinoids, which acted as a retrograde messenger to induce LTD. LTP required a large [Ca<sup>2+</sup>]<sub>i</sub> transient evoked by NMDA receptor activation. Blocking mGluRs abolished the induction of LTD and uncovered the Ca<sup>2+</sup>-dependent induction of LTP.

We conclude that the volume-averaged peak elevation of [Ca<sup>2+</sup>]<sub>i</sub> in spines of layer 2/3 pyramids determines the magnitude of long-term changes in synaptic efficacy. The direction of the change is controlled, however, via a mGluR-coupled signaling cascade. mGluRs act in conjunction with PLC as sequence-sensitive coincidence detectors when postsynaptic precede presynaptic action potentials to induce LTD. Thus presumably two different Ca<sup>2+</sup> sensors in spines control the induction of spike-timing-dependent synaptic plasticity.

*Key words:* LTP; LTD; synaptic plasticity; calcium; two-photon microscopy; spine; spike-timing-dependent plasticity; mGluR; NMDAR

## Introduction
Long-term changes in synaptic efficacy are thought to be the cellular basis of information storage and memory formation (Bliss and Collingridge, 1993; Whitlock et al., 2006). Modifications in the efficacy of transmission at synaptic contacts can be induced by coincident presynaptic and postsynaptic activity (Magee and Johnston, 1997; Markram et al., 1997; Debanne et al., 1998). The precise timing and the order of presynaptic and postsynaptic action potentials (APs) determine the magnitude and the direction of the change in synaptic strength (Markram et al., 1997; Bi and Poo, 1998; Feldman, 2000; Sjostrom et al., 2001; Froemke and Dan, 2002). A postsynaptic AP that follows a presynaptic AP within a time window of tens of milliseconds results in long-term potentiation (LTP), whereas the reverse order results in depression (LTD). Therefore, spike-timing-dependent plasticity (STDP) is one possible cellular model for the induction of local synaptic modifications, which could account for experience-driven changes of the connectivity in neuronal networks (Song and Abbott, 2001; Senn, 2002).

For STDP pairing protocols as well as for other plasticity induction protocols, like presynaptic tetanic stimulation (Lynch et al., 1983) or repetitive low-frequency synaptic stimulation (Mulkey and Malenka, 1992), the elevation of postsynaptic [Ca<sup>2+</sup>]<sub>i</sub> is essential. [Ca<sup>2+</sup>]<sub>i</sub> probably acts as one second messenger on downstream metabolic cascades that are responsible for the eventual modification of synaptic efficacy (Lisman, 1989). A long-standing hypothesis suggests that the peak amplitude of postsynaptic [Ca<sup>2+</sup>]<sub>i</sub> elevation determines the direction and the magnitude of such modifications (Bear et al., 1987; Artola and Singer, 1993; Hansel et al., 1997). Several studies have supported this hypothesis by measuring dendritic [Ca<sup>2+</sup>]<sub>i</sub> (Cormier et al., 2001; Ismailov et al., 2004; Gall et al., 2005). Nevertheless, the problem of how a single variable, the increase in global [Ca<sup>2+</sup>]<sub>i</sub>, can control the differential induction of changes in synaptic efficacy has not been explained conclusively. According to this "Ca<sup>2+</sup> control hypothesis" the level of [Ca<sup>2+</sup>]<sub>i</sub> acts differentially on downstream protein cascades, activating either kinases or phosphatases, which, respectively, phosphorylate or dephosphorylate postsynaptic AMPA receptors (Lisman, 1989; Lee et al., 2000). Modified versions of this hypothesis include veto mechanisms for LTD at moderate Ca<sup>2+</sup> levels (Rubin et al., 2005), differential microdomain Ca<sup>2+</sup> signaling (Franks and Sejnowski, 2002), or the proposition of a second coincidence detector (Karmarkar and Buonomano, 2002) to account for the differential induction of LTD and LTP.

Here we characterize the Ca<sup>2+</sup> signals in spines of basal dendrites of layer 2/3 (L2/3) pyramidal neurons in the somatosensory cortex during LTP- and LTD-inducing protocols. We found that postsynaptic elevation of [Ca<sup>2+</sup>]<sub>i</sub> in spines is necessary, but by itself it is not sufficient, to account for both potentiation and

Received April 25, 2006; revised Sept. 13, 2006; accepted Sept. 14, 2006.
We thank Matthew Larkum and Hans-Rudolf Lüscher for their support and helpful comments on this manuscript and Marlies Kaiser and Karl Schmidt for excellent technical support.
Correspondence should be addressed to Thomas Nevian at his present address: Institute for Physiology, Bern University, Bühlplatz 5, CH-3012 Bern, Switzerland. E-mail: nevian@pyl.unibe.ch.
DOI:10.1523/JNEUROSCI.1749-06.2006
Copyright © 2006 Society for Neuroscience 0270-6474/06/2611001-13$15.00/0

11002 • J. Neurosci., October 25, 2006 • 26(43):11001–11013
Nevian and Sakmann • Spine Ca<sup>2+</sup> and Spike-Timing-Dependent Plasticity

depression of synaptic strength. Induction of LTD requires, in addition, the activation of a metabotropic glutamate receptor-dependent (mGluR-dependent) signaling cascade resulting in the synthesis of endocannabinoids, which then act as a retrograde messenger. We conclude that the inductions of LTP and LTD both depend on a rise in [Ca<sup>2+</sup>]<sub>i</sub> but represent two separate processes triggered by two Ca<sup>2+</sup> sensors. Thus [Ca<sup>2+</sup>]<sub>i</sub> elevation and the coactivation of mGluRs are necessary to account for the bidirectional changes in synaptic efficacy observed with STDP protocols.

## Materials and Methods

*Electrophysiology.* Sagittal brain slices were prepared from postnatal day 13–15 Wistar rats. All experimental procedures were in accordance with the animal welfare guidelines of the Max Planck Society. Experiments were performed at physiological temperatures (32–36°C) in extracellular solution containing the following (in m<small>M</small>): 125 NaCl, 25 NaHCO<sub>3</sub>, 2.5 KCl, 1.25 NaH<sub>2</sub>PO<sub>4</sub>, 1 MgCl<sub>2</sub>, 25 glucose, 2 CaCl<sub>2</sub>, and 0.01 glycine, bubbled with 95% O<sub>2</sub>/5% CO<sub>2</sub>. Inhibitory inputs were blocked by the bath application of 10 $\mu$M bicuculline. Pyramidal neurons in L2/3 of the somatosensory cortex were visualized by using infrared (IR) gradient contrast video microscopy. Patch pipettes (5–7 M$\Omega$) were filled with a low chloride intracellular solution containing the following (in m<small>M</small>): 130 K-gluconate, 10 K-HEPES, 10 Na-phosphocreatine, 4 Mg-ATP, 0.3 Na-GTP, 4 NaCl, and 10 Na-gluconate. Whole-cell current-clamp recordings were made with an AxoClamp-2B (Molecular Devices, Union City, CA) patch-clamp amplifier. Voltage signals were filtered at 3 kHz and digitized at 10 kHz with an ITC-16 (InstruTech, Port Washington, NY). Access resistance (10–20 M$\Omega$) usually did not change during the course of the experiment. Input resistance was monitored constantly by a brief hyperpolarizing current pulse. Experiments were excluded if input resistance or membrane potential changed significantly over the time course of the experiment.

An extracellular stimulation pipette filled with extracellular solution was placed close to the basal dendrites of the L2/3 pyramidal neuron (50–150 $\mu$m from soma). Stimulation strength (3–10 $\mu$A; 100 $\mu$s duration) was adjusted to evoke baseline single component EPSPs with amplitudes between 1 and 3 mV. Baseline EPSPs were recorded for 10 min at 0.1 Hz stimulation. Then EPSPs were paired with one to three APs at different frequencies and variable onset times, $\Delta t'$. Pairings were repeated 60 times at 0.1 Hz stimulation. The time interval $\Delta t'$ was defined as the time between the onset of the AP burst (first AP in the burst) and the onset of the compound EPSP. The time interval $\Delta t$ was defined as the time between the AP closest in time to the EPSP and the onset of the EPSP. The change in EPSP amplitude was evaluated 20–40 min after the end of the pairing period and normalized to the baseline EPSP amplitude. Data are presented as the mean $\pm$ SEM. Paired Student's $t$ tests were applied as statistical tests if not indicated otherwise, and statistical significance was asserted for $p < 0.05$.

*Calcium imaging.* For two-photon excitation fluorescence microscopy an IR femtosecond-pulsed titanium sapphire laser (Mira 900, pumped by a 5 W Verdi; Coherent, Santa Clara, CA) was coupled directly to a confocal-scanning unit (LCS SP2RS, Leica Microsystems, Mannheim, Germany) attached to an upright microscope (DMLFS, Leica) equipped with a 63$\times$ objective (HCX APO W63$\times$ UVI; numerical aperture 0.9; Leica) (Rathenberg et al., 2003). Cells were filled with a combination of the Ca<sup>2+</sup>-insensitive dye Alexa 594 (50 $\mu$M; Invitrogen, Carlsbad, CA) and the Ca<sup>2+</sup> indicator Oregon Green BAPTA-6F (500 $\mu$M; Invitrogen) added to the low chloride intracellular solution. Thus all Ca<sup>2+</sup> imaging experiments represent a different set of experiments from the separately performed experiments addressing the change in synaptic strength. Dyes were excited at $\lambda = 820$ nm. Excitation IR laser light and fluorescence emission light were separated at 670 nm (excitation filter 670DCXXR, AHF Analysentechnik, Tübingen, Germany). The emission spectra were separated by a dichroic mirror at 560 nm (beam splitter 560DCXR, AHF) and corresponding bandpass (HQ525/50, HQ630/60, AHF) and IR-block filters (E700SP, BG39; AHF) and were detected by using non-descanned detection. In addition, the forward-scattered IR laser light was filtered spatially and imaged onto a photomultiplier tube to generate an IR-scanning gradient contrast image of the unstained brain slice (Wimmer et al., 2004). Line scans through a spine that responded with a [Ca<sup>2+</sup>]<sub>i</sub> transient to subthreshold synaptic stimulation were made to measure the [Ca<sup>2+</sup>]<sub>i</sub> transients evoked by the different pairing protocols. Relative changes in fluorescence are expressed as the following: $\Delta G/R = (G(t) - G_0)/R$, where $G(t)$ is the fluorescence signal integrated from the region of interest covering the spine in green channels, $G_0$ is the basal fluorescence averaged for 50 ms before electrical stimulation, and $R$ is the averaged fluorescence of the red channels (Oertner et al., 2002). Single exponential fits to the decay phase of the [Ca<sup>2+</sup>]<sub>i</sub> transients yielded the peak amplitude $(\Delta G/R)_{\text{max}}$. The nonlinearity factor of Ca<sup>2+</sup> signals was evaluated by normalizing the peak [Ca<sup>2+</sup>]<sub>i</sub> amplitude for a given pairing protocol to the expected linear sum of EPSP-evoked and AP-evoked [Ca<sup>2+</sup>]<sub>i</sub> transients (Nevian and Sakmann, 2004).

*Chemical compounds.* L-type voltage-dependent Ca<sup>2+</sup> channels (L-VDCCs) were blocked by the bath application of nimodipine (10 $\mu$M) and T-type VDCCs (T-VDCCs) by NiCl<sub>2</sub> (50 $\mu$M). NMDA receptors were blocked by the bath application of <small>D</small>-5-amino-phosphonopentanoic acid (<small>D</small>-APV) (50 $\mu$M) or intracellular application of (+)-5-methyl-10,11-dihydro-<sup>5</sup>$H$-dibenzo[$a,d$]cyclohepten-5,10-imine maleate (MK-801) (1 m<small>M</small>). mGluRs were blocked by the broadband mGluR antagonist ($S$)-$\alpha$-methyl-4-carboxyphenylglycine (MCPG) (500 $\mu$M), phospholipase C (PLC) was blocked by 1-(6-[(17$\beta$-methoxyestra-1,3,5 [10]-trien-17-yl) amino] hexyl)-<sup>1</sup>$H$-pyrrole-2,5-dione (U73122) (5 $\mu$M), and cannabinoid type 1 (CB1) receptors were blocked by $N$-(piperidin-1-yl)-5-(4-iodophenyl)-1-(2,4-dichlorophenyl)-4-methyl-<sup>1</sup>$H$-pyrazole-3-carboxamide (AM251) (2 $\mu$M). Inositol 1,4,5-triphosphate-mediated (IP<sub>3</sub>-mediated) Ca<sup>2+</sup> release from internal stores was blocked by intracellular application of heparin (400 U/ml). Cytosolic elevations of Ca<sup>2+</sup> were buffered by dialyzing the cells with EGTA and BAPTA at various concentrations (0.1–2 m<small>M</small>).

## Results

### Burst-timing-dependent potentiation and depression in L2/3 pyramids
Whole-cell *in vivo* recordings from morphologically identified barrel-related L2/3 pyramidal neurons indicate large compound EPSPs in response to tactile stimuli but sparse AP responses (Zhu and Connors, 1999; Brecht et al., 2003), with rare burst activity (Svoboda et al., 1997, 1999). The induction of long-lasting changes in synaptic efficacy has been reported to depend strongly on the timing between presynaptic and postsynaptic APs (Markram et al., 1997; Bi and Poo, 1998). Therefore, we first investigated in L2/3 pyramidal neurons whether such a timing dependence also existed for pairing a single compound EPSP that was elicited by focal extracellular stimulation with a burst of three APs occurring at 50 Hz in the postsynaptic cell. We found that the degree of both potentiation and depression depended on the time interval $\Delta t'$ between onset of the AP burst and onset of the EPSP (Figs. 1, 2). When the AP burst preceded the EPSP by $\Delta t' = -90$ ms or followed the EPSP by $\Delta t' = +50$ ms, the pairing had no effect ($\Delta t' = -90$ ms: $1.00 \pm 0.09, p > 0.5, n = 6$; $\Delta t' = +50$ ms: $0.92 \pm 0.11, p > 0.5, n = 4$). Pairing at $\Delta t' = -50$ ms resulted in a long-lasting depression of the EPSP amplitude, with a reduction of the EPSP amplitude by a factor of $0.68 \pm 0.05$ ($p < 0.01; n = 10$) as compared with control (Figs. 1*C*, 2*B*). Two APs before and one AP after the EPSP ($\Delta t' = -30$ ms) did not change the EPSP amplitude ($0.98 \pm 0.12; p > 0.5; n = 8$). In contrast, one AP before and two APs after the EPSP ($\Delta t' = -10$ ms) resulted in moderate potentiation ($1.42 \pm 0.19; p < 0.05; n = 12$). Strong potentiation by a factor of $2.01 \pm 0.22$ ($p < 0.01; n = 11$) was found if the three APs followed the EPSP by $\Delta t' = +10$ ms (Figs. 1*B*, 2).

Second, we tested the requirement for postsynaptic bursting to induce a change in synaptic strength (Fig. 3*A*, *B*). The time

Nevian and Sakmann • Spine Ca<sup>2+</sup> and Spike-Timing-Dependent Plasticity
J. Neurosci., October 25, 2006 • 26(43):11001–11013 • 11003

**A**
[The image shows a fluorescence micrograph of a Layer 2/3 pyramidal neuron with a recording pipette attached to the soma and an extracellular stimulation pipette placed near the basal dendrites. A scale bar indicates 20 μm.]

**B**
+10 ms
$\Delta t'$
pre [___]
post [___]
$\Delta t$
EPSP 3 APs
[Diagram showing the timing of an EPSP followed by a burst of 3 APs]
$\Delta t$ +10 ms

<table>
  <tbody>
    <tr>
        <td>time (min)</td>
        <td>A<sub>EPSP</sub> (mV)</td>
        <td>V<sub>m</sub> (mV)</td>
        <td>R<sub>i</sub> (MΩ)</td>
    </tr>
    <tr>
        <td>-10</td>
        <td>2</td>
        <td>-70</td>
        <td>200</td>
    </tr>
    <tr>
        <td>0</td>
        <td>2</td>
        <td>-70</td>
        <td>200</td>
    </tr>
    <tr>
        <td>10</td>
        <td>2.5</td>
        <td>-70</td>
        <td>200</td>
    </tr>
    <tr>
        <td>20</td>
        <td>3.2</td>
        <td>-70</td>
        <td>200</td>
    </tr>
    <tr>
        <td>30</td>
        <td>3.5</td>
        <td>-70</td>
        <td>200</td>
    </tr>
    <tr>
        <td>40</td>
        <td>3.5</td>
        <td>-70</td>
        <td>200</td>
    </tr>
    <tr>
        <td>50</td>
        <td>3.5</td>
        <td>-70</td>
        <td>200</td>
    </tr>
    <tr>
        <td>60</td>
        <td>3.5</td>
        <td>-70</td>
        <td>200</td>
    </tr>
  </tbody>
</table>
[Traces showing EPSP before and after pairing, with an increase in amplitude from 2 to 3.5 mV. Scale bars: 2 mV, 50 ms.]

**C**
-50 ms
$\Delta t'$
pre [___]
post [___]
$\Delta t$
3 APs EPSP
[Diagram showing the timing of a burst of 3 APs followed by an EPSP]
$\Delta t$ -10 ms

<table>
  <tbody>
    <tr>
        <td>time (min)</td>
        <td>A<sub>EPSP</sub> (mV)</td>
        <td>V<sub>m</sub> (mV)</td>
        <td>R<sub>i</sub> (MΩ)</td>
    </tr>
    <tr>
        <td>-10</td>
        <td>2</td>
        <td>-70</td>
        <td>200</td>
    </tr>
    <tr>
        <td>0</td>
        <td>2</td>
        <td>-70</td>
        <td>200</td>
    </tr>
    <tr>
        <td>10</td>
        <td>1.2</td>
        <td>-70</td>
        <td>200</td>
    </tr>
    <tr>
        <td>20</td>
        <td>1</td>
        <td>-70</td>
        <td>200</td>
    </tr>
    <tr>
        <td>30</td>
        <td>1</td>
        <td>-70</td>
        <td>200</td>
    </tr>
    <tr>
        <td>40</td>
        <td>1</td>
        <td>-70</td>
        <td>200</td>
    </tr>
    <tr>
        <td>50</td>
        <td>1</td>
        <td>-70</td>
        <td>200</td>
    </tr>
    <tr>
        <td>60</td>
        <td>1</td>
        <td>-70</td>
        <td>200</td>
    </tr>
  </tbody>
</table>
[Traces showing EPSP before and after pairing, with a decrease in amplitude from 2 to 1 mV. Scale bars: 2 mV, 50 ms.]

**Figure 1.** Timing-dependent induction of LTP and LTD by pairing an EPSP and a short burst of APs. **A**, Illustration of the recording configuration for the synaptic plasticity experiments. The extracellular stimulation pipette was placed close to the basal dendrites, ~50 μm from the soma of the L2/3 pyramidal neuron. **B**, The induction protocol for LTP is depicted on the left. An EPSP evoked by extracellular stimulation (pre) was paired with a short burst of three APs at 50 Hz elicited by current injections into the postsynaptic cell (post). The first AP followed the onset of the EPSP by $\Delta t = 10$ ms. The time interval $\Delta t$ is defined as the time between the EPSP and the AP closest in time to the EPSP. For bursts that follow the EPSP, $\Delta t$ is equivalent to $\Delta t'$, which is defined as the time interval between the EPSP and the first AP in the burst. Pairing was repeated 60 times, every 10 s. To the right, the EPSP amplitude, membrane potential, and input resistance are plotted over time during the experiment. The dashed line indicates the average EPSP amplitude before the pairing. The pairing protocol depicted to the left resulted in potentiation of the EPSP amplitude. The EPSPs averaged over the times indicated by the red bars are shown on the bottom. The average EPSP amplitude increased from 2 to 3.5 mV at 20–40 min after the pairing period. **C**, The protocol for the induction of LTD is depicted to the left. A short

**A**
[Six graphs showing normalized EPSP amplitude over time for different $\Delta t'$ intervals: -90 ms, -50 ms, -30 ms, -10 ms, +10 ms, and +50 ms. Insets show the relative timing of EPSP and the 3-AP burst.]

**B**
<table>
  <tbody>
    <tr>
        <td>Δt' (ms)</td>
        <td>change in EPSP amplitude</td>
    </tr>
    <tr>
        <td>-90</td>
        <td>1.0</td>
    </tr>
    <tr>
        <td>-50</td>
        <td>0.67</td>
    </tr>
    <tr>
        <td>-30</td>
        <td>1.0</td>
    </tr>
    <tr>
        <td>-10</td>
        <td>1.5</td>
    </tr>
    <tr>
        <td>10</td>
        <td>2.08</td>
    </tr>
    <tr>
        <td>50</td>
        <td>1.0</td>
    </tr>
  </tbody>
</table>
[The data is plotted as a burst-timing-dependent plasticity curve showing depression at -50 ms and potentiation at -10 ms and +10 ms.]

**Figure 2.** Burst-timing-dependent plasticity curve. **A**, Pooled and normalized EPSP amplitudes for different time intervals $\Delta t'$ between the onset of the EPSP and the first AP in the burst. The burst consisted of three APs (50 Hz). Insets show schematically the timing between the EPSP and the three APs. **B**, Summary of the change in the EPSP amplitude for different EPSP and AP burst-timing intervals $\Delta t'$. If the AP burst preceded the EPSP by more than -90 ms or followed by more than +150 ms, no change in EPSP amplitude was found ($p > 0.5; n = 4-6$). A timing interval of $\Delta t' = -50$ ms resulted in depression of the EPSP amplitude by a factor of $0.67 \pm 0.06$ (mean $\pm$ SEM; **$p < 0.01; n = 9$), whereas a timing interval of $\Delta t' = +10$ ms resulted in potentiation by a factor of $2.08 \pm 0.25$ (**$p < 0.01; n = 11$). If the EPSP was evoked within the AP burst, either no change ($\Delta t' = -30$ ms; $p > 0.5; n = 8$) or potentiation by a factor of $1.5 \pm 0.3$ ($\Delta t' = -10$ ms; *$p < 0.05; n = 6$) was observed. The horizontal dashed line represents no change in EPSP amplitude, and the vertical dashed line represents the onset of the EPSP.

interval $\Delta t$, defined as the time interval between the onset of the EPSP and the AP closest in time to the EPSP, was kept constant at $\Delta t = -10$ ms or $\Delta t = +10$ ms. Pairing an EPSP and one AP at $\Delta t = +10$ ms did not change the EPSP amplitude ($1.04 \pm 0.08; p > 0.5; n = 10$), whereas the reverse order resulted in a decrease in EPSP amplitude ($0.80 \pm 0.07; p < 0.05; n = 5$). A "minimal

$\leftarrow$
burst of three APs at 50 Hz was paired with a following EPSP. Here the first AP in the burst preceded the onset of the EPSP by $\Delta t' = -50$ ms, which is equivalent to the definition that the last AP in the burst preceded the onset of the EPSP by $\Delta t = -10$ ms. To the right, the EPSP amplitude, membrane potential, and input resistance are plotted over time for the experiment in which the burst preceded the EPSP. This pairing protocol resulted in depression of the EPSP amplitude. The voltage recordings averaged over the times indicated by the red bars are shown on the bottom. The average EPSP amplitude decreased from 2 to 1 mV at 20–40 min after the pairing period.

11004 • J. Neurosci., October 25, 2006 • 26(43):11001–11013
Nevian and Sakmann • Spine Ca<sup>2+</sup> and Spike-Timing-Dependent Plasticity

burst” of two APs (50 Hz) at $\Delta t = +10$ ms resulted in significant potentiation ($1.95 \pm 0.31$; $p < 0.05$; $n = 9$). The reverse order (two APs and EPSP; $\Delta t = -10$ ms) resulted in depression ($0.72 \pm 0.12$; $p < 0.05$; $n = 5$).

Synaptic stimulation without pairing with postsynaptic APs (EPSPs only) did not result in a change in EPSP amplitude ($0.97 \pm 0.08$; $p > 0.5$; $n = 5$). Similarly, a burst of three APs at 50 Hz without pairing with an EPSP (three APs only) had no effect ($1.05 \pm 0.2$; $p > 0.9$; $n = 6$).

Last, we varied the AP frequency in a burst of three APs between 20 and 100 Hz (Fig. 3C,D). A burst of three APs at 20 Hz that followed the EPSP by $\Delta t = +10$ ms had no effect on the EPSP amplitude after pairing ($1.09 \pm 0.27$; $p > 0.5$; $n = 5$). A burst at 20 Hz preceding the EPSP by $\Delta t = -10$ ms resulted in depression ($0.72 \pm 0.14$; $p < 0.05$; $n = 7$). A burst of three APs at 100 Hz that followed the EPSP by $\Delta t = +10$ ms induced strong potentiation ($2.29 \pm 0.48$; $p < 0.01$; $n = 7$), whereas a burst of three APs at 100 Hz preceding the EPSP by $\Delta t = -10$ ms resulted in depression ($0.52 \pm 0.12$; $p < 0.05$; $n = 4$).

From these results we conclude that LTP and LTD depend on the timing between the AP burst and onset of the compound EPSP. The induction of LTD is less sensitive to the properties of the burst, whereas LTP requires at least a burst of two APs at frequencies $>20$ Hz. These results indicate that small variations in postsynaptic AP firing with respect to a synaptically evoked EPSP can result in different directions of the change in synaptic strength.

### LTP and LTD are equally sensitive to fast and slow Ca<sup>2+</sup> buffers
Changes in synaptic strength depend on a transient [Ca<sup>2+</sup>]<sub>i</sub> elevation in the dendrites of the postsynaptic neuron. It is, however, unclear whether a putative "Ca<sup>2+</sup> sensor" that may trigger the induction of changes in synaptic strength is sensitive to global, volume-averaged [Ca<sup>2+</sup>]<sub>i</sub> or whether it is sensitive to local [Ca<sup>2+</sup>]<sub>i</sub> transients in Ca<sup>2+</sup> microdomains around sources of Ca<sup>2+</sup> influx (Franks and Sejnowski, 2002). One way to differentiate between global and local [Ca<sup>2+</sup>]<sub>i</sub> signaling is to load the postsynaptic cell with different Ca<sup>2+</sup> buffers, such as EGTA or BAPTA. These have a similar affinity ($K_D \approx 200$ nM) but faster (BAPTA) or slower (EGTA) equilibration kinetics.

Loading cells with EGTA or BAPTA at 2 mM was sufficient to block LTP induction (relative magnitude of LTP with EGTA: $-0.05 \pm 0.18$, $p > 0.7$, $n = 6$; BAPTA: $-0.16 \pm 0.14$, $p > 0.7$, $n = 5$) (Fig. 4A). Lower concentrations of EGTA and BAPTA blocked LTP induction less (0.5 mM EGTA: $0.64 \pm 0.20$, $p < 0.05$, $n = 6$; 0.5 mM BAPTA: $0.63 \pm 0.3$, $p < 0.05$, $n = 6$). A sigmoid fit to the concentration dependence indicated a half-effective concentration for blocking LTP of 0.6 mM for EGTA and 0.5 mM for BAPTA (Fig. 4B).

[The image contains four panels of graphs labeled A, B, C, and D.]

**Figure 3.** LTP and LTD dependence on the number of APs and AP frequency in the burst. **A**, Pooled and normalized EPSP amplitudes for different numbers of APs in the burst and a time interval between the AP closest to the EPSP and the onset of the EPSP of $\Delta t = -10$ ms and $\Delta t = +10$ ms. Insets depict the pairing protocol. Pairing an EPSP with a single AP at $\Delta t = +10$ ms had no effect on EPSP amplitude. The sequence of AP and EPSP at $\Delta t = -10$ ms resulted in depression. An EPSP paired with a burst of two APs (50 Hz) at $\Delta t = +10$ ms resulted in potentiation. Two APs (50 Hz) preceding the EPSP resulted in depression of the EPSP amplitude. **B**, LTP induction at $\Delta t = +10$ ms (open circles) depended on the number of APs that followed the EPSP, whereas LTD induction (open boxes) was correlated only weakly with the number of APs preceding the EPSP by $\Delta t = -10$ ms. The open diamond to the left corresponds to EPSP-only stimulation during the pairing period. Significant differences, *$p < 0.05$ and **$p < 0.01$. **C**, Pooled and normalized EPSP amplitudes for different AP frequencies in the burst and a time interval between the AP closest to the EPSP and the onset of the EPSP of $\Delta t = -10$ ms and $\Delta t = +10$ ms. Pairing an EPSP with a burst of three APs at 20 Hz did not result in significant potentiation. A burst of three APs at 20 Hz preceding an EPSP resulted in depression. Pairing an EPSP with a burst of three APs at 100 Hz that followed at $\Delta t = +10$ ms resulted in strong potentiation. A burst of three APs at 100 Hz preceding the EPSP resulted in depression. **D**, Summary plot of the change in EPSP amplitude versus the AP frequency in a burst of three APs. The induction of LTP is frequency-dependent (open circles), whereas LTD induction did not depend on the burst frequency. The dashed lines indicate no change in EPSP amplitude.

LTD induction was blocked completely by 1 mM EGTA (relative magnitude of LTD, $-0.2 \pm 0.25$; $p > 0.3$; $n = 3$) (Fig. 4C), whereas 0.25 mM EGTA blocked LTD induction correspondingly less ($-0.86 \pm 0.35$; $p < 0.05$; $n = 6$). BAPTA (1 mM) also blocked LTD induction completely ($-0.09 \pm 0.51$; $p > 0.5$; $n = 7$), and 0.25 mM BAPTA blocked LTD induction less ($0.92 \pm 0.2$; $p < 0.05$; $n = 6$). The half-effective concentration for the block of LTD was 0.39 mM by EGTA and 0.36 mM by BAPTA (Fig. 4D).

We conclude that postsynaptic elevation of [Ca<sup>2+</sup>]<sub>i</sub> is necessary for the induction of both LTP and LTD. The comparable concentration dependence of BAPTA and EGTA in blocking the induction or expression of potentiation or depression suggests that the putative Ca<sup>2+</sup> sensor or sensors that trigger long-lasting changes in synaptic efficacy respond to the volume-averaged increase in [Ca<sup>2+</sup>]<sub>i</sub>. Either they are separated by a relatively large distance from the site of Ca<sup>2+</sup> influx (Neher, 1998; Augustine et al., 2003) or they have slower Ca<sup>2+</sup>-binding kinetics than EGTA.

Nevian and Sakmann • Spine Ca<sup>2+</sup> and Spike-Timing-Dependent Plasticity J. Neurosci., October 25, 2006 • 26(43):11001–11013 • 11005

[The image contains four graphs labeled A, B, C, and D, illustrating the effects of Ca<sup>2+</sup> buffers on LTP and LTD.]

**Figure 4.** LTP and LTD are equally sensitive to fast and slow Ca<sup>2+</sup> buffers. **A**, Buffering of postsynaptic Ca<sup>2+</sup> with 2 mM EGTA (open circles) or 2 mM BAPTA (filled circles) blocked the induction of LTP by the sequence of an EPSP and three APs (50 Hz) at $\Delta t = +10$ ms. A lower concentration of EGTA or BAPTA of 0.5 mM blocked 40% of LTP. **B**, Titration curve for different concentrations of EGTA (open circles) and BAPTA (filled circles) and their effects on the induction of LTP. Data were normalized to the change in EPSP amplitude for no buffer added. The slow Ca<sup>2+</sup> buffer EGTA blocked LTP with the same concentration dependence as the fast Ca<sup>2+</sup> buffer BAPTA. A sigmoidal fit yielded a half-concentration of 0.6 and 0.5 mM for EGTA (gray line) and BAPTA (black line), respectively. **C**, A concentration of 1 mM EGTA (open squares) or BAPTA (filled squares) was sufficient to blocked the induction of LTD by the sequence of three APs (50 Hz) and an EPSP at $\Delta t = -10$ ms. A lower concentration of EGTA or BAPTA of 0.25 mM had no significant effect on the induction of LTD. **D**, Titration curve for different concentrations of EGTA (open boxes) and BAPTA (filled boxes) and their effect on the induction of LTD. Data were normalized to the change in EPSP amplitude for no buffer added. A sigmoidal fit yielded a half-concentration of 0.39 and 0.36 mM for EGTA (gray line) and BAPTA (black line), respectively. The dashed lines indicate no change in EPSP amplitude.

### Ca<sup>2+</sup> transients in spines evoked during STDP induction protocols
The results described above suggested that the induction of synaptic potentiation and depression depended on a rise in postsynaptic $[Ca^{2+}]_i$, presumably in dendritic spines. We used two-photon excitation fluorescence microscopy (Denk et al., 1990) to measure the $[Ca^{2+}]_i$ transients in single spines evoked during the stimulation patterns described above (Fig. 5A–C).

First, varying the time interval $\Delta t'$ between an EPSP and a burst of three APs (50 Hz) allowed us to map the peak $[Ca^{2+}]_i$ amplitude and the "nonlinearity" factor (peak $[Ca^{2+}]_i$ amplitude, normalized to the expected linear sum of EPSP-evoked and AP-evoked $[Ca^{2+}]_i$ transients) (Nevian and Sakmann, 2004) for the time intervals relevant for the induction of changes in synaptic strength (Fig. 5D). The plot of the peak $[Ca^{2+}]_i$ amplitude for different time intervals $\Delta t'$ revealed that pairing an EPSP with three APs at $\Delta t' = +10$ ms resulted in a significantly larger $[Ca^{2+}]_i$ transient than pairing an EPSP with three APs at the other time intervals that were tested ($p < 0.01$; $n = 4–39$; ANOVA; Newman–Keuls). These other pairing protocols evoked $[Ca^{2+}]_i$ transients with similar peak $[Ca^{2+}]_i$ amplitudes ($p > 0.2$; $n = 4–34$; ANOVA; Newman–Keuls). The summation of the $Ca^{2+}$ signals showed linear or supralinear summation, depending on the timing interval $\Delta t'$, similar to other types of cortical cells (Köster and Sakmann, 1998; Nevian and Sakmann, 2004). Pairing an EPSP with three APs at $\Delta t' = -10, +10,$ and $+50$ ms resulted in supralinear summation, with the largest degree of supralinearity for $\Delta t' = +10$ ms (nonlinearity factor of $1.8 \pm 0.1$; $p < 0.01$; $n = 39$). Pairing at $\Delta t' = -30, -50,$ and $-90$ ms did not deviate significantly from linear summation of the $Ca^{2+}$ signals ($p > 0.1$; $n = 4–21$).

Second, peak amplitude and nonlinearity of summation were analyzed as a function of the number of APs in a burst for a frequency of 50 Hz and a given time interval ($\Delta t = -10$ ms and $\Delta t = +10$ ms) between EPSP and AP burst. Peak amplitudes increased linearly with the number of APs for all of the stimulation protocols that were tested (linear regression lines, $r^2 > 0.98$) (Fig. 5E). The summation of the $Ca^{2+}$ signals was approximately twofold larger than the linear sum for stimulation protocols in which the APs followed the EPSP by $\Delta t = +10$ ms, independent of the number of APs in the burst ($p < 0.01$; $n = 7–39$). The summation of the $Ca^{2+}$ signals for stimulation protocols in which the APs preceded the EPSP was independent of the number of APs in the burst and not different from linear summation ($p > 0.5$; $n = 5–21$).

Third, the peak $[Ca^{2+}]_i$ amplitude increased linearly as a function of AP frequency in a burst of three APs for all of the stimulation protocols that were tested (Fig. 5F). The summation of the $Ca^{2+}$ signals was larger than the expected linear sum for stimulation protocols in which the burst followed the EPSP by $\Delta t = +10$ ms with burst frequencies $>20$ Hz. The summation of the $Ca^{2+}$ signals for stimulation protocols in which the burst preceded the EPSP was independent of the frequency of APs in the burst and not different from the expected linear sum ($p > 0.5$; $n = 3–21$).

### Dependence of LTP and LTD induction on Ca<sup>2+</sup> influx via NMDA receptors and VDCCs
The asymmetry in the peak $[Ca^{2+}]_i$ amplitude timing curve (Fig. 5D) presumably is caused by the coincident activation of the NMDA receptor (NMDAR) channel after binding of glutamate and relief of the $Mg^{2+}$ block caused by membrane depolarization attributable to backpropagating APs (Yuste and Denk, 1995; Köster and Sakmann, 1998; Nevian and Sakmann, 2004). We used the protocols consisting of an EPSP and three APs (50 Hz) at $\Delta t = -10$ ms and $\Delta t = +10$ ms to test the requirement of NMDAR activation for the timing-dependent difference in the $Ca^{2+}$ signals and for the induction of potentiation and depression (Fig. 6). Bath application of the NMDAR antagonist D-APV (50 $\mu$M) significantly reduced the peak $[Ca^{2+}]_i$ amplitude for both time intervals (reduction to $33 \pm 9\%$ for $\Delta t = +10$ ms, $p < 0.05$; $n = 4$ and reduction to $45 \pm 19\%$ for $\Delta t = -10$ ms, $p < 0.05$; $n = 4$). The $[Ca^{2+}]_i$ transients were reduced to the level of AP-evoked $Ca^{2+}$ influx, abolishing the timing dependence (Fig. 6B, C). Bath application of D-APV abolished potentiation and depression ($\Delta t = +10$ ms: $1.08 \pm 0.11, p > 0.1, n = 3; \Delta t = -10$ ms: $0.98 \pm 0.11, p > 0.5, n = 6$), indicating that NMDAR activation was necessary for the induction of LTP and LTD (Fig. 6D, E).

The induction of spike-timing-dependent LTD, but not LTP, might be mediated by the activation of presynaptic NMDARs (Sjostrom et al., 2003; Bender et al., 2006). Postsynaptic NMDARs can be blocked specifically by loading the postsynaptic cell with MK-801 (1 mM), an open channel blocker for NMDARs, through the patch pipette. The induction of LTP was abolished

11006 • J. Neurosci., October 25, 2006 • 26(43):11001–11013
Nevian and Sakmann • Spine Ca<sup>2+</sup> and Spike-Timing-Dependent Plasticity

($0.95 \pm 0.19$; $p > 0.5$; $n = 6$) with MK-801 present in the intracellular solution, similar to bath application of D-APV (Fig. 6$F,H$). In contrast, intracellular MK-801 had no effect on the induction of LTD ($0.60 \pm 0.11$; $p < 0.05$; $n = 7$) (Fig. 6$G,I$). Therefore, we conclude that the induction of LTP depends on a large postsynaptic Ca<sup>2+</sup> influx through NMDARs, whereas AP burst-pairing-induced LTD is independent of postsynaptic activation of NMDARs. In the latter case Ca<sup>2+</sup> influx through VDCCs is presumably sufficient to induce LTD (Fig. 6$B$).

L-VDCCs might control the induction of changes in synaptic strength (Magee and Johnston, 1997, 2005). Therefore, we tested the contribution of L-VDCCs on the [Ca<sup>2+</sup>]<sub>i</sub> transients evoked by the protocols consisting of an EPSP and three APs (50 Hz) at $\Delta t = -10$ ms and $\Delta t = +10$ ms (Fig. 7$A–C$). Bath application of the L-VDCC blocker nimodipine (10 $\mu$M) had no effect on the [Ca<sup>2+</sup>]<sub>i</sub> transient evoked by pairing an EPSP and three APs at $\Delta t = +10$ ms. For the sequence of three APs and an EPSP at $\Delta t = -10$ ms, nimodipine significantly reduced the peak [Ca<sup>2+</sup>]<sub>i</sub> amplitude to $69 \pm 1\%$ of control ($p < 0.05$; $n = 3$). Next we tested the effect of nimodipine on plasticity. Pairing an EPSP with three APs (50 Hz) at $\Delta t = +10$ ms had no effect on LTP induction in the presence of nimodipine (Fig. 7$D$). The potentiation of the EPSP amplitude was $1.92 \pm 0.31$ ($p < 0.01$; $n = 7$), similar to control. Nimodipine also had no effect on the induction of LTD ($0.70 \pm 0.10$; $p < 0.05$; $n = 8$), corresponding to the only minor reduction of the [Ca<sup>2+</sup>]<sub>i</sub> transient.

T-VDCCs have been suggested to be involved in spike-timing-dependent induction of LTD (Bender et al., 2006). Consistently, pairing one AP with an EPSP at $\Delta t = -10$ ms in the presence of Ni<sup>2+</sup> abolished the induction of LTD ($1.05 \pm 0.16$; $p > 0.5$; $n = 3$) (Fig. 7$F$). On the other hand, a burst of three APs (50 Hz) paired with an EPSP at $\Delta t = -10$ ms in the presence of Ni<sup>2+</sup> resulted in LTD ($0.70 \pm 0.03$; $p < 0.05$; $n = 4$) (Fig. 7$G,I$). Ni<sup>2+</sup> reduced the [Ca<sup>2+</sup>]<sub>i</sub> transient evoked in this case to $65 \pm 12\%$ ($n = 3$) of control (Fig. 7$H$). Bath application of nimodipine together with Ni<sup>2+</sup> resulted in a block of LTD by pairing three APs (50 Hz) with an EPSP at $\Delta t = -10$ ms ($1.00 \pm 0.05$; $p > 0.5$; $n = 3$) (Fig. 7$I$).

These experiments show that VDCCs are necessary for the induction of LTD. Single-spike LTD requires Ca<sup>2+</sup> influx through T-VDCCs; however, during pharmacological block of these channels, bursts of APs evoke a Ca<sup>2+</sup> influx via other subtypes, which is sufficiently large to induce LTD. Therefore, the requirement for a specific VDCC-subtype for the induction of LTD is masked by AP burst pairing.

**Figure 5.** [Ca<sup>2+</sup>]<sub>i</sub> transients in single spines evoked by pairing protocols. **A**, Two-photon fluorescence image of a L2/3 pyramidal neuron overlaid with the simultaneously acquired IR image. The region indicated by the dashed box is shown on an expanded scale in **B**. **B**, Fluorescence image of a spine responding with a [Ca<sup>2+</sup>]<sub>i</sub> transient to synaptic stimulation. The dashed line indicates the position of the line scan. An example line scan during a synaptically evoked EPSP paired with three APs (50 Hz) at $\Delta t = +10$ ms is shown in the bottom image. The green fluorescence of the Ca<sup>2+</sup>-sensitive indicator Oregon Green BAPTA-6F (500 $\mu$M) is overlaid with the red fluorescence from the Ca<sup>2+</sup>-insensitive dye Alexa 594 (50 $\mu$M). Stimulation begins 50 ms after the beginning of the recording. The bar at the top of the line scan indicates the region of interest over which fluorescence was averaged for each time point. **C**, Somatic voltage recordings (top traces) and the corresponding [Ca<sup>2+</sup>]<sub>i</sub> transients (bottom traces) for an EPSP, three APs, three APs (50 Hz) and an EPSP at $\Delta t' = -50$ ms, and an EPSP and three APs (50 Hz) at $\Delta t' = +10$ ms. The solid gray lines are exponential fits to the decay phase of the [Ca<sup>2+</sup>]<sub>i</sub> transient yielding the peak amplitude $(\Delta G/R)_{max}$. **D**, Peak amplitude (top graph) and nonlinearity factor (bottom graph) plotted for the time intervals $\Delta t'$. The dashed line indicates linear summation of the [Ca<sup>2+</sup>]<sub>i</sub> transients. Significant differences, *$p < 0.05$ and **$p < 0.01$. **E**, Peak amplitude (top graph) increases linearly with the number of APs in the burst for APs only (triangles), an EPSP and APs (50 Hz) at $\Delta t = +10$ ms (circles), and the sequence of APs (50 Hz) and an EPSP at $\Delta t = -10$ ms (boxes). The diamond to the left indicates the peak amplitude for the EPSP-evoked [Ca<sup>2+</sup>]<sub>i</sub> transient. The nonlinearity factor (bottom graph) does not depend on the number of APs in the burst, but it depends on the relative timing between the EPSP and APs. APs that follow the EPSP by $\Delta t = +10$ ms result in supralinear summation of the [Ca<sup>2+</sup>]<sub>i</sub> transients, whereas APs preceding the EPSP by $\Delta t = -10$ ms result in linear summation. **F**, Peak amplitude (top graph) increases linearly with the frequency in a burst of three APs for APs only (triangles), an EPSP and APs at $\Delta t = +10$ ms (circles), and the sequence of APs and an EPSP at $\Delta t = -10$ ms (boxes). The nonlinearity factor (bottom graph) indicates linear summation of the [Ca<sup>2+</sup>]<sub>i</sub> transients for APs preceding the EPSP, independent of AP frequency and supralinear summation of [Ca<sup>2+</sup>]<sub>i</sub> transients for APs that follow the EPSP for frequencies $>20$ Hz.

Nevian and Sakmann • Spine Ca<sup>2+</sup> and Spike-Timing-Dependent Plasticity
J. Neurosci., October 25, 2006 • 26(43):11001–11013 • 11007

![Figure 6](image_placeholder)
**Figure 6.** The contribution of NMDARs to Ca<sup>2+</sup> signaling and synaptic plasticity. **A**, Fluorescence image of a spine responding with a [Ca<sup>2+</sup>]<sub>i</sub> transient to synaptic stimulation. The dashed line indicates the position of the line scan. **B**, [Ca<sup>2+</sup>]<sub>i</sub> transients for the sequence of an EPSP and three APs (50 Hz) for $\Delta t = +10$ ms (top trace) and $\Delta t = -10$ ms (bottom trace) for control (gray trace) and after the bath application of the NMDAR blocker D-APV (50 $\mu$M; black trace). **C**, Average peak amplitudes for pairing an EPSP with three APs (50 Hz) at $\Delta t = -10$ ms and $\Delta t = +10$ ms for control (gray circles) and after the bath application of D-APV (black circles). D-APV significantly reduced the peak amplitude for both time intervals (*p < 0.05; n = 3). **D**, Normalized EPSP amplitude over time for the pairing protocol of an EPSP and three APs (50 Hz) at $\Delta t = +10$ ms in the presence of D-APV (bath application). LTP induction was abolished for $\Delta t = +10$ ms. Dashed lines indicate no change in EPSP amplitude. **E**, Bath application of D-APV also abolished the induction of LTD for $\Delta t = -10$ ms. **F**, Intracellular application of the open NMDAR channel blocker MK-801 (1 mM) blocked the induction of LTP. **G**, In contrast, intracellular application of MK-801 had no effect on the induction of LTD. **H**, Either bath application of D-APV or intracellular application of MK-801 blocked the induction of LTP by the pairing protocol of an EPSP and three APs (50 Hz) at $\Delta t = +10$ ms ($p > 0.3; n = 3–6$). **I**, Bath application of D-APV blocked the induction of LTD ($p > 0.5; n = 6$) by the pairing protocol of three APs (50 Hz) and an EPSP at $\Delta t = -10$ ms, whereas intracellular application of MK-801 had no effect on LTD (*p < 0.05; n = 7).

### Spine Ca<sup>2+</sup> transients cannot account for the direction of changes in synaptic efficacy
The analysis of the peak [Ca<sup>2+</sup>]<sub>i</sub> amplitude dependence on timing, number of APs, and frequency of APs suggested that pairing protocols, which induce either potentiation or depression, can raise [Ca<sup>2+</sup>]<sub>i</sub> to the same volume-averaged level. Therefore, we compared the average peak [Ca<sup>2+</sup>]<sub>i</sub> amplitude evoked by the different stimulation protocols to the average effect on long-term changes in synaptic efficacy. The resulting plot represented the average change in EPSP amplitude as a function of the average peak [Ca<sup>2+</sup>]<sub>i</sub> amplitude that a spine on a basal dendrite encounters during the pairing period (Fig. 8). The data showed no correlation between the peak [Ca<sup>2+</sup>]<sub>i</sub> amplitude and the direction of change in synaptic strength, because several stimulation protocols, which gave rise to similar peak [Ca<sup>2+</sup>]<sub>i</sub> amplitudes, induced either LTP or LTD (Fig. 8, shaded area). Thus although a postsynaptic increase in [Ca<sup>2+</sup>]<sub>i</sub> is necessary for the induction of changes in efficacy, the peak volume-averaged [Ca<sup>2+</sup>]<sub>i</sub> amplitude is not a unique determinant for the selective induction and expression of bidirectional synaptic plasticity. Additional factors, independent of the Ca<sup>2+</sup> signal, must contribute to the induction pathway. The observation that protocols in which the APs preceded the EPSP resulted in LTD whereas the reverse order resulted in LTP (Figs. 2G, 3 B, D) led to the clustering of the data into two groups according to the timing of the EPSP and APs. The clustering of stimulation protocols in which the APs either preceded or followed the EPSP indicated that the induction of LTD and LTP were two differential processes, both of which were controlled separately by the peak [Ca<sup>2+</sup>]<sub>i</sub> amplitude. Sigmoid fits to the two groups of data showed a high correlation between the induction of either LTD or LTP and the peak [Ca<sup>2+</sup>]<sub>i</sub> amplitude. The half-to-minimum value for the induction of LTD was $(\Delta G/R)_{max} = 0.07$, and the half-to-maximum value for the induction of LTP was approximately twofold larger [$(\Delta G/R)_{max} = 0.13$]. The stimulation patterns in which the EPSP was evoked within the burst of APs fell close to the fit of the Ca<sup>2+</sup> dependence for the LTP group, indicating that the induction of LTD can be blanked by APs that follow the AP–EPSP sequence.

We conclude that the induction of LTP and LTD is controlled differentially by the elevation of [Ca<sup>2+</sup>]<sub>i</sub> in the postsynaptic spine. The sequence of EPSP and postsynaptic APs seems to activate additional mechanisms, which determine the direction of the change in synaptic strength.

### LTD requires activation of mGluRs
In hippocampal as well as in cortical neurons the induction of some forms of LTD is sensitive to the activation of mGluRs (Otani and Connor, 1998; Anwyl, 1999, 2006; Egger et al., 1999; Cho et al., 2001). We tested whether LTD induced by pairing three APs (50 Hz) with an EPSP at $\Delta t = -10$ ms required the activation of mGluRs. Bath application of the broadband mGluR antagonist MCPG (500 $\mu$M) resulted in a complete block of LTD ($1.06 \pm 0.16; p > 0.5; n = 7$) (Fig. 9E). In contrast, MCPG had no effect on LTP induction by pairing an EPSP and three APs (50 Hz) at $\Delta t = +10$ ms ($1.79 \pm 0.33; p <$

11008 • J. Neurosci., October 25, 2006 • 26(43):11001–11013
Nevian and Sakmann • Spine Ca<sup>2+</sup> and Spike-Timing-Dependent Plasticity

$0.05; n = 9$ (Fig. 9D). The $[Ca^{2+}]_i$ transients evoked by pairing an EPSP with three APs (50 Hz) at $\Delta t = -10\text{ ms}$ and $\Delta t = +10\text{ ms}$ were not different from control after the bath application of MCPG ($p > 0.1; n = 4$) (Fig. 9A–C), indicating that mGluRs had no direct effect on the $[Ca^{2+}]_i$ transients (Brenowitz and Regehr, 2005) and that the block of LTD was attributable to inactivation of a G-protein-coupled signaling cascade (Piomelli, 2003).

The LTD-inducing protocol of three APs (100 Hz) paired with an EPSP at $\Delta t = -10\text{ ms}$ evoked $[Ca^{2+}]_i$ transients, which had a peak amplitude comparable to protocols that induced LTP. Bath application of MCPG also had no effect on the $[Ca^{2+}]_i$ transient in this case ($88 \pm 6\%$ of control; $p > 0.1; n = 4$). Surprisingly, pairing three APs (100 Hz) with an EPSP at $\Delta t = -10\text{ ms}$ in the presence of MCPG resulted in the induction of LTP ($1.53 \pm 0.24; p < 0.05; n = 11$) instead of inducing LTD, as under control conditions (Fig. 9F). We conclude that the activation of mGluRs is necessary for the induction of LTD. If the mGluR-dependent pathway is blocked, LTP can be induced, depending on the peak $[Ca^{2+}]_i$ amplitude (Fig. 9G).

### Signaling pathway for the induction of LTD
Next we investigated the downstream signaling cascade involved in the induction of LTD after mGluR activation. One product of this cascade is the synthesis of endocannabinoids, which can act as retrograde messengers, resulting in a decrease in synaptic efficacy (Sjöström et al., 2003, 2004). The endocannabinoid 2-arachidonoylglycerol (2-AG) is synthesized via a PLC and diacylglycerol (DAG) lipase-dependent pathway. Bath application of the CB1 receptor antagonist AM251 ($2\text{ }\mu\text{M}$) resulted in the block of LTD by the protocol consisting of three APs (50 Hz) and an EPSP at $\Delta t = -10\text{ ms}$ ($1.09 \pm 0.14; p > 0.5; n = 5$), suggesting the involvement of endocannabinoid signaling (Fig. 10E). The PLC inhibitor U73122 ($5\text{ }\mu\text{M}$) also blocked LTD ($1.11 \pm 0.07; p > 0.1; n = 4$) without any effect on the $[Ca^{2+}]_i$ transient evoked by the sequence of three APs (50 Hz) and an EPSP at $\Delta t = -10\text{ ms}$ ($99 \pm 11\%$ of control; $n = 4$) (Fig. 10A, B).

A side product of the synthesis of DAG by PLC is IP<sub>3</sub>, which could result in the release of $Ca^{2+}$ from internal stores. The observation that neither the block of mGluRs nor the block of PLC resulted in a significant change in the spine $[Ca^{2+}]_i$ transients suggested that release from internal stores did not contribute to the $[Ca^{2+}]_i$ transients that were required for the induction of LTD by a short burst of APs preceding synaptic activation. Another possibility might be that the concentration of IP<sub>3</sub> within the spine gradually increases during the 60 pairings, resulting in $Ca^{2+}$ release from internal stores only for later pairings. We recorded the $[Ca^{2+}]_i$ transients in a spine evoked by three APs (50 Hz) and an EPSP at $\Delta t = -10\text{ ms}$ during the pairing period. No correlation between the peak $[Ca^{2+}]_i$ amplitude and the number of pairings was found ($r^2 =$

[The page contains a composite figure labeled Figure 7 with panels A-I.]

**Figure 7.** The contribution of VDCCs to $Ca^{2+}$ signaling and synaptic plasticity. **A**, Fluorescence image of a spine responding with a $[Ca^{2+}]_i$ transient to synaptic stimulation. The dashed line indicates the position of the line scan. **B**, $[Ca^{2+}]_i$ transients for the sequence of an EPSP and three APs (50 Hz) for $\Delta t = +10\text{ ms}$ (top trace) and $\Delta t = -10\text{ ms}$ (bottom trace) for control (gray trace) and after the bath application of the L-VDCC blocker nimodipine ($10\text{ }\mu\text{M}$; black trace). **C**, Average peak amplitudes for pairing an EPSP with three APs (50 Hz) at $\Delta t = -10\text{ ms}$ and $\Delta t = +10\text{ ms}$ for control (gray circles) and after the bath application of nimodipine (black circles). Nimodipine significantly reduced the peak amplitude for $\Delta t = -10\text{ ms}$ ($*p < 0.05; n = 3$) but had no effect on the peak amplitude for $\Delta t = +10\text{ ms}$. **D**, The L-VDCC blocker nimodipine had no effect on the induction of LTP. **E**, Nimodipine also had no effect on the induction of LTD. **F**, Blocking T-VDCCs with $Ni^{2+}$ ($50\text{ }\mu\text{M}$) abolished the induction of LTD for the pairing protocol of one AP and an EPSP at $\Delta t = -10\text{ ms}$. **G**, In contrast, $Ni^{2+}$ had no effect on the induction of LTD for the pairing protocol of three APs (50 Hz) and an EPSP at $\Delta t = -10\text{ ms}$. **H**, The sequence of three APs (50 Hz) and an EPSP for $\Delta t = -10\text{ ms}$ evoked $[Ca^{2+}]_i$ transients in the spine, indicated in the fluorescence image by the dashed line under control conditions (gray trace) and a corresponding smaller $[Ca^{2+}]_i$ transient after the bath application of $Ni^{2+}$ (black trace). **I**, Summary of LTD induction with pharmacological block of VDCCs. $Ni^{2+}$ blocked the induction of LTD if one AP preceded the EPSP during pairing ($p > 0.5; n = 3$). Neither $Ni^{2+}$ nor nimodipine alone had an effect on LTD induction if a burst of three APs (50 Hz) preceded the EPSP ($*p < 0.05; n = 4-8$). Bath application of nimodipine and $Ni^{2+}$ together blocked LTD ($p > 0.5; n = 3$).

[Data from Figure 7I represented as a table]
<table>
  <thead>
    <tr>
        <th>Condition</th>
        <th>change in EPSP amplitude</th>
    </tr>
  </thead>
  <tbody>
    <tr>
        <td>3a1e, nimodipine</td>
        <td>0.7 [approx]</td>
    </tr>
    <tr>
        <td>3a1e, Ni2+</td>
        <td>0.7 [approx]</td>
    </tr>
    <tr>
        <td>3a1e, Ni2+ + nimodipine</td>
        <td>1.0 [approx]</td>
    </tr>
    <tr>
        <td>1a1e, Ni2+</td>
        <td>1.0 [approx]</td>
    </tr>
  </tbody>
</table>
*(Note: Values in the table are estimated from the bar chart in Figure 7I. Asterisks indicate statistical significance.)*

Nevian and Sakmann • Spine Ca<sup>2+</sup> and Spike-Timing-Dependent Plasticity
J. Neurosci., October 25, 2006 • 26(43):11001–11013 • 11009

[The image shows a scatter plot with sigmoid curve fits. The y-axis is "change in EPSP amplitude" ranging from 0 to 3. The x-axis is "$(\Delta G/R)_{max}$" ranging from 0.00 to 0.30. There are three data series: open circles (APs following EPSP), open boxes (APs preceding EPSP), and red circles (APs preceded and followed EPSP). A shaded gray area highlights the region where similar peak levels of $[Ca^{2+}]_i$ can result in either LTP or LTD.]

**Figure 8.** The peak $[Ca^{2+}]_i$ amplitude does not predict LTP or LTD. Summary plot of the average change in EPSP amplitude versus the average peak $[Ca^{2+}]_i$ amplitude expressed as $(\Delta G/R)_{max}$. The change in EPSP amplitude is not a unique function of the peak $[Ca^{2+}]_i$ amplitude. Similar peak levels of $[Ca^{2+}]_i$ can result in either LTP or LTD (data points in shaded area). The clustering of induction protocols in which the APs either precede the EPSP (open boxes) or follow the EPSP (open circles) indicate that the inductions of LTP and LTD are separate processes. Fitting sigmoid functions to the data sets (solid blue line, APs following; solid green line, APs preceding) shows that LTP and LTD both depend on the peak $[Ca^{2+}]_i$ amplitude. The protocols in which APs preceded and followed the EPSP (red circles) fall close to the curve for the $Ca^{2+}$ dependence of LTP.

$0.2; n = 7$), suggesting that the recorded spine $[Ca^{2+}]_i$ transients did not change significantly during the pairing period (Fig. 10$D$). Finally, blocking $IP_3$-dependent release of $Ca^{2+}$ from internal stores by intracellular application of heparin (400 U/ml) through the patch pipette had no effect on LTD ($0.68 \pm 0.12; p < 0.05; n = 9$) (Fig. 10$C$). The resulting pharmacological signature for LTD (Fig. 10$F$) suggests retrograde endocannabinoid signaling, mediated by G-coupled mGluR and PLC activation.

## Discussion
We tested several *in vivo*-like AP activity patterns in a local L2/3-to-L2/3 pyramid connection of somatosensory cortex for the induction of LTP or LTD and measured the corresponding $[Ca^{2+}]_i$ transients in spines. The results indicate that synaptic potentiation and depression are two processes separately controlled by a rise in postsynaptic $[Ca^{2+}]_i$. The magnitude of the long-term change in synaptic transmission depends on the volume-averaged elevation of $[Ca^{2+}]_i$ in a spine, but the direction of the change is controlled by the activation of a mGluR-dependent signaling cascade.

### Burst-timing-dependent plasticity
LTP induction required the pairing of an EPSP with a minimal burst of two APs at frequencies $>20$ Hz, whereas the induction of LTD was less sensitive to the properties of the AP burst. These results are comparable to findings in pairs of L5B pyramidal neurons in which low-frequency pairings of a unitary EPSP and a single AP failed to induce potentiation (Markram et al., 1997; Sjostrom et al., 2001). Pyramidal neurons in the CA1 region of the hippocampus show a developmental switch from single AP-induced LTP to burst LTP (Pike et al., 1999). In contrast, cultured hippocampal neurons (Bi and Poo, 1998), L2/3 pyramids in visual cortex (Froemke and Dan, 2002), and the L4-to-L2/3 spiny stellate-to-pyramid connection in somatosensory cortex (Feldman, 2000) are potentiated by a single AP that follows an EPSP. The $[Ca^{2+}]_i$ transients evoked by pairing an EPSP with a single AP in our experimental conditions indicated that the rise in $[Ca^{2+}]_i$ was not sufficient to induce LTP. This might be different for other connections, depending on the ion channel distribution in dendrites and spines, which influence the backpropagation of APs and $Ca^{2+}$ signaling (Magee and Johnston, 1995; Schiller et al., 1995; Spruston et al., 1995; Sabatini and Svoboda, 2000; Sabatini et al., 2001; Waters and Helmchen, 2004; Gasparini and Magee, 2006). The requirement for AP bursts to induce synaptic modifications might provide stability of connectivity (Lisman and Spruston, 2005), because *in vivo* cortical activity is sparse and single APs occur more frequently than bursts (Lee et al., 2006; Waters and Helmchen, 2006).

### Postsynaptic $Ca^{2+}$ transients and synaptic plasticity
The induction of both LTP and LTD was equally sensitive to loading of the postsynaptic cell with the rapidly equilibrating $Ca^{2+}$ buffer BAPTA or the slower buffer EGTA. Therefore, the volume-averaged elevation of spine $[Ca^{2+}]_i$ is an important factor for the induction of changes in synaptic strength. The putative $Ca^{2+}$ sensors that trigger the induction of LTP or LTD presumably are separated from the $Ca^{2+}$ entry site by several tens of nanometers (Neher, 1998). They might be mobile $Ca^{2+}$ buffers, which compete with $Ca^{2+}$ extrusion and additional fixed $Ca^{2+}$ buffers. Another possibility is that the $Ca^{2+}$ sensors have slower binding kinetics for $Ca^{2+}$ than EGTA. In this case the $Ca^{2+}$ sensors might be localized in close proximity to $Ca^{2+}$ channels, but still they would be sensitive mainly to the volume-averaged increases in $[Ca^{2+}]_i$. In the case of LTP, calmodulin is a possible candidate $Ca^{2+}$ sensor. It can activate $Ca^{2+}$/calmodulin-dependent protein kinase II (CaMKII), which translocates to the plasma membrane of spines (Bayer et al., 2001; Gleason et al., 2003; Otmakhov et al., 2004) and phosphorylates AMPA receptors (Malenka et al., 1989; Barria et al., 1997; Hayashi et al., 2000; Lisman et al., 2002). Its mobility and its presumably almost homogeneous distribution in the spine head could account for the similar effects of BAPTA and EGTA.

One hypothesis relating changes in synaptic efficacy and postsynaptic $[Ca^{2+}]_i$ transients suggests that the level of postsynaptic $[Ca^{2+}]_i$ elevation determines whether a synapse is potentiated or depressed (Lisman, 1989; Artola and Singer, 1993; Hansel et al., 1997). Below a threshold of $[Ca^{2+}]_i$ the synaptic efficacy remains unaffected. In an "intermediate" range of $[Ca^{2+}]_i$ elevations LTD is induced, whereas "large" increases in $[Ca^{2+}]_i$ induce LTP (Zucker, 1999). A number of reports are in accordance with this view (Cummings et al., 1996; Cho et al., 2001). Uncaging experiments of $Ca^{2+}$ showed a correlation between the peak level of volume-averaged $[Ca^{2+}]_i$ and the direction of changes in synaptic strength (Yang et al., 1999). $Ca^{2+}$ imaging in dendrites and somata also suggested a direct correlation (Hansel et al., 1996, 1997; Cormier et al., 2001; Ismailov et al., 2004; Gall et al., 2005). Simulations using the volume-averaged peak $[Ca^{2+}]_i$ amplitude as a readout to determine the direction of change in synaptic efficacy reproduced in part the experimentally determined STDP curves (Karmarkar et al., 2002; Shouval et al., 2002; Shouval and Kalantzis, 2005). In contrast, other reports failed to find a correlation between peak $[Ca^{2+}]_i$ and the direction of change in synaptic strength (Neveu and Zucker, 1996; Wang et al., 2005). Therefore, it was an unresolved issue whether and how a single variable, the global peak $[Ca^{2+}]_i$ amplitude, could induce LTD or LTP differentially. It was suggested that the time course of the $Ca^{2+}$

11010 • J. Neurosci., October 25, 2006 • 26(43):11001–11013
Nevian and Sakmann • Spine Ca<sup>2+</sup> and Spike-Timing-Dependent Plasticity

transient (Rubin et al., 2005) or that the stochastic opening of the NMDA receptor (Shouval and Kalantzis, 2005) can account for the differential induction of LTP and LTD.

Our results show clearly that induction protocols that increase volume-averaged [Ca<sup>2+</sup>]<sub>i</sub> to similar levels can induce either LTP or LTD. We conclude that the volume-averaged peak [Ca<sup>2+</sup>]<sub>i</sub> amplitude is not the only determinant for the direction of modification in synaptic efficacy. Thus the Ca<sup>2+</sup> dependence of STDP induction protocols does not conform to a Ca<sup>2+</sup> control hypothesis (Bear et al., 1987; Artola and Singer, 1993). However, we confirmed that both the inductions of LTD and LTP require the rise of [Ca<sup>2+</sup>]<sub>i</sub> above a threshold, which is approximately twofold higher for the induction of LTP than for LTD. The large Ca<sup>2+</sup> influx during coincident activation of NMDARs by glutamate and depolarization by backpropagating APs resulted in LTP. In contrast, the induction of LTD required a transient increase in [Ca<sup>2+</sup>]<sub>i</sub> through VDCCs and the coactivation of mGluRs. In accordance with these results it was suggested that a second coincidence detector, triggering only LTD induction, could resolve the problem of bidirectional plasticity (Karmarkar and Buonomano, 2002; Bender et al., 2006). mGluRs are a component in several signaling cascades resulting in changes of synaptic strength (Anwyl, 1999, 2006; Bortolotto et al., 1999). Block of mGluRs did not modulate the postsynaptic [Ca<sup>2+</sup>]<sub>i</sub> transients (Brenowitz and Regehr, 2005), excluding a contribution of Ca<sup>2+</sup> release from internal stores (Emptage et al., 2003). In agreement, blocking Ca<sup>2+</sup> release from internal stores had no effect on burst-pairing-induced LTD in contrast to other forms of LTD (Otani et al., 2002; Bender et al., 2006). We suggest that mGluR activation can be regarded as a postsynaptic switch, which sets the sign for the direction of the change in synaptic strength. Activation of this switch results in the induction of LTD, independent of the [Ca<sup>2+</sup>]<sub>i</sub> levels that have been reached. If the switch is not activated and [Ca<sup>2+</sup>]<sub>i</sub> levels are sufficiently large, LTP is induced. The induction of LTD is triggered by AP-evoked Ca<sup>2+</sup> influx through VDCCs and presumably the subsequent binding of Ca<sup>2+</sup> to proteins of the mGluR-coupled signaling cascade before its activation. In support of this view, several proteins of the G-protein-coupled signaling cascade require the binding of Ca<sup>2+</sup> and are thought to be involved in the induction of LTD (Artola and Singer, 1993; Daniel et al., 1998; Kemp and Bashir, 2001). The sequence for the induction of LTD is qualitatively different from the sequence of activation for the induction of LTP. In the case of LTP, synaptic activation and thus mGluR activation precede (or are simultaneous to) Ca<sup>2+</sup> influx, rendering the LTD pathway silent. The requirement for Ca<sup>2+</sup> influx preceding mGluR activation suggests the presence of a second coincidence detection mechanism in L2/3 pyramidal neurons. As a consequence, spines devoid of mGluRs might express only LTP or require a different mechanism for the induction of LTD, like dephosphorylation or removal of AMPA receptors (Carroll et al., 1999; Lee et al., 2000; Froemke et al., 2005; O’Connor et al., 2005).

### Signaling cascade for the induction of LTD
We identified the downstream signaling cascade for spike-timing-dependent LTD to involve PLC and endocannabinoids, similar to heterosynaptic depression of GABAergic transmission in hippocampus (Chevaleyre and Castillo, 2003). This suggested the signaling by group I mGluRs (Coutinho and Knopfel, 2002).

[The image shows a series of panels labeled A through G illustrating the effect of the mGluR blocker MCPG on calcium transients and synaptic plasticity.]

**Figure 9.** Effect of the mGluR blocker MCPG on [Ca<sup>2+</sup>]<sub>i</sub> transients and synaptic plasticity. **A**, Two-photon image of an active synaptic spine. The dashed line indicates the position of the line scan. **B**, [Ca<sup>2+</sup>]<sub>i</sub> transients for the sequences of an EPSP and three APs (50 Hz) at $\Delta t = +10$ ms (top traces), three APs (50 Hz) and an EPSP at $\Delta t = -10$ ms (middle traces), and three APs (100 Hz) and an EPSP at $\Delta t = -10$ ms (bottom traces) for control (gray traces) and after bath application of the mGluR blocker MCPG (500 $\mu$M; black traces). The fluorescence traces in the presence of MCPG do not differ from control. **C**, Summary plot of the peak [Ca<sup>2+</sup>]<sub>i</sub> amplitudes for control versus peak [Ca<sup>2+</sup>]<sub>i</sub> amplitude after the bath application of MCPG. The points for the stimulation protocols depicted in **B** fall close to the unity line (dashed line), indicating no significant effect of MCPG on the peak amplitude of the [Ca<sup>2+</sup>]<sub>i</sub> transients. **D**, Summary of normalized EPSP amplitudes for pairing an EPSP with three APs (50 Hz) at $\Delta t = +10$ ms under control conditions (gray trace) and for experiments in MCPG (black trace). LTP induction in the presence of MCPG was not different from control. Dashed lines indicate no change in EPSP amplitude. **E**, Summary of normalized EPSP amplitudes for pairing three APs (50 Hz) with an EPSP at $\Delta t = -10$ ms under control conditions (gray trace) and for experiments in MCPG (black trace). LTD was abolished in this case. **F**, Summary of normalized EPSP amplitudes for pairing three APs (100 Hz) with an EPSP at $\Delta t = -10$ ms under control conditions (gray trace) and for experiments in MCPG (black trace). The LTD under control conditions was reversed to LTP by MCPG. **G**, Summary plot of change in EPSP amplitude for the pairing protocols depicted in **B** as a function of peak [Ca<sup>2+</sup>]<sub>i</sub> amplitude under control conditions (gray) and in the presence of MCPG (black). In the presence of MCPG no LTD is induced, and the change in EPSP amplitude becomes a monotonically increasing function of the peak [Ca<sup>2+</sup>]<sub>i</sub> amplitude.

Nevian and Sakmann • Spine Ca<sup>2+</sup> and Spike-Timing-Dependent Plasticity
J. Neurosci., October 25, 2006 • 26(43):11001–11013 • 11011

# A
[The image shows a graph of normalized EPSP versus time (min). At time 0, a protocol of three APs at 50 Hz with an EPSP at $\Delta t = -10$ ms is applied in the presence of U73122. The normalized EPSP remains around 1.0 for 60 minutes.]

# B
[The image shows a two-photon image of an active spine with a scale bar of 1 µm. Below it is a line scan trace showing $\Delta G/R$ transients for control (gray) and U73122 (black) conditions, both showing similar peaks of approximately 0.2 $\Delta G/R$ over 50 ms.]

# C
[The image shows a graph of normalized EPSP versus time (min). At time 0, the pairing protocol is applied in the presence of heparin. The normalized EPSP decreases over 60 minutes to approximately 0.5.]

# D
[The image shows a graph of $(\Delta G/R)_{max}$ versus time (min) during the 10 min induction period. The peak amplitudes remain stable around 0.15.]

# E
[The image shows a graph of normalized EPSP versus time (min). At time 0, the pairing protocol is applied in the presence of AM251. The normalized EPSP remains around 1.0 for 60 minutes.]

# F
<table>
  <thead>
    <tr>
        <th>Condition</th>
        <th>change in EPSP amplitude</th>
    </tr>
  </thead>
  <tbody>
    <tr>
        <td>U73122</td>
        <td>1.1</td>
    </tr>
    <tr>
        <td>heparin</td>
        <td>0.7*</td>
    </tr>
    <tr>
        <td>AM251</td>
        <td>1.1</td>
    </tr>
  </tbody>
</table>

**Figure 10.** Signaling pathway for the induction of LTD. ***A***, Summary of normalized EPSP amplitudes for pairing three APs (50 Hz) with an EPSP at $\Delta t = -10$ ms in the presence of the PLC blocker U73122 (5 µM). U73122 blocked the induction of LTD. ***B***, Two-photon image of an active spine. The dashed line indicates the position of the line scan. U73122 had no effect on the $[Ca^{2+}]_i$ transients evoked by the sequence of three APs (50 Hz) and an EPSP at $\Delta t = -10$ ms (black trace) as compared with control (gray trace). ***C***, Heparin (400 U/ml), a blocker of $Ca^{2+}$ release from internal stores mediated by $IP_3$ receptors, had no effect on the induction of LTD. ***D***, Amplitude of averaged peak $[Ca^{2+}]_i$ transients measured during the induction of LTD with the protocol of three APs (50 Hz) and an EPSP at $\Delta t = -10$ ms under control conditions. Every sixth transient was recorded during the 10 min (0.1 Hz) induction period. No correlation between the peak amplitude and the number of stimuli was found (linear correlation, $r^2 = 0.2; n = 7$). ***E***, Bath application of the CB1 receptor antagonist AM251 (2 µM) blocked the induction of LTD. ***F***, Summary of the signaling pathway for the induction of LTD. Block of PLC and CB1 receptors abolished the induction of LTD ($p > 0.1; n = 4-5$), whereas the block of $Ca^{2+}$ release from internal stores had no effect on the induction of LTD ($*p < 0.05; n = 9$). Dashed lines represent no change in EPSP amplitude.

Endocannabinoids act as a retrograde messenger for the induction of LTD (Sjostrom et al., 2003; Duguid and Sjostrom, 2006). They are synthesized via a PLC-dependent pathway, the efficiency of which is modulated greatly by a rise of $[Ca^{2+}]_i$ (Hashimotodani et al., 2005; Maejima et al., 2005). This suggests that the occurrence of postsynaptic before presynaptic coincident activity is detected by PLC and finally results in the release of endocannabinoids to induce LTD. In this case the expression of LTD is presynaptic, and the changes in synaptic efficacy are attributed to a decrease in release probability (Sjostrom et al., 2003). Nevertheless, the trigger for these changes requires the postsynaptic elevation of $[Ca^{2+}]_i$ coupling postsynaptic induction to presynaptic expression of LTD.

We conclude that postsynaptic $Ca^{2+}$ influx is a necessary trigger for modifications of synaptic strength. However, the amplitude of postsynaptic $[Ca^{2+}]_i$ elevation evoked by naturally occurring neuronal activity patterns is not sufficient to decode the direction of the change. Therefore, LTP and LTD are induced by two separate $Ca^{2+}$ sensors. The selective activation of one of two coincidence detectors for the relative timing of presynaptic and postsynaptic activity can account for the capability of synapses for bidirectional modifications.

### References

Anwyl R (1999) Metabotropic glutamate receptors: electrophysiological properties and role in plasticity. Brain Res Brain Res Rev 29: 83–120.

Anwyl R (2006) Induction and expression mechanisms of postsynaptic NMDA receptor-independent homosynaptic long-term depression. Prog Neurobiol 78:17–37.

Artola A, Singer W (1993) Long-term depression of excitatory synaptic transmission and its relationship to long-term potentiation. Trends Neurosci 16:480–487.

Augustine GJ, Santamaria F, Tanaka K (2003) Local calcium signaling in neurons. Neuron 40:331–346.

Barria A, Muller D, Derkach V, Griffith LC, Soderling TR (1997) Regulatory phosphorylation of AMPA-type glutamate receptors by CaMKII during long-term potentiation. Science 276:2042–2045.

Bayer KU, De Koninck P, Leonard AS, Hell JW, Schulman H (2001) Interaction with the NMDA receptor locks CaMKII in an active conformation. Nature 411:801–805.

Bear MF, Cooper LN, Ebner FF (1987) A physiological basis for a theory of synapse modification. Science 237:42–48.

Bender VA, Bender KJ, Brasier DJ, Feldman DE (2006) Two coincidence detectors for spike timing-dependent plasticity in somatosensory cortex. J Neurosci 26:4166–4177.

Bi GQ, Poo MM (1998) Synaptic modifications in cultured hippocampal neurons: dependence on spike timing, synaptic strength, and postsynaptic cell type. J Neurosci 18:10464–10472.

Bliss TV, Collingridge GL (1993) A synaptic model of memory: long-term potentiation in the hippocampus. Nature 361:31–39.

Bortolotto ZA, Fitzjohn SM, Collingridge GL (1999) Roles of metabotropic glutamate receptors in LTP and LTD in the hippocampus. Curr Opin Neurobiol 9:299–304.

Brecht M, Roth A, Sakmann B (2003) Dynamic receptive fields of reconstructed pyramidal cells in layers 3 and 2 of rat somatosensory barrel cortex. J Physiol (Lond) 553:243–265.

Brenowitz SD, Regehr WG (2005) Associative short-term synaptic plasticity mediated by endocannabinoids. Neuron 45:419–431.

Carroll RC, Lissin DV, von Zastrow M, Nicoll RA, Malenka RC (1999) Rapid redistribution of glutamate receptors contributes to long-term depression in hippocampal cultures. Nat Neurosci 2:454–460.

Chevaleyre V, Castillo PE (2003) Heterosynaptic LTD of hippocampal GABAergic synapses: a novel role of endocannabinoids in regulating excitability. Neuron 38:461–472.

Cho K, Aggleton JP, Brown MW, Bashir ZI (2001) An experimental test of the role of postsynaptic calcium levels in determining synaptic strength using perirhinal cortex of rat. J Physiol (Lond) 532:459–466.

Cormier RJ, Greenwood AC, Connor JA (2001) Bidirectional synaptic plas-

11012 • J. Neurosci., October 25, 2006 • 26(43):11001–11013
Nevian and Sakmann • Spine Ca<sup>2+</sup> and Spike-Timing-Dependent Plasticity

ticity correlated with the magnitude of dendritic calcium transients above a threshold. J Neurophysiol 85:399–406.

Coutinho V, Knöpfel T (2002) Metabotropic glutamate receptors: electrical and chemical signaling properties. The Neuroscientist 8:551–561.

Cummings JA, Mulkey RM, Nicoll RA, Malenka RC (1996) Ca<sup>2+</sup> signaling requirements for long-term depression in the hippocampus. Neuron 16:825–833.

Daniel H, Levenes C, Crepel F (1998) Cellular mechanisms of cerebellar LTD. Trends Neurosci 21:401–407.

Debanne D, Gahwiler BH, Thompson SM (1998) Long-term synaptic plasticity between pairs of individual CA3 pyramidal cells in rat hippocampal slice cultures. J Physiol (Lond) 507:237–247.

Denk W, Strickler JH, Webb WW (1990) Two-photon laser scanning fluorescence microscopy. Science 248:73–76.

Duguid I, Sjostrom PJ (2006) Novel presynaptic mechanisms for coincidence detection in synaptic plasticity. Curr Opin Neurobiol 16:312–322.

Egger V, Feldmeyer D, Sakmann B (1999) Coincidence detection and changes of synaptic efficacy in spiny stellate neurons in rat barrel cortex. Nat Neurosci 2:1098–1105.

Emptage NJ, Reid CA, Fine A, Bliss TV (2003) Optical quantal analysis reveals a presynaptic component of LTP at hippocampal Schaffer–associational synapses. Neuron 38:797–804.

Feldman DE (2000) Timing-based LTP and LTD at vertical inputs to layer II/III pyramidal cells in rat barrel cortex. Neuron 27:45–56.

Franks KM, Sejnowski TJ (2002) Complexity of calcium signaling in synaptic spines. BioEssays 24:1130–1144.

Froemke RC, Dan Y (2002) Spike-timing-dependent synaptic modification induced by natural spike trains. Nature 416:433–438.

Froemke RC, Poo MM, Dan Y (2005) Spike-timing-dependent synaptic plasticity depends on dendritic location. Nature 434:221–225.

Gall D, Prestori F, Sola E, D’Errico A, Roussel C, Forti L, Rossi P, D’Angelo E (2005) Intracellular calcium regulation by burst discharge determines bidirectional long-term synaptic plasticity at the cerebellum input stage. J Neurosci 25:4813–4822.

Gasparini S, Magee JC (2006) State-dependent dendritic computation in hippocampal CA1 pyramidal neurons. J Neurosci 26:2088–2100.

Gleason MR, Higashijima S, Dallman J, Liu K, Mandel G, Fetcho JR (2003) Translocation of CaM kinase II to synaptic sites in vivo. Nat Neurosci 6:217–218.

Hansel C, Artola A, Singer W (1996) Different threshold levels of postsynaptic [Ca<sup>2+</sup>]<sub>i</sub> have to be reached to induce LTP and LTD in neocortical pyramidal cells. J Physiol (Paris) 90:317–319.

Hansel C, Artola A, Singer W (1997) Relation between dendritic Ca<sup>2+</sup> levels and the polarity of synaptic long-term modifications in rat visual cortex neurons. Eur J Neurosci 9:2309–2322.

Hashimotodani Y, Ohno-Shosaku T, Tsubokawa H, Ogata H, Emoto K, Maejima T, Araishi K, Shin HS, Kano M (2005) Phospholipase Cβ serves as a coincidence detector through its Ca<sup>2+</sup> dependency for triggering retrograde endocannabinoid signal. Neuron 45:257–268.

Hayashi Y, Shi SH, Esteban JA, Piccini A, Poncer JC, Malinow R (2000) Driving AMPA receptors into synapses by LTP and CaMKII: requirement for GluR1 and PDZ domain interaction. Science 287:2262–2267.

Ismailov I, Kalikulov D, Inoue T, Friedlander MJ (2004) The kinetic profile of intracellular calcium predicts long-term potentiation and long-term depression. J Neurosci 24:9847–9861.

Karmarkar UR, Buonomano DV (2002) A model of spike-timing-dependent plasticity: one or two coincidence detectors? J Neurophysiol 88:507–513.

Karmarkar UR, Najarian MT, Buonomano DV (2002) Mechanisms and significance of spike-timing-dependent plasticity. Biol Cybern 87:373–382.

Kemp N, Bashir ZI (2001) Long-term depression: a cascade of induction and expression mechanisms. Prog Neurobiol 65:339–365.

Köster HJ, Sakmann B (1998) Calcium dynamics in single spines during coincident pre- and postsynaptic activity depend on relative timing of back-propagating action potentials and subthreshold excitatory postsynaptic potentials. Proc Natl Acad Sci USA 95:9596–9601.

Lee AK, Manns ID, Sakmann B, Brecht M (2006) Whole-cell recordings in freely moving rats. Neuron 51:399–407.

Lee HK, Barbarosie M, Kameyama K, Bear MF, Huganir RL (2000) Regulation of distinct AMPA receptor phosphorylation sites during bidirectional synaptic plasticity. Nature 405:955–959.

Lisman J (1989) A mechanism for the Hebb and the anti-Hebb processes underlying learning and memory. Proc Natl Acad Sci USA 86:9574–9578.

Lisman J, Spruston N (2005) Postsynaptic depolarization requirements for LTP and LTD: a critique of spike timing-dependent plasticity. Nat Neurosci 8:839–841.

Lisman J, Schulman H, Cline H (2002) The molecular basis of CaMKII function in synaptic and behavioural memory. Nat Rev Neurosci 3:175–190.

Lynch G, Larson J, Kelso S, Barrionuevo G, Schottler F (1983) Intracellular injections of EGTA block induction of hippocampal long-term potentiation. Nature 305:719–721.

Maejima T, Oka S, Hashimotodani Y, Ohno-Shosaku T, Aiba A, Wu D, Waku K, Sugiura T, Kano M (2005) Synaptically driven endocannabinoid release requires Ca<sup>2+</sup>-assisted metabotropic glutamate receptor subtype 1 to phospholipase Cβ4 signaling cascade in the cerebellum. J Neurosci 25:6826–6835.

Magee JC, Johnston D (1995) Synaptic activation of voltage-gated channels in the dendrites of hippocampal pyramidal neurons. Science 268:301–304.

Magee JC, Johnston D (1997) A synaptically controlled, associative signal for Hebbian plasticity in hippocampal neurons. Science 275:209–213.

Magee JC, Johnston D (2005) Plasticity of dendritic function. Curr Opin Neurobiol 15:334–342.

Malenka RC, Kauer JA, Perkel DJ, Mauk MD, Kelly PT, Nicoll RA, Waxham MN (1989) An essential role for postsynaptic calmodulin and protein kinase activity in long-term potentiation. Nature 340:554–557.

Markram H, Lübke J, Frotscher M, Sakmann B (1997) Regulation of synaptic efficacy by coincidence of postsynaptic APs and EPSPs. Science 275:213–215.

Mulkey RM, Malenka RC (1992) Mechanisms underlying induction of homosynaptic long-term depression in area CA1 of the hippocampus. Neuron 9:967–975.

Neher E (1998) Usefulness and limitations of linear approximations to the understanding of Ca<sup>2+</sup> signals. Cell Calcium 24:345–357.

Neveu D, Zucker RS (1996) Postsynaptic levels of [Ca<sup>2+</sup>]<sub>i</sub> needed to trigger LTD and LTP. Neuron 16:619–629.

Nevian T, Sakmann B (2004) Single spine Ca<sup>2+</sup> signals evoked by coincident EPSPs and backpropagating action potentials in spiny stellate cells of layer 4 in the juvenile rat somatosensory barrel cortex. J Neurosci 24:1689–1699.

O’Connor DH, Wittenberg GM, Wang SS (2005) Dissection of bidirectional synaptic plasticity into saturable unidirectional processes. J Neurophysiol 94:1565–1573.

Oertner TG, Sabatini BL, Nimchinsky EA, Svoboda K (2002) Facilitation at single synapses probed with optical quantal analysis. Nat Neurosci 5:657–664.

Otani S, Connor JA (1998) Requirement of rapid Ca<sup>2+</sup> entry and synaptic activation of metabotropic glutamate receptors for the induction of long-term depression in adult rat hippocampus. J Physiol (Lond) 511:761–770.

Otani S, Daniel H, Takita M, Crepel F (2002) Long-term depression induced by postsynaptic group II metabotropic glutamate receptors linked to phospholipase C and intracellular calcium rises in rat prefrontal cortex. J Neurosci 22:3434–3444.

Otmakhov N, Tao-Cheng JH, Carpenter S, Asrican B, Dosemeci A, Reese TS, Lisman J (2004) Persistent accumulation of calcium/calmodulin-dependent protein kinase II in dendritic spines after induction of NMDA receptor-dependent chemical long-term potentiation. J Neurosci 24:9324–9331.

Pike FG, Meredith RM, Olding AW, Paulsen O (1999) Postsynaptic bursting is essential for “Hebbian” induction of associative long-term potentiation at excitatory synapses in rat hippocampus. J Physiol (Lond) 518:571–576.

Piomelli D (2003) The molecular logic of endocannabinoid signaling. Nat Rev Neurosci 4:873–884.

Rathenberg J, Nevian T, Witzemann V (2003) High-efficiency transfection of individual neurons using modified electrophysiology techniques. J Neurosci Methods 126:91–98.

Rubin JE, Gerkin RC, Bi GQ, Chow CC (2005) Calcium time course as a signal for spike-timing-dependent plasticity. J Neurophysiol 93:2600–2613.

Sabatini BL, Svoboda K (2000) Analysis of calcium channels in single spines using optical fluctuation analysis. Nature 408:589–593.

Nevian and Sakmann • Spine Ca<sup>2+</sup> and Spike-Timing-Dependent Plasticity J. Neurosci., October 25, 2006 • 26(43):11001–11013 • 11013

Sabatini BL, Maravall M, Svoboda K (2001) Ca<sup>2+</sup> signaling in dendritic spines. Curr Opin Neurobiol 11:349–356.

Schiller J, Helmchen F, Sakmann B (1995) Spatial profile of dendritic calcium transients evoked by action potentials in rat neocortical pyramidal neurones. J Physiol (Lond) 487:583–600.

Senn W (2002) Beyond spike timing: the role of nonlinear plasticity and unreliable synapses. Biol Cybern 87:344–355.

Shouval HZ, Kalantzis G (2005) Stochastic properties of synaptic transmission affect the shape of spike time-dependent plasticity curves. J Neurophysiol 93:1069–1073.

Shouval HZ, Bear MF, Cooper LN (2002) A unified model of NMDA receptor-dependent bidirectional synaptic plasticity. Proc Natl Acad Sci USA 99:10831–10836.

Sjöström PJ, Turrigiano GG, Nelson SB (2001) Rate, timing, and cooperativity jointly determine cortical synaptic plasticity. Neuron 32:1149–1164.

Sjöström PJ, Turrigiano GG, Nelson SB (2003) Neocortical LTD via coincident activation of presynaptic NMDA and cannabinoid receptors. Neuron 39:641–654.

Sjöström PJ, Turrigiano GG, Nelson SB (2004) Endocannabinoid-dependent neocortical layer-5 LTD in the absence of postsynaptic spiking. J Neurophysiol 92:3338–3343.

Song S, Abbott LF (2001) Cortical development and remapping through spike timing-dependent plasticity. Neuron 32:339–350.

Spruston N, Schiller Y, Stuart G, Sakmann B (1995) Activity-dependent action potential invasion and calcium influx into hippocampal CA1 dendrites. Science 268:297–300.

Svoboda K, Denk W, Kleinfeld D, Tank DW (1997) In vivo dendritic calcium dynamics in neocortical pyramidal neurons. Nature 385:161–165.

Svoboda K, Helmchen F, Denk W, Tank DW (1999) Spread of dendritic excitation in layer 2/3 pyramidal neurons in rat barrel cortex in vivo. Nat Neurosci 2:65–73.

Wang HX, Gerkin RC, Nauen DW, Bi GQ (2005) Coactivation and timing-dependent integration of synaptic potentiation and depression. Nat Neurosci 8:187–193.

Waters J, Helmchen F (2004) Boosting of action potential backpropagation by neocortical network activity in vivo. J Neurosci 24:11127–11136.

Waters J, Helmchen F (2006) Background synaptic activity is sparse in neocortex. J Neurosci 26:8267–8277.

Whitlock JR, Heynen AJ, Shuler MG, Bear MF (2006) Learning induces long-term potentiation in the hippocampus. Science 313:1093–1097.

Wimmer VC, Nevian T, Kuner T (2004) Targeted in vivo expression of proteins in the calyx of Held. Pflügers Arch 449:319–333.

Yang SN, Tang YG, Zucker RS (1999) Selective induction of LTP and LTD by postsynaptic [Ca<sup>2+</sup>]<sub>i</sub> elevation. J Neurophysiol 81:781–787.

Yuste R, Denk W (1995) Dendritic spines as basic functional units of neuronal integration. Nature 375:682–684.

Zhu JJ, Connors BW (1999) Intrinsic firing patterns and whisker-evoked synaptic responses of neurons in the rat barrel cortex. J Neurophysiol 81:1171–1183.

Zucker RS (1999) Calcium- and activity-dependent synaptic plasticity. Curr Opin Neurobiol 9:305–313.