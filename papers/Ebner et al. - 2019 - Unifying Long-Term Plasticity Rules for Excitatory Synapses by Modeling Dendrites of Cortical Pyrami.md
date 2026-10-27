Article

Unifying Long-Term Plasticity Rules for Excitatory
Synapses by Modeling Dendrites of Cortical
Pyramidal Neurons
Graphical Abstract

Authors
Christian Ebner, Claudia Clopath,
Peter Jedlicka, Hermann Cuntz

Correspondence
ebner@fias.uni-frankfurt.de

In Brief
Synaptic plasticity is shaped by local
dynamic processes within dendritic trees.
Ebner et al. present a biologically inspired
plasticity rule and study its implications in
a detailed model of a pyramidal cell. They
provide a unified description of rate,
timing, and location dependence and
predict cooperative plasticity in
dendrites.

Highlights
d

A phenomenological synaptic plasticity rule is applied to a
pyramidal neuron model

d

Model reproduces rate-, timing-, and location-dependent
plasticity results

d

Active dendrites allow plasticity via dendritic spikes and
subthreshold events

d

Cooperative plasticity exists across the dendritic tree and
within single branches

Ebner et al., 2019, Cell Reports 29, 4295–4307
December 24, 2019 ª 2019 The Authors.
https://doi.org/10.1016/j.celrep.2019.11.068

Cell Reports

Article
Unifying Long-Term Plasticity Rules
for Excitatory Synapses by Modeling
Dendrites of Cortical Pyramidal Neurons
Christian Ebner,1,2,3,4,9,* Claudia Clopath,5 Peter Jedlicka,1,6,7,8 and Hermann Cuntz1,2,8
1Frankfurt Institute for Advanced Studies, 60438 Frankfurt am Main, Germany
2Ernst Stru
€ngmann Institute (ESI) for Neuroscience in Cooperation with Max Planck Society, 60528 Frankfurt am Main, Germany
3NeuroCure Cluster of Excellence, Charité–Universitätsmedizin Berlin, 10117 Berlin, Germany
4Institute for Biology, Humboldt-Universität zu Berlin, 10117 Berlin, Germany
5Computational Neuroscience Laboratory, Bioengineering Department, Imperial College London, London SW7 2AZ, UK
6Institute of Clinical Neuroanatomy, Neuroscience Center, Goethe University Frankfurt, 60528 Frankfurt am Main, Germany
7ICAR3R–Interdisciplinary Centre for 3Rs in Animal Research, Faculty of Medicine, Justus-Liebig-University, 35392 Giessen, Germany
8These authors contributed equally
9Lead Contact
*Correspondence: ebner@fias.uni-frankfurt.de
https://doi.org/10.1016/j.celrep.2019.11.068

SUMMARY

A large number of experiments have indicated that
precise spike times, firing rates, and synapse locations crucially determine the dynamics of long-term
plasticity induction in excitatory synapses. However,
it remains unknown how plasticity mechanisms of
synapses distributed along dendritic trees cooperate
to produce the wide spectrum of outcomes for
various plasticity protocols. Here, we propose a
four-pathway plasticity framework that is well
grounded in experimental evidence and apply it to
a biophysically realistic cortical pyramidal neuron
model. We show in computer simulations that
several seemingly contradictory experimental landmark studies are consistent with one unifying set of
mechanisms when considering the effects of signal
propagation in dendritic trees with respect to synapse location. Our model identifies specific spatiotemporal contributions of dendritic and axo-somatic
spikes as well as of subthreshold activation of synaptic clusters, providing a unified parsimonious explanation not only for rate and timing dependence but
also for location dependence of synaptic changes.
INTRODUCTION
Adaptive behavior, guided by learning and memory processes,
can be seen as a macroscopic manifestation of microscopic
long-term changes in synaptic strength (Bliss and Collingridge,
1993). Such changes have been proposed to be related to the
causal contribution of a presynaptic (pre) cell to the excitation
of a postsynaptic (post) cell according to Hebbian theory
(Hebb, 1949). Thus, a number of studies exploring various ‘‘spike
timing-dependent plasticity’’ (STDP) (Abbott and Nelson, 2000)
protocols have investigated the relationship of the precise timing

between presynaptic and postsynaptic action potentials (APs)
on the efficacy of synapses. In the simplest arrangement, preAPs preceding post-APs (pre-post, positive timing) by a few milliseconds typically result in synaptic long-term potentiation
(LTP), whereas the opposite order (post-pre, negative timing)
leads to long-term depression (LTD) (Bi and Poo, 1998; Markram
et al., 1997). However, each newly tested plasticity protocol has
led to the discussion of new parameters. In particular, when
bursts of APs are considered, the frequency of these bursts
heavily influences the results of the simple STDP concept. Higher
frequencies tend to increase the strength of LTP at positive timings (Markram et al., 1997; Sjöström et al., 2001) and sometimes
even convert LTD at negative timings into LTP, bypassing the
pre-post timing requirement (Sjöström et al., 2001). Also, the
location of synapses along the dendritic tree was shown to
play an important role, with LTD often becoming more prominent
in distal synapses (Froemke et al., 2005; Sjöström and Häusser,
2006), a likely consequence of voltage attenuation of backpropagating action potentials (bAPs) in dendrites (Stuart et al., 1997),
where LTP was recovered by boosting bAPs through dendritic
current injection or cooperative synaptic inputs (Sjöström and
Häusser, 2006). In addition to these effects of frequency and
location on STDP, plasticity can also be induced by depolarization that originates from other sources besides bAPs in the postsynaptic neuron, e.g., dendritic Ca2+ spikes (Golding et al., 2002;
Kampa et al., 2006; Letzkus et al., 2006), N-methyl-D-aspartate
(NMDA) spikes (Brandalise et al., 2016; Gordon et al., 2006), or
excitatory postsynaptic potentials (EPSPs) alone (Sandler
et al., 2016; Weber et al., 2016). For all these reasons, the
concept of classical STDP as a self-contained mechanism has
been debated (Clopath and Gerstner, 2010; Clopath et al.,
2010; Goldberg et al., 2002; Lisman and Spruston, 2005;
Shouval et al., 2010). It stands to reason that the principle of
STDP is only one manifestation of an underlying general plasticity framework (Feldman, 2012; Shouval et al., 2010). In that
case, the question emerges as to which biophysical pathways
contribute to the results from classical STDP protocols and in
which ways they are related to other plasticity protocols. A large

Cell Reports 29, 4295–4307, December 24, 2019 ª 2019 The Authors. 4295
This is an open access article under the CC BY license (http://creativecommons.org/licenses/by/4.0/).

number of theories and models have been developed with both
phenomenological (Morrison et al., 2008) as well as biophysical
(Graupner and Brunel, 2010) backgrounds that explore these
questions, but only a few have recently proposed a unifying
concept of multiple pre- and postsynaptic plasticity pathways
(Costa et al., 2015) in neuron models with extended dendrites
(Bono and Clopath, 2017; Kastellakis et al., 2016; Krieg and
Triesch, 2014; Solinas et al., 2019; Urbanczik and Senn, 2014).
Although many biophysical details of excitatory long-term synaptic plasticity are still not fully understood, it is widely accepted
that postsynaptic Ca2+ plays a fundamental role. According to
some theories and experiments, low levels of Ca2+ lead to no
changes in synaptic strength, whereas intermediate levels cause
LTD and high levels lead to LTP (Artola and Singer, 1993; Artola
et al., 1990; Graupner and Brunel, 2012; Lisman, 1989; Shouval
et al., 2002). However, more recent experiments have indicated
that the levels of postsynaptic Ca2+ by themselves are not always good predictors for plasticity (Nevian and Sakmann,
2006), and increasing evidence suggests that multiple partly independent signaling routes that use Ca2+ exist (Bender et al.,
2006; Jedlicka and Deller, 2017; Oliet et al., 1997; Sjöström
et al., 2003, 2007), ultimately leading to a mixture of synaptic
changes both expressed at presynaptic and postsynaptic sites
(Sjöström et al., 2007). In our phenomenological plasticity model,
we incorporated four signaling routes that are loosely related to
signaling routes in long-term synaptic plasticity that have been
characterized previously. Our plasticity model is based on and
extends an existing phenomenological voltage-dependent
STDP rule (Clopath and Gerstner, 2010; Clopath et al., 2010).
We show in our simulations that a single, dendritic-location-independent plasticity mechanism is able to reconcile many of the
differences found in experiments, including plasticity measurements that previous models were not able to account for. We
propose, in line with previous suggestions (Feldman, 2012;
Shouval et al., 2010), that concepts such as the Ca2+ level hypothesis mentioned above and classical STDP rules could all
be consequences of the same pathways that strongly depend
on local interactions at the synapse.
RESULTS
A Plasticity Rule Based on Pre- and Postsynaptic
Pathways
In our plasticity model, we introduced four pathways that
contributed to changes in both pre- and postsynaptic weight
factors. Although implemented as a phenomenological rule, its
mechanisms were inspired by well-established biophysical
pathways described in a multitude of experimental studies on
long-term synaptic plasticity.
Briefly, presynaptically expressed LTD (pre-LTD; Figure 1A,
left) occurs when metabotropic glutamate receptors (mGluRs)
and postsynaptic voltage-gated Ca2+ channels (VGCCs) are
activated simultaneously (Heifets and Castillo, 2009). Phospholipase C (PLC) then integrates these two signals in the process
of synthesizing endocannabinoids (eCBs) (Hashimotodani
et al., 2005), which retrogradely act on presynaptic type 1 cannabinoid receptors (CB1Rs) to reduce transmitter release probability (Heifets and Castillo, 2009), causing pre-LTD. Presynaptically

4296 Cell Reports 29, 4295–4307, December 24, 2019

expressed LTP (pre-LTP) is thought to occur when postsynaptic
L-type VGCCs (L-VGCCs) are activated, presumably triggering
synthesis of nitric oxide (NO) (Padamsey et al., 2017; Pigott
and Garthwaite, 2016), possibly by calmodulin (CaM) at nitric oxide synthases (NOSs) (Abu-Soud et al., 1994). NO retrogradely
acts on presynaptic guanylyl cyclase (GC) (Koesling et al.,
2004), triggering a signaling chain that is combined with a presynaptic signal by a presynaptic coincidence detector that has
yet to be discovered (Padamsey et al., 2017). Postsynaptically
expressed LTD and LTP (post-LTD/-LTP; Figure 1A, right) are
described as both being driven by coincident binding of glutamate and depolarization of postsynaptic NMDA receptors
€scher and Malenka, 2012). Strong NMDAR-gated
(NMDARs) (Lu
Ca2+ influx activates protein kinases, such as Ca2+/CaM-dependent protein kinase II (CaMKII), whereas weak Ca2+ influx activates their counterpart molecules, protein phosphatases such
as protein phosphatase 1 and calcineurin (Lisman, 1989).
Kinases increase and phosphatases decrease synaptic efficacy
determined by a-amino-3-hydroxy-5-methyl-4-isoxazolepropionic acid receptors (AMPARs), essentially forming complementary mechanisms of postsynaptic LTD and LTP induction
€scher and Malenka, 2012). Studies indicate that CaMKII is
(Lu
able to phosphorylate itself (autophosphorylation) due to its specific subunit structure and that this process is more likely to take
effect if pulses of Ca2+ bound to CaM are applied rapidly (De Koninck and Schulman, 1998), suggesting that this mechanism
could play a role in the frequency dependence of plasticity.
In our model, pre-LTD (indicated by the variable E; Figure 1B,
left; see also Figure S1) was induced when the low-pass-filtered
postsynaptic voltage trace T coincided with the brief presynaptic
signal D. Due to the transient nature of D, pre-LTD was only
induced if the postsynaptic cell experienced depolarization
shortly before the presynaptic signal, e.g., if the stimulation
included a post-pre pair (Figure 1C, medium green color in left
versus right column). Consequently, this mechanism only detected post-pre timings and was insensitive to pre-post timings,
consistent with pre-LTD in experimental studies (Nevian and
Sakmann, 2006; Sjöström et al., 2003). Pre-LTP (indicated by
the variable X; Figure 1B, left) in our model required that two
consecutively filtered traces based on postsynaptic voltage Na
and Nb coincided to result in a trace N. Owing to Nb being filtered
from Na , N was sensitive to the frequency of postsynaptic events
during postsynaptic activity due to summation (Figure 1C, light
violet color). Only if N was sufficiently elevated during the occurrence of the slow presynaptic signal Z, pre-LTP was switched
on. Post-LTD (indicated by the variable P; Figure 1B, right) was
modeled by calculating the coincidence of the slow presynaptic
signal G and a portion of membrane voltage u. The resulting variable C was subjected to an activation function with two thresh+
olds, namely, q
C and qC (Figure 1D). This formalism was consistent with Ca2+-level-based rules (Artola and Singer, 1993;
Lisman, 1989; Shouval et al., 2002) that have previously been
used for modeling synaptic plasticity (see Discussion). No synaptic weight changes were induced below q
C . Whenever C
+
resided between q
C and qC , post-LTD was activated (Figure 1C,
orange color). To model post-LTP (indicated by the variable K;
Figure 1B, right), the portion of C above qC+ , named Ka , was
used to compute the two slower traces Kb and Kg , which were

A

C

B

D

E

Figure 1. A Plasticity Rule Based on Separate Pathways for Pre- and Postsynaptic Plasticity
(A) Simplified illustration of the biophysical pathways that inspired our model. Pre-LTD (left) is induced if postsynaptic Ca2+ influx through VGCCs coincides with
an mGluR-mediated signaling cascade, causing eCB release by PLC and subsequent downregulation of transmitter release probability (minus sign) by CB1Rs.
Pre-LTP (left) reportedly requires Ca2+ influx through VGCCs, triggering NO synthesis, which is then detected by GC and integrated with a presynaptic signal by
an unknown presynaptic coincidence detector to increase release probability (plus sign). Post-LTD and post-LTP (right) are driven by coincident depolarization
and activation of NMDARs. Lower amounts of NMDAR-gated Ca2+ (thin arrow) activate phosphatases (P), causing a reduction of AMPAR efficacy (minus sign).
Higher amounts of NMDAR-gated Ca2+ (bold arrow) activate kinases (K), causing an increase in AMPAR efficacy (plus sign). Coincidence detectors of pre- and
postsynaptic activity are indicated with orange color.
(B) Model interpretation and abstraction of the biophysical pathways in (A). Pre-LTD (left) is dependent on the coincidence between the presynaptic signal D and a
postsynaptic signal based on membrane voltage, reflected in the variable T. The resulting trace E indicates the amount of pre-LTD (minus sign). Pre-LTP (left)
requires coincidence of a presynaptic signal Z and a postsynaptic trace N (based on membrane voltage by Na and Nb ). The amount of pre-LTP (plus sign) is
indicated by X. Post-LTD and post-LTP (right) depend on coincidence between the presynaptic signal G and a portion of membrane voltage u. If the resulting
trace C reaches lower levels, P is activated, indicating post-LTD (minus sign), whereas higher levels activate K, indicating post-LTP (plus sign).
(C) Traces computed by the plasticity rule for two sample stimulation patterns (see top): post-pre-post pairing (left column) and pre-post-post pairing (right
column). The various rows include all fundamental variables of the plasticity rule with color code from (B). Overall synaptic weight W is shown in the bottom row in
red. Loose analogies of the model’s variables to biophysical processes are given in italic type to the right.
(D) Transfer function for post-LTD and post-LTP. Activation of either P or Ka is shown as a function of C.
(E) Voltage clamp simulation while a single presynaptic event is evoked. As a function of clamped voltage, absolute contribution of each of the four pathways is
plotted (black lines), as well as the overall relative weight change (red line).
See also Figure S1.

filtered versions of Ka . The variable r limited the sum of Ka and
Kb . Post-LTP was only switched on when Ka , Kb , and Kg were
nonzero. Thus, similarly to pre-LTP, this mechanism was frequency dependent. A pre-post-post protocol, therefore, evoked

considerably more post-LTP than a post–pre-post protocol (Figure 1C, turquoise color in right versus left column).
The rule’s voltage dependence is demonstrated in Figure 1E.
One single presynaptic event was evoked while the postsynaptic

Cell Reports 29, 4295–4307, December 24, 2019 4297

A

B

C

Figure 2. Effects of Synapse Location on Rate- and Timing-Dependent Plasticity
(A) Stimulation protocol (Sjöström et al., 2001) (inset) and postsynaptic voltage
profile of the neuron model with five somatic spikes (50 Hz) starting at time
10 ms. Voltage profile is shown as a function of time and location along one
specific path from the soma into the apical dendrite.
(B) Model results as relative weight changes (lines) when fit to experimental
data (filled circles, mean ± SEM) at a proximal location (90 mm in the model;
left) and at a distal location (669 mm in the model; right) for different intra-burst
frequencies. Pre- and postsynaptic bursts were shifted by either +10 ms (blue)
or 10 ms (red). Experimental data was recreated from Sjöström et al. (2001)
Figure 1D and Figure 7B and from Sjöström and Häusser (2006) Figure 3 (by
using the exponential fit of inset data).
(C) Spatiotemporal plasticity windows showing relative weight changes (color
coded) as a function of burst timing (x axis) and distance from the soma (y axis;
see cell morphology on the left). Synaptic weight changes were calculated at
41 different locations and 101 different timings.
See also Figure S2.

4298 Cell Reports 29, 4295–4307, December 24, 2019

cell was clamped to values between 75 mV and 15 mV. Voltages below 60 mV led to no change in weight, whereas
voltages between 60 mV and about 28 mV caused LTD and
voltages above that caused net LTP. Consistent with experiments on voltage dependence of LTD pathways (Oliet et al.,
1997), post-LTD more strongly depended on depolarization
than pre-LTD. To consider all the effects of realistic firing
behavior, active dendrites, and synapse location, we incorporated our plasticity model into a highly detailed cortical layer
5b (L5b) pyramidal cell model (Hay et al., 2011).
Effects of Synapse Location on Rate- and TimingDependent Plasticity
In the first stimulation protocol we used, regular bursts of five
pre- and five postsynaptic APs were evoked at either pre-post
(Dt = + 10 ms) or post-pre (Dt =  10 ms) timings (Figure 2A
inset; see STAR Methods). The frequency within the bursts
(intra-burst frequency) varied between 0.1 Hz and 50 Hz. Even
at 50 Hz, distal dendrites in the neuron model experienced
only weak depolarization due to bAP attenuation (Figure 2A).
We optimized the plasticity rule’s parameters (Table S1, set 1)
to match the experimental data (Sjöström and Häusser, 2006;
Sjöström et al., 2001). At proximal locations (90 mm from the
soma; Figure 2B, left panel), a pre-post timing at a frequency
of 0.1 Hz led to no change in weight, whereas at and above
10 Hz, in accordance with experiments (Sjöström et al., 2001),
LTP was induced. In the model, 0.1-Hz bursts were unable to
cause any relevant summation of postsynaptic traces in either
LTP pathway. In contrast, at 10 Hz and above, such summation
was achieved, leading to LTP. A post-pre timing caused LTD
below a frequency of about 30 Hz and LTP beyond 30 Hz both
in experiments (Sjöström et al., 2001) and in the model. Here,
mainly pre-LTD was initiated in the model at lower frequencies.
However, at higher frequencies, summation in both LTP pathways caused the overall switch. When pre-LTD was blocked in
the model, this caused even stronger LTP, whereas blockade
of pre-LTP substantially reduced the amount of LTP (Figure S2),
which is in line with experimental studies (Sjöström et al., 2007).
At distal locations (669 mm from the soma; Figure 2B, right panel),
LTP was absent for both timings and across all frequencies,
whereas frequencies above 20 Hz resulted in slight LTD. Except
for pre-post at 50 Hz (Sjöström and Häusser, 2006), no further
experimental data were available so all other conditions can be
regarded as predictions by the model. Due to the strong attenuation of bAPs (Figure 2A), the local voltage signal at this distance
was not strong enough to drive LTP pathways and rather caused
LTD.
To illustrate the overall interactions between spike timing, frequency, and location, we created spatiotemporal plasticity windows for this protocol (Figure 2C). At 10 Hz, the plasticity curve at
the shortest distance (about 30 mm from the soma) resembled
the typical relationship described by classical STDP paradigms
(Bi and Poo, 1998). There were local maxima of LTD and LTP
close to 0-ms timing, and with increasing delays, the amount
of plasticity exponentially decayed to zero in both directions.
With increasing distance to the soma, the LTD window (postpre) became wider, whereas the LTP window (pre-post) rapidly
became smaller, as reported in experiments on layer 2/3 (L2/3)

Figure 3. Burst-Timing-Dependent Plasticity
(A) STDP induction protocols pairing single presynaptic events with postsynaptic bursts. Comparison between experimental data (blue, mean ± SEM;
from Nevian and Sakmann, 2006) and our simulations (orange) for only pre- or
postsynaptic activity and for long delays at 50 Hz. Spike timing is defined as
the interval between the onset of the presynaptic event and the first step
current injection of the postsynaptic burst.
(B) As in (A) but with shorter delays between pre- and postsynaptic activity at
50 Hz.
(C) As in (A) but with either one or two instead of three postsynaptic events at
50 Hz.
(D) As in (A) but with three postsynaptic events at either 20 Hz or 100 Hz.

pyramidal cells (Froemke et al., 2005) and ultimately vanished at
far distal locations. Post–pre timings led to LTD across a wide
range of distances from the soma, whereas pre-post timings
led to LTP only below 200 mm, which is in accordance with
experimental results (Sjöström and Häusser, 2006). At higher frequencies, local peaks were visible in addition to the one at 0 ms.
These peaks repeated at multiples of the period (the inverse of
the intra-burst frequency). For example, the 20-Hz condition exhibited peaks at 50 ms, 0 ms, and +50 ms. Furthermore, at
40 Hz and 50 Hz, pre-post timings showed a tendency toward
overlapping LTP phases, overwriting proximal LTD entirely.
This result was partly due to the fact that the presynaptic signals
Z and G (which could loosely represent an unknown presynaptic
signal and glutamate activation of NMDARs, respectively) decayed slowly and summated at higher frequencies (see Figure 1C). These windows clearly demonstrate that no configuration in this stimulation protocol could evoke distal LTP based
solely on the backpropagation of axo-somatic APs, as seen in
experiments (Sjöström and Häusser, 2006).

Burst-Timing-Dependent Plasticity
In a second experimental study, the dynamics of burst-timingdependent plasticity were studied in basal dendrites of L2/3
pyramidal cells (Nevian and Sakmann, 2006). In this case, one
presynaptic spike was paired with a burst consisting of one to
three postsynaptic spikes at different frequencies. We implemented the corresponding protocols in a proximal basal dendrite
(55 mm from the soma) of the L5 pyramidal neuron model,
assuming that the electrophysiological properties of basal dendrites in L2/3 and L5 cells did not differ fundamentally with
respect to this stimulation protocol. Using our plasticity rule
with adjusted parameters (Table S1, set 2), we showed that all
experimental outcomes were reproduced (Figure 3). Plasticity
was not induced when there was only pre- or postsynaptic
activity or when the delays were too long (Figure 3A). Pure
post-pre and pre-post burst pairings at 50 Hz led to similar
amounts of LTD and LTP, respectively, whereas intermediate
timings caused intermediate effects (Figure 3B). In our model,
this was explained by fast saturation of T (which could be loosely
linked to postsynaptic VGCC-Ca2+; see Figure 1C). Pairings with
single postsynaptic events at 50 Hz led to either LTD (10 ms) or
weak LTP (+10 ms), and pairings with two postsynaptic events
caused similar results as with three events (Figure 3C). At
20 Hz, pairing one presynaptic spike with three postsynaptic
spikes resulted in similar outcomes compared to evoking only
one postsynaptic spike at 50 Hz (Figure 3D) because most
time constants in the model were too small to cause substantial
summation effects. When three postsynaptic spikes were generated at a frequency of 100 Hz, a pre-post timing led to very
strong LTP, whereas a post-pre timing led to LTD (Figure 3D).
Our plasticity model, thus, faithfully captured all of the outcomes
observed in the experimental study. The plasticity changes
induced through this set of stimulation protocols were mainly
due to pre-LTD and post-LTP in our plasticity rule. Considering
the loose biophysical analogies of the pathways, this matched
experimental results after pharmacological manipulation, which
indicated that LTD depended on activation of mGluRs and
Ca2+ influx through VGCCs, whereas LTP depended on Ca2+
influx through postsynaptic NMDARs (Nevian and Sakmann,
2006).
Rapid Bursts and Dendritic Ca2+ Spikes
Next, we tested whether active properties of dendrites combined with our plasticity rule could reproduce plastic changes
measured at distal synapses where a pre-post-post-post pairing
led to LTD and a post-post-post-pre pairing induced LTP (Letzkus et al., 2006). Here, single presynaptic spikes were paired
with a rapid (200 Hz) burst of three postsynaptic spikes (Figure 4A, inset). Such rapid bursts were found to sum up distally
and evoke dendritic Ca2+ spikes, which was reproduced by the
L5 pyramidal cell model (Hay et al., 2011) (Figure 4A). We then
found a set of parameters (Table S1, set 3) for the plasticity
model where simulation results matched the experimental outcomes of all four combinations of timings (+10/10 ms) and locations (proximal at 90 mm from the soma and distal at 669 mm
from the soma). At proximal locations, results of both the experiments and our simulations matched classical STDP, as pre-post
led to LTP and post-pre led to LTD (Figure 4B, left). Without any

Cell Reports 29, 4295–4307, December 24, 2019 4299

A

B

C

Figure 4. Rapid Bursts and Dendritic Ca2+ Spikes
(A) Stimulation protocol (Letzkus et al., 2006) (inset) and postsynaptic voltage
profile of the neuron model with a dendritic spike induced by high-frequency
somatic stimulation (200 Hz) starting at time 10 ms.
(B) Relative weight changes in the model (orange) and experimental data (blue,
mean ± SEM; from Letzkus et al., 2006, their Figure 5) are plotted as a function
of burst timing for proximal (left, 90 mm) and distal (right, 669 mm) locations
along the apical dendrite.
(C) Spatiotemporal plasticity windows showing relative weight change (color
coded) as a function of burst timing (x axis) and distance from the soma (y axis;
see cell morphology on the left). In the control condition (left), all pathways
were functional. Downward triangles indicate proximal locations, and upward
triangles indicate distal locations used for the model simulations in (B). In the
blockade condition (right), post-LTD was deactivated.
See also Figure S3.

distance-dependent changes to the plasticity rule, our simulations then also captured the experimentally observed plasticity
switch at distal locations, where pre-post timings caused LTD
and post-pre timings caused LTP (Figure 4B, right). We found
the crucial property of these results to be the dendritic spike. It
was delayed by about 20 ms compared to the first somatic AP

4300 Cell Reports 29, 4295–4307, December 24, 2019

and provided a long-lasting depolarization (Figure 4A). In the
model, with pre-post stimulation, the presynaptic signal G at
the distal synapse had already decayed substantially when the
dendritic spike occurred and, thus, caused only intermediate
elevation of C (which could be an abstraction of NMDAR-gated
Ca2+) over most of the duration (Figure S3, left). Consequently,
P was activated more strongly than K (loosely describing phosphatase-kinase competition), which caused post-LTD to surpass post-LTP. Conversely, with post-pre stimulation, the presynaptic event, although 10 ms late, strongly coincided with
the delayed dendritic spike, elevating C beyond the threshold
qC+ for a long duration and, thus, causing post-LTP to surpass
post-LTD (Figure S3, right). To investigate this in detail, we also
visualized the spatiotemporal plasticity window for this induction
protocol (Figure 4C, left panel). There was a switch at around
200 mm from the soma where the classical proximal timing requirements changed to more complex distal ones that were
shaped by the characteristic voltage curve of the dendritic
Ca2+ spike. Interestingly, we found that between 200 mm and
400 mm, even the combined depolarization of bAPs and the forward-propagating dendritic spike were below LTP requirements
so that only LTD was induced. The four conditions of the experiment (Figures 4C and 4D, triangles) matched the distinct areas in
the spatiotemporal plot. We removed post-LTD to illustrate its
contribution to these effects, which loosely corresponded to
pharmacological inhibition of phosphatases or any other crucial
component within this pathway (Figure 4C, right panel). The simulations predict that LTD is then entirely abolished at positive
timings, leading to larger LTP areas, whereas pre-LTD remains
present at negative timings. We conclude that our plasticity
rule, when implemented at synapses on biophysically realistic
dendrites with local dendritic Ca2+ electrogenesis, is able to
reproduce the counterintuitive STDP data for distal synapses
(Letzkus et al., 2006).
Subthreshold Activation of Small Synaptic Clusters
A recent study found that subthreshold activation of a small cluster of four synapses led to plasticity at thin dendritic branches
that depended on the relative location of the cluster on the
respective branch (Weber et al., 2016). Due to massive increases
in input resistance, such thin, distal dendrites are expected to
experience considerably more powerful voltage transients
caused by synaptic inputs than proximal dendrites (Williams
and Stuart, 2002), which could explain these results. We tested
this idea by implementing a stimulation protocol in which four
synapses located on the same branch segment were strongly
activated in rapid succession and repeated simulations for every
single segment of the neuron model. Such activation of a small
cluster of synapses did indeed result in plasticity, revealing a
gradient in relation to the dendritic location (Figure 5A, boxes
and center). At proximal, thick dendritic segments, LTD was predominant, which switched to LTP at more distal locations close
to the dendritic tips. This gradient was due to dramatic differences in local synaptic potentials that occurred even within single branches (Figure 5B, top left panel). The experimental study
focused mainly on oblique dendrites of hippocampal CA1 pyramidal cells (Weber et al., 2016). In our configuration, we found
weaker relative weight changes in oblique dendrites of the L5

Control

Tuft

Na+/Ca2+ channel block

40
20
0
-20
-40

200 μm

Oblique

Control

B

Relative weight change (%)

A

Na+/Ca2+ channel block

Voltage (mV)

5

-40

Relative weight change (%)

-85
0

25
Time (ms)

50 0
Weber et al. (2016)
Our model

50

25
Time (ms)

50

25

0

-25

Oblique

Tuft

Oblique

Tuft

Figure 5. Subthreshold Activation of Small Synaptic Clusters
(A) Color-coded dendritic maps of plasticity during subthreshold activation of a
cluster of four synapses. Simulations were repeated for all possible locations in
the neuron model (boxes on the left show zoomed-in view of example cluster
locations). Dendritic maps show average cluster weight changes for control
conditions (center) and during block of dendritic Na+ and Ca2+ channels (right).
(B) Voltage traces and average relative weight changes for the example
locations shown in the boxes in (A) for both control (left) and channel block
(right) conditions. Plasticity data recreated from Weber et al., 2016 Figure 7f,
shown as mean ± SEM).

pyramidal cell model, although our results in distal tuft dendrites
agree with the experimental data (Figure 5B, bottom left panel).
Generally, at branch points (light orange and green example
locations in Figure 5), no plasticity or LTD was common, whereas
LTP was typically induced close to the tips (dark orange and
green example locations in Figure 5). Due to the absence of a
clear post-pre timing in this protocol, pre-LTD was ineffective
and the majority of changes were caused by post-LTD and
both LTP pathways. We noticed that synaptic potentials in the

cell model were boosted by dendritic voltage-dependent Na+
and Ca2+ channels. Blocking both of these channel types resulted in overall less LTP and more LTD close to the dendritic
tips in our specific case (Figures 5A and 5B, right). However,
the experimental study showed no significant change in plasticity at oblique dendrites after blocking voltage-dependent
Na+ channels alone (Weber et al., 2016). We conclude that based
on strong input resistance increases, it is possible in thin dendrites to induce bidirectional location-dependent plasticity with
subthreshold synaptic inputs alone, as recently observed in
some experiments (Sandler et al., 2016; Weber et al., 2016).
Coincident Activation of Basal and Apical Tuft Inputs
After confirming the functionality of our plasticity model by reproducing experimental data of several different stimulation protocols, we predicted how synapses in different locations change
in response to more naturally occurring input patterns to pyramidal cells, for which no experimental data exist as of now. As L5
pyramidal cells possess the exclusive property of spanning all
six cortical layers, they could potentially act as integrating units
for different streams of information (Larkum, 2013). A possible
integration mechanism, called backpropagation-activated Ca2+
(BAC) firing, involves coincidence of strong proximal and distal
inputs that may lead to dendritic spikes and bursts of axosomatic APs, thereby changing the output mode of the neuron
(Larkum et al., 1999), which we expected to have an effect on
synaptic plasticity. We implemented a scenario where synapses
were placed randomly across basal and tuft regions of the L5b
neuron model in a similar way as done in a recent study (Shai
et al., 2015) (Figure 6A, shaded areas). In addition to these
‘‘background’’ input synapses, we placed a subset of 10 synapses close to each of the two main spiking zones (Figure 6A,
pipette symbols) and equipped them with our plasticity rule.
The stimulation protocol consisted of a 100-ms phase of random
synaptic activity. In one example, basal activity alone (Figure 6B,
left) evoked a few irregularly occurring axo-somatic APs with
weak impact on tuft dendrites. Apical synapses did not show
plasticity under these circumstances, whereas basal synapses
showed a tendency toward LTD. Apical activity alone (Figure 6B,
center) did not cause any plasticity in either group of synapses.
Coincident activation of both basal and apical synapses (Figure 6B, right) caused BAC firing in the neuron model, evoking
both dendritic spikes as well as bursts of APs. There was no absolute switch toward either LTP or LTD in any of the two groups
of synapses. Instead, synaptic weights diverged from baseline in
both directions. What determines whether a synapse potentiates
or depresses during BAC firing? We monitored the time course
of synaptic weight during the example simulation for both the
most potentiated and most depressed synapse of each group
(Figure 6C). This revealed that synapses underwent LTP if they
were active during BAC firing but experienced LTD if they were
active slightly before or after BAC firing. To assess the underlying
weight distributions, we ran 100 simulations with different
random seeds, resulting in a total of 1,000 plastic apical and
basal synapses in each condition (Figure 6D). We noticed that
dendritic spikes were much more common when both groups
were concurrently active (probability of about 10% in apical
only, 15% in basal only, and 95% in apical and basal conditions).

Cell Reports 29, 4295–4307, December 24, 2019 4301

A

B

Basal

Tuft

Basal + Tuft

Tuft
Voltage (mV)

15

-35

-85

Weight (arb. unit)

1.007

1.000

0.993
0

100

200
Time (ms)

Basal

300

400

0

100

200
Time (ms)

basal plasticity
Basal + Tuft, ba

C

400 0

300

100

200
Time (ms)

300

400

Basal + Tuft, ap
apical plasticity

Voltage (mV)

15

-35

-85

G (arb. unit)

1.0
0.5

Weight (arb. unit)

0
1.007
1.000
0.993
0

200

400
Time (ms)

600

800 0

200

400
Time (ms)

600

800

D

# of synapses

1000

Tuft

Basal

Basal + Tuft

500

0
0.98

1.00

1.02

0.98

Weight (arb. unit)

1.00

Weight (arb. unit)

1.02

0.98

1.00

1.02

Weight (arb. unit)

Figure 6. Coincident Activation of Basal and Apical Tuft Inputs
(A) Morphology of the neuron model highlighting locations of random background inputs at basal (light purple) and tuft (light green) dendrites, as well as locations
of synapses equipped with the plasticity rule at basal (purple pipette) and tuft (green pipette) dendrites.
(B) Example voltage traces and weight changes experienced by plastic synapses at basal (purple) and tuft (green) dendrites. Shaded areas indicate intervals of
active background inputs, which were either basal alone (left), tuft alone (center), or both (right).
(C) Example comparison of most potentiated (red) and most depressed (blue) synapses during coincident basal and tuft inputs, including voltage traces (from B;
top right), presynaptic activation G and weight changes over time.
(D) Weight distribution histograms of 100 randomly initialized simulations for all three conditions.

The plasticity results showed that basal activity alone shifted
basal synaptic weights toward the LTD regime (Figure 6D, left),
whereas apical activity alone did not cause much plasticity (Fig-

4302 Cell Reports 29, 4295–4307, December 24, 2019

ure 6D, center). Coincident basal and apical activity led to almost
the same distribution of basal weights but opened up LTP for
apical synapses (Figure 6D, right). Thus, our simulations predict

that BAC firing potentially gates profound bidirectional changes
in synaptic weights, especially at apical locations.
Local Heterosynaptic Effects and NMDA Spikes
Finally, we wanted to know what our model’s predictions would
be regarding heterosynaptic plasticity effects based on local
voltage differences. In principle, strong depolarization at one
point in the dendrite, e.g., by high-frequency synaptic input,
should have the potential to lead to depression in neighboring
synapses with low-frequency activity (Jedlicka et al., 2015; Jungenitz et al., 2018) if these synapses experienced depolarization
below LTP requirements due to attenuation. We tested this by
placing two clusters of 8 synapses each to a thin apical tuft
dendrite in the neuron model (location a: 1,077 mm from the
soma, location b: 950 mm from the soma; Figure 7A, left). Synapses in each cluster were then randomly activated in two modes,
either uniform (mimicking 8-Hz spontaneous Poisson activity) or
synchronized (mimicking 8-Hz stochastic oscillatory activity;
Figure 7A, right). Using this protocol, we found that when both
clusters were uniformly activated for 350 ms, levels of depolarization were moderate (below 40 mV), NMDA conductances
were relatively small, and weights barely changed (Figure 7B,
left column). However, when distal synapses at location a were
switched to synchronized activation, NMDA spikes were elicited,
causing strong local depolarization (up to about 12 mV) by substantial increases in NMDA conductance and leading to LTP on
average (Figure 7B, right column). In contrast, proximal synapses at location b then experienced moderate prolonged depolarization (up to about 35 mV) without major increases in NMDA
conductance, resulting mainly in LTD (Figure 7B, right column).
Here, the distance of synapses at location b from the origin of
NMDA spikes was far enough so that depolarization had already
been attenuated considerably, preventing LTP and, thus, leading
to LTD. These simulations show that local voltage-based heterosynaptic effects can be modeled using our plasticity rule. The
results suggest that NMDA spikes could serve as powerful triggers for LTP, but due to their spatially restricted profile might
cause opposing effects in weakly active neighboring synapses.
DISCUSSION
In this study, we developed a synaptic plasticity rule that accounts for a wide range of diverse plasticity experiments and
reconciles rate-, timing-, and location-dependent plasticity results with classical Ca2+-level-based rules. Our model has two
major advantages compared to previous plasticity models. First,
whereas most previous rules were developed for point neurons
and neglected dendrite morphology, our rule was implemented
in a realistic dendritic tree, revealing insights into the interaction
between local dynamics of dendritic voltage and plasticity
mechanisms. Second, it allows for reproduction of experimental
results regarding the dendritic spike-induced switch of LTD/LTP
windows at distal apical dendrites of cortical pyramidal cells by
using the same set of mechanisms that produce classical
STDP at proximal dendrites. In contrast to more traditional plasticity rules that specifically rely on spike timing, our approach accounts for spike timing, frequency, and dendritic location dependence of plasticity induction by computing local signals at the

synapse and provides loose analogies to underlying biophysical
mechanisms and pathways. In addition to stimulation protocols
that involved exclusively axo-somatic APs (Nevian and Sakmann, 2006; Sjöström and Häusser, 2006; Sjöström et al.,
2001), our plasticity rule also reproduced results of protocols
that led to more complex voltage curves (Letzkus et al., 2006)
and subthreshold activation (Weber et al., 2016), strengthening
the concept of a more general system of plasticity where STDP
is only one emergent property of many (Feldman, 2012; Shouval
et al., 2010).
The distance-dependent switch of LTP into LTD in AP-based
protocols (Sjöström and Häusser, 2006) (Figure 2) can be explained by voltage attenuation of bAPs, which in our simulations
were not powerful enough to activate LTP pathways distally. In
contrast, in the protocol with burst-induced dendritic spikes (Letzkus et al., 2006) (Figure 4), several factors contributed to the unique
distance-dependent switch of timing requirements. Proximal results in our model were dominated mainly by activation of preLTD, pre-LTP, and post-LTP, which have been proposed to be
the main pathways involved in STDP protocols (Sjöström et al.,
2007). Distally, due to the specific delay of the burst-evoked dendritic spike and the prolonged depolarization it produced, postLTD became much more prominent, being the only pathway that
could cause a decrease in synaptic weight for the pre-post timing
(i.e., presynaptic spike preceding postsynaptic burst; Figure 4B;
Figure S3). In this specific case, the unique property lay in the combination of a delayed pre-post timing (which enables post-LTD
instead of post-LTP due to decay of the presynaptic trace G)
with prolonged depolarization (giving post-LTD enough time to
have an effect despite its low amplitude). We conclude that intervals that are proximally designated as pre-post and post-pre in
this protocol convert to ‘‘more delayed pre-post’’ and ‘‘less delayed pre-post’’ at distal locations, respectively. As such, we think
that this behavior represents a mismatch between proximal and
distal definitions of spike timing. We further conclude that dendritic-spike-based plasticity is indeed timing dependent but uses
a different plasticity window that reflects the different temporal
properties of dendritic versus axo-somatic spikes (Figure 4C). All
of this was achieved without any distance-dependent mechanistic
changes to the rule, suggesting that active dendritic processes are
the main determinants of plasticity. At thin, far distal dendrites with
high input resistance, activation of small synaptic clusters could
cause depolarizing events strong enough to induce plasticity
even without local spikes (Weber et al., 2016) (Figure 5). This could
also be a possible explanation for the subthreshold plasticity
recently found in apical tuft dendrites of L5 pyramidal cells (Sandler
et al., 2016). In addition, our results show that the degree to which
synaptic cooperativity at a subthreshold level leads to plasticity
strongly depends on local dendritic excitability that varies with
location. They further indicate that active dendritic properties in
the form of voltage-dependent Na+ and Ca2+ channels could
potentially serve to boost these cooperativity-mediated signals,
possibly strengthening the plasticity gradient. The demonstration
of subthreshold cooperative plasticity is in accordance with modern dendrite-centered theories of memory (Kastellakis et al., 2015,
2016; Legenstein and Maass, 2011), and this opens up the question of how sub- and suprathreshold plasticity signals interact in
dendrites.

Cell Reports 29, 4295–4307, December 24, 2019 4303

Figure 7. Local Heterosynaptic Effects and NMDA Spikes
(A) Morphology of the cell model with zoomed-in view of synaptic cluster locations (left; blue and purple) and probability functions of synaptic events
(right), either uniform (orange) or synchronized (green). Probability is given as
the chance of a single synapse to activate per step of 1 ms.
(B) Simulation results of uniform activity within both clusters (left column) and
when cluster a was switched to synchronized activity (right column). Traces
show local voltage (black, with voltage trace of cluster a shown in gray for
reference in the panel of cluster b in right column), as well as AMPA (blue), and
NMDA (red) conductance values summed over all synapses in clusters and
weight values of all synapses (light blue and purple lines; dark lines represent
averages). Color code as according to (A).

Our simulations of randomly activated basal and tuft inputs
(Figure 6) predict that BAC firing could act as a gateway mechanism for considerable plasticity at both poles of the neuron by
generating both bursts of APs and dendritic spikes. They further
suggest that only synapses actively contributing to the initiation
of BAC firing undergo LTP, whereas those active during other
times tend to experience LTD, as expected from a rule with
Hebbian character. A possible extension to Hebb’s postulate
in this view could be that synapses that cooperate on their quest
to associate different inputs potentiate, whereas synapses that
do not cooperate and/or do not succeed to establish an associational signal depress. It could, thus, be speculated that BAC
firing, although possibly serving as a signal that couples feedfor-

4304 Cell Reports 29, 4295–4307, December 24, 2019

ward and feedback information in pyramidal cells (Larkum,
2013), also supports potentiation of those synapses that cause
it, thereby increasing the probability that this select subset of
synapses leads to BAC firing at the next time they are active.
Intriguingly, the plasticity of distal, feedback-associated synapses is a current topic in studies exploring the idea of deep
learning principles in the brain (Richards and Lillicrap, 2019).
Finally, we used our plasticity rule to explore heterosynaptic
effects within dendritic branches (Figure 7). Notably, we found
that NMDA spikes may serve as powerful triggers for LTP. This
is biologically plausible because a strong cooperation of synapses
is required to elicit NMDA spikes (Major et al., 2013), creating
a localized coincidence signal without the need of further
synaptic integration. Our results predict that synchronized synaptic activity may cause NMDA spikes and thereby strong LTP (Bono
and Clopath, 2017), but this could potentially depress neighboring
synapses with uncorrelated activity (Jedlicka et al., 2015; Jungenitz et al., 2018). Although such a process might contribute to synaptic homeostasis (Watt and Desai, 2010), this idea is based
purely on dendritic voltage differences in our simulations and
currently neglects other mechanisms of heterosynaptic signaling,
e.g., by astrocytes (Min et al., 2012) or by competition for resources (Triesch et al., 2018). From the perspective of Ca2+level-based rules, these results further highlight that plasticity gradients may exist not only in time (i.e., via STDP) but also in space
by means of localized dendritic potentials.
Our synaptic plasticity rule circumvents two issues that
occurred in previous Ca2+-amplitude-based models. First, it
did not rely on APs with an after-depolarizing tail component to
explain LTD at post-pre timings (Shouval et al., 2002), which
would be needed to achieve intermediate Ca2+ levels in such a
condition. In our rule, this interval was covered by pre-LTD,
which used low-pass-filtered voltage to detect post-pre timings.
Second, a pure amplitude-based model will always exhibit a second LTD window at more delayed pre-post timings (Rubin et al.,
2005), which most experimental results do not support. This second LTD window originates from the fact that the coincidence
signal has to pass the intermediate zone again every time it
decreases, which happens as a consequence of increased
pre-post delay. In our plasticity rule, because post-LTD was
set to a low amplitude compared to post-LTP in all simulations,
it was induced in negligible amounts during brief depolarizations,
such as from APs, even at delayed pre-post timings. However,
long-lasting depolarizations, such as dendritic spikes, could
lead to temporal integration of considerable amounts of postLTD, as seen in our simulations of the burst-induced dendritic
spike protocol. Thus, in addition to the amplitude hypothesis,
our rule also relates to the so-called duration hypothesis, which
states that long-lasting intermediate Ca2+ pulses are needed for
(post-)LTD, whereas short-duration, high-amplitude Ca2+ pulses
induce LTP (Evans and Blackwell, 2015). Another hypothesis that
our rule harmonizes well with is that there might be segregated
Ca2+ pools, which could come in the form of Ca2+ micro- or
nanodomains (Evans and Blackwell, 2015), originating
from distinct sources (Jedlicka and Deller, 2017) and feeding
strongly localized messenger chains. Such a system involving
multiple coincidence detectors is in accordance with previous
theoretical and experimental studies indicating that independent

VGCC- and NMDAR-activated pathways exist (Bender et al.,
2006; Karmarkar and Buonomano, 2002; Oliet et al., 1997; Padamsey et al., 2017; Pigott and Garthwaite, 2016; Sjöström
et al., 2003, 2007).
Because our model is based on a previously developed
voltage-dependent STDP rule (Clopath and Gerstner, 2010; Clopath et al., 2010), it might be worth highlighting some common
key elements. Notably, the implementation of pre-LTD in our
plasticity rule was almost identical to LTD in the Clopath et al.
(2010) rule, where the combination of a discrete presynaptic
event with low-pass-filtered postsynaptic voltage allowed for
precise post-pre coincidence detection. There were also similarities between LTP pathways in our plasticity rule and LTP in the
Clopath et al. (2010) rule. In both implementations, a slow presynaptic trace was multiplied with multiple factors based on
postsynaptic voltage to detect pre-post coincidence. However,
in contrast to our model presented here, there is no LTD mechanism in the Clopath et al. (2010) model that is activated at delayed pre-post timings, even during long-lasting depolarizations.
Numerous extensions could be made to our framework to increase precision and flexibility with regard to the broad range of
plasticity-induction protocols. As more detailed biophysical
models of receptors and proteins emerge, they could replace
the phenomenological components in our plasticity framework.
This could eventually lead to a realistic, fully biophysical model
of plasticity induction. For instance, our presynaptic signal D
could be replaced by realistic modeling of mGluR signaling,
and our postsynaptic pathway activation variables P and K could
be substituted by kinetic models of phosphatases, kinases, and
their binding agents. Furthermore, we had to readjust plasticity
amplitudes (representing the impact of each plasticity pathway)
as parameters to reproduce different experiments (see Table S1;
see also Figure S4). These adjustments could reflect differences
between these experiments, including methodological details
(e.g., animal age, recording temperature, and ionic composition
of solutions), differences between synapse types (Larsen and
Sjöström, 2015), and the lack of plasticity maintenance mechanisms in our rule. In addition, for simulations of longer time intervals and in networks, concepts such as short-term plasticity
(Zucker and Regehr, 2002), neuromodulation (Foncelle et al.,
2018; Gerstner et al., 2018), synaptic scaling (Turrigiano, 2008),
and metaplasticity (Jedlicka et al., 2015) would be needed.
In summary, our simulations indicate that a single general plasticity rule is sufficient to reproduce different outcomes of plasticity
experiments at various dendritic locations, providing a unification
of classical STDP and Ca2+-level-based rules. Our plasticity rule
can be readily combined with detailed neuron models to explore
STDP as well as plasticity mediated by dendritic Ca2+ and Na+
spikes, NMDA spikes, subthreshold activation of synaptic clusters, and any combination of these concepts.
STAR+METHODS
Detailed methods are provided in the online version of this paper
and include the following:
d KEY RESOURCES TABLE
d LEAD CONTACT AND MATERIALS AVAILABILITY

d METHOD DETAILS

B Plasticity Rule
B Synaptic currents
B Neuron Model
B Implementation
B Stimulation protocols
B Sensitivity Analysis
d QUANTIFICATION AND STATISTICAL ANALYSIS
d DATA AND CODE AVAILABILITY

SUPPLEMENTAL INFORMATION
Supplemental Information can be found online at https://doi.org/10.1016/j.
celrep.2019.11.068.
ACKNOWLEDGMENTS
The work was supported by a BMBF grant (no. 01GQ1406, Bernstein Award
2013 to H.C.); BBSRC BB/N013956/1 and BB/N019008/1, Wellcome Trust
200790/Z/16/Z, Simons Foundation 564408, and EPSRC EP/R035806/1 (to
C.C.); University Medical Center Giessen and Marburg (UKGM; to P.J.); and
by LOEWE CePTER-Center for Personalized Translational Epilepsy Research
(to P.J.). It was further supported by a NeuroCure PhD fellowship (DFG Exc 257
to C.E.).
AUTHOR CONTRIBUTIONS
C.E., C.C., P.J., and H.C. conceived the study and wrote the paper. C.E. designed the plasticity rule, performed the numerical simulations, and analyzed
the data.
DECLARATION OF INTERESTS
The authors declare that no competing interests exist.
Received: January 26, 2018
Revised: May 2, 2019
Accepted: November 15, 2019
Published: December 24, 2019
REFERENCES
Abbott, L.F., and Nelson, S.B. (2000). Synaptic plasticity: taming the beast.
Nat. Neurosci. 3 (Suppl), 1178–1183.
Abu-Soud, H.M., Yoho, L.L., and Stuehr, D.J. (1994). Calmodulin controls
neuronal nitric-oxide synthase by a dual mechanism. Activation of intra- and
interdomain electron transfer. J. Biol. Chem. 269, 32047–32050.
Artola, A., and Singer, W. (1993). Long-term depression of excitatory synaptic
transmission and its relationship to long-term potentiation. Trends Neurosci.
16, 480–487.
Artola, A., Bröcher, S., and Singer, W. (1990). Different voltage-dependent
thresholds for inducing long-term depression and long-term potentiation in
slices of rat visual cortex. Nature 347, 69–72.
Bender, V.A., Bender, K.J., Brasier, D.J., and Feldman, D.E. (2006). Two coincidence detectors for spike timing-dependent plasticity in somatosensory
cortex. J. Neurosci. 26, 4166–4177.
Bi, G.Q., and Poo, M.M. (1998). Synaptic modifications in cultured hippocampal neurons: dependence on spike timing, synaptic strength, and postsynaptic
cell type. J. Neurosci. 18, 10464–10472.
Bliss, T.V.P., and Collingridge, G.L. (1993). A synaptic model of memory: longterm potentiation in the hippocampus. Nature 361, 31–39.
Bono, J., and Clopath, C. (2017). Modeling somatic and dendritic spike mediated plasticity at the single neuron and network level. Nat. Commun. 8, 706.

Cell Reports 29, 4295–4307, December 24, 2019 4305

Brandalise, F., Carta, S., Helmchen, F., Lisman, J., and Gerber, U. (2016). Dendritic NMDA spikes are necessary for timing-dependent associative LTP in
CA3 pyramidal cells. Nat. Commun. 7, 13480.

Jedlicka, P., and Deller, T. (2017). Understanding the role of synaptopodin and
the spine apparatus in Hebbian synaptic plasticity—new perspectives and the
need for computational modeling. Neurobiol. Learn. Mem. 138, 21–30.

Carnevale, N.T., and Hines, M.L. (2006). The NEURON Book (Cambridge University Press).

Jedlicka, P., Benuskova, L., and Abraham, W.C. (2015). A voltage-based
STDP rule combined with fast BCM-like metaplasticity accounts for LTP and
concurrent ‘‘heterosynaptic’’ LTD in the dentate gyrus in vivo. PLoS Comput.
Biol. 11, e1004588.

Chao, L.H., Pellicena, P., Deindl, S., Barclay, L.A., Schulman, H., and Kuriyan,
J. (2010). Intersubunit capture of regulatory segments is a component of cooperative CaMKII activation. Nat. Struct. Mol. Biol. 17, 264–272.
Clopath, C., and Gerstner, W. (2010). Voltage and spike timing interact in
STDP—a unified model. Front. Synaptic Neurosci. 2, 25.
€sing, L., Vasilaki, E., and Gerstner, W. (2010). Connectivity reClopath, C., Bu
flects coding: a model of voltage-based STDP with homeostasis. Nat. Neurosci. 13, 344–352.
Costa, R.P., Froemke, R.C., Sjöström, P.J., and van Rossum, M.C.W. (2015).
Unified pre- and postsynaptic long-term plasticity enables reliable and flexible
learning. eLife 4, e09457.
De Koninck, P., and Schulman, H. (1998). Sensitivity of CaM kinase II to the frequency of Ca2+ oscillations. Science 279, 227–230.
Evans, R.C., and Blackwell, K.T. (2015). Calcium: amplitude, duration, or location? Biol. Bull. 228, 75–83.
Feldman, D.E. (2012). The spike-timing dependence of plasticity. Neuron 75,
556–571.
Foncelle, A., Mendes, A., Je˛drzejewska-Szmek, J., Valtcheva, S., Berry, H.,
Blackwell, K.T., and Venance, L. (2018). Modulation of spike-timing dependent
plasticity: Towards the inclusion of a third factor in computational models.
Front. Comput. Neurosci. 12, 49.
Froemke, R.C., Poo, M.-M., and Dan, Y. (2005). Spike-timing-dependent synaptic plasticity depends on dendritic location. Nature 434, 221–225.
Gerstner, W., Lehmann, M., Liakoni, V., Corneil, D., and Brea, J. (2018). Eligibility traces and plasticity on behavioral time scales: Experimental support of
NeoHebbian three-factor learning rules. Front. Neural Circuits 12, 53.
Goldberg, J., Holthoff, K., and Yuste, R. (2002). A problem with Hebb and local
spikes. Trends Neurosci. 25, 433–435.
Golding, N.L., Staff, N.P., and Spruston, N. (2002). Dendritic spikes as a mechanism for cooperative long-term potentiation. Nature 418, 326–331.

Jungenitz, T., Beining, M., Radic, T., Deller, T., Cuntz, H., Jedlicka, P., and
Schwarzacher, S.W. (2018). Structural homo- and heterosynaptic plasticity
in mature and adult newborn rat hippocampal granule cells. Proc. Natl.
Acad. Sci. USA 115, E4670–E4679.
Kampa, B.M., Letzkus, J.J., and Stuart, G.J. (2006). Requirement of dendritic
calcium spikes for induction of spike-timing-dependent synaptic plasticity.
J. Physiol. 574, 283–290.
Karmarkar, U.R., and Buonomano, D.V. (2002). A model of spike-timing
dependent plasticity: one or two coincidence detectors? J. Neurophysiol.
88, 507–513.
Kastellakis, G., Cai, D.J., Mednick, S.C., Silva, A.J., and Poirazi, P. (2015). Synaptic clustering within dendrites: an emerging theory of memory formation.
Prog. Neurobiol. 126, 19–35.
Kastellakis, G., Silva, A.J., and Poirazi, P. (2016). Linking memories across time
via neuronal and dendritic overlaps in model neurons with active dendrites.
Cell Rep. 17, 1491–1504.
Koesling, D., Russwurm, M., Mergia, E., Mullershausen, F., and Friebe, A.
(2004). Nitric oxide-sensitive guanylyl cyclase: structure and regulation. Neurochem. Int. 45, 813–819.
Krieg, D., and Triesch, J. (2014). A unifying theory of synaptic long-term
plasticity based on a sparse distribution of synaptic strength. Front. Synaptic
Neurosci. 6, 3.
Larkum, M. (2013). A cellular mechanism for cortical associations: an organizing principle for the cerebral cortex. Trends Neurosci. 36, 141–151.
Larkum, M.E., Zhu, J.J., and Sakmann, B. (1999). A new cellular mechanism for
coupling inputs arriving at different cortical layers. Nature 398, 338–341.
Larsen, R.S., and Sjöström, P.J. (2015). Synapse-type-specific plasticity in
local circuits. Curr. Opin. Neurobiol. 35, 127–135.

Gordon, U., Polsky, A., and Schiller, J. (2006). Plasticity compartments in basal
dendrites of neocortical pyramidal neurons. J. Neurosci. 26, 12717–12726.

Legenstein, R., and Maass, W. (2011). Branch-specific plasticity enables selforganization of nonlinear computation in single neurons. J. Neurosci. 31,
10787–10802.

Graupner, M., and Brunel, N. (2010). Mechanisms of induction and maintenance of spike-timing dependent plasticity in biophysical synapse models.
Front. Comput. Neurosci. 4, 136.

Letzkus, J.J., Kampa, B.M., and Stuart, G.J. (2006). Learning rules for
spike timing-dependent plasticity depend on dendritic synapse location.
J. Neurosci. 26, 10420–10429.

Graupner, M., and Brunel, N. (2012). Calcium-based plasticity model explains
sensitivity of synaptic changes to spike pattern, rate, and dendritic location.
Proc. Natl. Acad. Sci. USA 109, 3991–3996.

Lisman, J. (1989). A mechanism for the Hebb and the anti-Hebb processes underlying learning and memory. Proc. Natl. Acad. Sci. USA 86, 9574–9578.

Hashimotodani, Y., Ohno-Shosaku, T., Tsubokawa, H., Ogata, H., Emoto, K.,
Maejima, T., Araishi, K., Shin, H.-S., and Kano, M. (2005). Phospholipase
Cbeta serves as a coincidence detector through its Ca2+ dependency for triggering retrograde endocannabinoid signal. Neuron 45, 257–268.

Lisman, J., and Spruston, N. (2005). Postsynaptic depolarization requirements
for LTP and LTD: a critique of spike timing-dependent plasticity. Nat. Neurosci.
8, 839–841.
Lisman, J., Yasuda, R., and Raghavachari, S. (2012). Mechanisms of CaMKII
action in long-term potentiation. Nat. Rev. Neurosci. 13, 169–182.

€rmann, F., Markram, H., and Segev, I. (2011). Models of
Hay, E., Hill, S., Schu
neocortical layer 5b pyramidal cells capturing a wide range of dendritic and
perisomatic active properties. PLoS Comput. Biol. 7, e1002107.

€scher, C., and Malenka, R.C. (2012). NMDA receptor-dependent long-term
Lu
potentiation and long-term depression (LTP/LTD). Cold Spring Harb.
Perspect. Biol. 4, 1–16.

Hebb, D.O. (1949). The Organization of Behaviour: A Neuropsychological
Theory (John Wiley & Sons).

Mainen, Z.F., Malinow, R., and Svoboda, K. (1999). Synaptic calcium transients in single spines indicate that NMDA receptors are not saturated. Nature
399, 151–155.

Heifets, B.D., and Castillo, P.E. (2009). Endocannabinoid signaling and longterm synaptic plasticity. Annu. Rev. Physiol. 71, 283–306.
Ishikawa, T., Sahara, Y., and Takahashi, T. (2002). A single packet of transmitter does not saturate postsynaptic glutamate receptors. Neuron 34,
613–621.
Jahr, C.E., and Stevens, C.F. (1990). Voltage dependence of NMDA-activated
macroscopic conductances predicted by single-channel kinetics. J. Neurosci.
10, 3178–3182.

4306 Cell Reports 29, 4295–4307, December 24, 2019

Major, G., Larkum, M.E., and Schiller, J. (2013). Active properties of neocortical pyramidal neuron dendrites. Annu. Rev. Neurosci. 36, 1–24.
€bke, J., Frotscher, M., and Sakmann, B. (1997). Regulation of
Markram, H., Lu
synaptic efficacy by coincidence of postsynaptic APs and EPSPs. Science
275, 213–215.
McDougal, R.A., Morse, T.M., Carnevale, T., Marenco, L., Wang, R., Migliore,
M., Miller, P.L., Shepherd, G.M., and Hines, M.L. (2017). Twenty years of

ModelDB and beyond: building essential modeling tools for the future of
neuroscience. J. Comput. Neurosci. 42, 1–10.
Meyer, T., Hanson, P.I., Stryer, L., and Schulman, H. (1992). Calmodulin trapping by calcium-calmodulin-dependent protein kinase. Science 256, 1199–
1202.
Min, R., Santello, M., and Nevian, T. (2012). The computational power of astrocyte mediated synaptic plasticity. Front. Comput. Neurosci. 6, 93.
Morrison, A., Diesmann, M., and Gerstner, W. (2008). Phenomenological
models of synaptic plasticity based on spike timing. Biol. Cybern. 98, 459–478.
Nevian, T., and Sakmann, B. (2006). Spine Ca2+ signaling in spike-timingdependent plasticity. J. Neurosci. 26, 11001–11013.
Oliet, S.H., Malenka, R.C., and Nicoll, R.A. (1997). Two distinct forms of longterm depression coexist in CA1 hippocampal pyramidal cells. Neuron 18,
969–982.
Padamsey, Z., and Emptage, N. (2013). Two sides to long-term potentiation: a
view towards reconciliation. Philos. Trans. R. Soc. Lond. B Biol. Sci. 369,
20130154.
Padamsey, Z., Tong, R., and Emptage, N. (2017). Glutamate is required for
depression but not potentiation of long-term presynaptic function. eLife 6,
e29688.
Pepke, S., Kinzer-Ursem, T., Mihalas, S., and Kennedy, M.B. (2010). A dynamic model of interactions of Ca2+, calmodulin, and catalytic subunits of
Ca2+/calmodulin-dependent protein kinase II. PLoS Comput. Biol. 6,
e1000675.
Persechini, A., and Stemmer, P.M. (2002). Calmodulin is a limiting factor in the
cell. Trends Cardiovasc. Med. 12, 32–37.
Pigott, B.M., and Garthwaite, J. (2016). Nitric oxide is required for L-type
Ca(2+) channel-dependent long-term potentiation in the hippocampus. Front.
Synaptic Neurosci. 8, 17.
Poleg-Polsky, A. (2015). Effects of neural morphology and input distribution on
synaptic processing by global and focal NMDA-spikes. PLoS One 10,
e0140254.
Rhodes, P. (2006). The properties and implications of NMDA spikes in neocortical pyramidal cells. J. Neurosci. 26, 6704–6715.
Richards, B.A., and Lillicrap, T.P. (2019). Dendritic solutions to the credit
assignment problem. Curr. Opin. Neurobiol. 54, 28–36.
Rubin, J.E., Gerkin, R.C., Bi, G.-Q., and Chow, C.C. (2005). Calcium time
course as a signal for spike-timing-dependent plasticity. J. Neurophysiol. 93,
2600–2613.
Sandler, M., Shulman, Y., and Schiller, J. (2016). A novel form of local plasticity
in tuft dendrites of neocortical somatosensory layer 5 pyramidal neurons.
Neuron 90, 1028–1042.

Shai, A.S., Anastassiou, C.A., Larkum, M.E., and Koch, C. (2015). Physiology
of layer 5 pyramidal neurons in mouse primary visual cortex: coincidence
detection through bursting. PLoS Comput. Biol. 11, e1004090.
Shouval, H.Z., Bear, M.F., and Cooper, L.N. (2002). A unified model of NMDA
receptor-dependent bidirectional synaptic plasticity. Proc. Natl. Acad. Sci.
USA 99, 10831–10836.
Shouval, H.Z., Wang, S.S.-H., and Wittenberg, G.M. (2010). Spike timing
dependent plasticity: a consequence of more fundamental learning rules.
Front. Comput. Neurosci. 4, 19.
Sjöström, P.J., and Häusser, M. (2006). A cooperative switch determines the
sign of synaptic plasticity in distal dendrites of neocortical pyramidal neurons.
Neuron 51, 227–238.
Sjöström, P.J., Turrigiano, G.G., and Nelson, S.B. (2001). Rate, timing, and
cooperativity jointly determine cortical synaptic plasticity. Neuron 32, 1149–
1164.
Sjöström, P.J., Turrigiano, G.G., and Nelson, S.B. (2003). Neocortical LTD via
coincident activation of presynaptic NMDA and cannabinoid receptors.
Neuron 39, 641–654.
Sjöström, P.J., Turrigiano, G.G., and Nelson, S.B. (2007). Multiple forms of
long-term plasticity at unitary neocortical layer 5 synapses. Neuropharmacology 52, 176–184.
Solinas, S.M.G., Edelmann, E., Leßmann, V., and Migliore, M. (2019). A kinetic
model for Brain-Derived Neurotrophic Factor mediated spike timing-dependent LTP. PLoS Comput. Biol. 15, e1006975.
Stuart, G., Spruston, N., Sakmann, B., and Häusser, M. (1997). Action potential
initiation and backpropagation in neurons of the mammalian CNS. Trends
Neurosci. 20, 125–131.
Triesch, J., Vo, A.D., and Hafner, A.S. (2018). Competition for synaptic building
blocks shapes synaptic plasticity. eLife 7, e37836.
Turrigiano, G.G. (2008). The self-tuning neuron: synaptic scaling of excitatory
synapses. Cell 135, 422–435.
Urbanczik, R., and Senn, W. (2014). Learning by the dendritic prediction of
somatic spiking. Neuron 81, 521–528.
Watt, A.J., and Desai, N.S. (2010). Homeostatic plasticity and STDP: Keeping a
neuron’s cool in a fluctuating world. Front. Synaptic Neurosci. 2, 5.
Weber, J.P., Andrásfalvy, B.K., Polito, M., Magó, Á., Ujfalussy, B.B., and
Makara, J.K. (2016). Location-dependent synaptic plasticity rules by dendritic
spine cooperativity. Nat. Commun. 7, 11380.
Williams, S.R., and Stuart, G.J. (2002). Dependence of EPSP efficacy on synapse location in neocortical pyramidal neurons. Science 295, 1907–1910.
Zucker, R.S., and Regehr, W.G. (2002). Short-term synaptic plasticity. Annu.
Rev. Physiol. 64, 355–405.

Cell Reports 29, 4295–4307, December 24, 2019 4307

STAR+METHODS
KEY RESOURCES TABLE

REAGENT or RESOURCE

SOURCE

IDENTIFIER

Cortical L5b pyramidal cell model

Hay et al., 2011

ModelDB: 139653

Four-pathway phenomenological synaptic plasticity model

This paper

ModelDB: 251493

NEURON

Carnevale and Hines, 2006

RRID: SCR_005393

ModelDB

McDougal et al., 2017

RRID: SCR_007271

Software and Algorithms

LEAD CONTACT AND MATERIALS AVAILABILITY
Further information and requests for resources should be directed to and will be fulfilled by the Lead Contact, Christian Ebner (ebner@
fias.uni-frankfurt.de). This study did not generate new unique reagents.
METHOD DETAILS
Plasticity Rule
Our plasticity rule quantifies the activation of its four separate pathways (Figures 1B, 1C, and S1) directly from the timing of a presynaptic event and the local postsynaptic membrane voltage.
Presynaptic LTD
Presynaptic LTD in this model was inspired by mGluR-CB1R-LTD (Heifets and Castillo, 2009) (Figures 1A and 1B, left panels). For
simplification, local postsynaptic membrane potentials u were dimensionless quantities. T was then a low-pass filtered version of
the portion of u that was above a threshold qTu with a time constant tT
tT $



d
$TðtÞ =  TðtÞ + uðtÞ  qTu + ;
dt

(1)

where for any given value x, the notation ½x + indicated a rectifier, defined as being x for positive values of x and 0 in all other cases.
Using T, we calculated the trace T


TðtÞ = tanh bT $ TðtÞ
(2)

bT =

lnðmT Þ
:
2

(3)

The hyperbolic tangent was used as a sigmoid saturation function in multiple instances below to provide a soft boundary for variables of the model and to loosely relate to binding kinetics of the agents involved. The characteristic saturation in the case of T was
determined by the specific slope mT (Table S1) and different slopes according to Equation 3 were used to compute other traces (see
below). We further defined the presynaptic variable D as a series of delta pulses
X
DðtÞ =
dðt  ti Þ;
(4)
i

with ti representing times of presynaptic events. The coincidence of pre- and postsynaptic signals E was therefore given by
EðtÞ = DðtÞ$TðtÞ;

(5)

which was used as a direct indicator of pre-LTD. A possible link to biophysical processes could be the following: The postsynaptic trace T could loosely represent the amount of VGCC-gated Ca2+ that was bound to PLC at a given time. A minimal depolarization qTu would then be required to open the VGCCs (simplified with a linear increase in permeability), while the binding/unbinding
rate of Ca2+ from PLC would be determined by tT . The delta pulses in D could be related to the signaling cascades evoked by
mGluRs upon glutamate binding and E could loosely represent the amount of synthesized eCBs by PLC (Hashimotodani et al.,
2005).
Presynaptic LTP
Presynaptic LTP was inspired by NO-LTP (Padamsey et al., 2017; Pigott and Garthwaite, 2016; Sjöström et al., 2007) (Figures 1A and
1B, left panels). In a similar way as with T, we defined another trace Na based on low-pass filtered postsynaptic voltage

e1 Cell Reports 29, 4295–4307.e1–e6, December 24, 2019

taN $



d
$Na ðtÞ =  Na ðtÞ + uðtÞ  qNu + ;
dt


Na ðtÞ = tanh baN $ Na ðtÞ ;

(6)

(7)

using a slope maN (Table S1) via baN (as according to Equation 3). We then defined a second postsynaptic trace Nb from
tbN $

d
$Nb ðtÞ =  Nb ðtÞ + Na ðtÞ;
dt



Nb ðtÞ = tanh bbN $ Nb ðtÞ ;

(8)

(9)

using a slope mbN via bbN . The product of Na and Nb surpassing a threshold qN was defined as N
NðtÞ = ½Na ðtÞ$Nb ðtÞ  qN  + :

(10)

In addition, a presynaptic activity trace Z was shaped by the difference of two exponentials and application of the hyperbolic
tangent
ZðtÞ = tanhðbZ $ ðZb ðtÞ  Za ðtÞÞÞ;

(11)

taZ $

d
$Za ðtÞ =  Za ðtÞ + εZ $DðtÞ;
dt

(12)

tbZ $

d
$Zb ðtÞ =  Zb ðtÞ + εZ $DðtÞ;
dt

(13)

1

εZ =
e
uZ =

u
 aZ
tZ

u ;
 bZ

(14)

+ e tZ

 b
taZ $tbZ
t
$ln Za :
b
a
tZ
tZ  tZ

(15)

Here, the triggering of a presynaptic event via the event times in D (Equation 4) also elevated Z in a time course characterized
by the time constants taZ and tbZ . The sole purpose of εZ was to normalize the peak of the trace to 1. Coincidence of presynaptic
signals Z and postsynaptic signals N yielded X
XðtÞ = ZðtÞ$NðtÞ;

(16)

which was used as the indicator for pre-LTP. A possible link of our implementation of pre-LTP to biophysical mechanisms could be
a recently described presynaptic form of LTP (Padamsey and Emptage, 2013). In this view, Na could be related to the influx of Ca2+ via
L-VGCCs (Pigott and Garthwaite, 2016) with a relatively high voltage threshold qN
u (see Table S1), while Nb could be related to a
slower process based on Na such as CaM binding. Based on this perspective, N could be a loose analogy to NOS activation and
NO synthesis via CaM (Abu-Soud et al., 1994), while X could refer to a yet unknown presynaptic coincidence detector based on a
presynaptic signal Z (Padamsey et al., 2017).
Postsynaptic LTD
€scher and Malenka, 2012) (Figures 1A and 1B, right panels).
Postsynaptic LTD in our model was loosely based on NMDAR-LTD (Lu
Here, a third presynaptic trace G was shaped by the difference of two exponential functions in the same way as Z (Equations 11–15),
but with time constants t aG and tbG and a saturation slope mG . We then computed C as the coincidence of G and u above qC
u


C
CðtÞ = GðtÞ$ uðtÞ  qu + :
(17)
Based on C, a trace P was calculated as


  +
1
PðtÞ = CðtÞ  q
2 :
C + $ qC  CðtÞ + $ +
qC  q
C
2

(18)

+
P was chosen to be a quadratic function of C between the thresholds q
C and qC . The peak of the quadratic function was normalized

+
to 1 by removing its dependence on qC and qC , which we found to be useful for optimizing the plasticity amplitude. The amount of
post-LTD was correspondingly directly dependent on P. This implementation of post-LTD could be interpreted as a loose analogy to
the following biophysical processes: G could represent the total amount of activated NMDARs following glutamate binding, where

Cell Reports 29, 4295–4307.e1–e6, December 24, 2019 e2

binding and unbinding kinetics could be determined by t aG and t bG , respectively. In our implementation, multiple events summed up in
G but were limited to a maximum of 1 via the saturating process. This could correspond to the existing proposal that NMDARs of a
synapse are not fully saturated upon a single release event (Ishikawa et al., 2002; Mainen et al., 1999). In this view, C could loosely
2+
block.
represent the total fraction of open NMDARs, where qC
u could be related to the minimal voltage required to release the Mg
2+
could
mark
the
minimal
amount
of
Ca
needed
to
activate
phosphatases
and
therefore
the
start
of
LTD
induction
The threshold q
C
along the Ca2+ continuum. The threshold qC+ on the other hand could then designate the amount of Ca2+ where the competition between phosphatases and kinases reaches an equilibrium and therefore would mark the start of LTP along the continuum. The
2+
+
where phosmaximum of the function, located at the center between q
C and qC , could then loosely represent the amount of Ca
phatases are most active (Lisman, 1989).
Postsynaptic LTP
€scher and Malenka, 2012) (Figures 1A and 1B, right panels). In our model, postPostsynaptic LTP was inspired by NMDAR-LTP (Lu
LTP depended on the coincidence of three traces, denoted Ka , Kb and Kg . Ka from



Ka ðtÞ = tanh baK $ CðtÞ  qC+ + $rðtÞ;
(19)
was limited to a maximum of 1 via the hyperbolic tangent using a slope maK via baK (Table S1; see Equation 3). Ka was additionally
limited by a variable r
rðtÞ = 1  Kb ðtÞ:

(20)

Since Kb was dependent on Ka (see below), r served as a negative feedback signal, thus ensuring that the sum of Ka and Kb could
not be greater than 1. Kb itself was a low-pass filtered version of Ka
d
$K b ðtÞ =  K b ðtÞ + Ka ðtÞ;
dt

tbK $



Kb ðtÞ = tanh bbK $ sbK $ K b ðtÞ ;

(21)

(22)

where sbK was a factor simply used to scale up K b into the range of the saturation function with slope mbK via bbK . It was then further
low-pass filtered with a time constant t gK to compute Kg
d
$Kg ðtÞ =  Kg ðtÞ + Kb ðtÞ:
dt

(23)

KðtÞ = Ka ðtÞ$Kb ðtÞ$Kg ðtÞ;

(24)

tgK $
K was then the product of all three traces

which was directly used as the indicator of post-LTP. A loose analogy to biophysical processes could be the following: Ka could be
related to the instantaneous activation of CaMKII by Ca2+-CaM after enough Ca2+ passed NMDARs so that kinase activation surpassed phosphatase activation, illustrated by the amount of C above qC+ (Lisman, 1989). In this view, the negative feedback trace
r could represent competition among proteins in different states regarding the limited amounts of free CaM in dendritic spines (Persechini and Stemmer, 2002). Kb could then be related to a slower process, such as the amount of a certain configuration of Ca2+-CaM
bound to CaMKII (Pepke et al., 2010) or possibly trapped (Meyer et al., 1992), decaying with a time constant t b . Kg could be linked to a
process based on Kb , such as another configuration of Ca2+-CaM with even slower kinetics (Pepke et al., 2010) or a slow conformational change which might be required for autophosphorylation of CaMKII (Chao et al., 2010). Finally, K could loosely illustrate the
amount of CaMKII that reaches the autonomous state and/or binds to NR2B subunits at any given time, which both have been
proposed to be crucial for LTP (Lisman et al., 2012).
Synaptic Weight
Synaptic weight was the product of both pre- and postsynaptic weight factors
wðtÞ = wpre ðtÞ$wpost ðtÞ:

(25)

The factors were each updated by the sum of their respective pathway indicators
d
LTP
$ wpre ðtÞ =  ALTD
pre $EðtÞ + Apre $XðtÞ$h;
dt

(26)

d
LTP
$ wpost ðtÞ =  ALTD
post $PðtÞ$h + Apost $KðtÞ$h;
dt

(27)

LTP
LTD
LTP
where ALTD
pre , Apre , Apost and Apost were the respective pathway amplitudes (Table S1). Due to its calculation via delta pulses, the preLTD pathway was inherently invariant to changes in integration step size. In contrast, the other three pathways were integrated over
time and continuously (i.e., at each step) updated. We thus introduced a learning rate h (h = 0:025 ms1 for all of our simulations) to

e3 Cell Reports 29, 4295–4307.e1–e6, December 24, 2019

make them independent of changes in integration step size. The presynaptic weight wpre could be related to transmitter release probability, although our model did not explicitly calculate probabilities and synaptic responses should all be regarded as averages. We
thus set the hard bounds to 0%wpre %1. The postsynaptic weight wpost , which would be interpreted as a factor contributing to postsynaptic current, was limited via 0%wpost %5. These bounds prevented both the occurrence of negative weights and excessively
init
strong synapses. In the beginning of each simulation, weight factors were initialized to winit
pre = 0:5 and wpost = 2, leading to a total
weight winit = 1.
Synaptic currents
Synaptic currents were computed as sums of both AMPAR- and NMDAR-mediated components. The AMPAR component was
calculated via

(28)
gAMPA ðtÞ = sAMPA $wpost ðtÞ$ gbAMPA ðtÞ  gaAMPA ðtÞ ;
t aAMPA $

d a
$g
ðtÞ =  gaAMPA ðtÞ + εAMPA $DðtÞ$gmax ;
dt AMPA

(29)

t bAMPA $

d b
$g
ðtÞ =  gbAMPA ðtÞ + εAMPA $DðtÞ$gmax ;
dt AMPA

(30)

where D indicated event times (see Equation 4). The time constants describing glutamate kinetics were t aAMPA = 0:2 ms and
b
tAMPA = 2 ms and εAMPA was calculated from these time constants according to (Equations 14–15). gmax was the maximum synaptic
conductance, which in addition to NMDA/AMPA ratio constants sAMPA and sNMDA was set individually for each stimulation protocol
(see further below). The NMDAR component was calculated via
gNMDA ðtÞ = sNMDA $winit
post $GðtÞ$gmax $

1
;
1 + e0:08$uðtÞ $3:571

(31)

where glutamate kinetics were modeled via the difference of exponentials in G, using taG = 2 ms and t bG = 50 ms (Poleg-Polsky,
2015) (Table S1) and the last factor described the Mg2+ block depending on local voltage u (Jahr and Stevens, 1990; Rhodes,
2006). Finally, synaptic currents were computed from
gsyn ðtÞ = wpre ðtÞ$ðgAMPA ðtÞ + gNMDA ðtÞÞ;

(32)


Isyn ðtÞ = gsyn ðtÞ$ VðtÞ  Esyn ;

(33)

where V was the local membrane voltage in the cell model and Esyn was the reversal potential, which we set to 0 mV.
Neuron Model
In our simulations, the plasticity rule was applied to a pyramidal neuron model developed by Hay et al. (2011). It represents a cortical
L5b pyramidal cell of the rat whose morphology was reconstructed in 3D via light microscopy. Electrophysiological properties
were acquired by current injection protocols combined with whole-cell recording techniques in vitro and subsequently reproduced
in the cell model using a multi-objective genetic algorithm. Of the four biophysical ion channel configurations provided by the authors,
we selected the fourth one (see their supplementary materials, Hay et al., 2011), where the APs are generated in the axon initial
segment.
Implementation
Simulations were run using the NEURON 7.4 environment (Carnevale and Hines, 2006) using a constant integration time step of
0.025 ms. Presynaptic signals were sent directly to the synapse without explicit modeling of a presynaptic cell. Fitting of the plasticity
model parameters was done via manual search in a two-stage process. First, all time constants, thresholds and saturation slopes
were adjusted so that all simulations would describe their respective experimental data qualitatively. Then, plasticity pathway amplitudes were fine-tuned for the three stimulation protocols used in experiments, aiming at quantitative matches wherever possible
(Table S1). To save simulation time of repetitive stimulation protocols, we ran one sweep at a time and approximated the final
outcome via
  
 



sweep
sweep
init
init
(34)
 winit
wfinal = winit
pre + wpre
pre $ n $ wpost + wpost  wpost $ n ;
where wsweep
and wsweep
pre
post were weight factors after one sweep and n was the total number of sweeps in the protocol.

Cell Reports 29, 4295–4307.e1–e6, December 24, 2019 e4

Stimulation protocols
Voltage Clamp
The membrane potential of the neuron model was clamped to voltages in the range of –75 to –15 mV (Figure 1E). Single presynaptic
events were directly sent to the plasticity rule (using parameter set 1; see Table S1) and ten of these sweeps were taken into account
to calculate final weight changes.
Pre- and Postsynaptic Bursts
The stimulation procedure was implemented according to previous experiments (Sjöström and Häusser, 2006; Sjöström et al., 2001)
(Figure 2A). Stimulation was performed by injecting step currents (5 ms at 2.7 nA) into the soma of the cell model to evoke bursts of
five axo-somatic APs. Simulations included one pre- and one postsynaptic burst at a time, shifted by either +10 ms or –10 ms. Ten of
these sweeps were considered for an intra-burst frequency of 0.1 Hz (representing 50 spikes in total) and 15 sweeps were used for
intra-burst frequencies of 10, 20, 40 and 50 Hz (representing 75 spikes in total), matching the experimental procedure (Sjöström et al.,
2001). To assess location differences, one proximal (90 mm from soma) and one distal (669 mm from soma) location along the apical
dendrite were chosen, at which the plasticity rule (using parameter set 1; see Table S1) was placed. Proximal locations mimicked
L5/L5 connections, while distal locations mimicked L2/3/L5 connections (Sjöström and Häusser, 2006). We set gmax = 3:5 nS
and sAMPA = sNMDA = 0:5.
Single Presynaptic Events and Postsynaptic Bursts
In this stimulation protocol, several combinations of timings, burst frequencies and numbers of spikes were tested (Nevian and Sakmann, 2006) (Figure 3). To match the experiments, we applied the plasticity rule (using parameter set 2; see Table S1) to a proximal
basal dendrite (55 mm from the soma). In our simulations, postsynaptic APs were evoked via somatic step current injection (5 ms at
2.1 nA). Single sweeps were simulated and 60 sweeps used for calculation of weights (corresponding to 60 low-frequency repetitions). Exact configurations of all the different spike patterns are visualized in Figure 3. We set gmax = 3:5 nS and sAMPA = sNMDA = 0:5.
Burst-Induced Dendritic Spikes
Following the corresponding experiments (Letzkus et al., 2006), one presynaptic event was either followed (+10 ms) or preceded
(–10 ms) by a burst of three postsynaptic APs at 200 Hz (Figure 4A). APs in our simulations were evoked via somatic step current
injection (2 ms at 5.5 nA) to reliably produce distal dendritic Ca2+ spikes (see voltage traces in Figures 4A and S3). Proximal and distal
locations for the plasticity rule (using parameter set 3; Table S1) in these simulations were 90 mm and 669 mm from the soma, respectively. Single sweeps were simulated and 100 sweeps used for calculations, representing 100 low-frequency repetitions. We set
gmax = 3:5 nS and sAMPA = sNMDA = 0:5.
Voltage Profiles and Plasticity Windows
To generate the spatiotemporal voltage plots (Figures 2C and 4C), we selected one termination point of an apical dendrite in the
model, calculated the exact path down to the soma and selected locations roughly every 30 mm (limited by the compartmental resolution of the cell model), leading to a set of 41 more or less evenly distributed locations along the path. We then applied the exact
current injection protocol of each given experiment and measured voltages at all selected locations over the entire duration. Voltage
curves did not change considerably with respect to distance for alternative paths (i.e., where another termination point was chosen).
For the spatiotemporal plasticity windows, we used the set of 41 locations along one specific path to apply the plasticity rule to and
simulated each stimulation protocol with different timings in an interval of [–50 50] ms at steps of 1 ms. Each data point of the image
thus corresponded to one of the resulting 4,141 single simulations.
Subthreshold Activation of Synapse Clusters
This stimulation protocol involved activation of a cluster of four synapses in rapid succession (0.1 ms interval), imitating two-photon
glutamate uncaging experiments (Weber et al., 2016). We repeated the protocol, each time placing the cluster (using parameter set 3;
see Table S1) at a different segment of the neuron model for all possible segments and then mapped plasticity outcomes onto the
morphology (Figure 5A). We used 50 sweeps, representing 50 low-frequency repetitions. In the channel block condition, we simply
set the conductance of all apical Na+ and Ca2+ channels of the neuron model to zero. Example locations (Figure 5A, boxes) were
chosen to be at 20% and 90% of the total branch length, respectively, in accordance with experiments (Weber et al., 2016). For
each of the four synapses, we set gmax = 2:5 nS, sAMPA = 0:8 and sNMDA = 0:2 to prevent excessive NMDA currents, as these were
not reported in the study (Weber et al., 2016).
Random Basal and Tuft Inputs
We randomly distributed non-plastic input synapses across parts of the dendritic tree, amounting to 50 basal and 300 tuft synapses
with an AMPAR-exclusive conductance of 2.5 nS each. In addition, we placed ten plastic synapses (using parameter set 3; see Table
S1) close to each of the two main spiking zones of the cell (basal: 32 mm from the soma; apical: 672 mm from the soma; Figure 6A). A
single sweep in the simulations consisted of a phase of synaptic activity (100 ms). During active phases, either basal, apical or all
synapses were randomly activated independently using a Poisson distribution at an average frequency of 10 Hz. For the weight distribution histograms (Figure 6D), we ran 100 simulations with different random seeds per condition, leading to a total of 1,000 plastic
synapses per location and condition. For each of the plastic synapses, we set gmax = 2:5 nS, sAMPA = 0:8 and sNMDA = 0:2.
Heterosynaptic Effects and NMDA Spikes
This stimulation protocol involved two clusters of eight synapses each, which we placed on a far distal apical tuft dendrite of the cell
model (Figure 7A). The distances to the soma were 950 and 1,077 mm, respectively. Synaptic activation patterns were available in two
modes, uniform and synchronized. Uniform activation was modeled using a Poisson distribution with an average frequency of 8 Hz.

e5 Cell Reports 29, 4295–4307.e1–e6, December 24, 2019

Synchronized activation was modeled using a sinusoid with a frequency of 8 Hz and amplitude of 0.05 oscillating around 0.005, where
positive function values gave probabilities of synapses being activated per 1 ms. The protocol consisted of one phase of
synaptic activity at both clusters with a duration of 350 ms. For each of the synapses, we set gmax = 2:5 nS and
sAMPA = sNMDA = 0:5.
Sensitivity Analysis
We performed sensitivity analysis by varying each single parameter in the model by four different factors for each set of amplitudes
(Figure S4). The model is relatively sensitive especially to changes in threshold parameters, which in most extreme cases can lead to
pathways being activated even at resting potential.
QUANTIFICATION AND STATISTICAL ANALYSIS
For experimental results reproduced by our model, original data is always given as mean ± SEM (see Figures 2B, 3, 4B, and 5B).
Simulation data presented in the histograms of Figure 6D was acquired using 100 different random seeds for generating Poissondistributed event sequences in 10 plastic synapses each, leading to n = 1,000 plastic synapses per location and condition. Identical
seeds were used between conditions.
DATA AND CODE AVAILABILITY
The accession number for the stimulation procedures and the plasticity model reported in this paper is ModelDB: 251493.

Cell Reports 29, 4295–4307.e1–e6, December 24, 2019 e6

