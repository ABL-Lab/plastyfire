Article

Heterosynaptic cross-talk of pre- and postsynaptic
strengths along segments of dendrites
Graphical Abstract

Authors
Rudi Tong, Thomas Edward Chater,
Nigel John Emptage, Yukiko Goda

Correspondence
rudi.tong@mcgill.ca (R.T.),
nigel.emptage@pharm.ox.ac.uk (N.J.E.),
yukiko.goda@riken.jp (Y.G.)

In Brief
Tong et al. show that eliciting plasticity at
groups of spatially clustered synapses
triggers changes at neighboring
synapses along the dendrite where spine
plasticity is accompanied also by
changes in presynaptic strength. The
polarity of heterosynaptic plasticity is
distance dependent, and CaMKII,
calcineurin, and nitric oxide (NO)
signaling play distinct roles.

Highlights
d

Potentiation of spatially clustered synapses induces
heterosynaptic plasticity

d

Heterosynaptic plasticity is expressed at both pre- and
postsynaptic terminals

d

The polarity of neighboring pre- and postsynaptic plasticity is
distance dependent

d

Multiple signaling pathways differentially contribute to
heterosynaptic plasticity

Tong et al., 2021, Cell Reports 34, 108693
January 26, 2021 ª 2021 The Authors.
https://doi.org/10.1016/j.celrep.2021.108693

ll

ll
OPEN ACCESS

Article

Heterosynaptic cross-talk
of pre- and postsynaptic strengths
along segments of dendrites
Rudi Tong,1,4,* Thomas Edward Chater,2 Nigel John Emptage,1,3,* and Yukiko Goda2,3,5,*
1University of Oxford, Department of Pharmacology, Mansfield Rd, Oxford OX1 3QT, UK
2RIKEN Center for Brain Science, 2-1 Hirosawa, Wako-shi, Saitama 351-0198, Japan
3Senior author
4Present address: Montreal Neurological Institute & Hospital, 3801 University Street, Montreal, QC H3A 2B4, Canada
5Lead contact

*Correspondence: rudi.tong@mcgill.ca (R.T.), nigel.emptage@pharm.ox.ac.uk (N.J.E.), yukiko.goda@riken.jp (Y.G.)
https://doi.org/10.1016/j.celrep.2021.108693

SUMMARY

Dendrites are crucial for integrating incoming synaptic information. Individual dendritic branches are thought
to constitute a signal processing unit, yet how neighboring synapses shape the boundaries of functional dendritic units is not well understood. Here, we address the cellular basis underlying the organization of the
strengths of neighboring Schaffer collateral-CA1 synapses by optical quantal analysis and spine size measurements. Inducing potentiation at clusters of spines produces NMDA-receptor-dependent heterosynaptic
plasticity. The direction of postsynaptic strength change shows distance dependency to the stimulated synapses where proximal synapses predominantly depress, whereas distal synapses potentiate; potentiation
and depression are regulated by CaMKII and calcineurin, respectively. In contrast, heterosynaptic presynaptic plasticity is confined to weakening of presynaptic strength of nearby synapses, which requires CaMKII
and the retrograde messenger nitric oxide. Our findings highlight the parallel engagement of multiple
signaling pathways, each with characteristic spatial dynamics in shaping the local pattern of synaptic
strengths.

INTRODUCTION
The maintenance of complex neural network dynamics imposes
constraints on synaptic strength regulation. Cellular principles
underlying neural network function, such as local dendritic integration or the minimization of metabolic cost, are thought to give
rise to an optimal range of synaptic strength distributions, which,
in turn, best support network function (Poirazi and Mel, 2001;
Winnubst and Lohmann, 2012; Govindarajan et al., 2006; Harris
et al., 2012). Indeed, synapses display several stereotypical
organizational features of their strengths, including an apparent
stable distribution across a synapse population that is skewed
toward low magnitudes and a non-uniform spatial distribution
(Buzsáki and Mizuseki, 2014; Barbour et al., 2007). How these
synaptic patterns emerge from elementary molecular mechanisms is however not fully understood.
Previous studies have reported functional and anatomical
similarities between neighboring synapses (Liang et al., 2018;
McBride et al., 2008; Branco et al., 2008; Lee et al., 2019; Winnubst and Lohmann, 2012; Fu et al., 2012; Letellier et al., 2019;
but see also Jia et al., 2010; Chen et al., 2011) that were mostly
confined to dendritic branches. For instance, pre- and postsynaptic strengths were found to be more similar along short
stretches of dendrites than across different branches of den-

drites (Branco et al., 2008; Letellier et al., 2019). Similarly, synapses in close proximity were more likely to be synchronously
active (Takahashi et al,. 2012; Winnubst et al., 2015; Kleindienst
et al., 2011) and functionally correlated (Lee et al., 2019; Scholl
et al., 2017). These synaptic patterns could be either explained
by hardwired developmental mechanisms (Druckmann et al.,
2014; Deguchi et al., 2011) or suggest additional activitydependent mechanisms that facilitate coordination of synaptic
strengths along short dendritic segments. In support of the
contribution of activity-dependent mechanisms, pharmacological inhibition of NMDA receptor (NMDAR) activity abolished
the spatial clustering of temporally correlated inputs (Winnubst
et al., 2015; Kleindienst et al., 2011), and local dendritic activity
was shown to adjust pre- and postsynaptic strengths (Branco
et al., 2008; Letellier et al., 2019). Furthermore, functional synaptic clusters were more pronounced after environmental
enrichment (Nikonenko et al., 2013) and prism rearing in barn
owls (McBride et al., 2008). Alternatively, these local synaptic
patterns may also arise as a consequence of local dendritic
spiking (Poirazi et al., 2003; Branco and Häusser, 2010; Mel,
1992) or local homeostatic constraints on resources that promote the normalization of synaptic strengths in a homeostatic
manner (Bourne and Harris, 2011; Sutton et al., 2006; Ju
et al., 2004).

Cell Reports 34, 108693, January 26, 2021 ª 2021 The Authors. 1
This is an open access article under the CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/).

ll
OPEN ACCESS

A

Article
B

C

D

E

Figure 1. Experimental setup to study pre- and postsynaptic heterosynaptic plasticity
(A) Schematic of the experimental setup (top), cl-LTP induction protocol (middle), and timeline (bottom).
(B) Example neuron filled with 1 mM OGB-1 and 0.5 mM AF594 for 1 min. The boxed area is shown in (E). SE, stimulation electrode. Scale bar: 20 mm.

(legend continued on next page)

2 Cell Reports 34, 108693, January 26, 2021

ll
OPEN ACCESS

Article
Collectively, these observations suggest that neighboring synapses interact with each other and that synaptic plasticity plays a
role in implementing this coordination and possibly also competition between synapses. In agreement with such a view of intersynaptic interactions, a large body of evidence shows that activity-dependent changes of synaptic strength occur at both active
(homosynaptic plasticity) and inactive synapses (heterosynaptic
plasticity) (Royer and Paré, 2003; Schuman and Madison, 1994;
Scanziani et al., 1996; Lynch et al., 1977). For instance, changes
in synaptic strength at active synapses can either spread to adjacent synapses (Engert and Bonhoeffer, 1997; Wiegert and Oertner, 2013; Hayama et al., 2013; Letellier et al., 2019), induce a
compensatory response at adjacent synapses (Bian et al.,
2015; Oh et al., 2015; El-Boustani et al., 2018; Lee et al., 2013),
or affect metaplastic properties to bias future modifications (Harvey and Svoboda, 2007; Harvey et al., 2008; Fonseca et al.,
2004; Govindarajan et al., 2011). A better understanding of
when, where, and how neighboring synapses interact will help
provide insights into the emergence of higher level organizational
principles of synapses, such as the minimal computational unit,
and how its bound is set (Govindarajan et al., 2006).
Although a variety of phenomena and mechanisms have been
reported for synaptic interactions, a common framework to
explain the spatial coordination of synaptic strengths has been
difficult to formulate (Kruijssen and Wierenga, 2019). For
instance, it has been suggested that low-level organizational
patterns might emerge from the unique spatiotemporal properties of signaling molecules that interact during synaptic plasticity
(Nishiyama and Yasuda, 2015). However, the interplay of parallel
plasticity pathways during synaptic plasticity has been difficult to
observe and dissociate experimentally. In addition, even though
it is generally accepted that synaptic plasticity occurs at both
pre- and postsynaptic loci (Bliss and Collingridge, 1993) and heterosynaptic signaling between postsynaptic terminals along
dendrites has been explored (Chater and Goda, 2020), little to
nothing is known about heterosynaptic plasticity at the presynaptic terminal. This knowledge gap can mainly be attributed to
the technical difficulties of measuring presynaptic strength,
especially in more intact systems, such as brain slices, or in vivo
(O’Rourke et al., 2012; Burette et al., 2015). However, investigating synaptic plasticity at both loci is particularly important
because pre- and postsynaptic terminals might represent functionally distinct and independently regulated compartments (Padamsey et al., 2017; Costa et al., 2017; Xu et al., 2013; Sáez and
Friedlander, 2009; Letellier et al., 2019) and are therefore likely
to also differ in their spatial regulation with consequences on
the property of information transfer (Costa et al., 2017; Tong
et al., 2020).

Here, we have addressed the cellular basis underlying the organization of neighboring synaptic strengths at the Schaffer
collateral inputs to the dendrites of CA1 pyramidal neurons.
Our findings reveal a spread of plasticity along dendrite segments with a pattern distinct between the pre- and the postsynaptic sides, and our observations point to the engagement of
independent signaling pathways that operate in parallel.
RESULTS
Local structural potentiation of groups of spines
induces bi-directional postsynaptic heterosynaptic
plasticity
We first investigated the synaptic cross-talk that might
contribute to the formation of stable functional clusters. We hypothesized that the local coordination of synaptic strength is implemented by a tight balance of cooperation and competition
between synapses. This predicts that experimentally imposing
strong coordinated activity locally at a subset of synapses will
induce a reorganization of synaptic strengths at all nearby synapses due to underlying competitive processes. We therefore
emulated the activation of groups of spatially clustered synapses
along segments of dendrites by using two-photon glutamate uncaging. Quasi-synchronous stimulation of groups of synapses
was paired with cell-wide depolarization of the postsynaptic
neuron to promote dendritic signaling events (Figure 1A). Postsynaptic strength was monitored by estimating the changes in
spine volume, a parameter known to tightly correlate with the
number of functional postsynaptic a-amino-3-hydroxy-5methyl-4-isoxazolepropionic acid (AMPA) receptors (Matsuzaki
et al., 2001, 2004; Smith et al., 2003; Béı̈que et al., 2006). As a
readout of presynaptic strength, we used optical quantal analysis to determine release probability (Pr) at selected synapses
(Padamsey et al., 2019). Briefly, CA1 pyramidal neurons in organotypic hippocampal slices were loaded (45–60 s) by the patch
pipette with the Ca2+ indicator OGB-1 (1 mM) and a fluorescent
dye (AF594; 0.5 mM) to visualize fine dendrites. Dendritic segments of interest were locally stimulated to measure Pr of select
synapses (Figure 1; Pr was estimated as the fraction of trials that
resulted in successful release to the first pulse of pairs of stimulation pulses, as monitored by an excitatory postsynaptic Ca2+
transient [EPSCaT]; see STAR methods). For glutamate uncaging, 4-methoxy-7-nitroindolinyl (MNI)-caged glutamate (10 mM)
was locally perfused, and the uncaging laser power was
adjusted to match the magnitude and temporal dynamics of
Ca2+ responses evoked by electrical stimulation (Figure 1D).
To elicit long-term potentiation (LTP), the same neuron was repatched, and five to seven adjacent dendritic spines were

(C) Example of Pr measurements by optical quantal analysis before (left) and after (right) cl-LTP induction. Raw fluorescence recording traces (top) sorted into two
groups showing either successful (red) or unsuccessful release (gray) to the first action potential. The plots below show the peak intensity of the response to the
first (middle) and second pulse (bottom). The thick gray bar corresponds to two standard deviations of the baseline used as detection criteria for automated
analysis. Pr was estimated as the fraction of trials that triggered successful release (top right corner). Facilitation to the second pulse of the pairwise stimuli is
evident as an increase in Pr.
(D) Glutamate photolysis was adjusted to elicit Ca2+ transients similar in strength compared with EPSCaTs (laser power varied from 2.4–3.5 mW measured at the
sample). The gray trace corresponds to the nearby dendrite.
(E) cl-LTP was induced by quasi-synchronous photolysis of glutamate at five to seven spines paired with postsynaptic depolarization, 303 at 0.5 Hz. An increase
in spine size can be seen immediately after photolysis. For optical quantal analysis, line scans were taken through the spine of interest (s) and the dendritic branch
(d). See STAR methods for details. Scale bar: 10 mm.

Cell Reports 34, 108693, January 26, 2021 3

ll
OPEN ACCESS

Article

A

B

Ci

Cii

Ciii

Di

Dii

Diii

Figure 2. Induction of structural plasticity of groups of spines leads to heterosynaptic bi-directional changes in spine size
(A) Example of cl-LTP induction at six spines by glutamate photolysis. Spine structural changes were assessed for 30–45 min. Scale bar: 5 mm.
(B) The distance metric used for analysis of the spine cluster was defined as the distance along the dendritic branch to the center of mass of the stimulated cluster
(Dx).
(Ci) Spine size changes of stimulated and unstimulated spines. cl-LTP was induced at t = 0 min. Stimulated spines showed robust structural LTP, whereas
unstimulated spines enlarged transiently (n = 17). (Cii) Spine size changes of unstimulated spines, grouped by distance (see also Figure S4). Spines located close
to the stimulated cluster decreased in size. More distal spines, however, increased in size. (Ciii) Persistent spine size change at 30–45 min depends on the
distance to the stimulated cluster. Spines belonging to within 3 mm and 5 to 7 mm groups were significantly different from control experiments summarized in (Diii)
(p < 0.05, Mann-Whitney U test). The polarity change associated with distance dependency of the mean spine size change is highlighted in red.
(D) Same as (C), but postsynaptic depolarization was prevented by holding the neuron at 70 mV during cl-LTP induction (n = 9). Neither stimulated nor unstimulated spines showed substantial structural changes. Data points are shown as mean ± SEM. See also Figures S1 and S2.

targeted for glutamate photolysis (303; 0.5 Hz) paired with depolarization of the patch-clamped cell to 0 mV (Figure 1E). LTP was
induced within 10 s of establishing the whole-cell patch mode,
and the neuron was immediately patched off after induction.
The re-patch led to a slight drop in the concentration of the fluorescent dye delivered during the first patch due to intracellular
dialysis, which could confound the estimation of spine volume.
We therefore normalized all spine fluorescence measurements
to the average fluorescence intensity of the underlying dendritic
branch (see STAR methods), which was assumed to remain
largely structurally stable (Nimchinsky et al., 2004). We re-examined presynaptic strength after 30 min by repeating optical
quantal analysis; spine structure was monitored for up to

4 Cell Reports 34, 108693, January 26, 2021

45 min post-induction. We term LTP induced by this form of clustered stimulation as ‘‘cl-LTP.’’
Clustered stimulation led to structural enlargement of spine
head volumes of a large fraction of stimulated spines, which
was stable for at least 45 min (cl-LTP; Figure 2Ci, 59.5% ±
4.1% of spines with >20% increase in spine size at 30–45 min;
see Figures S1 and S2 for spine Ca2+ changes induced by stimulation). The average increase in spine size after 30–45 min was
DV30-45min = 38.7% ± 4.5% of baseline (n = 17 experiments).
Neighboring unstimulated spines also transiently increased in
size shortly after cl-LTP induction (Figure 2Ci, DV10min = 10.8%
± 4.0%) before slowly decaying back to baseline by 45 min
(DV30-45min = 1.1% ± 2.2%). When the postsynaptic cell was

ll
Article
held at 70 mV during clustered stimulation, stimulated spines
did not increase in size, neither transiently nor over the long
term (Figure 2Di; DV30-45min = 4.9% ± 5.9%; n = 9). Without
postsynaptic depolarization, neighboring unstimulated spines
also remained unchanged (Figure 2Di; DV10min = 2.8% ±
4.2%; DV30-45min = 3.1% ± 2.9%). Local application of AP5
(500 mM; locally perfused) also abolished both homo- and heterosynaptic spine structural changes (Figure S3Ai).
Spatial proximity plays an important role in heterosynaptic
plasticity because spatially delimited diffusion of signaling molecules between synapses is thought to mediate the local coordination of synaptic strengths (Nishiyama and Yasuda, 2015). We
therefore grouped unstimulated spines based on their distance
along the dendrite to the center of mass of the stimulated cluster
(Figure 2B; see STAR methods). Unstimulated spines within 7 mm
of the stimulated cluster showed significant bi-directional structural changes (Figure S4). The distance dependence of the bidirectional changes was most apparent when we compared
spines within and spines beyond 3.8 mm of the stimulated cluster
(Figure S4). We found that spines located within 3.8 mm showed a
consistent decrease in size, which lasted for at least 45 min (Figure 2Cii; DV30-45min = 23.0% ± 5.3%). In contrast, spines at a
distance of 3.8–7 mm exhibited weak but significant potentiation,
which lasted for the duration of the experiment (DV30-45min =
18.1% ± 4.0%). Spines located further away showed a slow
but transient increase in size, reaching a maximum at 10 min after stimulation and decaying back to baseline by the end of the
experiment (DV30-45min = 2.3% ± 2.9%). Figure 2Ciii summarizes the distance-dependent bi-directional regulation of postsynaptic strength, which is absent without postsynaptic depolarization (Figure 2Diii; cl-LTP versus no depolarization, p < 0.05 for
spines located <3 mm and 5–7 mm; Mann-Whitney U test) and
when NMDAR signaling was blocked with AP5 (Figure S3).
Local structural potentiation of groups of spines
induces presynaptic heterosynaptic weakening
Next, we examined presynaptic changes following cl-LTP. Due
to the comparatively low throughput nature of optical quantal
analysis in which Pr measurement in general was feasible for
only a single synapse per experiment, we focused mainly on synapses located close to the stimulated cluster (synapse located
<4 mm to the stimulated cluster in 13/17 experiments). cl-LTP induction resulted in a decrease in Pr of such proximal synapses
(Figure 3C; DPr = 0.27 ± 0.05, n = 13, <4 mm), whereas Pr of
synapses located more distally did not change (Figure 3C;
DPr = 0.02 ± 0.06, n = 4, >4 mm). Given that the activation of
presynaptic NMDARs can lead to the expression of presynaptic
long-term depression (LTD; Padamsey et al., 2017), the large
amounts of glutamate released during photolysis could have resulted in local accumulation and spillover of glutamate to cause a
decrease in Pr at neighboring unstimulated synapses. However,
holding the neuron at 70 mV during cl-LTP induction abolished
presynaptic changes (Figure 3D; no depolarization: DPr = 0.02
± 0.06, n = 11, p < 0.01, Kruskal-Wallis H-test with post hoc
Dunn’s test). Thus, glutamate uncaging by itself was not sufficient to induce heterosynaptic LTD (hetLTD) of Pr, and concurrent postsynaptic depolarization was required. Notably, postsynaptic depolarization alone is known to influence Pr (Branco

OPEN ACCESS

et al., 2008; Volgushev et al., 1997). Therefore, to control for
this possibility, we repeated the experiment without glutamate
photolysis by either locally perfusing artificial cerebrospinal fluid
(aCSF) without caged glutamate or by omitting the photolysis
pulse while applying postsynaptic depolarization. Under these
conditions, presynaptic weakening was not observed (Figure 3D;
no photolysis: DPr = 0.01 ± 0.03, n = 18, p < 0.01, KruskalWallis H-test with post hoc Dunn’s test). Presynaptic hetLTD
was also abolished by the local application of AP5 during the
induction of cl-LTP (Figure S3C). Thus, presynaptic hetLTD requires both postsynaptic depolarization and synaptic stimulation
of nearby synapses.
The decrease in Pr at synapses neighboring the stimulated
spines was correlated with the initial Pr, which we did not
observe for control conditions (Figure 3F). This dependence of
the extent hetLTD of Pr on the initial Pr, which in itself was highly
heterogeneous (mean initial Pr = 0.58, SD = 0.21), could explain
the variance of the data for the change in Pr as a function of distance. Specifically, the magnitude of presynaptic changes was
not strongly related to the distance to the stimulated cluster (Figure 3G). Instead, the observation of presynaptic hetLTD was
limited to synapses located within ~4 mm of the stimulated cluster, unlike the magnitude of spine size changes that were correlated with the distance to the stimulated cluster for up to ~7 mm.
Moreover, despite exhibiting some trend, the magnitudes of preand postsynaptic changes were not significantly correlated
(Pearson r = 0.487, n = 13, p = 0.09; Figure 3H).
Altogether, strong stimulation of groups of synapses led to a
restructuring of synaptic strengths along the dendrite. The postsynaptic strength of unstimulated synapses changed bi-directionally, which was dependent on the distance to the stimulated
synapses. In contrast, the presynaptic terminal weakened only at
synapses close to the stimulated cluster. The different patterns
of heterosynaptic pre- and postsynaptic changes and the lack
of correlation between pre- and postsynaptic changes at individual synapses suggest the parallel engagement of independent
signaling pathways.
Presynaptic hetLTD requires CaMKII activation,
whereas postsynaptic changes are driven by parallel
signaling pathways involving calcineurin and CaMKII
Next, we sought to further confirm the suggested involvement of
distinct signaling pathways by gaining insights into the underlying
molecular pathways involved in the heterosynaptic cross-talk. The
spatial reorganization of synaptic strengths might reflect the
biochemical computations driven by the interaction of signaling
pathways (Nishiyama and Yasuda, 2015). Strong evidence suggests that heterosynaptic postsynaptic weakening requires the
activation and diffusion of the phosphatase calcineurin (Oh
et al., 2015; Hayama et al., 2013). In contrast, the activation of
the kinase CaMKII and its downstream targets including the small
GTPases was reported to facilitate LTP at nearby synapses (Harvey and Svoboda, 2007; Harvey et al., 2008; Hedrick et al., 2016).
Therefore, calcineurin and CaMKII are potential candidates for
plasticity factors mediating the depression and potentiation,
respectively. Our observation of a bi-directional regulation of postsynaptic strength could be explained by the orchestrated spatial
diffusion of these two plasticity factors, with CaMKII causing

Cell Reports 34, 108693, January 26, 2021 5

ll
OPEN ACCESS

Article

A

B

C

D

E

F

G

H

Figure 3. Induction of structural plasticity of groups of spines leads to heterosynaptic presynaptic weakening
(A) Pr was measured using optical quantal analysis. Top: Line scans were taken through the spine of interest (red arrow, s) and the dendrite (d) while the synapse
was stimulated. Scale bar: 2 mm. Bottom: Glutamate release was detected as Ca2+ transients in the spine head (ESPCaTs).
(B) Top: Raw Ca2+ transients show all-or-none response profile to repeated stimulation, before (left) and 30 min after (right) the induction of cl-LTP. Vertical gray
lines indicate time of stimulation. Responses to the first stimulation pulse were used to estimate Pr. Bottom: Peak amplitude of the EPSCaT of the first pulse is
shown as function of trial number (inter-trial interval = 15 s). The gray bar represents two standard deviations from the baseline. An increase in release failures is
observed 30 min after cl-LTP induction (right panel).
(C) Pr 30 min after cl-LTP induction (post) versus baseline Pr (pre). Pr decreased for synapses located close to the stimulated cluster (<4 mm; unstimulated
proximal) but not in distal synapses (>4 mm; unstimulated distal).
(D) No Pr change was observed when omitting glutamate photolysis (no photolysis, triangles) or postsynaptic depolarization (no depolarization, circles).
(E) Summary of Pr changes (cl-LTP unstim, n = 13, all synapses within 4 mm of the stimulated cluster; CTRL no depol., n = 11; CTRL no photolysis, n = 18).
Significance is shown for comparisons with label cl-LTP unstim (p < 0.01, Kruskal-Wallis H-test with post hoc Dunn’s test).
(F) The change in Pr as function of the initial Pr.
(G) Pr changes relative to distance to the stimulated cluster.
(H) Correlation between pre- and postsynaptic heterosynaptic changes.

6 Cell Reports 34, 108693, January 26, 2021

ll
Article
potentiation and calcineurin causing depression. We thus predicted that (1) selective inhibition of calcineurin should abolish
the postsynaptic weakening of proximal synapses without
affecting the potentiation at more distal synapses and (2) inhibition
of CaMKII should lead to the opposite effect. To test these predictions, we repeated the experiment under the pharmacological
blockade of either calcineurin or CaMKII.
The vehicle control recapitulated our previous results (Figure 2)
with the exception of a slight shift in the spatial bounds that
best explained the bi-directional heterosynaptic changes (Figure 4A; < 2.8 mm, DV30-45min = 15.5% ± 4.7%; 2.8–6.5 mm,
DV30-45min = 11.4% ± 5.7%; >6.5 mm, 5.9% ± 3.1%, n = 19; see
also Figure S4), which could have been attributed to culture variability. In the presence of bath-applied 10 mM KN62, an inhibitor
of CaMKII, the cl-LTP of stimulated spines was blocked (Figure 4B;
DV30-45min = 6.4% ± 9.2%, n = 10, p < 0.01), consistent with previous reports (Matsuzaki et al., 2004). In the presence of KN62, unstimulated spines located close to the stimulated
spines decreased in size, as observed above. However,
similarly to the stimulated spines, the enlargement of spines
located more distally was completely abolished (Figure 4Bii;
<2.8 mm, DV30-45min = 10.6% ± 8.2%; 2.8–6.5 mm, DV30-45min =
6.9% ± 5.6%; >6.5 mm, DV30-45min = 3.1% ± 3.2%). In contrast,
inhibition of calcineurin by bath application of 2 mM FK506
completely abolished heterosynaptic postsynaptic weakening of
proximal spines (Figure 4C). Instead, upon blocking calcineurin activity, we observed potentiation of both proximal and distal spines
(<2.8 mm, DV30-45min = 24.7% ± 8.8%; 2.8–6.5 mm, DV30-45min =
17.7% ± 8.3%; >6.5 mm, DV30-45min = 6.6% ± 4.9%, n = 8).
FK506 did not affect the magnitude of cl-LTP at stimulated spines
(DV30-45min = 40.4% ± 4.3%). Data summary is shown in Figure 4D.
Inhibition of CaMKII but not calcineurin also abolished presynaptic
heterosynaptic weakening (Figure 4E).
Altogether, these experiments demonstrate a clear pharmacological dissociation of heterosynaptic postsynaptic potentiation
and depression. Consistent with our hypothesis, our observations
suggest that postsynaptic heterosynaptic plasticity is regulated
by at least two opposing factors, which depend on calcineurin
and CaMKII activity. Moreover, selective block of presynaptic
but not postsynaptic weakening by KN62 also points to a dissociation between pre- and postsynaptic plasticity. Notably, upon inhibiting calcineurin, presynaptic weakening was accompanied by
postsynaptic potentiation. Presynaptic hetLTD is therefore not
likely caused by a matching between pre- and postsynaptic
strengths (see e.g., Tang et al., 2016; Meyer et al., 2014; Bayazitov
et al., 2007; Kay et al., 2011; Tokuoka and Goda, 2008) but is subject to regulation by an independent signaling pathway that acts in
parallel to the postsynaptic regulatory mechanism.
NO is the retrograde messenger for presynaptic hetLTD
The heterosynaptic presynaptic weakening in our experimental
paradigm necessitated some form of retrograde signaling, as
the effects of stimulation targeted to the postsynaptic spines
were prevented in the absence of depolarization of the postsynaptic neuron. The findings so far suggested that the retrograde
messenger was likely produced downstream of NMDAR activation, involved CaMKII signaling either upstream or downstream,
and was expected to rapidly diffuse over short distances from

OPEN ACCESS

the stimulated synapse cluster. Any one of a range of known
candidate molecules fulfills these requirements. In particular,
strong local dendritic depolarization can trigger the production
and release of nitric oxide (NO) (Padamsey et al., 2017), a highly
diffusible signaling molecule. Furthermore, NO signaling is
known to interact with NMDAR and CaMKII signaling (Hardingham et al., 2013). NO is also implicated in (homosynaptic) presynaptic LTP (Padamsey et al., 2017) and LTD (Stanton et al.,
2003). We therefore hypothesized that strong local activity triggers the production and release of NO, which, due to its gaseous
nature, diffuses over a confined distance to neighboring presynaptic terminals to elicit presynaptic weakening. We tested this
hypothesis by applying the NO synthase inhibitor L-NAME.
Bath application of 100 mM L-NAME completely abolished
presynaptic hetLTD (Figures 5B and 5D; DPr = 0.02 ± 0.06,
n = 8, p < 0.01). In contrast, L-NAME did not impair cl-LTP induction and expression at stimulated spines (DV30-45min = 53.6% ±
9.5%, n = 8), the transient increase in spine size immediately after induction (DV10min = 11.7% ± 1.4%), or the distance-dependent bi-directional heterosynaptic regulation of spine size
(Figure 5A). In fact, the fraction of stimulated spines
expressing stable structural LTP was significantly larger when
NO production was inhibited (77.4% ± 4.5% of spines with
DV30-45min >20%, p < 0.05 versus vehicle).
We next sought to further confirm the suggested role of NO
in a more physiological experimental paradigm. To do so, we
devised an experimental scheme by using acute brain slices
that closely followed our optical setup in slice cultures. We assessed the plasticity of Schaffer-collateral-CA1 synapses within
relatively close spatial proximity of stimulated synapses and yet
at a distance at which basal responses remained independent
(Figure 6; see also STAR Methods). Along individual stretches
of dendrites, which were identified by an intracellular fluorescent
dye (AF488) introduced to the patched CA1 neuron, we positioned two stimulation electrodes (SEs) 10–40 mm apart to increase the probability of stimulating sets of synapses located
close to each other. We used a cross-facilitation test to ensure
the stimulation of independent Schaffer-collateral pathways.
We induced LTP in one pathway by pairing presynaptic stimulation with strong postsynaptic depolarization by current injection,
which elicited bursts of action potentials 60 times at 5 Hz (Padamsey et al., 2017). We then examined changes in the unstimulated pathway.
Induction of LTP led to an increase in the excitatory postsynaptic potential (EPSP) slope of the stimulated pathway (Figure 6B; DEPSP = 94.9% ± 13.1%, n = 7) and a decrease in the
paired-pulse ratio (PPR) of EPSP slopes (Figure 6C;
DPPR = 0.203 ± 0.032), which was stable for at least 30 min,
consistent with previous findings (Padamsey et al., 2017). This
result supported a likely presynaptic expression locus.
Conversely, the unstimulated pathway showed mild depression
of the EPSP slope with an increase in PPR (DEPSP = 19.1% ±
15.4%, DPPR = 0.307 ± 0.092). The heterosynaptic PPR
changes were strongly correlated with the fractional change of
the EPSP slope (Figure S5; r = 0.842, p < 0.05), which suggested that electrically induced heterosynaptic plasticity was expressed as a decrease in Pr, a finding that was in agreement
with heterosynaptic presynaptic weakening that accompanied

Cell Reports 34, 108693, January 26, 2021 7

ll
OPEN ACCESS

Article
vehicle

Ai

Aii

< 2.8 μm
2.8-6.5 μm
> 6.5 μm

∆Spine size
(% baseline)

∆Spine size
(% baseline)

∆Spine size
(% baseline)

Stimulated
Unstimulated

t (min)

Distance from
LTP cluster (μm)

t (min)

KN62

Bi

Bii

< 2.8 μm
2.8-6.5 μm
> 6.5 μm

Biii

∆Spine size
(% baseline)

∆Spine size
(% baseline)

∆Spine size
(% baseline)

Stimulated
Unstimulated

t (min)

Distance from
LTP cluster (μm)

t (min)

FK-506

Ci

Aiii

Cii

Ciii

< 2.8 μm
2.8-6.5 μm
> 6.5 μm

t (min)

< 2.8 μm
n.s.

2.8-6.5 μm

**

Distance from
LTP cluster (μm)

t (min)

*

n.s.

E

> 6.5 μm
n.s.

n.s.

**

n.s.

KN62

FK-506

∆Pr

∆Spine size (% baseline)

D

∆Spine size
(% baseline)

∆Spine size
(% baseline)

∆Spine size
(% baseline)

Stimulated
Unstimulated

vehicle KN62 FK-506

vehicle KN62 FK-506

vehicle KN62 FK-506

vehicle

Figure 4. CaMKII is required for postsynaptic heterosynaptic LTP, whereas calcineurin is required for postsynaptic heterosynaptic LTD
(Ai) Spine size changes of stimulated (orange) and unstimulated (gray) spines in vehicle control group. Spine enlargement was robust and persistent at stimulated
spines (DV30-45min = 42.5% ± 5.1%) and transient at unstimulated spines (DV10min = 14.6% ± 3.3%) (n = 19). (Aii) Spine size changes of unstimulated spines,
grouped by distance. (Aiii) Spine size changes of unstimulated spine, grouped in 2-mm bins.
(Bi–Biii) Effects of bath-applied CaMKII inhibitor KN62 (10 mM) on spine size changes. The expression of homo- and heterosynaptic LTP but not LTD was
abolished (n = 7).
(Ci–Ciii) Effects of bath-applied calcineurin inhibitor FK-506 (2 mM). Heterosynaptic shrinkage of proximal spines was completely abolished, and spine
enlargement was observed instead (n = 8). Homosynaptic LTP remained unchanged (DV30-45min = 40.4% ± 4.3%, n = 8). Error bars represent SEM.
(D) Comparison of heterosynaptic spine size changes for proximal (<2.8 mm), distal (2.8–6.5 mm), and more distant spines (>6.5 mm). Statistical comparison
includes L-NAME experiments shown in Figure 5 for the multiple comparison.
(E) Comparison of heterosynaptic Pr changes for spines located within 4 mm of the stimulated spine cluster (Kruskal-Wallis H-test with post hoc Dunn’s test).
Significance indicates comparison with vehicle control group. Data points are shown as mean ± SEM. See also Figures S3 and S4.

8 Cell Reports 34, 108693, January 26, 2021

ll
OPEN ACCESS

Article
Ai

L-NAME

Aii

< 2.8 μm
2.8-6.5 μm
> 6.5 μm

∆Spine size
(% baseline)

∆Spine size
(% baseline)

t (min)

Distance from
LTP cluster (μm)

t (min)

C

< 2.8 μm

2.8-6.5 μm

n.s.

n.s.

Pr pre

D

**

∆Pr

∆Spine size
(% baseline)

L-NAME

Pr post

B

Aiii

∆Spine size
(% baseline)

Stimulated
Unstimulated

vehicle L-NAME

vehicle

L-NAME

vehicle

L-NAME

Figure 5. NO is the retrograde messenger for heterosynaptic presynaptic weakening
Spine size changes in the presence of bath-applied NO synthase inhibitor L-NAME (100 mM). (Ai) Spine size changes of stimulated and unstimulated spines (n = 8).
(Aii) Spine size changes of unstimulated spines, grouped by distance (<2.8 mm: DV30-45min = 14.3% ± 4.7%; 2.8–6.5 mm: DV30-45min = 20.3% ± 6.9%; >6.5 mm:
11.7% ± 4.0%). (Aiii) Spine size changes of unstimulated spines, grouped in 2-mm bins. Inhibition of NO production by L-NAME did not affect spine size changes
of either stimulated or unstimulated spines.
(B) Pr before and 30 min after the induction cl-LTP.
(C) Summary of spine size changes.
(D) Summary of presynaptic changes (n = 8, p < 0.01). Statistical significance was tested using Kruskal-Wallis H-test and post hoc Dunn’s test, and comparisons
were made with the vehicle control group. Data points are shown as mean ± SEM.

cl-LTP in slice cultures. In contrast, PPR changes of the homosynaptic pathway were not correlated with the EPSP slope change
(r = 0.373, not significant), which suggested the occurrence of a
postsynaptic change that was concomitant to the expression of
presynaptic plasticity. Next, we inhibited NO synthase by bath
application of 100 mM L-NAME. As observed previously (Padamsey et al., 2017), this abolished the homosynaptic presynaptic
LTP (Figures 6D and 6E; DPPR = 0.029 ± 0.052, n = 7, p <
0.05). Furthermore, the heterosynaptic increase in PPR was
also blocked and a stable increase of the EPSP slope was unmasked (DEPSP = 56.0% ± 12.3%, DPPR = 0.020 ± 0.036, p <
0.05). Altogether, our data suggest that NO signaling orchestrates both potentiation and depression of presynaptic strength
in a context-dependent manner, that is, potentiation of release at
stimulated synapses and depression of release at neighboring
non-stimulated synapses.
DISCUSSION
Studies of heterosynaptic plasticity to date have mostly focused
on local structural plasticity of spines along dendritic branches
(Oh et al., 2015; Hayama et al., 2013; Wiegert and Oertner,
2013; Govindarajan et al., 2011; El-Boustani et al., 2018; Bian
et al., 2015; ; Harvey and Svoboda, 2007; but see also Letellier
et al., 2019). Both compensatory and non-compensatory forms

of heterosynaptic changes relative to the direction of plasticity
induced at the stimulated spines have been reported, whereas
the underlying basis for when, where, and what form of heterosynaptic plasticity is expressed have remained unclear. Here, we
found that heterosynaptic spine changes that likely represented
both compensatory and non-compensatory processes could be
induced at individual synapses by the same cl-LTP stimulus
although with different spatial profiles along the dendritic
segment. HetLTD was expressed at spines close to the stimulated synapses, whereas heterosynaptic LTP was observed at
spines further away, which altogether resulted in a spatial profile
that resembled a Mexican hat wavelet. What might constitute the
emergence of this spatial profile? In light of previous findings, an
attractive explanation for our observations could involve the differential properties of signaling pathways triggered by calcineurin and CaMKII and how they influence the neighboring
spines. Calcineurin, which promotes structural depression of
spines (Oh et al., 2015), is itself reported to diffuse into the
dendrite and to reach neighboring spines over short distances
(Fujii et al., 2013). Furthermore, small GTPases such as H-Ras,
Rac, or Rho, which are downstream effectors of the potentiating
factor CaMKII, also diffuse along the dendrite and can reach
substantial distances (Harvey and Svoboda, 2007; Harvey
et al., 2008; Hedrick et al., 2016) to mediate heterosynaptic
metaplastic changes that manifest as a facilitation to subsequent

Cell Reports 34, 108693, January 26, 2021 9

ll
OPEN ACCESS

Article

A

B

C

D

E

F

G

Figure 6. Heterosynaptic LTD in acute hippocampal slices is expressed presynaptically and requires NO signaling
(A) Experimental scheme. CA1 pyramidal neurons in acute hippocampal slices were filled by patching with an electrode containing 200 mM AF488 for 5 min. Two
tungsten stimulation electrodes (SE1/SE2) were placed along a segment of dendrite, 10–40 mm apart. Spike-timing-dependent LTP was induced in one pathway
by pairing presynaptic stimulation with postsynaptic current injection that elicited complex spikes (603 at 5 Hz). Scale bar: 10 mm.

(legend continued on next page)

10 Cell Reports 34, 108693, January 26, 2021

ll
Article
induction of potentiation (Harvey et al., 2008; Hedrick et al.,
2016). Therefore, in one possibility, differential diffusion along
the dendrite of calcineurin and active small GTPases could in effect produce the Mexican hat wavelet profile of spine structural
plasticity that we have observed. Such a scenario stipulates
that the effective diffusion of calcineurin must be weaker than
that of CaMKII-dependent signaling events, and this possibility
remains to be tested.
The overall spatial pattern of heterosynaptic spine plasticity
that we found was mostly consistent with previous reports of
heterosynaptic plasticity apart from some notable differences,
such as the presence of the Mexican hat profile discussed
above. The heterosynaptic spine shrinkage that we observed
was in accord with the study by Oh et al. (2015), in which the
sequential potentiation of spines in close spatial proximity produced shrinkage of unstimulated spines within ~3.4 mm of the
stimulated cluster in a process requiring the activation of calcineurin but not CaMKII. However, unlike our study, Oh et al.
(2015) did not report an enlargement of spines located further
away. Furthermore, several studies of single spine LTP showed
heterosynaptic facilitation to subsequent induction of LTP (Harvey and Svoboda, 2007; Harvey et al., 2008; Hedrick et al.,
2016). These studies, however, did not observe a shrinkage
of nearby spines, which is in contrast to our study and the
study by Oh et al. (2015). What might be the cause of these differences in heterosynaptic spine changes accompanying LTP
induction at target spines? Previously, CaMKII activation was
shown to display a highly supralinear response to stimulation
frequency and number (Fujii et al., 2013). In contrast, calcineurin activation was sublinear and required a longer integration time window (Fujii et al., 2013) from multiple sources (Oh
et al., 2015), suggesting that calcineurin activation might be
strongly thresholded (Oh et al., 2015). Based on these findings,
the qualitative differences in heterosynaptic spine dynamics
might be explained not only by the possible contribution of
the differential diffusion properties of signals associated with
activation of calcineurin and CaMKII at stimulated spines (see
above) but also by the differential activation of calcineurin
and CaMKII at stimulated and neighboring spines. In other
words, different stimulation used to trigger LTP would produce
distinct spatiotemporal pattern of calcium increases along the
dendritic segment (Figures S1 and S2), and this could in turn
promote distinct activation patterns of calcineurin and CaMKII
along the dendrite to bias the plasticity outcome at individual
spines. Accordingly, the heterosynaptic spine shrinkage in the
absence of potentiation of distal spines was triggered with

OPEN ACCESS

sequential low-frequency (1 Hz) glutamate uncaging of individual spines of groups of targeted synapses (Oh et al., 2015).
Such a paradigm involving a slow integration time window between stimulation would favor activation of calcineurin over that
of CaMKII (Fujii et al., 2013), and hetLTD would dominate. In
contrast, stimulation to elicit LTP that is confined to a single
synapse might not reach the activation threshold of calcineurin
while favoring the supralinear activation of CaMKII; without
hetLTD, heterosynaptic LTP or metaplasticity for subsequent
potentiation would be promoted. In our study, we stimulated
groups of synapses quasi-synchronously by using glutamate
photolysis, which resulted in the activation of both CaMKII
and calcineurin, as supported by the sensitivity to inhibitors.
Although the postsynaptic neuron was depolarized during
glutamate photolysis, under such LTP triggering stimulation, a
trend for a higher rise in calcium was observed for proximal
spines located within 3 mm of the stimulated cluster than
more distal spines (Figure S2). Locally boosted calcineurin activation by the spatially clustered stimulation, therefore, could
have promoted the shrinkage of spines close to the stimulated
spines, whereas a potentially more efficient spread of CaMKII
signaling than those of calcineurin and/or a less strongly
thresholded CaMKII activation at distal spines could have promoted their potentiation. Taken together, the expression of heterosynaptic plasticity at individual synapses is strongly biased
by ongoing spatiotemporal patterns of synaptic activity. Here,
we discussed the simplest possibility in an attempt to explain
the observed simultaneous engagement of heterosynaptic
potentiation and depression, based on the presumed differential activation thresholds and diffusional properties of calcineurin and CaMKII signaling. Nevertheless, interactions of synaptic scaffolding proteins and regulatory proteins with
calcineurin and CaMKII at individual spines (Dudilot et al.,
2020; Woolfrey and Dell’Acqua, 2015) are also likely to
contribute to the landscape of heterosynaptic responses. A
better understanding of the coordinated responses of signaling
networks to patterned activity across neighboring synapses is
warranted.
Here, somewhat unexpectedly, our observations suggest a
requirement for NO in presynaptic hetLTD. Although NO has
been traditionally studied in the context of presynaptic LTP
(e.g., Arancio et al., 1996; Stanton et al., 2005; Ratnayaka
et al., 2012; Padamsey et al., 2017), several studies have implicated NO in homosynaptic LTD (Zhuo et al., 1994; Reyes-Harde
et al., 1999; Stanton et al., 2003; Gage et al., 1997). These
studies collectively support a putative model for NO-dependent

(B) Induction of spike-timing-dependent LTP led to long-lasting increase in the EPSP slope of the stimulated pathway (homosynaptic) and decrease in the EPSP
slope of the unstimulated pathway (heterosynaptic). Example traces illustrate the average EPSP waveform before LTP induction and during the last 5 min of the
recording of the stimulated (black) and unstimulated (red) pathway.
(C) PPR before and after the induction of LTP. The stimulated pathway showed a significant decrease of the PPR (Wilcoxon rank-sum test). Conversely, the PPR
increased in the unstimulated pathway.
(D) Bath application of 100 mM L-NAME abolished heterosynaptic LTD but did not affect the magnitude of homosynaptic LTP. Spread of LTP was unmasked in the
unstimulated pathway.
(E) PPR changes of both stimulated and unstimulated pathways were abolished.
(F) Summary of EPSP slope changes.
(G) Summary of PPR changes (see also Figure S3). Statistical significance was assessed using Mann-Whitney U test. Data points are shown as mean ± SEM. See
also Figure S5.

Cell Reports 34, 108693, January 26, 2021 11

ll
OPEN ACCESS

LTD at Schaffer-collateral-CA1 synapses that involves the activation of presynaptic internal Ca2+ stores by cyclic guanosine
monophosphate (cGMP)-dependent protein kinase (PKG), a
known downstream effector of NO signaling, along with ADP-ribosyl cyclase and ryanodine receptor, and also requires the concurrent activation of presynaptic CaMKII (Reyes-Harde et al.,
1999; Stanton and Gage, 1996; Gage et al., 1997). The involvement of CaMKII in NO-mediated presynaptic depression is in
line with our observation that inhibiting CaMKII abolished presynaptic hetLTD, although we did not discern the exact locus
of CaMKII activation in our experiments. Whether or not postsynaptically generated NO at the stimulated spines diffuses short
distances to neighboring presynaptic terminals remains to be
confirmed, although the spatially confined heterosynaptic
spread of presynaptic LTD (within 4 mm of stimulated cluster) is
consistent with the direct involvement of short-lived NO as the
extracellularly diffusible messenger.
Previously, a crucial role for NO in presynaptic LTP was
demonstrated by the Emptage group (Padamsey et al., 2017).
Specifically, NO synthesis driven by activation of L-type
voltage-gated Ca2+ channels (L-VGCCs) coupled to presynaptic
activity promoted an increase in Pr. The NO in the current study
was unlikely to have originated by the same mechanism due to
the fast inactivation of L-VGCCs at depolarized membrane potential used during the induction of cl-LTP. Such differences between the two studies suggest distinct functional roles for NO
that depend on the context in which NO is generated. In the previous study, local photolysis of NO at single synapses on its own
was not sufficient to induce LTP or LTD in the absence of presynaptic stimulation, whereas after dialyzing caged NO into the
postsynaptic neuron, photolyzing a larger area could elicit LTD
(Figure 5 Suppl. 2 in Padamsey et al., 2017). Therefore, induction
of presynaptic hetLTD might require a critical amount of NO in
conjunction with synaptic activity that is dependent on NMDARs
but not L-VGCCs. In contrast, L-VGCC is needed for NO-dependent presynaptic potentiation.
Lastly, our study identified a slight trend for pre- and postsynapses to undergo correlated changes. Nevertheless, the correlation was not strong, and moreover, an apparent lack of correlation was noted when we inhibited calcineurin and observed
presynaptic LTD in combination with postsynaptic LTP. The
observed uncoupling of the direction of plasticity and dissociation of the pharmacological sensitivity of pre- and postsynaptic
changes suggest that the pre- and the postsynaptic plasticity
might represent in part independent parallel processes that
follow unique dynamics. That the pre- and the postsynaptic
sides of the same synapse do not necessarily change their
strengths in the same direction may seem counterintuitive
based on the overall matching of the pre- and the postsynaptic
strengths at basal state across a synapse population (Tokuoka
and Goda, 2008; Kay et al., 2011). Nevertheless, a similar disjunctional plasticity of presynaptic potentiation associated with
postsynaptic depression was observed in cultured hippocampal neurons (Xu et al., 2013; Letellier et al., 2019) and between
layer 4 neurons in visual cortical slices (Sáez and Friedlander,
2009). Such a mechanism could confer homeostatic stability
to synapses, although the precise properties of disjunctional
plasticity warrant further studies. Notably, the independent

12 Cell Reports 34, 108693, January 26, 2021

Article
regulation of pre- and postsynaptic strengths has recently
been explained by the difference in functional impact on the
postsynaptic neuron (Costa et al., 2017, Tong et al., 2020).
For future studies on the regulation of synaptic strength, it is
therefore crucial to consider both, pre- and postsynaptic
strengths for fully interpreting the consequences on neural
network activity.
STAR+METHODS
Detailed methods are provided in the online version of this paper
and include the following:
d KEY RESOURCES TABLE
d RESOURCE AVAILABILITY

B Lead contact
B Materials availability
B Data and code availability
d EXPERIMENTAL MODEL AND SUBJECT DETAILS
d METHOD DETAILS

B Preparation of organotypic hippocampal slice culture
B Preparation of acute hippocampal slices
B 2-photon imaging
B Optical quantal analysis
B Glutamate photolysis
B Induction of structural LTP
B Analysis of spine structural changes
B Spine calcium imaging
B Electrophysiological recordings in acute slices
B Pharmacology
d QUANTIFICATION AND STATISTICAL ANALYSIS

SUPPLEMENTAL INFORMATION
Supplemental Information can be found online at https://doi.org/10.1016/j.
celrep.2021.108693.
ACKNOWLEDGMENTS
We would like to thank Tiago Branco and Simon Butt for helpful discussion and
feedback and Karen Zito for her generous technical guidance. R.T. received
the RIKEN-SdV (‘‘Studienstiftung des Deutschen Volkes’’) fellowship and
was jointly funded by the Clarendon Fund, University of Oxford, and the Medical Research Council UK. The study was also supported by funds from the
RIKEN Brain Science Institute, RIKEN Center for Brain Science, and JSPS
Core-to-Core Program A Advanced Research Networks (Y.G.).
AUTHOR CONTRIBUTIONS
R.T., Y.G., and N.J.E. conceived the study and designed the experiments. R.T.
conducted and analyzed the experiments. T.E.C. performed and analyzed the
calcium imaging experiments. The manuscript was drafted by R.T., and all authors edited and approved the final manuscript.
DECLARATION OF INTERESTS
The authors declare no competing interests.
Received: May 5, 2020
Revised: August 13, 2020
Accepted: January 5, 2021
Published: January 26, 2021

ll
Article

OPEN ACCESS

REFERENCES

Fujii, H., Inoue, M., Okuno, H., Sano, Y., Takemoto-Kimura, S., Kitamura, K.,
Kano, M., and Bito, H. (2013). Nonlinear decoding and asymmetric representation of neuronal input information by CaMKIIa and calcineurin. Cell Rep. 3,
978–987.

Arancio, O., Kiebler, M., Lee, C.J., Lev-Ram, V., Tsien, R.Y., Kandel, E.R., and
Hawkins, R.D. (1996). Nitric oxide acts directly in the presynaptic neuron to
produce long-term potentiation in cultured hippocampal neurons. Cell 87,
1025–1035.
Barbour, B., Brunel, N., Hakim, V., and Nadal, J.P. (2007). What can we learn
from synaptic weight distributions? Trends Neurosci. 30, 622–629.
Bayazitov, I.T., Richardson, R.J., Fricke, R.G., and Zakharenko, S.S. (2007).
Slow presynaptic and fast postsynaptic components of compound longterm potentiation. J. Neurosci. 27, 11510–11521.
Béı̈que, J.-C., Lin, D.-T., Kang, M.-G., Aizawa, H., Takamiya, K., and Huganir,
R.L. (2006). Synapse-specific regulation of AMPA receptor function by PSD95. Proc. Natl. Acad. Sci. USA 103, 19535–19540.
Bian, W.-J., Miao, W.-Y., He, S.-J., Qiu, Z., and Yu, X. (2015). Coordinated
Spine Pruning and Maturation Mediated by Inter-Spine Competition for Cadherin/Catenin Complexes. Cell 162, 808–822.
Bliss, T.V.P., and Collingridge, G.L. (1993). A synaptic model of memory: longterm potentiation in the hippocampus. Nature 361, 31–39.
Bourne, J.N., and Harris, K.M. (2011). Coordination of size and number of
excitatory and inhibitory synapses results in a balanced structural plasticity
along mature hippocampal CA1 dendrites during LTP. Hippocampus 21,
354–373.
Branco, T., and Häusser, M. (2010). The single dendritic branch as a fundamental functional unit in the nervous system. Curr. Opin. Neurobiol. 20,
494–502.

Gage, A.T., Reyes, M., and Stanton, P.K. (1997). Nitric-oxide-guanylylcyclase-dependent and -independent components of multiple forms of longterm synaptic depression. Hippocampus 7, 286–295.
Govindarajan, A., Kelleher, R.J., and Tonegawa, S. (2006). A clustered plasticity model of long-term memory engrams. Nat. Rev. Neurosci. 7, 575–583.
Govindarajan, A., Israely, I., Huang, S.-Y., and Tonegawa, S. (2011). The dendritic branch is the preferred integrative unit for protein synthesis-dependent
LTP. Neuron 69, 132–146.
Hardingham, N., Dachtler, J., and Fox, K. (2013). The role of nitric oxide in presynaptic plasticity and homeostasis. Front. Cell. Neurosci. 7, 190.
Harris, J.J., Jolivet, R., and Attwell, D. (2012). Synaptic energy use and supply.
Neuron 75, 762–777.
Harvey, C.D., and Svoboda, K. (2007). Locally dynamic synaptic learning rules
in pyramidal neuron dendrites. Nature 450, 1195–1200.
Harvey, C.D., Yasuda, R., Zhong, H., and Svoboda, K. (2008). The spread of
Ras activity triggered by activation of a single dendritic spine. Science 321,
136–140.
Hayama, T., Noguchi, J., Watanabe, S., Takahashi, N., Hayashi-Takagi, A., Ellis-Davies, G.C.R., Matsuzaki, M., and Kasai, H. (2013). GABA promotes the
competitive selection of dendritic spines by controlling local Ca2+ signaling.
Nat. Neurosci. 16, 1409–1416.

Branco, T., Staras, K., Darcy, K.J., and Goda, Y. (2008). Local dendritic activity
sets release probability at hippocampal synapses. Neuron 59, 475–485.

Hedrick, N.G., Harward, S.C., Hall, C.E., Murakoshi, H., McNamara, J.O., and
Yasuda, R. (2016). Rho GTPase complementation underlies BDNF-dependent
homo- and heterosynaptic plasticity. Nature 538, 104–108.

Burette, A., Collman, F., Micheva, K.D., Smith, S.J., and Weinberg, R.J. (2015).
Knowing a synapse when you see one. Front. Neuroanat. 9, 100.

Jia, H., Rochefort, N.L., Chen, X., and Konnerth, A. (2010). Dendritic organization of sensory input to cortical neurons in vivo. Nature 464, 1307–1312.

Buzsáki, G., and Mizuseki, K. (2014). The log-dynamic brain: how skewed distributions affect network operations. Nat. Rev. Neurosci. 15, 264–278.

Ju, W., Morishita, W., Tsui, J., Gaietta, G., Deerinck, T.J., Adams, S.R., Garner,
C.C., Tsien, R.Y., Ellisman, M.H., and Malenka, R.C. (2004). Activity-dependent regulation of dendritic synthesis and trafficking of AMPA receptors.
Nat. Neurosci. 7, 244–253.

Chater, T.E., and Goda, Y. (2020). My Neighbour Hetero-deconstructing the
mechanisms underlying heterosynaptic plasticity. Curr. Opin. Neurobiol. 67,
106–114.
Chen, X., Leischner, U., Rochefort, N.L., Nelken, I., and Konnerth, A. (2011).
Functional mapping of single spines in cortical neurons in vivo. Nature 475,
501–505.
Costa, R.P., Padamsey, Z., D’Amour, J.A., Emptage, N.J., Froemke, R.C., and
Vogels, T.P. (2017). Synaptic Transmission Optimization Predicts Expression
Loci of Long-Term Plasticity. Neuron 96, 177–189.e7.
Deguchi, Y., Donato, F., Galimberti, I., Cabuy, E., and Caroni, P. (2011).
Temporally matched subpopulations of selectively interconnected principal
neurons in the hippocampus. Nat. Neurosci. 14, 495–504.
Druckmann, S., Feng, L., Lee, B., Yook, C., Zhao, T., Magee, J.C., and Kim, J.
(2014). Structured synaptic connectivity between hippocampal regions.
Neuron 81, 629–640.
Dudilot, A., Trillaud-Doppia, E., and Boehm, J. (2020). RCAN1 Regulates Bidirectional Synaptic Plasticity. Curr. Biol. 30, 1167–1176.e2.
El-Boustani, S., Ip, J.P.K., Breton-Provencher, V., Knott, G.W., Okuno, H.,
Bito, H., and Sur, M. (2018). Locally coordinated synaptic plasticity of visual
cortex neurons in vivo. Science 360, 1349–1354.
Emptage, N., Bliss, T.V., and Fine, A. (1999). Single synaptic events evoke
NMDA receptor-mediated release of calcium from internal stores in hippocampal dendritic spines. Neuron 22, 115–124.
Engert, F., and Bonhoeffer, T. (1997). Synapse specificity of long-term potentiation breaks down at short distances. Nature 388, 279–284.
Fonseca, R., Nägerl, U.V., Morris, R.G.M., and Bonhoeffer, T. (2004).
Competing for memory: hippocampal LTP under regimes of reduced protein
synthesis. Neuron 44, 1011–1020.
Fu, M., Yu, X., Lu, J., and Zuo, Y. (2012). Repetitive motor learning induces coordinated formation of clustered dendritic spines in vivo. Nature 483, 92–95.

Kay, L., Humphreys, L., Eickholt, B.J., and Burrone, J. (2011). Neuronal activity
drives matching of pre- and postsynaptic function during synapse maturation.
Nat. Neurosci. 14, 688–690.
Kleindienst, T., Winnubst, J., Roth-Alpermann, C., Bonhoeffer, T., and Lohmann, C. (2011). Activity-dependent clustering of functional synaptic inputs
on developing hippocampal dendrites. Neuron 72, 1012–1024.
Kruijssen, D.L.H., and Wierenga, C.J. (2019). Single Synapse LTP: A Matter of
Context? Front. Cell. Neurosci. 13, 496.
Lee, K.J., Park, I.S., Kim, H., Greenough, W.T., Pak, D.T.S., and Rhyu, I.J.
(2013). Motor skill training induces coordinated strengthening and weakening
between neighboring synapses. J. Neurosci. 33, 9794–9799.
Lee, K.-S., Vandemark, K., Mezey, D., Shultz, N., and Fitzpatrick, D. (2019).
Functional Synaptic Architecture of Callosal Inputs in Mouse Primary Visual
Cortex. Neuron 101, 421–428.e5.
Letellier, M., Levet, F., Thoumine, O., and Goda, Y. (2019). Differential role of
pre- and postsynaptic neurons in the activity-dependent control of synaptic
strengths across dendrites. PLoS Biol. 17, e2006223.
Liang, L., Fratzl, A., Goldey, G., Ramesh, R.N., Sugden, A.U., Morgan, J.L.,
Chen, C., and Andermann, M.L. (2018). A Fine-Scale Functional Logic to
Convergence from Retina to Thalamus. Cell 173, 1343–1355.e24.
Lynch, G.S., Dunwiddie, T., and Gribkoff, V. (1977). Heterosynaptic depression: a postsynaptic correlate of long-term potentiation. Nature 266, 737–739.
Matsuzaki, M., Ellis-Davies, G.C., Nemoto, T., Miyashita, Y., Iino, M., and Kasai, H. (2001). Dendritic spine geometry is critical for AMPA receptor expression in hippocampal CA1 pyramidal neurons. Nat. Neurosci. 4, 1086–1092.
Matsuzaki, M., Honkura, N., Ellis-Davies, G.C.R., and Kasai, H. (2004). Structural basis of long-term potentiation in single dendritic spines. Nature 429,
761–766.

Cell Reports 34, 108693, January 26, 2021 13

ll
OPEN ACCESS

McBride, T.J., Rodriguez-Contreras, A., Trinh, A., Bailey, R., and Debello,
W.M. (2008). Learning drives differential clustering of axodendritic contacts
in the barn owl auditory system. J. Neurosci. 28, 6960–6973.
Mel, B.W. (1992). The Clusteron: Toward a Simple Abstraction for a Complex
Neuron. Advances in Neural Information Processing Systems (Morgan Kaufmann), pp. 35–42.
Meyer, D., Bonhoeffer, T., and Scheuss, V. (2014). Balance and stability of synaptic structures during synaptic plasticity. Neuron 82, 430–443.
Nikonenko, I., Nikonenko, A., Mendez, P., Michurina, T.V., Enikolopov, G., and
Muller, D. (2013). Nitric oxide mediates local activity-dependent excitatory
synapse development. Proc. Natl. Acad. Sci. USA 110, E4142–E4151.
Nimchinsky, E.A., Yasuda, R., Oertner, T.G., and Svoboda, K. (2004). The
number of glutamate receptors opened by synaptic stimulation in single hippocampal spines. J. Neurosci. 24, 2054–2064.

Article
Schuman, E.M., and Madison, D.V. (1994). Locally distributed synaptic potentiation in the hippocampus. Science 263, 532–536.
Smith, M.A., Ellis-Davies, G.C.R., and Magee, J.C. (2003). Mechanism of the
distance-dependent scaling of Schaffer collateral synapses in rat CA1 pyramidal neurons. J. Physiol. 548, 245–258.
Stanton, P.K., and Gage, A.T. (1996). Distinct synaptic loci of Ca2+/calmodulin-dependent protein kinase II necessary for long-term potentiation and
depression. J. Neurophysiol. 76, 2097–2101.
Stanton, P.K., Winterer, J., Bailey, C.P., Kyrozis, A., Raginov, I., Laube, G.,
€ller, W. (2003). Long-term depression of preVeh, R.W., Nguyen, C.Q., and Mu
synaptic release from the readily releasable vesicle pool induced by NMDA receptor-dependent retrograde nitric oxide. J. Neurosci. 23, 5936–5944.

Nishiyama, J., and Yasuda, R. (2015). Biochemical Computation for Spine
Structural Plasticity. Neuron 87, 63–75.

€ller, W. (2005). Imaging LTP of
Stanton, P.K., Winterer, J., Zhang, X.-L., and Mu
presynaptic release of FM1-43 from the rapidly recycling vesicle pool of
Schaffer collateral-CA1 synapses in rat hippocampal slices. Eur. J. Neurosci.
22, 2451–2461.

O’Rourke, N.A., Weiler, N.C., Micheva, K.D., and Smith, S.J. (2012). Deep molecular diversity of mammalian synapses: why it matters and how to measure
it. Nat. Rev. Neurosci. 13, 365–379.

Sutton, M.A., Ito, H.T., Cressy, P., Kempf, C., Woo, J.C., and Schuman, E.M.
(2006). Miniature neurotransmission stabilizes synaptic function via tonic suppression of local dendritic protein synthesis. Cell 125, 785–799.

Oh, W.C., Parajuli, L.K., and Zito, K. (2015). Heterosynaptic structural plasticity
on local dendritic segments of hippocampal CA1 neurons. Cell Rep. 10,
162–169.
Padamsey, Z., Tong, R., and Emptage, N. (2017). Glutamate is required for
depression but not potentiation of long-term presynaptic function. eLife 6, 839.
Padamsey, Z., Tong, R., and Emptage, N. (2019). Optical Quantal Analysis Using Ca2+ Indicators: A Robust Method for Assessing Transmitter Release Probability at Excitatory Synapses by Imaging Single Glutamate Release Events.
Front. Synaptic Neurosci. 11, 5.
Poirazi, P., and Mel, B.W. (2001). Impact of active dendrites and structural
plasticity on the memory capacity of neural tissue. Neuron 29, 779–796.
Poirazi, P., Brannon, T., and Mel, B.W. (2003). Pyramidal neuron as two-layer
neural network. Neuron 37, 989–999.
Ratnayaka, A., Marra, V., Bush, D., Burden, J.J., Branco, T., and Staras, K.
(2012). Recruitment of resting vesicles into recycling pools supports NMDA receptor-dependent synaptic potentiation in cultured hippocampal neurons.
J. Physiol. 590, 1585–1597.
Reyes-Harde, M., Potter, B.V., Galione, A., and Stanton, P.K. (1999). Induction
of hippocampal LTD requires nitric-oxide-stimulated PKG activity and Ca2+
release from cyclic ADP-ribose-sensitive stores. J. Neurophysiol. 82, 1569–
1576.
Royer, S., and Paré, D. (2003). Conservation of total synaptic weight through
balanced synaptic depression and potentiation. Nature 422, 518–522.
Sáez, I., and Friedlander, M.J. (2009). Plasticity between neuronal pairs in layer
4 of visual cortex varies with synapse state. J. Neurosci. 29, 15286–15298.
Scanziani, M., Malenka, R.C., and Nicoll, R.A. (1996). Role of intercellular interactions in heterosynaptic long-term depression. Nature 380, 446–450.
Schindelin, J., Arganda-Carreras, I., Frise, E., Kaynig, V., Longair, M., Pietzsch,
T., Preibisch, S., Rueden, C., Saalfeld, S., Schmid, B., et al. (2012). Fiji: an
open-source platform for biological-image analysis. Nat. Methods 9, 676–682.
Schneider, C.A., Rasband, W.S., and Eliceiri, K.W. (2012). NIH Image to ImageJ: 25 years of image analysis. Nat. Methods 9, 671–675.
Scholl, B., Wilson, D.E., and Fitzpatrick, D. (2017). Local Order within Global
Disorder: Synaptic Architecture of Visual Space. Neuron 96, 1127–1138.e4.

14 Cell Reports 34, 108693, January 26, 2021

Takahashi, N., Kitamura, K., Matsuo, N., Mayford, M., Kano, M., Matsuki, N.,
and Ikegaya, Y. (2012). Locally synchronized synaptic inputs. Science 335,
353–356.
Tang, A.-H., Chen, H., Li, T.P., Metzbower, S.R., MacGillavry, H.D., and Blanpied, T.A. (2016). A trans-synaptic nanocolumn aligns neurotransmitter release
to receptors. Nature 536, 210–214.
Tokuoka, H., and Goda, Y. (2008). Activity-dependent coordination of presynaptic release probability and postsynaptic GluR2 abundance at single synapses. Proc. Natl. Acad. Sci. USA 105, 14656–14661.
Tong, R., Emptage, N.J., and Padamsey, Z. (2020). A two-compartment model
of synaptic computation and plasticity. Mol. Brain 13, 79.
Volgushev, M., Voronin, L.L., Chistiakova, M., and Singer, W. (1997). Relations
between long-term synaptic modifications and paired-pulse interactions in the
rat neocortex. Eur. J. Neurosci. 9, 1656–1665.
Wiegert, J.S., and Oertner, T.G. (2013). Long-term depression triggers the selective elimination of weakly integrated synapses. Proc. Natl. Acad. Sci. USA
110, E4510–E4519.
Winnubst, J., and Lohmann, C. (2012). Synaptic clustering during development and learning: the why, when, and how. Front. Mol. Neurosci. 5, 70.
Winnubst, J., Cheyne, J.E., Niculescu, D., and Lohmann, C. (2015). Spontaneous Activity Drives Local Synaptic Plasticity In Vivo. Neuron 87, 399–410.
Woolfrey, K.M., and Dell’Acqua, M.L. (2015). Coordination of Protein Phosphorylation and Dephosphorylation in Synaptic Plasticity. J. Biol. Chem.
290, 28604–28612.
Xu, W., Tse, Y.C., Dobie, F.A., Baudry, M., Craig, A.M., Wong, T.P., and Wang,
Y.T. (2013). Simultaneous monitoring of presynaptic transmitter release and
postsynaptic receptor trafficking reveals an enhancement of presynaptic activity in metabotropic glutamate receptor-mediated long-term depression.
J. Neurosci. 33, 5867–5877.
Zhuo, M., Kandel, E.R., and Hawkins, R.D. (1994). Nitric oxide and cGMP
can produce either synaptic depression or potentiation depending on the frequency of presynaptic stimulation in the hippocampus. Neuroreport 5, 1033–
1036.

ll
OPEN ACCESS

Article
STAR+METHODS
KEY RESOURCES TABLE

REAGENT or RESOURCE

SOURCE

IDENTIFIER

Oregon Green 488 BAPTA-1
(hexapotassium salt)

ThermoFisher Scientific

Cat#O6806

Fluo-5F

ThermoFisher Scientific

Cat#F14221

AlexaFluor 488

ThermoFisher Scientific

Cat#A10436

AlexaFluor 594

ThermoFisher Scientific

Cat# A10438

MNI-glutamate

Tocris

Cat#1490; CAS: 295325-62-1

AP5

Abcam

Cat#ab120003; CAS: 79055-68-8

KN62

Tocris

Cat#1277; CAS: 127191-97-3

FK-506

Tocris

Cat#3631; CAS: 104987-11-3

L-NAME

Cayman Chemical

Cat#Cay80210; CAS: 51298-62-5

Chemicals, peptides, and recombinant proteins

Picrotoxin

Sigma Aldrich

Cat#P1675; CAS: 124-87-8

Trolox

Sigma Aldrich

Cat#238813; CAS: 53188-07-1

Ascorbic acid

Sigma Aldrich

Cat#A92902; CAS: 50-81-7

Harlan UK, Nihon SLC

RRID: RGD_13508588

ImageJ

NIH

https://imagej.nih.gov/ij/

pClamp

Molecular Devices

RRID: SCR_011323

Zen 2009

Zeiss Microscopy

https://www.zeiss.com/corporate/int/
home.html

Zeiss LSM780

Zeiss Microscopy

https://www.zeiss.com/corporate/int/
home.html

Leica DMLFSA

Leica

https://www.leica-microsystems.com

Experimental models: organisms/strains
Rat: WI Wistar rat (rattus norvegicus)
Software and algorithms

Other

RESOURCE AVAILABILITY
Lead contact
Further information and requests for resources and reagents should be directed to and will be fulfilled by the Lead Contact, Yukiko
Goda (yukiko.goda@riken.jp).
Materials availability
This study did not generate new unique reagents.
Data and code availability
The data generated in the current study and analysis code have not been deposited in a public repository but are available from corresponding authors upon request.
EXPERIMENTAL MODEL AND SUBJECT DETAILS
Male Wistar rat pups (RRID: RGD_13508588, age: postnatal day 6-7 and 14-21) were used for making organotypic and acute hippocampal slices. Animals were supplied by Harlan UK and Nihon SLC. All animal experiments were either approved by the RIKEN
Animal Experiments Committee and performed in accordance with the RIKEN rules and guidelines or were in accordance with the
Animals (Scientific Procedures) Act, 1986 (UK).

Cell Reports 34, 108693, January 26, 2021 e1

ll
OPEN ACCESS

Article

METHOD DETAILS
Preparation of organotypic hippocampal slice culture
Organotypic hippocampal slices were prepared as previously reported (Padamsey et al., 2017). Hippocampi of postnatal day P6-7
Wistar rat pups (Harlan UK, Nihon SLC) were isolated and cut into 350 mm-thick transverse slices on a McIlwain tissue chopper
(Mickle Laboratory Engineering Co. Ltd. and Cavey Laboratory Engineering Co. Ltd.). The brain was dissected in ice-cold Earle’s
Balanced Salt solution (EBSS) containing 35 mM glucose, 20-25 mM HEPES and adjusted to pH 7.35 with 5 mM NaOH. Slices
were transferred onto cell culture inserts (0.4 mm pore size, Merk Millipore) and placed in a 6-well plate filled with 1 ml/well of culture
media containing 50% Minimum Essential Medium (MEM, Thermo Fisher Scientific), 23% EBSS, 25% horse serum (Thermo Fisher
Scientific), and 36 mM glucose. Some slices (used for experiments in Figures 2 and 3) were cultured in media that included 2% B-27
Supplement (Thermo Fisher Scientific). Slices were maintained at 37 C and 5% CO2 and used for experiments at DIV10-15.
During experiments, slices were constantly perfused (1-2 ml/min) with artificial cerebrospinal fluid (aCSF) containing (in mM) 145
NaCl, 2.5 KCl, 16 NaHCO3, 1.2 NaH2PO4, 11 glucose, 3 CaCl2, and 2 MgCl2. For experiments with intensive imaging, 1 mM Trolox
(Sigma Aldrich) and 0.2 mM ascorbic acid were additionally included. The aCSF was bubbled with 95% O2 and 5% CO2 and heated to
near-physiological temperature (33-35 C) using an in-line heater (Warner Instruments).
Preparation of acute hippocampal slices
Acute hippocampal slices were prepared from P14-21 male Wistar rats (Harlan UK). Rats were sacrificed and brains were extracted
and immediately submerged in ice-cold dissection media saturated with 95% O2 and 5% CO2. Dissection media consisted of (in
mM): 65 sucrose, 85 NaCl, 2.5 KCl, 25 NaHCO3, 1.25 NaH2PO4, 10 glucose, 7 MgCl2, and 0.2 CaCl2. The cerebellum was manually
removed by a coronal cut using a razor blade to create a flat surface, and the brain was glued onto a platform of the vibratome (Microm HM 650V, Thermo Scientific). Brains were further stabilized by positioning a block of 2% agar to contact the dorsal surface of
the brain. 400 mm coronal slices were cut, and the hemispheres were separated and transferred into a recovery chamber containing
standard aCSF (in mM: 120 NaCl, 2.5 KCl, 26 NaHCO3, 1.2 NaH2PO4, 11 glucose, 1 MgCl2, and 2 CaCl2), which was bubbled with
95% O2 and 5% CO2. Slices were allowed to recover for 8 min at 37 C and 60 min at room temperature prior to use. Slices were
maintained at room temperature up to 4-5 h.
For experiments, slices were transferred to the recording chamber and secured using a ‘‘harp’’ (Warner Instruments). The
recording chamber was constantly perfused (2-3 ml/min) with aCSF, and 100 mM picrotoxin (Sigma Aldrich) was added to block
inhibitory synaptic transmission. The aCSF was bubbled with 95% O2 and 5% CO2 and heated to near-physiological temperature
(33-35 C) using a custom made in-line heater.
2-photon imaging
CA1 pyramidal cells were bolus-loaded with 0.5 mM AF594 (ThermoFisher Scientific) and 1 mM Oregon Green BAPTA-1 (OGB-1,
ThermoFisher Scientific) diluted in standard internal solution (in mM: 135 K-Gluconate, 10 KCl, 10 HEPES, 1 MgCl2, 2 Na2ATP,
and 0.4 Na3GTP) for 45-60 s. OGB-1 was included for optical quantal analysis (see Optical quantal analysis below). The dyes
were given at least 10 min and up to 45 min to allow for sufficient diffusion in proximal dendrites such that no obvious change in
mean fluorescence intensity over time was observed (Padamsey et al., 2017). After identifying a segment of dendrite with linear geometry that was largely contained within a single optical plane and hence suited for imaging, presynaptic strength was measured via
optical quantal analysis. Afterward, z stacks (0.5 mm steps) capturing ~35 mm of the dendrite centered on the previously measured
spine were taken every 5 min. Images were acquired on a Zeiss LSM780 confocal laser scanning microscope with a 63x 1.0 NA
objective (Plan-APOCHROMAT, Zeiss) using commercial software provided by Zeiss (Zen 2009, Version 6.0.0.303). Images were
taken at 4x zoom, which led to a lateral pixel size of 65.9 nm. Both fluorescent dyes were simultaneously excited using an
800 nm 2-photon laser source (Coherent) and emission was separated using bandpass emission filters. At least three baseline images were taken before applying clustered stimulation to elicit LTP (cl-LTP). Images were then taken at 1, 5, 10, 20, 30, and 45 min
after the cl-LTP induction. Experiments showing elevated Ca2+ levels or outgrowth of filopodia were discarded as these features
correlated to subsequent signs of cell death, such as fragmentation of the dendrite. At the end of the experiment, a z stack was taken
(2 mm steps, 1x zoom), which included the dendrite of interest, the cell body, and the position of the stimulation and puffing
electrodes.
Optical quantal analysis
We performed optical quantal analysis as previously reported (Padamsey et al., 2019). A stimulation electrode (tungsten monopolar in
glass pipette) was positioned within 5-10 mm of the dendrite of interest. xt-scans (line scans) covering large fractions of spines were
taken at 500 Hz while synapses were stimulated with a paired pulse given at 70 ms inter-pulse interval. The imaging software was
synchronized with Clampex software (Molecular Devices) via TTL pulses. Stimulation and the initiation of image acquisition were
controlled by Clampex. Responsive spines were detected as time-locked OGB-1 Ca2+ responses (EPSCaTs). In order to prevent potential confounds of AP failures, which will inflate the estimate of release failures, the stimulation strength was increased until the
observed Pr stabilized, such that the average failure rate remained constant over time (Emptage et al., 1999). For data acquisition,
line scans were oriented to capture both the spine of interest and the dendrite, which allowed also the detection of local dendritic

e2 Cell Reports 34, 108693, January 26, 2021

ll
Article

OPEN ACCESS

spikes or back-propagating APs. Before the induction of cl-LTP, optical quantal analysis was limited to 20-25 trials (inter-trial interval:
15 s) to prevent photodynamic damage. Optical quantal analysis was repeated 30 min after the induction of LTP and 20-30 trials were
taken. The stimulation strength was further increased at the end of the experiment and a paired-pulse at 5 ms inter-pulse interval was
given to ensure that release events could still be detected. This also ruled out the possibility that an observed decrease in Pr is due to
spatial drift of the stimulation electrode, a change in action potential initiation, or a decrease in the EPSCaT amplitude below detection threshold. Experiments were excluded if no release events could be detected. EPSCaTs were analyzed in ImageJ (Schindelin
et al., 2012; Schneider et al., 2012) and a custom written Python script. Fluorescence signals were background-subtracted, averaged
within the spatial dimension, and normalized to baseline (defined as the average fluorescence before the first stimulation pulse). Successful release events were counted manually. All results were cross-validated using an automated analysis script, which detected
successful release as an increase in the average Ca2+ signal within 6 ms after stimulation that exceeded two standard deviations
above the mean of the baseline. No qualitative difference was found between the two methods; however, the automated method
was not robust against imaging artifacts or unstable baselines. To distinguish successful release events from dendritic or somatic
spikes, the peak intensity was required to be at least 50% larger in the spine compared to the neighboring dendrite. Alternatively,
the identification of EPSCaTs was aided by the signal onset that preceded dendritic/somatic spikes. Pr was calculated as the ratio
between the number of successful release events to the first stimulation pulse and number of trials. The sampling error is given by the
binomial theorem
rﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃ
pð1  pÞ
s=
N
where p is the release probability and N the sample number. For an average sample number of N = 25, the maximum expected error is
therefore ~0.1.
Glutamate photolysis
Caged glutamate was locally applied during photolysis through a glass pipette coupled to a picospritzer (5-10 psi, Parker Instrumentation). 10 mM MNI-glutamate (Tocris) dissolved in Tyrode’s solution (in mM: 58.44 NaCl, 2.5 KCl, 20 HEPES, 30 glucose, 3 CaCl2,
2 MgCl2, pH-adjusted to 7.2 using 5 mM NaOH) was filtered (0.45 mm pore syringe filter, Merck Millipore) and loaded into a glass
pipette (3-4 MU). The pipette was positioned ~10-20 mm from the dendrite of interest. A 2-photon laser source (720 nm) was used
for focal photolysis of glutamate. Photolysis was controlled by a custom written script for the imaging software and was
synchronized with the electrophysiology via Clampex. Photolysis consisted of 4 ms pulses, and the laser power was adjusted for
each experiment to obtain Ca2+ transients similar in size compared to EPSCaTs measured in the same cell (laser power varied
from 2.4-3.5 mW measured at the sample). For the induction of clustered LTP, laser power was reduced slightly (to ~80%) to prevent
over-excitation, which would frequently cause excitotoxicity observed as sustained elevated Ca2+ in the dendrite. For quasi-synchronous photolysis, the inter-pulse interval was ~2.2 ms.
Induction of structural LTP
Clustered structural LTP (cl-LTP) was induced by quasi-synchronous glutamate photolysis at 5-7 spines paired with postsynaptic
depolarization, 303 at 0.5 Hz. The postsynaptic neuron was re-patched with internal solution containing 100 mM OGB-1 and
50 mM AF594 to reduce the osmotic pressure of the dye. Photolysis was set up after the GU-seal was established by focusing on
the focal plane containing at least 5 neighboring spine heads, and cl-LTP was initiated within 10 s of establishing whole-cell
mode in order to prevent dialysis. Cells were held at 0 mV in voltage clamp mode. The average access conductance was
Ga = 45.86 ± 2.56 nS. Cells were allowed to reseal after the induction by gently retracting the patch pipette. The success rate of
re-patch and reseal was ~80%. The maximum distance and mean pairwise distance between stimulated spines did not vary across
experimental conditions.
Analysis of spine structural changes
Spine size was analyzed in ImageJ. The integrated fluorescence of a spine contained in a rectangular regions of interest (ROI) was
used as an estimate of spine size. Integrated fluorescence was measured on the z-image slice that resulted in the maximal value. The
background fluorescence equivalent to the size of the ROI was subtracted. The re-patch during clustered LTP induction caused a
substantial global reduction in the fluorescence intensity. Fluorescence signals were therefore normalized to the mean fluorescence
intensity of the underlying dendrite, which was assumed to remain structurally stable (Nimchinsky et al., 2004). At least three dendritic
segments were averaged to calculate the mean fluorescence intensity of the dendrite. Spine size changes were calculated as fractional changes with respect to the average size prior to cl-LTP induction. Final spine size changes (at the end of the experiment) was
defined as the average change observed at 30 and 45 min.
All visible spines along the dendrite of interest were analyzed. Spines that were only partially visible or those located within 2 mm of
the photolysis spot were excluded from the analysis. Spines that showed substantial fluctuations (> 25%) during the baseline were
also excluded, although these were rare (< 5% of all spines analyzed). For subgroup analysis of spines, such as distance bins, spine
structural changes were averaged within experiments before averaging across experiments.

Cell Reports 34, 108693, January 26, 2021 e3

ll
OPEN ACCESS

Article

Spine calcium imaging
Hippocampal organotypic slices were continually perfused with aCSF as described above, with the addition of 4 mM MNI-glutamate.
Individual CA1 pyramidal neurons were patched with intracellular solution containing (in mM) 100 K-gluconate, 17 KCl, 5 NaCl,
5 MgCl2, 10 HEPES, 0.5 EGTA, 4 K2ATP, and 0.5 Na3GTP (pH 7.3, osmolarity adjusted to 280 mOsm). On the day of the experiment
200 mM Fluo-5F and 30 mM Alexa 594 were added to the intracellular solution. After allowing cells to fill for 15 min, dendritic segments
were chosen for imaging. Dendrites were imaged at 2 Hz using 800 nm laser light during which dendrites received one of three stimuli:
depolarization alone, MNI-glutamate photolysis alone, or a combination of both depolarization and MNI-glutamate photolysis. For
depolarization, the cell was held for 60 s at 0 mV. For glutamate photolysis, 303 4 ms pulses of 720 nm laser light were delivered
0.5 mm from spine heads at 0.5 Hz. For the depolarization alone condition, sham homosynapses were selected prior to imaging.
Images were analyzed in ImageJ. Individual spines were overlaid with ROIs, background subtracted, and Fluo-5F intensity was
measured for every time point throughout the experiment. Spine calcium levels were then normalized to the pre-stimulation baseline.
Peak calcium was defined as the average of the first 10 images immediately following stimulation. Plateau calcium was defined
as the average of the stable calcium transient that followed the peak (70 images, beginning 20 s after stimulation onset). Cells
were discarded if calcium levels did not return to pre-stimulation levels.
Electrophysiological recordings in acute slices
CA1 pyramidal cells located 50-100 mm below the slice surface were targeted for patch clamp recordings. For visualization of the
dendritic arbour, 200 mM AlexaFluor 488 (AF488) was included in the internal solution and cells were bolus-loaded for 5 min and
patched-off by gently retracting the patch pipette. The fluorescence signal was captured on a Leica DMLFSA confocal microscope
with a 63x water-immersion objective (NA = 0.9, HCX APO L 63x/0.9 W U-V-I, Leica) using the Leica Confocal Software (Version 2.61).
Neurons were stimulated using two monopolar tungsten electrodes placed directly into the tissue. In order to stimulate sets of synapses in close spatial proximity along individual segments of dendrite, stimulation electrodes were positioned within 10-20 mm of a
visually-identified segment of dendrite, 10-40 mm apart from each other. The position of the stimulation electrode could either be visualized via autofluorescence or as shadows in cases where tissue autofluorescence was strong. After positioning, a 10 min wait was
included for the stimulation electrodes to reach thermal equilibrium in order to minimize spatial drift. Cells were subsequently repatched with internal solution containing 200 mM AF488. The probability of successful re-patch was > 90%. Re-patch of the correct
cell was confirmed by fluorescence imaging.
A cross-facilitation test was used to ensure stimulation of two independent pathways. The two pathways were stimulated sequentially with 70 ms inter-pulse interval followed by stimulation in the reverse order, repeated 10 times. A lack of detection of short-term
facilitation of synaptic responses to the second pulse was used as a criterion for independence of the two pathways. The magnitude
of overlap was estimated by comparing the fractional difference in EPSP slope of a given pathway when stimulated as first or second
pulse. Significant overlap was rarely observed, which was also supported by the lack of post-tetanic potentiation in the unstimulated
pathway.
EPSPs were evoked at 0.06 Hz by 100 ms current pulses (20-80 mA), alternated between the stimulation electrodes. PPR was
measured every fifth sweep. A stable baseline was recorded for a maximum of 8 min in order to reduce potential confounds due
to dialysis of cytosolic contents. Experiments showing consistent run-up or run-down were discarded. LTP was induced in one
pathway (randomly chosen) by pairing stimulation with strong postsynaptic depolarization via current injection, 603 at 5 Hz, as previously reported (Padamsey et al., 2017). The unstimulated pathway was kept silent during the induction. EPSPs were recorded for at
least 30 min after the induction and PPRs were sampled throughout the recording period.
Pharmacology
CaMKII was blocked by bath application of 10 mM of the non-competitive inhibitor KN62 (Tocris). Due to its low solubility in water,
KN62 was first dissolved in DMSO and subsequently diluted in aCSF heated to 37 C. Calcineurin was blocked by bath application
of 2 mM FK506 (Tocris). Both, KN62 and FK506 were maintained in the aCSF throughout the experiment. NO signaling was blocked by
bath application of 100 mM L-NAME (Cayman Chemical), a nitric oxide synthase inhibitor. L-NAME was pre-incubated for at least 1 h
before the experiment to abolish tonic levels of NO. NMDA receptor signaling was blocked by local perfusion of 500 mM AP5 (Abcam)
diluted in Tyrode’s solution (additionally containing 10 mM MNI-glutamate for glutamate photolysis).
QUANTIFICATION AND STATISTICAL ANALYSIS
Non-parametric tests were used for all statistical comparisons. Mann-Whitney U test was used for comparisons of independent
means; Wilcoxon signed-rank test was used for dependent means. For multiple comparisons, Kruskal-Wallis H test was used followed by post hoc Dunn’s test for pairwise comparisons. We report the coefficient of determination (r2) for linear regression and
further evaluated linear correlations using the Pearson correlation coefficient (r). Sample sizes (n) are indicated in the corresponding
sections of the main text and refer to the number of experiments (note that each slice was used for a single experiment only). All measurements are stated as mean ± standard error of the mean (SEM). Significance is denoted as follows: *p < 0.05, **p < 0.01.

e4 Cell Reports 34, 108693, January 26, 2021

