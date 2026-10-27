7RZDUGVDXQLÀHGXQGHUVWDQGLQJRIV\QDSWLFSODVWLFLW\
     SDUVLPRQLRXVPRGHOLQJDQGVLPXODWLRQRI
        WKHJOXWDPDWHUJLFV\QDSVHOLIHF\FOH



                        THÈSE NO 8186 (2018)
                       PRÉSENTÉE LE 12 JANVIER 2018
                   À LA FACULTÉ DES SCIENCES DE LA VIE
                              PROJET BLUEBRAIN
                PROGRAMME DOCTORAL EN NEUROSCIENCES



       ÉCOLE POLYTECHNIQUE FÉDÉRALE DE LAUSANNE

           POUR L'OBTENTION DU GRADE DE DOCTEUR ÈS SCIENCES




                                       PAR


                       Giuseppe CHINDEMI




                         acceptée sur proposition du jury:

                     Prof. K. Hess Bellwald, présidente du jury
               Prof. H. Markram, Dr E. B. Muller, directeurs de thèse
                            Prof. A. Holtmaat, rapporteur
                            Prof. M. Tsodyks, rapporteur
                            Prof. C. Petersen, rapporteur




                                      Suisse
                                       2018
Acknowledgements
To Henry Markram and Eilif Muller, for supervising this work with passion and enthu-
siasm. Without your contribution nothing of this could have been done. Henry, it was
an honor for me to be a student of yours. Your vision has inspired this work more than
anything else. Nothing makes me prouder than to have contributed, for what I could,
to what you started. Eilif, you are the best co-supervisor I could have dreamed about.
I am grateful for every moment spent with you talking about science, for every little or
big piece of advice you gave me, for your support during the most diﬃcult moments of
this experience. Words will never be enough to thank you.
  To Michael Graupner, for all his help during the design of the functional plasticity
model. To Christian Rossert, for designing, parameterizing and implementing the Volt-
age Dependent Calcium Channels (VDCC) and calcium dynamics models. To James
King, for designing and implementing most of the simulator extensions required by this
work and for his constant support at all stages of this project. To Michael Reimann, for
designing the thalamocortical stimulus used for network simulations, his help to analyze
the results and precious comments to this manuscript. To Werner Van Geit, for helping
with model ﬁtting in BluePyOpt, fruitful discussion on optimization techniques and all
the good laughs. To Srikanth Ramaswamy for his help with the extension of the synaptic
transmission model. To Max Nolte and Taylor Newton, my companions in misfortune,
for their support. To Rodrigo Perin, for being my guide in the world of in vitro exper-
iments. To Pramod Kumbhar, for all his help with the simulation infrastructure and
the resignation threats. To Fabien Delalondre, for his full support to this project. To
Michael Hines, for his priceless help with NEURON and the stimulating discussions. To
Francesco Cremonesi and Francesco Casalegno, for helping with the simpliﬁcation and
numerical integration of several synapse model components. To Alexandre Beche, for
his precious support over the last years. To Juan Hernando, for making the Viz world
a better place. To Caroline Ingoglia for fruitful discussion and for her help with the
manuscript. To all my colleagues, you have always been the best part of this job.
  To Cani Malati, the worst pack of animals I have ever seen: Tallero, Palermo, Scheg-
gia, Awara, Porco, Crotalo, Manolo, Marcella and Marcellina. To Jota and MP, for
being more than brothers to me. To Chupi and Lopo, for our friendship and all the
unforgettable moments. To Perazza and Phil, for sharing with me these beautiful years.
To all my good friends.
  To my family, for always getting the worst of me and giving back their best anyway.
There was no accomplishment in my life that did not come from your sacriﬁces. From
the bottom of my heart, thank you.
  To Laura, for everything...



                                                                                      3
Summary
Understanding how learning and memory formation work in the brain is a major chal-
lenge in neuroscience, with important implications for many other ﬁelds, including
medicine and industry. It is nowadays widely accepted that synaptic plasticity is the
biological foundation of these higher order brain functions. So far, many diﬀerent plas-
tic behaviors have been intensively studied and characterized, leading to the deﬁnition
of several forms of plasticity: structural, functional, homeostatic, inhibitory, and many
others. Unfortunately, despite all the interest and eﬀorts of the scientiﬁc community, a
complete and consistent understanding of synaptic plasticity is still lacking.
   The main goal of this study is to unify data and theories on synaptic plasticity in a com-
prehensive model, suitable for studying learning and memory down to the synapse level.
To reach our objective, we identiﬁed a minimal set of biological mechanisms responsible
for plastic dynamics and integrated them into a single synapse model, relying whenever
possible on well accepted sub-models from literature. We designed a data-driven ﬁtting
and generalization strategy to parameterize all excitatory-to-excitatory synapses in a
large scale reconstruction of neocortical tissue [Markram et al., 2015]. Finally, we tested
the eﬀects of functional synaptic plasticity on neural circuits in simulations.
   Our model was able to capture not only the outcome of Spike Timing Dependent Plas-
ticity (STDP) protocols used in the training phase, but also the one of all others available
in the same experimental dataset [Markram et al., 1997b]. Moreover, it correctly repro-
duced results on distance-dependent synaptic plasticity [Sjöström and Häusser, 2006],
even though this dataset was never used during ﬁtting. In network simulations, we ob-
served the emergence of a self-regulatory homeostatic mechanism, preventing runaway
excitation. Furthermore, we noticed the strengthening of connections between neurons
that are similarly innervated, as previously shown in vitro [Perin et al., 2011].

Keywords: synaptic plasticity, rewiring, Long Term Potentiation (LTP), Long Term
Depression (LTD), neural networks, learning, memory, neocortex, simulation, neuro-
science, Blue Brain Project (BBP).




                                                                                           5
Sommario
Capire come l’apprendimento e la formazione della memoria lavorino nel cervello è uno
dei principali obiettivi della neuroscienza, con importanti implicazioni per molti altri
settori, tra cui la medicina e l’industria. È oggi ampiamente accettato che la plasticità
sinaptica sia la base biologica di queste alte funzioni cerebrali. Finora sono stati studiati e
caratterizzati diversi comportamenti plastici, portando alla deﬁnizione di svariate forme
di plasticità: strutturale, funzionale, omeostatica, inibitoria e molti altri. Purtroppo,
nonostante il grande interesse e gli sforzi della comunità scientiﬁca, non esiste ancora
una visione completa e coerente della plasticità sinaptica.
   L’obiettivo principale di questo studio è quello di uniﬁcare i dati e le teorie sulla pla-
sticità sinaptica in un modello completo, idoneo allo studio dell’apprendimento e della
memoria ﬁno al livello delle singole sinapsi. Per raggiungere il nostro obiettivo, abbia-
mo individuato un insieme minimo di meccanismi biologici responsabili delle dinamiche
plastiche e li abbiamo integrati in un unico modello di sinapsi, aﬃdandoci, ove possibile,
a risorse ben note in letteratura. Abbiamo progettato una strategia di adattamento e
generalizzazione per parametrizzare tutte le sinapsi eccitatorio-eccitatorie in una rico-
struzione di tessuto neocorticale su larga scala [Markram et al., 2015]. Inﬁne abbiamo
testato gli eﬀetti della plasticità sinaptica funzionale sui circuiti neurali in simulazione.
   Il nostro modello è stato in grado di catturare non solo il risultato dei protocolli Spike
Timing Dependent Plasticity (STDP) utilizzati nella fase di parametrizzazione, ma anche
di tutti gli altri disponibili nello stesso set di dati sperimentali [Markram et al., 1997b].
Inoltre, ha riprodotto correttamente i risultati sulla plasticità sinaptica dipendente dalla
distanza [Sjöström and Häusser, 2006], nonostante questo set di dati non fosse mai stato
utilizzato durante la parametrizzazione. Nelle simulazioni di rete, abbiamo osservato
l’emergere di un meccanismo di autoregolazione omeostatico, utile ad impedire il sovra-
eccitamento. Inoltre, abbiamo notato il raﬀorzamento delle connessioni tra neuroni
innervati in modo simile, come precedentemente mostrato in vitro [Perin et al., 2011].

Parole chiave: plasticità sinaptica, rewiring, Long Term Potentiation (LTP), Long
Term Depression (LTD), reti neurali, apprendimento, memoria, neocorteccia, simulazio-
ne, neuroscienze, Blue Brain Project (BBP).




                                                                                             7
Contents
1 Introduction                                                                            15

2 Physiology, anatomy and life cycle of synaptic connections                             19
  2.1 Structure and function of excitatory synapses . . . . . . . . . . . . . . . . 19
  2.2 Synapse life cycle . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 21

3 Modeling and parameterization of excitatory synapses                                    29
  3.1 The GluSynapse model . . . . . . . . . . . . . . . . . . . . . . . . . . . . .      30
      3.1.1 AMPAR and NMDAR . . . . . . . . . . . . . . . . . . . . . . . . .             31
      3.1.2 Vesicle release . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .   32
      3.1.3 VDCC . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .      34
      3.1.4 Spine . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .     35
      3.1.5 Postsynaptic calcium . . . . . . . . . . . . . . . . . . . . . . . . . .      35
      3.1.6 Long-term plasticity . . . . . . . . . . . . . . . . . . . . . . . . . .      38
      3.1.7 Rewiring . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .      43
  3.2 Fitting the long-term plasticity model . . . . . . . . . . . . . . . . . . . .      47
      3.2.1 Problem statement and initial conditions . . . . . . . . . . . . . .          48
      3.2.2 Pair selection . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .    52
      3.2.3 Parameter ﬁtting . . . . . . . . . . . . . . . . . . . . . . . . . . . .      54
      3.2.4 Pathway generalization . . . . . . . . . . . . . . . . . . . . . . . .        57
      3.2.5 Circuit generalization . . . . . . . . . . . . . . . . . . . . . . . . .      60
  3.3 Validation and predictions . . . . . . . . . . . . . . . . . . . . . . . . . . .    60
      3.3.1 Distance-dependent synaptic plasticity . . . . . . . . . . . . . . . .        61
      3.3.2 STDP map . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .        64
  3.4 NMODL implementation and HPC considerations . . . . . . . . . . . . .               67
      3.4.1 Spine as virtual compartment . . . . . . . . . . . . . . . . . . . . .        67
      3.4.2 GluSynapse integration methods . . . . . . . . . . . . . . . . . . .          67
      3.4.3 Inhomogeneous Poisson elimination . . . . . . . . . . . . . . . . . .         68

4 Eﬀects of functional synaptic plasticity on neural circuits                             71
  4.1 Fundamentals of microcircuit reconstruction and simulation . . . . . . . .          71
      4.1.1 Thalamocortical stimulation . . . . . . . . . . . . . . . . . . . . . .       72
  4.2 Preliminary analyses . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .    74
      4.2.1 Network balance . . . . . . . . . . . . . . . . . . . . . . . . . . . .       74
      4.2.2 Redistribution of synaptic eﬃcacy . . . . . . . . . . . . . . . . . .         77
      4.2.3 Layer speciﬁc changes . . . . . . . . . . . . . . . . . . . . . . . . .       80
  4.3 Common neighbors strength bias . . . . . . . . . . . . . . . . . . . . . . .        81




                                                                                           9
Contents

5 Conclusion                                                                            85
  5.1 Limitations . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 87
  5.2 Future work . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 87

Bibliography                                                                           90

Acronyms                                                                              101

Glossary                                                                              103




10
List of Figures
 2.1   The synapse life cycle . . . . . . . . . . . . . . . . . . . . . . . . . . . . .     22
 2.2   Ex novo synapse formation . . . . . . . . . . . . . . . . . . . . . . . . . .        23
 2.3   Multi synapse bouton formation . . . . . . . . . . . . . . . . . . . . . . .         24
 2.4   Long-term functional changes . . . . . . . . . . . . . . . . . . . . . . . . .       25
 2.5   Synapse elimination . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .      27

 3.1 GluSynapse logical structure. . . . . . . . . . . . . . . . . . . . . . . . . .        31
 3.2 Fractional NMDAR-mediated calcium current . . . . . . . . . . . . . . . .              37
 3.3 Calcium thresholding problem . . . . . . . . . . . . . . . . . . . . . . . . .         41
 3.4 Calcium thresholding solution . . . . . . . . . . . . . . . . . . . . . . . . .        42
 3.5 Rewiring model. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .        44
 3.6 Potential synaptic locations in the BBP’s reconstruction. . . . . . . . . . .          45
 3.7 Rewiring ﬁt. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .     47
                                    (d)       (p)
 3.8 Set of valid solutions for USE and USE . . . . . . . . . . . . . . . . . . . .         50
                 (initial)
 3.9 Bimodal USE           validation . . . . . . . . . . . . . . . . . . . . . . . . . .   51
 3.10 Exemplar postsynaptic cell used for ﬁtting. . . . . . . . . . . . . . . . . .         52
 3.11 Geometrical constraints on cell selection. . . . . . . . . . . . . . . . . . . .      53
 3.12 Eﬀects of STDP protocols on a single synapse in silico . . . . . . . . . . .          55
 3.13 LTP/LTD model ﬁtting results . . . . . . . . . . . . . . . . . . . . . . . .          56
 3.14 Pathway generalization results . . . . . . . . . . . . . . . . . . . . . . . .        58
 3.15 Classical STDP experiment in silico . . . . . . . . . . . . . . . . . . . . .         59
 3.16 Example of distance-dependent synaptic plasticity . . . . . . . . . . . . .           62
 3.17 Summary of distance-dependent synaptic plasticity . . . . . . . . . . . . .           63
 3.18 Distance-dependent synaptic changes. . . . . . . . . . . . . . . . . . . . .          63
 3.19 STDP map . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .        65
 3.20 LTP frequency dependence for uncharacterized pathways . . . . . . . . . .             66

 4.1 Circuit layout . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .     72
 4.2 Exemplar stimulus . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .        73
 4.3 Summary of all stimuli . . . . . . . . . . . . . . . . . . . . . . . . . . . . .       74
 4.4 Firing rate evolution during spontaneous network activity . . . . . . . . .            75
 4.5 Firing rate evolution during stimulated network activity . . . . . . . . . .           76
 4.6 Burst incidence per stimulus . . . . . . . . . . . . . . . . . . . . . . . . . .       77
 4.7 Redistribution of synaptic eﬃcacy during spontaneous activity . . . . . .              78
 4.8 Redistribution of synaptic eﬃcacy during thalamocortical stimulation . .               79
 4.9 Synaptic eﬃcacy distribution of a typical layer 5 TTPC . . . . . . . . . .             80
 4.10 Average release probability change after 10 minutes of simulation. . . . .            81



                                                                                            11
List of Figures

     4.11 Emergence of common neighbors strength bias (talamocortical stimulation) 83
     4.12 Emergence of common neighbors strength bias (spontaneous activity) . . 84




12
List of Tables
 3.1   Parameters of the long term plasticity model . . . . . . . . . . . . . . . . 48
 3.2   Fixed synaptic parameters. . . . . . . . . . . . . . . . . . . . . . . . . . . 51




                                                                                     13
1 Introduction
                                                       “Learning never exhausts the mind.”
                                                                       Leonardo Da Vinci

Learning and memory are undoubtedly among the most critical and fascinating brain
functions. Indeed complex organisms such as mammals would not survive in the external
world without the capability to acquire new knowledge and to retrieve it whenever needed
(i.e. to avoid predators). While this simple argument surely explains why learning and
memory are critical, the fascinating aspect deserves deeper analysis. Let us consider
the brain as a dynamical learning system, an engineering artifact. To begin with, our
system is constantly bombarded by very heterogeneous stimuli from the external world:
images, sounds, odors and so on. Of all the information received, very little matters.
The important bits need to be identiﬁed, processed, stored and, most importantly, linked
to existing knowledge, possibly of diﬀerent nature (i.e. the sound of rain and the need
of an umbrella). Furthermore, memory capacity is arguably limited and so it has to be
managed in an eﬀective way. That is, there must be space for new memories, but impor-
tant ones have to be kept for the entire lifetime of the organism. Each step of this chain
poses fundamental challenges to the design of such a system. Most of these challenges
are nowadays unsolved and they represent active areas of research in engineering and
science. If we now consider that not only the brain actually solves all of them, but also
it does that in an eﬃcient manner, we are left with one obvious question: how?
   Understanding how learning and memory formation work in the brain would have
important implications in many ﬁelds, including Medicine. Research on memory-related
pathologies, such as the Alzheimer disease, could greatly beneﬁt from a deeper under-
standing of information storage principles. Furthermore, the subject has also captured
the interest of Industry, due to the recent successes of Artiﬁcial Neural Networks (ANN)
and Deep Learning in many commercial applications, such as image and speech recog-
nition [LeCun et al., 2015]. Although ANN are only loosely inspired to biological neural
networks, certain features of biological learning represent the holy grail of modern Ar-
tiﬁcial Intelligence (AI). In particular, no AI-based system has ever been able to show
strong generalization capabilities. That is, artiﬁcial systems have limited ability to adapt
to unforeseen situations, as they are usually trained to perform a given set of tasks and
any variation to the set requires re-training. Brain research is not the only way to look
into these issues, but might represent the missing viewpoint to allow the design of next
generation AIs.
   The ﬁrst step towards a better understanding of how learning and memory work
is obviously to identify what are their biological building blocks. It is nowadays widely
accepted that synaptic plasticity, the main subject of this work, is the biological substrate



                                                                                          15
1 Introduction

of these higher order brain functions. A (chemical) synapse is a contact point between
two neurons allowing the exchange of information1 . Synapses are dynamic structures:
they can form and dissolve, but also adapt their responses to neural activity and other
environmental conditions. This set of capabilities is collectively referred to as synaptic
plasticity. As synapses are one of the most important factors determining the response
of neurons to stimuli and their behavior can change over time, they represent the ideal
foundation for learning and memory.
   So far many diﬀerent plastic behaviors have been intensively studied and characterized,
leading to the deﬁnition of several “forms” of plasticity, such as:

     • Functional plasticity, the capability of neurons to modify the strength of existing
       connections. This form of adaptation have been the subject of many modeling
       studies in recent years.

     • Structural plasticity, identiﬁed by a group of major changes of the brain mi-
       crocircuit, ranging from creation or elimination of synapses to axonal or dendritic
       remodeling.

     • Homeostatic plasticity, referring to a set of slow mechanisms regulating and
       maintaining the operating regimes of the neurons, despite ongoing plastic changes
       and turnover of the molecular substrate.

Regardless of this taxonomy, devised for technical and historical reasons, synaptic plas-
ticity is nowadays considered an integrated phenomenon. Already in 1999, it was shown
that persistent strengthening of synaptic connection, namely Long Term Potentiation
(LTP), and formation of new synapses are strongly correlated [Toni et al., 1999]. Years
later, also weakening of synaptic connections, or Long Term Depression (LTD), was
found to be accompanied by important structural changes, synapse elimination in this
case [Bastrikova et al., 2008]. The continuous accumulation of such evidences challenges
the idea of a sharp distinction between forms of plasticity. The plethora of plastic changes
observed in the brain may be the product of a single biological machinery, operating over
diﬀerent time scales.
   Unfortunately, despite all the interest and eﬀorts of the scientiﬁc community, a com-
plete and consistent understanding of synaptic plasticity is lacking. In particular, there
is still no consensus on what the most relevant plastic mechanisms are and how they are
orchestrated together to actually enable learning.
   The main diﬃculty to establish such a unifying theory of synaptic plasticity is the
sparsity of available data. Lab experiments involving synaptic plasticity are inherently
very complex. First of all, individual synapses are diﬃcult to observe and even more
diﬃcult to control. For this reason, establishing any cause-eﬀect link between network
dynamics and synaptic variables is hardly feasible. Second, plastic mechanisms span

 1
     In the classic view of synaptic transmission, information ﬂows in one direction only, from the pre-
     synaptic neuron (the sender) to the post-synaptic one (the receiver). However, it is now known that
     synapses are capable of retrograde signaling and indeed these mechanisms are crucial for synaptic
     plasticity.




16
multiple time scales, from a millisecond to several days. Such a large range makes it ex-
tremely hard to understand how plastic mechanisms inﬂuence each other. Furthermore,
experimental design is particularity diﬃcult for the slowest time scales. This is indeed
the main reason why plastic changes happening on diﬀerent time scales were tradition-
ally described as separate forms of plasticity. Last but not least, the heterogeneity of
synaptic pathways suggests that the number of experiments required to quantitatively
map the whole spectrum of plastic changes is currently beyond our reach.
   The modeling approach has been proposed as a viable solution to overcome all the
experimental limitations. Unfortunately, the lack of data has also limited our capability
to build global theories on synaptic plasticity with a solid empirical foundation. As a
consequence, most available models focus on very speciﬁc aspects only, such as LTP or
homeostatic plasticity.
   The main goal of this study was to unify data and theories on synaptic plasticity
in a comprehensive model, suitable for studying learning and memory in cortical net-
works, down to the neuronal and synaptic level. We decided to restrict our scope to
connections between excitatory neurons in the somatosensory cortex of the rat brain,
as this system oﬀers the largest datasets. Furthermore, connections between excitatory
neurons accounts for 80% of synapses and almost %100 percent of spiny ones, which are
thought to encapsulate the rich molecular machinery for plasticity. Particular attention
has been payed to the identiﬁcation of a minimal set of key plastic mechanisms and to
the characterization of their interactions. The simulation framework, developed to com-
plement this work, allows the observation of learning dynamics at the synapse level (i.e.
LTP/LTD or synaptogenesis), neuron level (i.e. N-methyl-D-aspartate receptor (NM-
DAR) spikes or synapse cooperativity) and network level (i.e. emergence of connectivity
motifs and structure-function relationships). For this reason, we propose this resource
as a new tool for investigating learning and memory dynamics in the brain. Our hope
is that this work will trigger more studies attempting to integrate data and models in
the ﬁeld of plasticity research. Moreover, we would value as a major achievement of this
thesis the adoption of our simulation framework by other scientists to help the design of
targeted in vitro and in vivo experiments.
   To reach our objective, we identiﬁed a minimal set of biological mechanisms respon-
sible for plastic dynamics and integrated them into a single synapse model, building
from well accepted sub-models in literature. We then designed a two-step procedure for
parameterizing the synapse model for any synaptic pathway. In the ﬁrst step we use
a custom Genetic Algorithm (GA) for ﬁtting the parameters of experimentally charac-
terized pathways to data. The second step generalizes the results of the previous, pro-
ducing a predicted parametrization for uncharacterized pathways. The procedure was
ﬁnally used to parameterize the plastic behavior of all 4 billion excitatory synapses in
a large-scale computer reconstruction of neocortical tissue [Markram et al., 2015]. The
reconstructed microcircuit was than simulated on the DOE ANL MIRA and BBPIV
IBM BlueGene/Q supercomputers. To allow the parametrization and simulation of such
a large piece of digital brain tissue, several technological challenges have been faced.
The methods herein developed to solve them have, in principle, broad applicability in
the ﬁelds of optimization, High Performance Computing (HPC) and simulation-based



                                                                                      17
1 Introduction

science.
   Over the course of this study, several artifacts have been produced. In particular, we
designed and implemented:

      i High-level model of the synapse lifecycle.

     ii NEURON synapse model, and the corresponding parameterization, suitable to
        investigate synaptic plasticity in large scale simulations

     iii Method to ﬁt and generalize synaptic parameters from sparse data

     iv Framework for synaptic plasticity modeling and simulation

   The rest of this dissertation follows a bottom-up approach: we build from the experi-
mental evidence on plastic mechanisms to the implications of synaptic plasticity on the
structure and function of digitally reconstructed neural circuits. In Chapter 2 we sum-
marize the key experimental results on synaptic plasticity and propose a parsimonious
model of the synapse life-cycle. The goal of this chapter is to provide a short review of
some key experimental results on synaptic dynamics, later used for model design. Chap-
ter 3 describes how our synaptic model was developed, parametrized on available data
and generalized to cover uncharacterized synaptic pathways. Chapter 4 presents the re-
sults of our ﬁrst large-scale simulations of plastic neocortical networks. Here we describe
the main eﬀects of functional synaptic plasticity on neural circuits. We investigated the
behavior of the network under two conditions: spontaneous activity and thalamocortical
stimulation. Finally in chapter 5 we draw our conclusions, identifying the limitations of
the present study and proposing future directions of research on the subject.




18
2 Physiology, anatomy and life cycle of
  synaptic connections
                                                     “In the case of all things which have
                                                     several parts and in which the totality
                                                     is not, as it were, a mere heap, but
                                                     the whole is something beside the
                                                     parts, there is a cause.”
                                                                  Aristotle, Metaphysics

The main objective of this thesis was to unify data and theories on synaptic plasticity in
a comprehensive synapse model. A systematic understanding of the diﬀerent synaptic
states and transformations was then a fundamental prerequisite for this work. Unfortu-
nately, such a complete overview has not been proposed yet. To move the ﬁrst steps in
this direction, we derived a summarizing description of what we can call the synapse life
cycle from the vast body of literature on synaptic plasticity.
   In this chapter we present our view of the synapse life cycle and set the founda-
tions for the synapse model proposed in chapter 3. We start by brieﬂy introducing the
main “components” of excitatory synapses and how they are orchestrated together to
implement the synaptic functions (Sec. 2.1). We then focus on the characterization
of the diﬀerent synaptic states and their relationships, as emerging from our review of
experimental literature (Sec. 2.2).


2.1 Structure and function of excitatory synapses
Synapses are biological electro-chemical transducers, required for information transmis-
sion and decoding in the brain. That is, they not only allow the propagation of electric
signals from a neuron to another (transmission), but also aﬀect the information content
of the signal (decoding), as we will discuss later in this section.
   Excitatory synapses consist of two parts, a presynaptic terminal (or bouton) and a
postsynaptic specialization, bound together by a dense network of extracellular proteins
in the so called synaptic cleft. The basic functioning of an excitatory synapse can be
roughly summarized as follows. When an action potential (AP) is ﬁred by the presynap-
tic neuron, the bouton might release one or more vesicles, containing neurotransmitter.
Neurotransmitter molecules diﬀuse then in the synaptic cleft and activate the appropriate
receptors, sitting on the postsynaptic side of the synapse. The activation of postsynap-
tic receptors causes an inward current which depolarizes the cell membrane, eﬀectively
propagating the electric signal form the presynaptic neuron to the postsynaptic one.



                                                                                         19
2 Physiology, anatomy and life cycle of synaptic connections

   The magnitude and time course of the postsynaptic depolarization depend on several
pre- and post-synaptic factors, such as: AP timing, location of the synapse on the
dendritic arborization, extracellular ionic conditions, number of vesicles released and
receptors available. For this reason, a synapse is not just propagating an electric pulse
from a neuron to another, but it is actually “translating” a relatively poor signal, the
AP, into a much richer representation, the postsynaptic depolarization. Considering
that many of the factors responsible of signal decoding (i.e. the number of released
vesicles) are actively adjusted by plastic mechanisms, we can reformulate the deﬁnition
of synaptic plasticity as the capability to dynamically adjust the readout of pre- and
postsynaptic state.
   Let us now focus on the anatomy of a synapse, starting from the presynaptic side. As
mentioned earlier, the bouton contains one or more vesicles, loaded with neurotransmit-
ter. The most prominent neurotransmitter of excitatory cells is glutamate 1 . Vesicles
ready to be released (or “docked”) accumulate at the synaptic surface of the bouton,
called the Active Zone (AZ). The machinery responsible for vesicle release is quite com-
plex and still not completely understood, but in ﬁrst approximation it consists of a set of
Voltage Dependent Calcium Channels (VDCC), coupled with SNAP (Soluble NSF At-
tachment Protein) REceptor (SNARE) proteins. When an AP invades the bouton, the
calcium inﬂux through the VDCC activate the SNARE complex, causing vesicle fusion
with the cell membrane and the consequent release of neurotransmitter. Vesicle release
is a stochastic process. An AP may or may not trigger a release event, depending on
several conditions, such as the availability of docked vesicles. The number of vesicles
in the AZ and their probability of release are subject to short- and long-term forms of
synaptic plasticity.
   The vast majority of excitatory synapses are located on dendritic spines. Spines are
small membrane protrusions, usually shorter than 2 µm and having a volume of approx-
imatively 1 fL. They come in diﬀerent sizes and shapes, ranging from thin cylindric
formations to “mushroom” like structures, with a clear neck and head. With respect
to synaptic plasticity, spines are extremely important for two reasons. First, they al-
low speciﬁcity of plastic changes through biochemical compartmentalization of signaling
molecules, such as calcium ions. That is, several ion species cannot easily escape the
spine, limiting the signal they convey to the synapses in situ. Second, they allow to
bridge the distance between a presynaptic axon and a postsynaptic dendrite, tremen-
dously extending the number of possible synaptic partners for a given neuron.
   Most of the postsynaptic receptors are located on the spine head, right in front
of the presynaptic AZ, packed in the so called Postsynaptic Density (PSD). The two
most important receptors at excitatory synapses are alpha-amino-3-hydroxy-5-methyl-
4-isoxazole propionate receptors (AMPARs) and N-methyl-D-aspartate receptors (NM-
DARs). AMPARs are the largest population of receptors and contribute most of the
 1
     Strictly speaking no neurotransmitter is excitatory nor inhibitory. The eﬀect of neurotransmitter re-
      lease on the postsynaptic neuron depends on the type of receptors binding to it. However, considering
      that a cell can only produce a single type of neurotransmitter and that the number of diﬀerent neu-
      rotransmitters and binding receptors in the brain is limited, the terms excitatory cell and inhibitory
      cell are quite common.




20
                                                                    2.2 Synapse life cycle

synaptic current under normal conditions. They are permeable to sodium and potassium
ions, but not calcium. Insertion and removal of AMPARs is one of the most powerful, and
noticeable, plastic change happening at excitatory synapses. NMDAR are undoubtedly
the single most important protein for synaptic plasticity. At resting membrane poten-
tial, the NMDAR is blocked by magnesium ions. That is, even if bound to glutamate,
the receptor does not carry any current. However, at depolarized potentials, the block
is removed allowing the ﬂux of several ion species, calcium in particular. Considering
that postsynaptic depolarization can be due to both a backpropagating action potential
(bAP) or a local electrical event such as the activation of a cluster of nearby synapses,
the NMDAR is considered a coincidence detector of pre- and post-synaptic activity. In
other words, two conditions are required to obtain an NMDAR-mediated current: the
release of glutamate from the presynaptic terminal and the depolarization of the postsy-
naptic membrane. Furthermore, the order of these two events is non-commutative. Due
to the latencies of the chemical and electrical dynamics of synaptic transmission, the
release of neurotransmitter shortly followed by postsynaptic depolarization produces a
much larger current than the reverse order of events. This property is extremely inter-
esting from the computational point of view because, thanks to the NMDAR, a synapse
can distinguish a causal relationship between the activity of the pre- and post-synaptic
neuron form an anti-causal one.
   Receptors are not the only channels inserted in the spine membrane. Voltage de-
pendent ion channels are a major class of proteins related to several aspects of neurons
physiology, from the electrical behavior to synaptic plasticity [Hille et al., 2001]. VDCCs
are a large family of ion channels, responsible for a large fraction of the total calcium
current in spines. Depending on their kinetic, VDCC can be classiﬁed as High Volt-
age Activated (HVA) and Low Voltage Activated (LVA) channels. HVA channels are
particularly relevant for several forms of synaptic plasticity, given their activation as a
consequence of large depolarizing events, such as a bAP.
   In this brief overview of the anatomy and physiology of synapses we neglected several
important components, such as metabotropic receptors or internal calcium stores. These
elements are surely important for synapse functioning, but their role with respect to
synaptic plasticity is either poorly understood or thought to be of secondary relevance.
As mentioned earlier in this section, our goal is to distill the most important aspects of
synaptic structure and function to produce a model suitable to study plasticity in large
scale simulations.


2.2 Synapse life cycle
Based on our analysis of experimental literature, we distilled a simpliﬁed description of
the synapse life cycle, shown in ﬁgure 2.1. We start by considering the process of synapse
formation. This subject was extensively reviewed in two seminal publications [Holtmaat
and Svoboda, 2009, Ziv and Garner, 2004] that guided our initial steps and provided the
foundation for our analysis. We can identify two main paths for synaptogenesis: ex novo
synapse formation and Multi Synapse Bouton (MSB) formation [Holtmaat and Svoboda,




                                                                                        21
2 Physiology, anatomy and life cycle of synaptic connections


                            $[RQ
                            'HQGULWH
                                                                        ,,
                                                                    QLQJ                 )LO
                                                                   X                        RS
                                                                 3U                           RG




                                                                                                 LX
                                                                                                   P
                                                                                                    S
                                                                                                      UR
                                                                                                        MHF
                                                                                                           WLR
                     ,




                                               ,
                   J




                                            QJ




                                                                                                              Q
                QLQ




                                       3UXQL
             3UX




                                                              /7
                                                                 3
                                                                           7KH6\QDSVH
                                                                    /7'


                                                                            /LIH&\FOH
                         3UXQLQJ,,




                                                                                                          LR Q
                                             6S
                                               LQ H




                                                                                                          DW
                                                                                                         UP
                                                    I R




                                                                                                     IR
                                                      UP




                                                           DW


                                                                                                Q
                                                                LR                                   R
                                                                     Q                          XW
                                                                                           %R



                                            )LO
                                                R   SRG                             RQ
                                                           L X P  S UR M H F W L


     Figure 2.1 – Simpliﬁed model of the synapse life cycle. We consider six possible axo-
                  dendritic conﬁguration. Each solid arrow represent a plastic transition from
                  one state to another. The dashed arrow identify functional changes at mature
                  synapses, without major structural adaptation.


2009]. Let us consider the former case (Fig. 2.2). The initial requirement for any de novo
synapse is a close axo-dendritic apposition that can be bridged by a spine [Fares and
Stepanyants, 2009]. Spines are preceded by the elongation of ﬁlopodia, thin actin protru-
sions originating from the dendrites [Ziv and Smith, 1996, Kasai et al., 2010]. Filopodia
have an extremely high turnover rate and often they get in contact with nearby axons
before being eliminated [Ziv and Garner, 2004, Ziv and Smith, 1996]. Under particular
circumstances 2 , these axo-dendritic contacts can lead to the ordered maturation of the
presynaptic (ﬁrst) and postsynaptic (later) structures [Ziv and Garner, 2004, Friedman
et al., 2000], and eventually to the formation of the so called “silent synapses”.       A
(postsynaptically) silent synapse, is a synapse lacking a PSD despite a well formed AZ
 2
     The mechanisms behind the stabilization of axo-dendritic contacts are mostly unknown, however, an
      interesting possibility is the involvement of mobile vesicles exocytosis [Ziv and Garner, 2004].




22
                                                                                         2.2 Synapse life cycle




                                                                  ,,
                                                             QLQJ                 )LO
                                                            X                        RS
                                                          3U                           RG




                                                                                          LX
                                                                                            P
                                                                                             S
                                                                                               UR
                                                                                                 MHF
                                                                                                    WLR
              ,




                                        ,
               J




                                                                                                       Q
                                     QJ
            QLQ




                                3UXQL
         3UX




                                                      /7
                                                          3
                                                              /7'      ([1RYR
                                                                       6\QDSVH
                                                                      )RUPDWLRQ
                   3UXQLQJ ,,




                                                                                                   LR Q
                                      6S
                                        LQ H




                                                                                                   DW
                                                                                                  UP
                                             I R




                                                                                              IR
                                               UP




                                                    DW


                                                                                         Q
                                                         LR                                   R
                                                              Q                          XW
                                                                                    %R



                                     )LO
                                         R   SRG                           RQ
                                                    L X P S UR M H F W L


  Figure 2.2 – Ex novo synapse formation. Filopodia can contact axons several microns
               away from their source dendrite. Occasionally, these events can lead to the
               formation of a presynaptic bouton and postsynaptic spine, followed by synapse
               maturation.


and, for this reason, not able to induce an Excitatory Post-Synaptic Potential (EPSP)
[Kullmann, 2011]. However, silent synapses usually possess postsynaptic NMDARs and
this is the key of their activity-dependent maturation, or “unsilencing”[Kullmann, 2011].
The ordered and close activation of the pre and postsynaptic neurons releases the mag-
nesium block from the NMDARs, allowing the ﬂux of ions, including calcium [Kullmann,
2011]. The calcium inﬂux causes the insertion of AMPARs in the postsynaptic mem-
brane and so the formation of a PSD at the new born functional synapse [Kullmann,
2011].
  Another possible pathway for synaptogenesis is what can be deﬁned as “bouton re-
cycling” (Fig. 2.3). It was shown that ﬁlopodia, and spines, tend to project towards
existing boutons and form synapses [Holtmaat and Svoboda, 2009, Knott et al., 2006].
The steps of this model of synaptogenesis are similar to those described above, with



                                                                                                            23
2 Physiology, anatomy and life cycle of synaptic connections




                                                                  ,,
                                                             QLQJ                   )LO
                                                            X                          RS
                                                          3U                             RG




                                                                                            LX
                                                                                              P
                                                                                                    SU
                                                                                                      RMH
                                                                                                         FWL
              ,




                                                                                                            RQ
                                        ,
               J




                                     QJ
            QLQ




                                3UXQL
         3UX




                                                       /7
                                                          3
                                                                           06%
                                                              /7'


                                                                        )RUPDWLRQ
                   3UXQLQJ ,,




                                                                                                       LR Q
                                      6S
                                        LQ H




                                                                                                     DW
                                                                                                    UP
                                             I R




                                                                                                IR
                                               UP




                                                    DW
                                                                                                R

                                                                                           Q
                                                         LR
                                                              Q                            XW
                                                                                      %R



                                     )LO
                                         R   SRG                             RQ
                                                    L X P  S UR M H F W L


  Figure 2.3 – Multi synapse bouton formation. An existing bouton can be contacted by a
               newly formed spine. Weather this conﬁguration can be maintained indeﬁnitely
               or not is currently unknown.


the main diﬀerence that the process ends with the creation of a MSB [Holtmaat and
Svoboda, 2009]. Nevertheless, some interesting insights deserve to be mentioned. First,
glutamate might be responsible of the attraction of ﬁlopodia and spines to existing bou-
tons [Holtmaat and Svoboda, 2009, Richards et al., 2005], supporting the idea that the
exocytosis of mobile vesicles might be crucial for the de novo synapse formation. Sec-
ond, the coexistence of more than one synapse on the same bouton in not necessarily a
stable condition: synapses on the same bouton might compete for presynaptic resources,
causing the pruning of one or more of them [Holtmaat and Svoboda, 2009, Toni et al.,
2007]. Last, MSBs potentially resolve the discrepancies found in the turnover rates of
spines and boutons [Holtmaat and Svoboda, 2009, De Paola et al., 2006, Holtmaat et al.,
2005].
  After the analysis of synaptogenesis, we now focus on functional synapses (Fig. 2.4).
Mature synapses show, among all the other synaptic states, the richest spectrum of



24
                                                                                         2.2 Synapse life cycle




                                                                 ,,
                                                            QLQJ                  )LO
                                                           X                         RS
                                                         3U                            RG




                                                                                          LX
                                                                                            P
                                                                                                  SU
                                                                                                    RMH
                                                                                                       FWL
              ,




                                                                                                          RQ
                                        ,
               J




                                     QJ
            QLQ




                                3UXQL
         3UX




                                                     /7
                                                         3
                                                            /7'     /RQJ7HUP
                                                                     )XQFWLRQDO
                                                                      &KDQJHV
                   3UXQLQJ ,,




                                                                                                     LR Q
                                      6S
                                        LQ H




                                                                                                   DW
                                                                                                  UP
                                             IR




                                                                                              IR
                                              UP




                                                   DW
                                                                                              R

                                                                                         Q
                                                        LR
                                                             Q                           XW
                                                                                    %R



                                     )LO
                                         R   SRG                          RQ
                                                   L X P S UR M H F W L


  Figure 2.4 – Long-term functional changes. Mature synapses can transiently or persis-
               tently adapt their response to presynaptic activation. Persistent changes,
               such as LTP (LTD), do not require substantial structural adaptation, al-
               though they are often accompanied by the enlargement (shrinkage) of spine
               and bouton.


plastic dynamics, spanning from morphological to functional changes. Long Term Po-
tentiation (LTP) and Long Term Depression (LTD) are among the most striking of these
dynamics. The term LTP refers to a persistent increase of the postsynaptic response due
to the activation of the presynaptic ﬁbers [Bliss and Lømo, 1973]. Without going into
the details of the expression and maintenance of LTP (for a review see [Malenka and
Bear, 2004]), we brieﬂy point out that LTP can have a presynaptic locus of expression
[Markram and Tsodyks, 1996] or a postsynaptic one or both [Sjöström et al., 2007].
Postsynaptic forms of LTP usually involve the insertion of AMPARs in the postsynap-
tic membrane (but see [Bliss and Collingridge, 1993]). Remarkably enough, LTP plays
also an important structural role. It was indeed shown that LTP and synaptogenesis are
strongly related [Holtmaat and Svoboda, 2009, Toni et al., 1999], raising the attention of




                                                                                                               25
2 Physiology, anatomy and life cycle of synaptic connections

the scientiﬁc community about the connections between structural and functional plas-
ticity. LTD is the specular process of LTP and can be oversimpliﬁed as a drop in release
probability [Sjöström et al., 2007] or the removal of AMPAR from the membrane, with
the consequent reduction of the synaptic response [Malenka and Bear, 2004]. One of the
most interesting and curious things about LTD is its relationship with synapse elimina-
tion. Since at least one decade, LTD is speculated to be correlated with spine pruning
(for a review see [Segal, 2005]) and, ultimately, to be one of the causes of synapse elim-
ination. However, as remarked by M. Segal, the conjecture of a link between LTD and
spine pruning is mainly due to theoretical reasoning rather than experimental evidence:
the well established connection between spinogenesis and LTP, plus the complementarity
of this last process with LTD, have suggested the logical conclusion that LTD and spine
elimination must be correlated, but nobody has proved it yet [Segal, 2005].
   Despite the disruptive power of plastic mechanisms, morphological variables at func-
tional synapses appear to be tightly correlated. Harris and Stevens [1989] showed that
PSD area is strongly correlated with several pre- and post-synaptic variables, such as
spine head volume, bouton volume and number of vesicles. Furthermore, spine head
volume is also correlated with the number of vesicles [Harris and Stevens, 1989] and
strong correlation between PSD area and spine head volume has been reported across
several cortical and non-cortical areas [Harris and Stevens, 1989, Schikorski and Stevens,
1999, Arellano et al., 2007]. Such ﬁndings suggest the presence of remarkable regulatory
mechanisms at the level of single synapses. As many forms of long-term plasticity have
well deﬁned morphological correlates, we can speculate that eventually any punctual
change will be matched over time by the whole synapse. For example, we could reason-
ably expect that even a purely postsynaptic form of LTP, expressed through the rapid
insertion of AMPAR, will be followed by presynaptic adaptation (i.e. an increase in size
of the Readily-Releasable Pool (RRP)).
   Spine (bouton) pruning is the last structural change needed to produce a consistent
model of the synapse life cycle (Fig. 2.5). In recent years, a large eﬀort has been
devoted to identify the link between LTD and synapse elimination. Several experiments
have shown that there is indeed a correlation between spine retraction and LTD [Nägerl
et al., 2004, Segal, 2005]. Although spine retraction does not necessary imply synapse
elimination [Segal, 2005], the discovery of a connection with LTD prompts us to wonder
about the relationships between functional and structural changes at depressed synapses.
In particular, the exact order of the structural changes leading to synapse elimination is
still a matter of debate, being that a consequence of LTD or not. While for many years
the trend was to consider the retraction of spines the hallmark of synapse elimination,
recent papers demonstrated a substantial involvement of presynaptic and extracellular
structures, leading to the idea that synapse elimination might start with bouton loss and
only later be followed by spine vanishing [Bastrikova et al., 2008, Becker et al., 2008].
For the moment we still cannot completely rule out any of these hypotheses and perhaps
the two dynamics might coexist in the brain.




26
                                                                                        2.2 Synapse life cycle




                                                                ,,
                                                             LQJ                 )LO
                                                          UXQ                       RS
                                                         3                            RG




                                                                                         LX
                                                                                           P
                                                                                                 SU
                                                                                                   RMH
                                                                                                      FWL
               ,




                                                                                                         RQ
                                        ,
             J




                                     QJ
          QLQ




                                3UXQL
       3UX




                                                     /7
                                                         3




                                                                    6\QDSVH
                                                             /7'




                                                                   (OLPLQDWLRQ
                   3UXQLQJ,,




                                                                                                    LR Q
                                      6S
                                        LQ H




                                                                                                  DW
                                                                                                 UP
                                             IR




                                                                                             IR
                                              UP




                                                   DW
                                                                                             R
                                                                                        Q




                                                        LR
                                                             Q                          XW
                                                                                   %R



                                     )LO
                                         R   SRG                          RQ
                                                   L X P S UR M H F W L


Figure 2.5 – Synapse elimination. Two possible paths for pruning: presynaptically initi-
             ated or postsynaptically initiated. In the former case, bouton loss (pruning
             I) is followed by spine retraction (pruning II). In the latter case, the order of
             events is reversed. The two processes are likely to be independent.




                                                                                                              27
3 Modeling and parameterization of
  excitatory synapses

                                                      “The devil is in the details.”
                                                                                Anonymous

The most tangible contribution of this work to the ﬁeld of Neuroscience was the de-
velopment of a novel synapse model, suitable for studying the eﬀects and implications
of synaptic plasticity in the context of biophysically detailed simulations of neocortical
tissue [Markram et al., 2015]. In the future, such a model could be used to investigate
higher order brain functions, such as learning, in a multiscale manner from networks
down to single neurons and synapses.
   Synapse modeling is an active research topic bridging diﬀerent ﬁelds, from biology to
mathematics and computer science. In recent years many aspects of synaptic physiology
have been successfully captured in a variety of theoretical frameworks. However, a single
coherent description of the most relevant synaptic dynamics for learning and memory
is yet not available. To ﬁll this gap, we designed a novel synapse model to meet the
following requirements:

    i Accounting for all the main physiological features of excitatory synapses;

   ii Constrained to experimental data whenever available;

  iii Computationally parsimonious;

  iv Coupled with a constrained network and a ﬁtting infrastructure;

   v Built reusing existent and consolidated results on synapse modeling whenever avail-
     able;

  vi Easily extensible as new data on neuronal and synaptic physiology and plasticity
     emerges.

The result is a simple but comprehensive implementation of the basic synapse life cycle,
as presented in Chapter 2, for the NEURON simulator [Carnevale and Hines, 2006].
  In this chapter we introduce the synapse model, which we will refer to as the GluSy-
napse, and describe all the steps and data sources required to parametrize it. In Section
3.1 we present all the sub-components of GluSynapse and how they are orchestrated to-
gether. Section 3.2 is dedicated to the ﬁtting strategy speciﬁcally designed to parametrize
the long-term plastic dynamics. Validation and predictions of the parametrized model



                                                                                        29
3 Modeling and parameterization of excitatory synapses

are then discussed in Section 3.3. Finally, in Section 3.4 we deal with some technical
aspects of the NMODL implementation (NEURON) and simulation.


3.1 The GluSynapse model
GluSynapse is a collections of sub-models, as illustrated in Figure 3.1. Brieﬂy, the model
has two structural states: functional and potential. The former represents a regular
excitatory synapse, while the latter describes a close apposition between an axon and a
dendrite, suitable for hosting a synapse but not yet fully realized. Transitions from one
structural state to the other are regulated by the rewiring mechanism (see Sec. 3.1.7).
   In the functional state, presynaptic vesicle release (see Sec. 3.1.2) induces a postsy-
naptic current by opening the alpha-amino-3-hydroxy-5-methyl-4-isoxazole propionate
receptors (AMPARs) and N-methyl-D-aspartate receptors (NMDARs) (see Sec. 3.1.1).
Furthermore, Voltage Dependent Calcium Channelss (VDCCs) can also allow current
ﬂuxes if the cell membrane is suﬃciently depolarized (see Sec. 3.1.3). We assumed that
all excitatory synapses are located on spines (see Sec. 3.1.4), although our results can be
easily generalized to shaft synapses. Postsynaptic calcium dynamics are explicitly taken
into account (see Sec. 3.1.5), as they drive long-term changes of synaptic eﬃcacy (see
Sec. 3.1.6). Calcium ions enters the spine mostly through NMDAR and VDCC con-
ductances at these synapses [Koester and Sakmann, 1998, Sabatini and Svoboda, 2000,
Sabatini et al., 2002]. Release from internal stores and diﬀusion from dendrites was not
taken into account in the model as their contribution to the total calcium current is still
debated [Bloodgood and Sabatini, 2007]. Finally, synapse elimination is described as a
stochastic process, dependent on synaptic eﬃcacy (see Sec. 3.1.7).
   In the potential state, the synapse model cannot induce any current into the parent
dendrite. The only enabled component is the rewiring model to allow dynamic synapse
formation (see Sec. 3.1.7).
   The basic components of the functional state, excluding synaptic plasticity, were de-
signed, parameterized and implemented by several scientists at the Blue Brain Project
(BBP) over the last 12 years. The following sections provide a global overview of the
whole model, emphasizing the original contributions of this work.




30
                                                               3.1 The GluSynapse model



                                                                /73/7'



                                          9'&&                 &DOFLXP



                                         10'$5                   6SLQH



          5HOHDVH                        $03$5



                                                                     5HZLULQJ



         'HQGULWH                              J


  Figure 3.1 – GluSynapse logical structure. AMPAR and NMDAR with double exponential
               conductance proﬁle and magnesium block dynamics (NMDAR only). Vesi-
               cle release mechanisms based on the Tsodyks-Markram (TM) model. VDCC
               of High Voltage Activated (HVA) family, based on R-type channel. Spine
               as an independent biochemical compartment. Calcium dynamics account-
               ing for spine volume and simpliﬁed buﬀering/extrusion processes, bounded
               calcium dynamics as leaky-integration of free calcium. Long-term plasticity
               featuring both NMDAR-dependent Long Term Potentiation (LTP) and Long
               Term Depression (LTD). Rewiring mechanism supporting synaptogenesis and
               eﬃcacy-dependent pruning.


3.1.1 AMPAR and NMDAR
AMPARs and NMDARs are among the most fundamental components of an excitatory
synapse, as their activation depolarizes the postsynaptic neuron, and indeed they were
already part of the BBP’s model of neocortical tissue before the beginning of this work
[Ramaswamy et al., 2012, Markram et al., 2015, Ramaswamy et al., 2015]. We did not
alter the available models nor their parameterization in any way, but for the sake of




                                                                                         31
3 Modeling and parameterization of excitatory synapses

completeness we present them here.
  The AMPAR and NMDAR conductances are commonly described using a double
exponential proﬁle:
                                   g = gmax (b − a)                 (3.1)


                                    1
                             ȧ = −    a + jump · δ(t − trel )                     (3.2a)
                                    τr
                                    1
                              ḃ = − b + jump · δ(t − trel )                      (3.2b)
                                    τd


                              jump = −e−tp/τr + e−tp/τd                            (3.3a)
                                              τr τd
                                 tp =                                             (3.3b)
                                      (τd − τr )log(τd /τr )

where gmax (nS) is the peak conductance of the receptor population; a is the rising
component of the conductance, with time constant τr (ms); b is the decaying component
of the conductance, with time constant τd (ms); trel is the time of a successful release
event (ms); tp (ms) is the time to peak of the conductance; jump is a normalization
factor such that when t = trel + tp, g = gmax .
   The current generated by the synaptic population of AMPAR is then simply computed
as
                                 iAMPAR = g · (Vm − E)                             (3.4)
while the one produced by the NMDARs depends also on the magnesium block dynamics

                                             1
                                   mg =                                            (3.5a)
                                       1 + e−c1 Vm c2
                              iNMDAR = mg · g · (Vm − E)                          (3.5b)

where E (mV) is the reversal potential of the receptor; Vm (mV) is the membrane
potential; c1 and c2 are appropriate constants (see Jahr and Stevens [1990]).

3.1.2 Vesicle release
Vesicle release is a complex electrochemical mechanism allowing information to ﬂow
from the presynaptic neuron to the postsynaptic one. Moreover, it is responsible for the
distinctive features of synaptic transmission, such as short-term facilitation and short-
term depression.
   These dynamics are often described using the so called TM model [Tsodyks and
Markram, 1997]. In the original formulation, the TM model consists of a set of three
Ordinary Diﬀerential Equations (ODEs), designed and parameterized to describe the
short-term dynamics of depressing synaptic connections, such as those between Thick-
tufted pyramidal cells (TTPC) in layer 5. Over the course of the years, the model



32
                                                                           3.1 The GluSynapse model

has been reﬁned and extended, for example to capture other short-term eﬀects such
as facilitation or to model independent release sites rather than whole connections (see
[Fuhrmann et al., 2002, Loebel et al., 2009]). The latter formulation is particularly in-
teresting because it reconciles the TM model with the quantal model of vesicle release
[Purves et al., 2014].
   As part of the BBP’s eﬀort to accurately reproduce synaptic transmission, Muller
et al. [2011] further elaborated on the previous works to design a model of single release
site dynamics featuring:

     1. No resource consumption in case of transmission failure,

     2. Explicit vesicle recovery dynamics,

while retaining convergence to canonical TM model for number of trials → ∞. In this
model, a release site R is assumed to be in one of two states: ready or recovering.
Whenever R is in the ready state, the arrival of a presynaptic spike at time t triggers
vesicle release with probability U (t) given by:

               U (t) = U (tsyn ) · e−(t−tsyn )/τfac + USE · (1 − U (tsyn ) · e−(t−tsyn )/τfac )      (3.6)

where USE is the absolute release probability, that is, the release probability of a synapse
not activated in a very long time; tsyn (ms) is the time of the last presynaptic spike; where
τfac is the facilitation time constant. In case of successful release event, the postsynaptic
AMPARs and NMDARs are activated as shown in Eq. 3.2 for trel = t) and the R state is
set to recovering. Vesicle recovery is modeled as a survival process, where the probability
of remaining in state recovering at time t Psurv is given by:

                                         Psurv (t) = e−(t−tsyn )/τdep                                (3.7)

where τdep is the depression time constant.
  We extended the model in Muller et al. [2011] by adding support for Multi Vesicular
Release (MVR) [Rahmon et al., in preparation], a key feature of synaptic transmission
with major implications on plasticity1 . In our model, a pool of release sites is assumed
to be available at every synapse. The dynamics of each release site are described as in
Muller et al. [2011]. The activation functions of AMPAR and NMDAR, presented in
Eq. 3.2), were adapted to ensure that the peak conductance is obtained only when all
available vesicles are released:
                                      1              Nrel
                                ȧ = −   a + jump ·       · δ(t − trel )                           (3.8a)
                                      τr            NRRP
                                      1              Nrel
                                ḃ = − b + jump ·         · δ(t − trel )                           (3.8b)
                                      τd            NRRP
 1
     To give an example, LTP requires calcium inﬂux through the NMDAR [Markram et al., 1997b] (see
     also Sec. 3.1.5 and Sec. 3.1.6). As vesicle release is required to activate theNMDAR, the number of
     vesicles ready to be released has a strong impact on the capability of individual synaptic contacts to
     potentiate.




                                                                                                        33
3 Modeling and parameterization of excitatory synapses

where Nrel is the number of released vesicles at time trel and NRRP is the size of the
readily releasable pool.
   As of the other sub-components of GluSynapse designed before the beginning of this
work, parameters for the release model are available online through the Neocortical
Microcircuit Collaboration (NMC) portal [Ramaswamy et al., 2015]. The sole excep-
tion is the number of readily releasable vesicles per synapse NRRP , ﬁxed to 2 from our
preliminary results on MVR at excitatory connections [Rahmon et al., in preparation].
Although this value is not yet ﬁnalized at the time of writing, recent experimental ﬁnd-
ings on pyramidal to interneuron connections [Molnár et al., 2016] seem to support our
choice.

3.1.3 VDCC
VDCCs are the most important source of calcium at dendritic spines after the NMDARs
[Sabatini and Svoboda, 2000]. The model of VDCCs used in this work was designed
and parametrized by Christian Rössert [Chindemi et al., in preparation], using data and
components from several sources [Magee and Johnston, 1995, Sabatini and Svoboda,
2000]. In brief, a simple but inactivating R-type VDCC was designed in the Hodgkin-
Huxley (HH) formalism. The population conductance is given by

                                      g = gmax m2 h                                  (3.9)


                                           m∞ − m
                                     ṁ =                                         (3.10a)
                                             τm
                                           h∞ − h
                                      ḣ =                                        (3.10b)
                                             τh


                                                1
                                  m∞ =          Vhm −Vm                           (3.11a)
                                          1+e     km

                                                1
                                   h∞ =         Vhh −Vm                           (3.11b)
                                          1+e      kh



where gmax (nS) is the peak conductance of the VDCCs population; m is the activation
variable, with time constant τm (ms) and steady state m∞ ; h is the inactivation variable,
with time constant τh (ms) and steady state h∞ ; Vhm (mV) is the half-maximum activa-
tion voltage and km (mV) is the slope factor; Vhh (mV) is the half-maximum inactivation
voltage and kh (mV) is the slope factor; Vm (mV) is the membrane potential, as usual.
   The current induced by the population of VDCCs is then computed as

                                    i = g(Vm − Eca )                               (3.12)

where Eca is the reversal potential for calcium.



34
                                                               3.1 The GluSynapse model

3.1.4 Spine
The vast majority of excitatory synapses are located on dendritic spines. With respect
to synaptic plasticity, spines are particularly important because they act as biochemical
compartments, granting speciﬁcity of synaptic changes. Furthermore, spine volume is
a key determinant of intracellular calcium concentration, the main driver of synaptic
plasticity.
   In this work we do not explicitly model spines for computational reasons (see Sec.
3.4), but we do account for biochemical compartmentalization and volume-related ef-
fects on calcium concentration. Speciﬁcally, a spherical spine head is assumed for all
computations and calcium ions are not allowed to diﬀuse into the parent dendrite.
   Spine volume parameterization was carried out by Christian Rössert [Chindemi et al.,
in preparation], using data from Ruth Benavides-Piccione (unpublished at the time
of writing). A crucial aspect of the parameterization was the correlation introduced
between spine volume and NMDAR conductance [Arellano et al., 2007]. This feature
produced normalized calcium transients across synapses. Consequently, we could easily
generalize the long-term plasticity model parameters obtained for a prototypical synapse
to all the others, as later explained in Sec. 3.2. The algorithm used to generate correlated
pairs of spine volume and conductance is described in elsewhere [Chindemi et al., in
preparation].

3.1.5 Postsynaptic calcium
The postsynaptic calcium model was designed and implemented by Christian Rössert
using data from several sources [Sabatini et al., 2002, Cornelisse et al., 2007]. In brief,
calcium current can ﬂow through cell membrane via two paths: NMDARs and VDCCs.
After entering the spine, calcium ions quickly bind to endogenous buﬀers and only a
small fraction remains free. Slower mechanisms, such as calcium pumps and diﬀusion,
are in place to re-establish the intracellular calcium concentration. These dynamics are
all summarized in a single ODE, modeling free calcium at the postsynaptic membrane
surface:
                                                                    (0)
                    d[Ca]i                         γ   ([Ca]i − [Ca]i )
                            = (iNMDAR + iVDCC ) −                                  (3.13)
                      dt                          V           τCa
where iNMDAR is the calcium component of the NMDAR-mediated current, as later
described; iVDCC is the VDCC-mediated current; γ is the fraction of free (non buﬀered)
                                    (0)
calcium; V is the spine volume; [Ca]i is the intracellular calcium concentration at rest;
τCa is the time constant of calcium transients.
  While VDCC current can all be considered calcium, NMDAR are permeable to several
other ions, including sodium and potassium. The fractional calcium current Pf through
the NMDAR can be calculated from the Goldman-Hodgkin-Katz (GHK) ﬂux equation
[Schneggenburger et al., 1993]:

                               iCa               4[Ca]o
                     Pf =            =                        Vm F                   (3.14)
                            iCa + iM   4[Ca]o + pM [M ](1 − e2 RT )
                                                  pCa




                                                                                         35
3 Modeling and parameterization of excitatory synapses

where iCa and iM are respectively the current due to calcium ions and the one due to all
monovalent ions; [Ca]o is the extracellular calcium concentration; [M ] is the concentra-
tion of monovalent ions, assumed to be identical inside and outside the cell membrane;
pCa and pM are respectively the permeability to calcium ions and the one to monovalent
ions; Vm is the membrane voltage; T is the temperature; F is the Faraday’s constant
and R is the ideal gas constant.
   Unfortunately Eq. 3.14 cannot be directly used in our model as a scaling factor for
the total NMDAR current because it loses physical meaning in the limit of Vm → 0
mV, a threshold often crossed during action potential (AP) (see Fig. 3.2). That is,
Pf → 1 for Vm → 0, but as also iNMDAR → 0, the total calcium current would be
incorrectly estimated. Based on the observation that the calcium current approaches 0
for Vm ≥ 40 mV [Schneggenburger et al., 1993], we approximated the calcium component
of the NMDAR conductance gNMDAR       as a separate ion channel with reversal potential
Erev = 40 mV:
                                          (−∞)
                            gNMDAR = αPf          gNMDAR                         (3.15a)
                             iNMDAR = gNMDAR
                                        
                                              (Vm − Erev )                       (3.15b)

where α is an appropriate scaling factor, determined by comparing model simulation
                                         (−∞)
results to Eq. 3.14 (see Fig. 3.2); Pf         is the fractional calcium current at very
negative potentials (Eq. 3.14, all other parameters as in Schneggenburger et al. [1993]);
Erev was assumed to be independent of extracellular calcium concentration, at least for
values around 1 mM (see [Jahr and Stevens, 1993]). The proposed model was then
tested against the ﬁndings of Schneggenburger et al. [1993] and accepted as a good
approximation in the physiological voltage range and at relevant extracellular calcium
contention (see Fig. 3.2).
   Intracellular calcium dynamics are key for the induction of long-term synaptic changes
and so they have been heavily studied in the context of Spike Timing Dependent Plastic-
ity (STDP) [Nevian and Sakmann, 2006]. Even tough it is nowadays clear that calcium
is the main driver for such changes, how exactly synapses could read this signal is still
debated. One of the most (ab)used theories is the so called calcium control hypothesis:
LTP is induced by large calcium transients, while LTD is tirggered by prolonged medium-
level transients. Although this simple explanation is very elegant and appealing from
the theoretical point of view, its experimental validation has proven to be particularity
tricky and many open questions remain. In particular, it was shown that high-frequency
stimulation triggers LTP while low-frequency does not [Markram et al., 1997b], but
whether the former causes calcium accumulation at the synapse is unknown [Manni-
nen et al., 2010]. Some experimental evidence seems to suggest that small-capacity and
fast-binding calcium reaction partners need to be saturated before the mechanisms re-
sponsible of LTP and LTD could be activated (for a review see [Manninen et al., 2010]).
Indeed preliminary analysis showed that calcium does not accumulate in our model at
high-frequency stimulation (data not shown). For this reason, we speculate that the
frequency dependence of LTP could not be explained based on free-calcium only and
requires to consider the diﬀerent calcium reaction partners. To account for such dy-



36
                                                                            3.1 The GluSynapse model

       $                                                    %
                                                                      
                                

                      9P P9
                                                                     

                               ಝ
                                                                      



                                                                      
                               

                               
                      L Q$




                                                                      



                                                                3I
                               

                                                                     ಝ


                       
                                                                     ಝ
                      ಝ
           LFD Q$




                      ಝ                                        ಝ

                      ಝ

                                                                     ಝ
                                                                    ಝ       
                                          WLPH PV                           9P P9


  Figure 3.2 – Fractional NMDAR-mediated calcium current. A top, a series of voltage
               steps has been applied to a simulated patch of membrane around our model
               of NMDAR; the total NMDAR current (middle) and its calcium component
               (bottom) have been measured during voltage steps; persistent activation of
               the NMDAR was simulated by setting the decay time constant τd in Eq. 3.2
               to an extremely large value. B, comparison of simulated fractional calcium
               current Pf (triangles) and predicted values from the GHK equation (solid
               line) [Schneggenburger et al., 1993]; simulated Pf was computed by dividing
               the value of calcium current by total current at 750 ms as reported in A.


namics without explicitly modeling the endogenous buﬀers, we introduced the calcium
integrator [Ca]i :
                               d[Ca]i    [Ca]i               (0)
                                       =−        + [Ca]i − [Ca]i                      (3.16)
                                 dt         τ
where τ is an appropriate time constant for calcium integration, here ﬁxed to 200 ms.
This value was chosen considering both the average duration of an NMDAR-mediated
spike [Schiller et al., 2000], considered the key signal for triggering LTP, and the typical
stimulation window of LTP-inducing protocols [Markram et al., 1997b]. The idea behind



                                                                                                 37
3 Modeling and parameterization of excitatory synapses

the introduction of [Ca]i is to prevent a single and isolated calcium transient from
inducing any long-term change of synaptic eﬃcacy, while allowing several pulses at high
frequency, such as those used in STDP protocols, or large cooperative events, such as
glsnmdar-mediated spikes, to produce enough activation to trigger a long-term synaptic
change.

3.1.6 Long-term plasticity
As discussed in chapter 2, calcium-dependent forms of synaptic plasticity are the most
studied and well characterized, both experimentally and theoretically. For this reason,
we decided to adapt a calcium-based model of STDP by Graupner and Brunel [2012]
to account for more realistic synaptic dynamics, such as those described above. Before
discussing the details of the extended model, we brieﬂy introduce the main features of
the original one [Graupner and Brunel, 2012], from now on simply the Graupner-Brunel
(GB) model.
  The GB model assumes bistable synapses [Petersen et al., 1998, Bagal et al., 2005,
O’Connor et al., 2005]. That is, the “eﬃcacy” of each synapse is a dynamical variable
with only two stable states: potentiated and depressed. Synaptic plasticity is then a
transition from one state to the other 2 . State transitions are governed by postsynaptic
calcium transients, through a thresholding mechanism. Whenever calcium crosses the
depression (potentiation) threshold, the synaptic eﬃcacy steps toward the depressed
(potentiated) state. The magnitude of the steps depends on the time spent above the
threshold.
  The synaptic eﬃcacy is an abstract concept of the model, to be mapped into an
actual synaptic quantity (i.e. density of AMPARs, release probability or number of
docked vesicles). In the original model, the eﬃcacy is linearly translated to a generic
“weight”, a correlate of the peak synaptic conductance.
  In equations, the Graupner and Brunel [2012] model is formulated as:
                        dρ
                    τ      = −ρ(1 − ρ)(ρ − ρ) + γp (1 − ρ)Θ[c(t) − θp ]
                        dt                                                                       (3.17a)
                              − γd ρΘ[c(t) − θd ] + Noise(t)
                        dc     c                                  
                           =−     + Cpre    δ(t − ti − D) + Cpost      δ(t − tj )                (3.17b)
                        dt    τca
                                             i                            j

                        w = w0 + ρ(w1 − w0 )                                                     (3.17c)

where ρ represents the synaptic “eﬃcacy”, with time constant τ ; ρ delimits the basin of
attraction of the two stable states; Θ is the Heaviside function; θd and θp are respectively
the depression and potentiation thresholds; γd and γp the depression and potentiation
rate; Noise(t) and activity dependent noise; c(t) is the postsynaptic calcium concentra-
tion, with time constant tauca ; Cpre is the calcium transient caused by a presynaptic
 2
     Bistable synapses should not be confused with binary synapses (i.e. having only two possible values).
     In a bistable model, the synaptic eﬃcacy is a continuous variable having two stable ﬁxed points.
     That is, in the absence of synaptic activity the eﬃcacy moves toward the closest ﬁxed point.




38
                                                                      3.1 The GluSynapse model

spike at time ti , the delay D account for the slow activation of the NMDAR; Cpost is the
calcium transient caused by a postsynaptic spike at time tj ; w is the synaptic weight;
w0 and w1 are the weight of the depressed and potentiated state respectively. Another
important variable of the GB model is the ratio of depressed synapses β, controlling the
maximum potentiation (depression) obtainable:
                                                w0
                               ΔEPSPd =                                                 (3.18a)
                                        β0 w0 + (1 − β0 )w1
                                                w1
                               ΔEPSPp =                                                 (3.18b)
                                        β0 w0 + (1 − β0 )w1

where ΔEPSPd and ΔEPSPp are respectively the maximum amount of depression and
the maximum amount of potentiation; β0 is the initial ratio of depressed synapses; w0
and w1 have the usual meaning. For a more precise description of the GB model and its
parameters, we invite the reader to refer to the original article by Graupner and Brunel
[2012].
   Our ﬁrst modiﬁcation of the GB model was dropping the noise term in the synap-
tic eﬃcacy (Eq. 3.17a), as we already account for several stochastic aspect of synaptic
transmission and membrane potential ﬂuctuations. Second, we replaced the linear cal-
cium component (Eq. 3.17b) with the more detailed description described in Sec 3.1.5,
as NMDARs and VDCCs are explicitly modeled and so their contribution to postsynap-
tic calcium concentration. For more details on the postsynaptic calcium model please
refer to Sec. 3.1.5. Finally, the abstract synaptic eﬃcacy (Eq. 3.17a) is mapped into
presynaptic release probability rather than a generic weight. We decided to restrict our-
selves to presynaptic changes as they seem to be the predominant plastic variable in the
somatosensory cortex [Tsodyks and Markram, 1997]. However, postsynaptic changes
were indeed observed in the somatosensory cortex and other neocortical areas [Sjöström
et al., 2007, Perin and Chindemi, in preparation] and they could be easily added to
the model in the future. To account for the slow implementation of synaptic changes
often observed in vitro, where the synaptic “strength” right after the induction protocol
remains basically unaltered (see Fig. 3C in Markram et al. [1997b]), we introduced a
delay factor in the eﬃcacy to release probability mapping:
                                               (target)
                                 dUSE  U        − USE
                                      = SE                                              (3.19a)
                                  dt       τchange
                                (target)      (d)         (p)   (d)
                              USE          = USE + ρ[USE − USE ]                        (3.19b)

where USE is the absolute release probability, as previously discussed (see Sec. 3.1.2);
                                            (target)
τchange is the time constant of the delay; USE       is the target absolute release probability;
  (d)                                                                 (p)
USE is the absolute release probability in the depressed state; USE is the absolute release
probability in the potentiated state.
   As we mentioned before, we included a new term, the calcium integrator [Ca]i , to
explain the frequency dependence of LTP. Now that we introduced the fundamentals of
the GB model, we can clarify why this term was needed. It is known that pairing pre- and



                                                                                             39
3 Modeling and parameterization of excitatory synapses

post-synaptic activity at high frequency induces LTP, while low frequency is ineﬀective
at triggering any change [Markram et al., 1997b]. In the GB model, potentiation and
depression depend on the time spent by the postsynaptic calcium concentration about
the respective thresholds θp and θd . However, as

     1. layer 5 TTPC connections have strongly depressing short-term dynamics [Markram
        et al., 1997a],

     2. we do not expect them to have a large Readily-Releasable Pool (RRP) of vesicles
        (see [Molnár et al., 2016]),

     3. the [Ca]i time constant τCa is fast (see Sec. 3.1.5),

high-frequency stimulation could not produce signiﬁcant calcium accumulation at these
synapses and low-frequency stimulation is much more eﬀective at triggering calcium
events (Fig. 3.3 and 3.4). For this reason, comparing the time spent by the [Ca]i
above any threshold leads us to the conclusion that a diﬀerent “read-out” would be
needed to explain why high-frequency stimulation is more eﬀective than low-frequency
to induce LTP. The results of previous studies indicate that saturation of fast calcium
buﬀers could be a requirement for LTP’s frequency-dependence [Holmes and Levy, 1990,
Markram et al., 1998]. The calcium integrator was the chosen read-out as it rewards a
series of calcium events in a short time, somehow mimicking the interplay of a fast and
a slow calcium buﬀer. Indeed using [Ca]i rather than [Ca]i sorts times above threshold
in the expected order (Fig. 3.4) and allowed us to ﬁt the synapse model to LTP and
LTD data (see Sec. 3.2).
  Fitting our extended model of long-term plastic changes is a non trivial task. The main
diﬃculty is the inaccessibility of the relevant synaptic variables, not directly observable
during standard plasticity experiments, neither in vivo nor in vitro. For this reason, we
designed of a special procedure to ﬁt single synapses on the basis of whole connection
data from STDP experiments in vitro, described in Sec. 3.2.




40
                                                                        3.1 The GluSynapse model




     $                     %                              K]PV67'3SURWRFRO
                                                 
                               VRPD9P P9


                                                  

                                             ಝ


                                                          K]PV67'3SURWRFRO
                                                 
                               VRPD9P P9




                                                  

                                             ಝ




                                                 
                                     >&D@L ͑0




                                                 
                                                                                     ͍
                                                  
                                                                     
                                                                 WLPH PV


Figure 3.3 – Calcium thresholding problem. A, 2D representation of a typical layer 5
             TTPC used in our experiments, green dot marks the position of a synapse used
             for the analysis in the next panel. B, spiking activity of cell in panel A when
             stimulated at low frequency (top) and high frequency (middle), black dots
             mark the spikes of the presynaptic neuron innervating the cell at the analyzed
             synapse; postsynaptic calcium transients (bottom) for the two stimulation
             conditions; black dashed line marks a hypothetical threshold θ.




                                                                                             41
3 Modeling and parameterization of excitatory synapses




      $                                                 %
                                                                            

                                                                              
                       




                                                            WLPHDERYH͍ V
                                                                              
          >&D@L ͑0




                                                                              
                       
                                                                              

                                             ͍                                
                        
                                                                              

                                                                             
                                                                                                  
                                 WLPH PV                                                WKUHVKROG͍ ͑0


      &                                                 '
                                                                           

                                                                              
                       
                                                            WLPHDERYH͍ V




                                                                              
        >&D@L൑  ͑0V




                                                                           
                                             ͍
                                                                              
                       
                                                                              
                       
                                                                              

                                                                           
                                                                                     
                                 WLPH PV                                                WKUHVKROG͍ ͑0V




                                                                                   
                                        VWLPXODWLRQIUHTXHQF\ +]


 Figure 3.4 – Calcium thresholding solution. A, exemplar calcium transient for a layer 5
              TTPC synapse stimulated at various frequencies during a single repetition of
              a typical STDP inducing protocol; black dashed line marks an hypothetical
              threshold θ. B, total time spent by postsynaptic calcium above threshold
              θ during the full STDP inducing protocol. C, same as A, but showing the
              time course of the calcium integrator instead. D, same as B, but time above
              thresholds is computed on the calcium integrator instead.




42
                                                              3.1 The GluSynapse model

3.1.7 Rewiring
The rewiring model was adapted from Fauth et al. [2015]. In this model the probability
of creating a new synapse pbuild at a potential location is assumed to be constant, while
the probability of eliminating an existing one pdel is a function of the current synap-
tic strength, reﬂecting the experimental ﬁndings of [Holtmaat et al., 2005] and Le Bé
and Markram [2006], who showed that weak spines/connections are more likely to be
eliminated:
                                                    2 4/3
                                  pdel = pkbuild e−a wij                           (3.20)
where k is a scaling factor controlling the generation to elimination rate; a is the gain
of the synaptic weight wij .
   The main modiﬁcations with respect to Fauth et al. [2015] were the use of ρ from the
GB model instead of the generic weights w and the replacement of Eq. 3.20 with a better
behaved logistic function. The latter has the advantage of being monotonic and deﬁned
for negative values, allowing a simpler handling of some limit cases (i.e. extremely weak
synapses where the eﬃcacy variable ρ is transiently smaller than zero):
                                                  L
                                   pdel =                                           (3.21)
                                            1 + e−k(ρ−ρ )
where L is the maximum value, that is, the elimination probability of a fully depressed
synapse; k controls is the steepness of the curve; ρ and ρ have their usual meaning (see
Sec. 3.1.6).
  Parametrization of the rewiring model is based on both experimental evidence and
theoretical considerations on the nature of potential synaptic contacts:
   • Total number of synapses decreases with age already in a time window of
     24 hours (from P16 to P17), about 30% of the imaged spines disappeared while
     only 25% new spines were created [Holtmaat et al., 2005].

   • Newly formed connections are weak and have low release probability [Le Bé
     and Markram, 2006].

   • Weak connections have an higher chance of being eliminated, although
     newly formed ones seem to be protected from immediate pruning [Le Bé and
     Markram, 2006].

   • Potential synapses greatly outnumber functional ones as noticed by Stepa-
     nyants et al. [2002].
Based on these observations, we assumed a transition schema for rewiring as the one
in Fig. 3.5. A synapse can only be in one of the two structural states: functional
or potential [Reimann et al., 2015]. We further distinguish between potentiated and
depressed synapses to account for the diﬀerent elimination probability. A transition
from state i to state j is a stochastic event, happening with probability pi,j .
   Fitting the rewiring model corresponds then to estimate the probability of the nine
possible state transitions. To solve the problem, we ﬁrst considered all the constraints we



                                                                                        43
3 Modeling and parameterization of excitatory synapses




                 S                                     S                       S

                       )XQFWLRQDO      S                            )XQFWLRQDO
                       GHSUHVVHG                                     SRWHQWLDWHG
                                                                         



                          S                                   S




                                       S                     S



                                             3RWHQWLDO
                                                




                                                S


  Figure 3.5 – Rewiring model. We consider three possible synaptic states [0, 1, 2], namely
               potential, functional depressed and functional potentiated. Transition proba-
               bilities are indicated by pij labels.


could derive from experimental evidence and from the BBP’s connectivity rules [Reimann
et al., 2015]. The probability of a new spine appearing from P16 to P17 is on average 23%
[Holtmaat et al., 2005]. Furthermore, potential synapses exceed functional synapses by
a factor 16 [Markram et al., 2015, Reimann et al., 2015] (see Fig. 3.6). In mathematical
terms:
                                              0.23
                                  p01 + p02 =       0.014                          (3.22)
                                               16
 The probability of a spine disappearing from P16 to P17 is on average 0.32 [Holtmaat
et al., 2005]. In our circuit 50% of the synapses are potentiated and 50% are depressed
(see Sec. 3.2.1). Leading to:

                                    0.5p10 + 0.5p20 = 0.32                               (3.23)

Newly formed synapses are small [Holtmaat et al., 2005] and whole new connections are
weak [Le Bé and Markram, 2006]. For this reason the probability to immediately create
a very strong synapse is negligible:

                                             p02 = 0                                     (3.24)



44
                                                                  3.1 The GluSynapse model

       $                                        %




       &
                 H                                            H
                                                         

                                                         

           
                                                            
                                                             3RWHQWLDO )XQFWLRQDO
                                                                6\QDSVHV 6\QDSVHV
           

           

           
                       7RXFK    5HGXFH          &XW          $FWLYH     )XQFWLRQDO
                       5XOHV                                )UDFWLRQ

  Figure 3.6 – Potential synaptic locations in the BBP’s reconstruction. A, all potential
               synaptic locations around the soma of a single layer 5 TTPC, as determined
               by the “touch detection” algorithm [Markram et al., 2015]; yellow spheres
               represent potential synapses, the round soma is an artifact of the detection
               algorithm, using a ﬁxed radius rather than the actual morphology for com-
               putational reasons. B, functional synapses on the same cell in A; the pruning
               process is carried on by a multi-step ﬁltering algorithm enforcing several bi-
               ological constraints (i.e. bouton density and number of synapses per connec-
               tion) [Reimann et al., 2015]. C, excitatory-to-excitatory synapses pruned at
               every step of the algorithm; blue bars represent the pool of potential synapses
               used in this work, green the initial count of functional synapses; the inset
               compares the overall number of potential synapses against functional ones.


The turnover rate for small spines is much larger than the one of big spines [Holtmaat
et al., 2005]. Although spine size and synaptic “strength” are well correlated, it is very
diﬃcult to quantitatively map synapse properties to elimination rate in our datasets. For
this reason, we simply assumed an order of magnitude diﬀerence between the elimination




                                                                                             45
3 Modeling and parameterization of excitatory synapses

probability of potentiated and depressed synapses:

                                     p10 − 10p20 = 0                                 (3.25)

We complemented the four equations above with other two, based on theoretical consid-
erations on the potentiation/depression processes. In particular, we assumed that the
probability of a depressed synapse to remain depressed or to potentiate are identical.
Indeed we don’t see why a newly formed synapse, depressed by deﬁnition as explained
before, would not potentiate in the 24h observation window. We than added the equa-
tion:
                                     p11 − p12 = 0                             (3.26)
Furthermore, we want the ratio of depressed and potentiated synapses to be maintained,
until further evidence demonstrate a diﬀerent trend. We enforce this constraint by
imposing the net ﬂux of the depressed state to be identical to the one of the potentiated
state:
                   p01 − p02 − p10 + p11 − 2p12 + p20 + 2p21 − p22 = 0             (3.27)
Finally, we included three constraints to preserve the consistency of states and transi-
tions:

                                    p00 + p01 + p02 = 1                              (3.28)
                                    p10 + p11 + p12 = 1                              (3.29)
                                    p20 + p21 + p22 = 1                              (3.30)

The resulting linear system can be easily solved analytically. The ﬁnal values for the pi,j
parameters are:
                                   ⎡ ⎤ ⎡            ⎤
                                     p00       0.986
                                   ⎢p01 ⎥ ⎢0.014⎥
                                   ⎢ ⎥ ⎢            ⎥
                                   ⎢p02 ⎥ ⎢ 0 ⎥
                                   ⎢ ⎥ ⎢            ⎥
                                   ⎢p10 ⎥ ⎢0.582⎥
                                   ⎢ ⎥ ⎢            ⎥
                                   ⎢p11 ⎥  ⎢0.209⎥                                 (3.31)
                                   ⎢ ⎥ ⎢            ⎥
                                   ⎢p12 ⎥ ⎢0.209⎥
                                   ⎢ ⎥ ⎢            ⎥
                                   ⎢p20 ⎥ ⎢0.058⎥
                                   ⎢ ⎥ ⎢            ⎥
                                   ⎣p21 ⎦ ⎣0.553⎦
                                     p22       0.388
  Fitting the logistic function in Eq. 3.21 to p10 = 0.582 and p20 = 0.058 leads us to the
ﬁnal parameters for the elimination function: L = 0.640 and k = −4.605, completing
the rewiring model (Fig. 3.7A). We then re-ﬁtted the curve proposed by [Fauth et al.,
2015] (see Eq.3.20) and compared it with our formulation (Fig. 3.7B). We noticed that
our model shows a more marked saturation for small values of synaptic eﬃcacy than the
one from [Fauth et al., 2015]. Furthermore, as mentioned before, the logistic formulation
here used (see Eq.3.21) has the advantage of being deﬁned and well behaved for eﬃcacy
values outside the (0, 1) range. In addition, we considered an alternative formulation of
Eq. 3.20, used by Fauth et al. [2015] to test the sensitivity of their results to the chosen



46
                                                                  3.2 Fitting the long-term plasticity model

parameters. They found that indeed their conclusions hold also for the alternative
formulation [Fauth et al., 2015]. As this second version of the elimination curve is
very similar to our sigmoid, we are conﬁdent that it was a valid replacement. Synapse

       $                                                          %
                                                                                             
                                               HOLPLQDWLRQ                                                         %%3




                                                                      HOLPLQDWLRQSUREDELOLW\
                                                                                             
                                               IRUPDWLRQ                                                           )DXWK
                                                                                                             )DXWK
           SUREDELOLW\




                                                                                             

                                                                                             

                                                                                             

                                                                                             

                                                                                             
                                                                                             
                                      V\QDSWLFHIILFDF\                                                V\QDSWLFHIILFDF\


  Figure 3.7 – Rewiring ﬁt. A, synapse formation and elimination probabilities; synapses
               emerge at a constant rate, while their elimination probability depends on
               synaptic eﬃcacy. B, comparison of the elimination curve used in this work
               (BBP, see Eq. 3.21) with the one used in Fauth et al. [2015] (Fauth 4/3, see
               Eq. 3.20); an alternative formulation where the wij exponent in Eq. 3.20 is
               changed from 4/3 to 2 is also shown (Fauth 2) [Fauth et al., 2015]; the x axis
               is purposely include points smaller than 0 and greater than 1 to show the
               behavior of the three functions for limit and pathological cases.

formation and elimination were then implemented as an Poisson processes of parameter
pbuild = p01 and pdel respectively (see 3.4.3 for technical details on the implementation
of the elimination process).


3.2 Fitting the long-term plasticity model
Observing the activity of single synapses undergoing functional plastic changes is nowa-
days generally not possible, with notable exceptions [Nevian and Sakmann, 2006]. For
this reason, the vast majority of available datasets contain only information at the level
of whole (multi-synaptic) connections [Markram et al., 1997b, Sjöström et al., 2001].
This lack complicates quite dramatically the task of constraining models of synaptic
plasticity to experimental data. To overcome the problem, we designed a ﬁtting proce-
dure capable of extracting a suitable parameterization for single synapses from whole
connection data.
  In brief, a pool of in silico connections between layer 5 TTPC is randomly chosen from
our circuit reconstruction [Markram et al., 2015]. The connections are then stimulated
by pairing the activity of the pre- and post-synaptic cells to trigger functional synaptic
plasticity, as in standard in vitro STDP experiments [Markram et al., 1997b, Sjöström



                                                                                                                               47
3 Modeling and parameterization of excitatory synapses

et al., 2001]. A genetic algorithm identiﬁes the best set of synaptic parameters by
comparing the in silico results with the in vitro dataset in Markram et al. [1997b]. The
obtained parameter set is then generalized to all other synaptic pathways.
  In the following sections we describe the details of the ﬁtting procedure. We start by
stating the optimization problem in Sec. 3.2.1 and how neuron pairs were selected in
Sec. 3.2.2. Then we characterize the actual ﬁtting method in Sec. 3.2.3. Finally, the
generalization steps are discussed in Sec. 3.2.4 and 3.2.5.

3.2.1 Problem statement and initial conditions
To ﬁt the long-term plasticity model to experimental data we need to ﬁnd appropriate
values for several parameters:

  τ               time constant of the synaptic eﬃcacy;

  θd , θp         threshold for depression and potentiation respectively;

  γd , γp         rate for depression and potentiation respectively;

  τchange         time constant of long-term release probability changes;
      (d)   (p)
  USE , USE       absolute release probability in the depressed and potentiated state re-
                  spectively;

  β0              initial ratio of depressed synapses.

We decided to ﬁx some parameters and to limit the search space of others to our best
guess to simplify the ﬁtting procedure (see table 3.1).

                              Parameter       Value        Unit
                                  τ             70          s
                                  θd        (0.25, 0.5)    µM s
                                  θp        (0.25, 0.5)    µM s
                                  γd         (1, 350)     scalar
                                  γp         (1, 350)     scalar
                               τchange          100         s
                                   (d)
                                 USE            0.2       scalar
                                   (p)
                                 USE            0.8       scalar
                                  β0            0.5       scalar

  Table 3.1 – Parameters of the long term plasticity model. Values between parentheses rep-
              resent search boundaries (lower, upper), while the others are ﬁxed parameters.

   The two time constants τ and τchange were estimated from preliminary ﬁtting attempts
and are compatible with previously used parameter sets (see Supplementary Information
in Graupner and Brunel [2012]). The initial ratio of depressed synapses β0 was inher-
ited from Graupner and Brunel [2012], who successfully used this value to ﬁt another



48
                                                           3.2 Fitting the long-term plasticity model

neocortical STDP dataset by Sjöström et al. [2001]. The release probability of the de-
                (d)                                      (p)
pressed state USE and the one of the potentiated state USE were chosen to satisfy both
constraints on the maximum potentiation and depression attainable and to preserve the
mean release probability of the pathway. Speciﬁcally:

     Mean release probability constraint as in Markram et al. [2015]
                                  (initial)         (d)                  (p)
                                USE           = β0 USE + (1 − β0 )USE = 0.5                       (3.32)

     Maximum depression as in Markram et al. [1997b]
                                                             (d)
                                                           USE
                                        ΔEPSPd =           (initial)
                                                                       ≤ 0.8                      (3.33)
                                                          USE

     Maximum potentiation as in Markram et al. [1997b]
                                                             (p)
                                                           USE
                                        ΔEPSPp =           (initial)
                                                                       ≥ 1.5                      (3.34)
                                                          USE

All previous formulas were derived from results in Graupner and Brunel [2012] and led
                                            (d)     (p)
us to the following set of constraints for USE and USE :
                                    ⎧     (d)         (p)
                                    ⎪
                                    ⎪ β U + (1 − β0 )USE = 0.5
                                    ⎨ 0 SE
                                            (d)
                                      0 ≤ USE ≤ 0.4                                               (3.35)
                                    ⎪
                                    ⎪
                                    ⎩           (p)
                                      0.75 ≤ USE ≤ 1

The solution of this system of equalities and inequalities is a linear set, graphically shown
                                                                   (d)        (p)
in Fig. 3.8. We validated the proposed parametrization for USE and USE by reproduc-
                                                                                        (initial)
ing in silico a typical in vitro connectivity test. Our goal was to verify that the USE
extracted at the connection level was still compatible with what reported experimentally
[Markram et al., 1997a]. We considered the standard “H20S8” stimulation protocol, as
used in our lab to probe connection properties: given a connected pair of neurons, the
presynaptic one is forced to spike 8 times at a frequency of 20 Hz; a delayed supplemen-
tary spike (Δt = 550ms) is then delivered to show recovery from depression dynamics;
the protocol is repeated 60 times and the average somatic Excitatory Post-Synaptic Po-
tential (EPSP) computed. We selected a pool of 200 in silico connections satisfying the
typical geometrical constraints holding in vitro (see Sec. 3.2.2). The synapse model
used for this test does not include any forms of long-term plasticity nor VDCC3 . After
                (initial)
simulation, USE           was calculated using a Genetic Algorithm (GA) [Van Geit et al.,
 3
     Long-term plasticity and the VDCC considered in this work are not likely to aﬀect the results of our
     test. The former operates at much longer time scales. The latter would have very little inﬂuence on
     the somatic EPSP, as only HVA channels are considered here.




                                                                                                      49
3 Modeling and parameterization of excitatory synapses


                       




                       




                       
                 86(
                S




                       




                                       G               S
                                ͆86(   ಜ ͆ 86(        
                                         G
                                     86( ೮ 
                                            S
                                    ೮ 86(  
                       
                                                             
                                                         G
                                                       86(


                                                (d)          (p)
  Figure 3.8 – Set of valid solutions for USE and USE . Red, green and blue sets represent
               data points satisfying the prescribed constraints in Eq. 3.35. The intersection
               of the three is a line segment, corresponding to the set of valid solutions.
                                                          (d)        (p)
               White circle shows the chosen value for USE and USE .


2016] to ﬁt the TM model to the somatic traces (see. Fig. 3.9). Our results show that,
despite bimodal synapses, the release probability extracted at the connection level is still
unimodal. It should be noticed that the mean value we estimated is slightly higher than
                                                        (initial)
expected. There is indeed a tendency to overestimate USE          for connections composed
by a mixture of potentiated and depresses synapses. This issue will be addressed in the
next iteration of the synapse model (see Sec. 5.2).
   After ﬁxing boundaries and values of several parameters of the long-term plasticity
model, we focused on the other synaptic parameters. Our ﬁrst goal was to obtain a mean
model of calcium-dependent plastic changes for layer 5 TTPC connections, neglecting
for a moment all the complications due to variability of other synaptic parameters (i.e.
peak conductance or release probability). For this reason we ﬁxed them to their mean re-
ported values (see table 3.2 or Ramaswamy et al. [2015]). The variability of transmission
parameters was later re-introduced and accounted for in the corresponding generaliza-
tion step (see Sec. 3.2.4). The ﬁtting problem is then ﬁnally reduced to the estimation
of just four parameters: the two thresholds θd and θp , plus the corresponding rates γd
and γp .




50
                                                        3.2 Fitting the long-term plasticity model

     $




     %                                                  &
                                                                      
         
                                                                      
         
                                                                      
                                                            8VHILW



         

                                                                   
         
                                                                      
         

                                                                   

                                                                     

                             8VHILW                                                 8VHLQVLOLFR


                              (initial)
Figure 3.9 – Bimodal USE       validation. A, sample connection and corresponding ﬁt; blue
             line is the mean response of 60 repetitions of the stimulation protocol, ﬁlled
             blue area represents the 95% conﬁdence interval; green line is the ﬁt produced
                                    (initial)
             by the algorithm. B, USE         distribution as extracted by the ﬁtting algorithm
             for every connection in the in silico dataset; blue line represents the kernel
                                                         (initial)                   (initial)
             density estimate. C, prescribed mean USE              versus extracted USE        ; the
                                                           (initial)
             ﬁtting algorithm tends to overestimate USE              .

                         Parameter                                      NMODL             Value     Unit
                     Peak conductance                                    gmax               1.5      nS
            Conductance decay time constant                              tau d             1.74      ms
                 Depression time constant                                 Dep              670       ms
                 Facilitation time constant                               Fac               17       ms
           Size of readily releasable vesicle pool                       Nrrp                2     scalar

                              Table 3.2 – Fixed synaptic parameters.




                                                                                                                51
3 Modeling and parameterization of excitatory synapses

3.2.2 Pair selection
We selected in silico pairs of connected layer 5 TTPC to match the key features of
standard STDP experiments in vitro [Markram et al., 1997b, Sjöström et al., 2001].
In particular, we considered two type of constraints, namely electrical and geometrical,
next described.
  Layer 5 TTPC were ﬁrst tested for stimulation protocol compliance. That is, each
cell was stimulated as in the STDP induction protocols and discarded in case of missing
or redundant spikes (see Fig. 3.10)4 .     Cells properly responding to all stimulation
protocols were then added to the pool of valid postsynaptic cells. The same constraint
was not enforced on presynaptic cells as they are not actually simulated in our ﬁtting
procedure. The presynaptic spikes are instead simply delivered to the synapse at the
appropriate time. This choice reduce the computational cost of ﬁtting and does not
impact on its quality in any way. That is, in our network model, and more generally
in NEURON, spikes are anyway simply represented as discrete events delivered to the
synapse, located on the postsynaptic neuron. The actual waveform of the spike has no
eﬀect on the response of the postsynaptic cell. Axonal delay is accounted for in the spike
delivery system and we added an additional 5 ms delay to compensate for the membrane
depolarization delay (see Fig. 3.10C).

          $                                %

                                               +]                                         P9


                                                                                            Q$

                                           &                                       PV




                                               +]                                        P9


                                                                                            Q$

                                                                                    PV

     Figure 3.10 – Exemplar postsynaptic cell used for ﬁtting. A, close-up of a typical layer 5
                   TTPC, selected as a candidate postsynaptic cell. B, spiking activity of cell
                   in A (top) at 5 Hz stimulation (bottom). C, spiking activity of cell in A
                   (top) at 40 Hz stimulation (bottom).

 4
     A relaxation of the electrical model for pyramidal cells was required to allow proper spiking at high-
     frequency stimulation. Speciﬁcally, the decay time constant of calcium dynamics at the soma was
     reduced to eliminate excessive hyper-polarization after a rapid series of spikes (data not shown).
     Although this modiﬁcation alters the electrical properties of the cell, its application is conﬁned to
     model ﬁtting. Furthermore, we would like to remark that the typical STDP stimulation protocols are
     quite artiﬁcial in terms of spiking activity and do not represent a typical network scenario. For this
     reason, we decided to avoid a proper re-ﬁtting of the electrical model and accepted the relaxation of
     a single parameter as a valid approximation.




52
                                                3.2 Fitting the long-term plasticity model

   After electrical screening, we focused on the geometrical constraints of typical STDP
experiments. Experimentalists tend to patch cells close to each other to increase the
chances of ﬁnding a synaptic connection, as connection probability decays rapidly with
cells distance [Perin et al., 2011]. Furthermore, to simplify the patching process, the
selected cells are often on the same focal plane (same depth in the tissue). To reproduce
these technical biases, each pair was selected by extracting a random (postsynaptic)
cell A from the pool of candidate previously deﬁned and building a bounding box of
100 m × 100 m × 10 m (x, y, z coordinates) around it. All other (presynaptic) cells in
the bounding box were identiﬁed and tested for connection with cell A. The ﬁrst match,
if any, was chosen to complete the pair. Following these constraints, we produced two
sets of 50 connections each, namely the training set and the testing set (Fig. 3.11). The

       $                                      %
                                                                     

                                                  FRQQHFWLRQFRXQW   

                                                                     

                                                                     

                                                                      

                                                                      
                                                                                           
                                                                          LQWHUVRPDWLFGLVWDQFH ͑P


                                              &
                                                                     

                                                                     
                                                  FRQQHFWLRQFRXQW




                                                                     

                                                                     

                                                                      

                                                                      
                                                                                              
                                                                           V\QDSVHVSHUFRQQHFWLRQ


  Figure 3.11 – Geometrical constraints on cell selection. A, exemplar layer 5 TTPC cells
                used for ﬁtting. B, intersomatic distance distribution for all pairs in the
                training and testing set. C, number of synapses per connection distribution
                for same pairs in B.

former was used for parameter ﬁtting via Multi-Objective Optimization (MOO), as later



                                                                                                       53
3 Modeling and parameterization of excitatory synapses

discussed, while the latter was used to asses the generalization power of the solution.

3.2.3 Parameter ﬁtting
To ﬁt the extended GB model we considered in vitro procedures in Markram et al.
[1997b], Sjöström et al. [2001] and reproduced them in silico. We focused on three
stimulation protocols5 and their outcome:
     H10minus an LTD-inducing protocol. Both the pre- and post-synaptic cells are forced
     to spike at a frequency of 10 Hz, but presynaptic current injections are delayed by 10
     ms with respect to the postsynaptic ones.
     H2plus a baseline protocol. Both the pre- and post-synaptic cells are forced to spike
     at a frequency of 2 Hz, but presynaptic current injections are anticipated by 10 ms
     with respect to the postsynaptic ones.
     H2plus a mild LTP-inducing protocol. Both the pre- and post-synaptic cells are forced
     to spike at a frequency of 10 Hz, but presynaptic current injections are anticipated by
     10 ms with respect to the postsynaptic ones.
The eﬀects of the three stimulation protocols on a typical synapse are shown in Fig.
3.12.
  A typical STDP protocol consists of three phases: (1) assessment of initial Post-
Synaptic Potential (PSP) amplitude, (2) plasticity induction, (3) assessment of ﬁnal
PSP amplitude. Assessing PSP amplitude is particularly time consuming because it
requires many trials (usually more than 30). Furthermore, in case of after-protocol PSP
assessment (i.e. phase 3 above), the eﬀects need to be evaluated after a suﬃciently long
time (up to 30 min) to be sure that the any change could be really considered long-term.
Unfortunately, simulating such a long protocol is an extremely expensive task. For this
reason, parameter ﬁtting could not be performed by exactly reproducing all the steps of
an STDP experiment. However, given that in our model the mean ΔEPSP only depends
on the ratio of depressed synapses β:
                                   (d)            (p)        (d)             (p)
                               βUSE + (1 − β)USE           βUSE + (1 − β)USE
                  ΔEPSP =          (d)             (p)
                                                       =                                  (3.36)
                               β0 USE + (1 − β0 )USE             USE
and that θd and θp thresholds are likely to be larger than calcium transients during
connectivity test, where the postsynaptic cell is not spiking, we could ﬁt the model by
simulating only the induction phase and computing β before and after it (β-strategy).
  We used a multi-objective GA Van Geit et al. [2016] to ﬁnd model parameters best
matching in vitro data (see 3.13). The ﬁnal solution found by the GA was then validated
by fully simulating all available STDP protocols (including connectivity test) on both the
training and the testing set. Our results show not only that the β-strategy was a valid
approach but also that the model well capture the typical LTP frequency-dependence
and the switch from LTD to LTP at high-frequency [Sjöström et al., 2001].
 5
     To standardize the induction procedure and to comply with the methods in more recent studies
     [Sjöström et al., 2001], we ﬁxed abs(Δt) = 10 ms, n = 15 and T = 10 s for all protocols.




54
                                                   3.2 Fitting the long-term plasticity model




                                                                           +PLQXV
       9P P9




                                                                          +SOXV
                                                                           +SOXV
                      ಜ



                                                                                 ͌S
        >&D@L൐  ͐0V




                                                                                 ͌G
                          


                          


                          


                          
                                                                                 ͕
         ͕




                          

                          
          UHOHDVHSURE




                          


                          
                                                              
                                                  WLPH V


Figure 3.12 – Eﬀects of STDP protocols on a single synapse in silico. From top to bottom:
              soma voltage of postsynaptic cell during a single repetition of the protocol,
              presynaptic spikes not shown for clarity; eﬀective calcium concentration at
              the synapse under analysis, an appropriate choice of the thresholds allows
              to correctly recognize the stimulation protocols; the LTP-inducing proto-
              col triggers an increase of synaptic eﬃcacy, while the other protocols do
              not; release probability is not immediately aﬀected by changes of synaptic
              eﬃcacy.




                                                                                          55
3 Modeling and parameterization of excitatory synapses




     A                                          B




     C                                          D




 Figure 3.13 – LTP/LTD model ﬁtting results. A, GA error evolution during ﬁtting; blue
               line represents the mean error for all candidate solutions in the population;
               gray area is the standard deviation of the population error; red line the
               error of the best solution in the population. B, Pareto front obtained after
               ﬁtting, marked in red the overall best solution. C, performance of the best
               solution on the training (blue) and testing (green) sets, compared to in vitro
               data (red), for the LTP inducing protocols; error bars represent SEM. D,
               performance of the best solution on the training (blue) and testing (green)
               sets, compared to in vitro data (red), for the LTD inducing protocols; error
               bars represent SEM.




56
                                                3.2 Fitting the long-term plasticity model

3.2.4 Pathway generalization
As previously mentioned, in the model ﬁtting step we intentionally neglected the vari-
ability of synaptic parameter, such as the density of glutamate receptors (Section 3.2.1).
   Synapses can be located on large or small spines and bear a large or small population
of calcium permeable channels. These two factors are, to ﬁrst approximation, the main
drivers of postsynaptic calcium concentration. Spine volume and Postsynaptic Density
(PSD) area are strongly correlated [Arellano et al., 2007]. This means that generally
“strong” (“weak”) synapses are located on large (small) spines. This fact provide support
for an intrinsic “calcium normalization” at the level of individual spines. We explicitly
modeled this relationship by correlating spine volume and conductance at each synapse,
as described in Section 3.1.4. However, this basic “calcium normalization” was not
suﬃcient to reproduce the in vitro STDP results (data not shown).
   To further investigate the relationship between synaptic parameters and calcium con-
centration, we simulated a long connectivity test protocol (40 minutes) and recorded
the peak of the eﬀective calcium concentration. We discovered that this value is well
correlated with gmax/vol and the relationship can be safely modeled with a linear func-
tion, as shown in Fig. 3.14A. We then hypothesized that conditioning the LTP and LTD
thresholds to gmax/vol could allow us to reproduce the in vitro results, such as
                                              gmax
                                   θ x = mx        + bx                             (3.37)
                                               vol
in particular we imposed the following two constraints:

   • mx = mpeak , same slope of peak and gmax/vol linear function

   • θx (defult) = mx default + b, as ﬁtted in step 2 (Sec. 3.2.3)

which in turn lead to the following:
                                           gmax
                                θd = 0.015      − 0.003
                                            vol                                     (3.38)
                                           gmax
                                θp = 0.015      + 0.051
                                            vol
Fig. 3.14B shows how the new threshold compare to the calcium peaks during con-
nectivity test. This generalization strategy ﬁnally allowed us to reproduce the in vitro
STDP results and recover the same performance of the mean model ﬁtted in Sec. 3.2.3,
as shown in Fig. 3.14C and Fig. 3.14D. After scaling all the thresholds using the
conductance-volume ratio, we can reproduce in silico a classical STDP experiment as in
Markram et al. [1997b] (Fig. 3.15).




                                                                                       57
3 Modeling and parameterization of excitatory synapses




     A                                          B




     C                                          D




 Figure 3.14 – Pathway generalization results. A, correlation between peaks of the eﬀective
               calcium transients during presynaptic-only stimulation and the conductance-
               volume ratio; blue line represent a linear ﬁt. B, threshold scaling as a func-
               tion of the conductance-volume ratio; red and green square represent the
               ﬁtted potentiation and depression threshold respectively. C, performance of
               the best solution after threshold scaling on the training (blue) and testing
               (green) sets, compared to in vitro data (red), for the LTP inducing proto-
               cols; error bars represent SEM. C, performance of the best solution after
               threshold scaling on the training (blue) and testing (green) sets, compared
               to in vitro data (red), for the LTD inducing protocols; error bars represent
               SEM.




58
                                                                 3.2 Fitting the long-term plasticity model




     $                                  %                                        &
                                                                                          EHIRUH
                                                                                          DIWHU
                                                                         P9                 P9



                                                                 PV               PV
                                        '
                                                               
                                          (363DPSOLWXGH P9




                                                               

                                                               

                                                               

                                                               
                                                                                                
                                                                              WLPH PLQ

Figure 3.15 – Classical STDP experiment in silico. A, pair of layer 5 TTPC used for
              the in silico experiment; presynaptic cell in blue, postsynaptic cell in green,
              synapses in red. B, spiking activity of the pre- and postsynaptic cell dur-
              ing stimulation at 40 Hz; the green trace has only illustration purposes, as
              the presynaptic cell is not simulated in our experiments. C, mean EPSP
              amplitude before (black) and after (red) the pairing. D, evolution of EPSP
              amplitude over time; gray dots represent individual trials, black dots the
              mean over 5 minutes of the individual trials; error bars represent SEM.




                                                                                                        59
3 Modeling and parameterization of excitatory synapses

3.2.5 Circuit generalization
The ﬁnal step towards a comprehensive model of calcium-based plastic changes is the
generalization to arbitrary pathways. Unfortunately the experimental validation of this
step is quite problematic, given that only a few pathways have been experimentally
characterized.
  We then assumed the main discriminants between pathways, with respect to synaptic
plasticity, to be the innervation proﬁle and the morphology types rather than a diﬀerent
set of plastic mechanisms. We conﬁrmed the validity of this assumption by reproducing
an in vitro experiment on soma distance-dependent plasticity, as later shown in Section
3.3.1. For this reason, a minimal variation of our ﬁtted solution, accounting only for the
diﬀerences in synaptic physiology, would be suﬃcient to parametrize any uncharacterized
pathways.
  The generalization procedure computes β, Used and Usep for all pathways, knowing
only the mean release probability Use. We used the Sequential Least SQuares Program-
ming (SLSQP) method [Kraft, 1988], as implemented in the Scipy package, to ﬁnd a
suitable solution to the minimization problem:

                         minimize Use − [βUsed + (1 − β)Usep ]                      (3.39)
                         β,Used ,Usep

The solution to this problem is expected to have several local minima. In the absence
of other constraints, we imposed the same β, Used and Usep used for layer 5 TTPC, as
the initial point for the optimization.


3.3 Validation and predictions
We previously showed how the long-term plasticity model was ﬁtted to available data
and how it was then generalized to parametrize a whole neocortical circuit (see Sec. 3.2).
We now focus on the validation of this strategy and on the predictions that can be made
about uncharacterized pathways.
   Validating our generalization strategy is a diﬃcult task, mostly because of the sparsity
of experimental data. In Section 3.3.1 we prove the quality of our parameterization by re-
producing in silico the results on distance-dependent synaptic plasticity by Sjöström and
Häusser [2006]. The idea of distance-dependent plasticity can be summarized as follows:
the direction and magnitude of plastic changes induced by a given stimulation protocol
depend on both the protocol itself and the innervation pattern of the connection. By
showing that distance-dependent plasticity is a consequence of our parametrization, even
though the dataset by Sjöström and Häusser [2006] was never considered in the ﬁtting
procedure, we provide direct evidence of the validity of our approach. Furthermore, we
show how synaptic dynamics, neuronal morphologies and innervation patterns jointly
determine the nature of plastic changes, demonstrating once more the need of such a
level of modeling detail to even consider the study of higher brain functions, such as
learning.




60
                                                             3.3 Validation and predictions

   In Section 3.3.2 we extend our STDP analysis to totally uncharacterized pathways.
In particular, we simulated the outcome of typical STDP protocols on every possible
excitatory to excitatory pathway in our circuit. Such a dataset do not exist in experi-
mental literature and contains interesting indications on the plastic regimes observable
in a network. Our goal is to show how the proposed model and parametrization can be
used to obtain insights on network functioning and to guide future in vitro and in vivo
experiments.

3.3.1 Distance-dependent synaptic plasticity
It has been shown that a standard 50 Hz STDP protocol, triggering potentiation of
layer 5 TTPC connections, results in depression of distal synapses, i.e. layer 2/3 to
layer 5 pyramidal connections [Sjöström and Häusser, 2006]. This ﬁnding conﬁrms the
hypothesis that excitatory pathways might exhibit diﬀerent forms of synaptic plasticity
as the result of diﬀerent innervation patterns rather than diﬀerent plasticity rules.
   To test this hypothesis and validate our ﬁtting generalization, we reproduced in silico
the experiments in Sjöström and Häusser [2006]. We randomly selected 50 connections
between layer 5 TTPC and 50 connections between layer 2/3 pyramidal cells (PCs). To
maximize the mean synaptic distance span we searched for cells in the half-top part
of layer 2/3 and the half-bottom part of layer 5. Plastic changes were triggered as in
previously described [Markram et al., 1997b, Sjöström et al., 2001, Sjöström and Häusser,
2006].
   Our results conﬁrm the LTD switch of distal connections, while the proximal ones
still exhibited LTP (Fig. 3.16 and Fig. 3.17). Furthermore, we extended our analysis to
single synapses and observed a clear transition between LTP and LTD as a function of
synaptic distance from soma (Fig. 3.18).




                                                                                           61
3 Modeling and parameterization of excitatory synapses




     A                                          B




     C                                          D




 Figure 3.16 – Example of distance-dependent synaptic plasticity. A, evolution of EPSP
               amplitude over time for a layer 5 TTPC to layer 5 TTPC connection, stim-
               ulated at 50 Hz; empty circles represent individual trials, ﬁlled circles the
               mean over 5 minutes od the individual trials; error bars represent SEM. B,
               mean EPSP amplitude before (blue) and after (red) the pairing for the con-
               nection in A. C, evolution of EPSP amplitude over time for a layer 2/3 PC
               to layer 5 TTPC connection, stimulated at 50 Hz; empty circles represent
               individual trials, ﬁlled circles the mean over 5 minutes od the individual tri-
               als; error bars represent SEM. D, mean EPSP amplitude before (blue) and
               after (red) the pairing for the connection in C.




62
                                                           3.3 Validation and predictions




Figure 3.17 – Summary of distance-dependent synaptic plasticity. Potentiation switches to
              depression for distally connected cells. Triangles show layer 5 TTPC to layer
              5 TTPC connections, ﬁlled circles layer 23 PC to layer 5 TTPC connections;
              dashed red line is a polynomial ﬁt to the data.




Figure 3.18 – Distance-dependent synaptic changes. Distal synapses are more prone to
              depression and almost never potentiate. Single data points from a random
              sample of 150 synapses; layer 5 TTPC to layer 5 TTPC connections repre-
              sented as triangles; layer 2/3 PC to layer 5 TTPC represented as circles.




                                                                                          63
3 Modeling and parameterization of excitatory synapses

3.3.2 STDP map
A systematic characterization of STDP in the neocortex is nowadays missing. Long-term
synaptic plasticity has been experimentally studied only at layer 5 TTPC connections
and a few other neocortical pathways [Markram et al., 1997b, Sjöström et al., 2001,
Froemke et al., 2005, Sjöström and Häusser, 2006] 6 , due to the many technical challenges:


      • Connections are rare or hard to ﬁnd for most pathways, especially if the pre- and
        post-synaptic cells belong to diﬀerent cortical layers;

      • Experiments are inherently long and have high failure rate;

      • (Sadly enough) only a full characterization of all pathways or the discovery of
        a computationally relevant diﬀerence between some of them would be worth the
        eﬀort.

However, it is known that the direction and magnitude of synaptic changes induced by
paired neural activity profoundly depend on the morphological features of the neurons
and the innervation patterns. For this reason, we could expect great variability of STDP
dynamics across pathways (see also Sec. 3.3.1).
   In this section we present our attempt to ﬁll the gap by generating a simulated STDP
map. The map not only provides a comprehensive view of STDP in the neocortex, but
also it could be used as a reference to guide future experiments. We simulated the
outcome of three STDP protocols (2, 10 and 50 Hz) for all 169 excitatory to excitatory
pathways in our circuit reconstruction [Markram et al., 2015]. The protocols have been
chosen to be representative of typical low-, medium- and high-frequency stimulation in
standard STDP experiments [Markram et al., 1997b, Sjöström et al., 2001].
   The STDP maps for the three protocols are shown in Fig. 3.19. As expected only few
pathways showed any change after the 2 Hz protocol, while most of them potentiated at
50 Hz stimulation. The 10 Hz protocol produced mixed outcomes, with a tendency for
depression. Furthermore, we observed some pathways which did not potentiate under
any circumstance, i.e. layer 2/3 PC to layer 4 SP connection.
   The variability across pathways is even more evident when comparing the classical
frequency-dependence curves (see Fig. 3.20). As previously mentioned, the 10 Hz pro-
tocol is particularly interesting because it diﬀerentially induces LTP, LTD or no change
depending on the identity of the pre- and post-synaptic neurons.




 6
     We do not consider studies lacking a precise identiﬁcation of the pre- and post-synaptic cells, i.e.
      studies employing extracellular electrodes for axonal stimulation.




64
                                                       3.3 Validation and predictions


                 A




                 B




                 C




Figure 3.19 – STDP map. A Stimulation at 2 Hz. B Stimulation at 10 Hz. C Stimulation
              at 50 Hz.




                                                                                   65
3 Modeling and parameterization of excitatory synapses




         Figure 3.20 – LTP frequency dependence for uncharacterized pathways.




66
                                  3.4 NMODL implementation and HPC considerations

3.4 NMODL implementation and HPC considerations
As mentioned at the beginning of this chapter, GluSynapse was implemented in NEU-
RON model description language (NMODL) for the NEURON simulator [Carnevale and
Hines, 2006]. Here we discuss some of the technical details of model implementation
and simulation. In particular, in Sec. 3.4.1 we describe how to account for dendritic
spines without aﬀecting simulation performance. In Sec. 3.4.2 we review integration
methods in NEURON and our strategy to reduce simulation time. Finally in Sec. 3.4.3
we propose an eﬃcient implementation of the synapse elimination mechanism.
  This section often relies on some basic knowledge of the NEURON simulator’s in-
ternals. For this reason, we encourage the reader to consult “The NEURON Book”
[Carnevale and Hines, 2006] whenever necessary.

3.4.1 Spine as virtual compartment
The most natural way to model dendritic spines in NEURON is to use one or more
dedicated compartments. A popular approach consists of one compartment for the
spine neck and one for the head. This choice allow to control in great details several
aspects of spines’ internal dynamics, such as diﬀusion or neck resistance. However, in
large network simulations it is practically unfeasible to explicitly model every spine.
A typical pyramidal neuron consists of a few hundreds of compartments but features
several thousands of spines. NEURON heavily rely on the assumption that the matrix
of compartments is relatively small and easy to invert. The explosion of the number of
compartment increases the complexity of the problem quite dramatically.
   For this reason, we did not implement spines as actual compartments bur rather as
virtual ones. That is, the synapse model itself takes care of every aspect related to the
spine, without the need to add more compartments to the simulation: Ionic diﬀusion is
neglected and calcium dynamics are explicitly managed (see 3.1.5); VDCC are imple-
mented as a point mechanism and current density estimation is carried on internally, as
opposed to the usual density mechanisms inserted into the membrane and directly man-
aged by NEURON; calcium compartmentalization do not rely on NEURON’s internal
features for ionic currents to avoid aﬀecting the parent dendrite; Spine neck is neglected
as it seems to have a small impact in our simulations (data not shown).
   Due to the extended “duties” of the synapse model, Glusynapse can be considered
more of a synaptic simulator within NEURON, rather than just a regular model. Overall
this approach allows our simulations to be as eﬃcient as those without spines and to
capture most of the spine’s aspects relevant for synaptic plasticity.

3.4.2 GluSynapse integration methods
Before digging into the integration methods used in this work, a premise on what is
available in NEURON seems required. For an exhaustive review of NEURON’s integra-
tion methods and a deeper analysis of numerical error, please refer to “The NEURON
Book” [Carnevale and Hines, 2006].




                                                                                       67
3 Modeling and parameterization of excitatory synapses

  NEURON oﬀers several ﬁxed time step integrators, speciﬁed in NMODL at model
design time:
     cnexp An implicit method for solving linear, non coupled, ODE of type ẋ = f (v, x).
     This method is particularly eﬃcient and relies on the existence of an analytical solu-
     tion, i.e. x = e− τ + c;
                       t



     derivimplicit An implicit method for solving general ODE. This method trades per-
     formance for robustness and range of applicability;

     euler An explicit method based on forward Euler. This method is extremely fast
     but its usage is discouraged, as it requires very small time steps to prevent numerical
     instability and errors accumulate over time.
Other choices are available, but they are either less common (i.e. runge method) or
designed for special purposes, such as solving kinetic schemes (i.e. sparse method).
  NEURON provides also two adaptive integrators: CVODE by SUNDIALS [Hindmarsh
et al., 2005] and DASPK [Brown et al., 1994]. The former is particularly interesting
because it supports local variable time step, a feature that can dramatically improve the
performance in certain contexts, such as simulations involving few AP to propagate7 .
  Network simulations performed at the BBP involve a large number of spike exchanges
and for this reason they are not suited for adaptive integration, at least for now. Stan-
dard, non-plastic, synapses are then integrated using the fast cnexp method, as their
equations satisfy the constraints imposed by this solver. GluSynapse instead is a much
more complex model, involving non-linear and/or coupled diﬀerential equations, and
requires the derivimplicit integration method. Unfortunately, the introduction of
GluSyanapse in the network model caused a major performance drop which could not be
tolerated, as it increased the computational requirements beyond our current resources.
To solve the problem, we decided to temporarily switch to the euler integrator for
GluSynapse. We assumed that the risk of numerical instability would be mitigated by
the small integration step used in all our network simulations (dt = 0.025 ms) . Further-
more, we explicitly veriﬁed that the integration error is small compared to the standard
precision of our simulations (data not shown).
  The euler integrator was only used for network simulations. Model ﬁtting (Sec. 3.2)
instead was a perfect use case for the fast variable time step implicit integrator CVODE.

3.4.3 Inhomogeneous Poisson elimination
Synapse elimination was modeled as an inhomogeneous Poisson process, as the pruning
probability is a function of synaptic eﬃcacy (see Sec. 3.1.7). For this reason, we cannot
assume an exponential distribution for the event intervals. How to calculate the intervals
for an inhomogeneous Poisson process is an active area of research and we do not have
 7
     Spike exchange is a global synchronization barrier: after spike delivery at the postsynaptic side, the
      system needs to be re-initialized and every “jump” performed by the adaptive solver to be undone.
      For this reason, network simulations might not beneﬁt from local variable time step integration if
      too many spikes are generated.




68
                                 3.4 NMODL implementation and HPC considerations

a general formula yet. However, under certain assumptions, we can approximate the
inter-arrival interval as an exponential distribution.
   The ﬁrst obvious approach would be to draw a random sample every time step, as the
synaptic eﬃcacy variable ρ evolves on much longer time scales and it could be considered
constant in between sampling. Unfortunately drawing so many random numbers would
be too expensive to be practical.
   A better solution can be designed if we consider that (a) the average ﬁring frequency
in vivo is approximatively 1 Hz and (b) abrupt changes in ρ require the NMDAR to be
involved. That is, in between spikes the ρ variable can be assumed to be constant, as the
inter-spike interval is approximatively 1 s and any extreme change in ρ happens almost
instantaneously after an incoming spike. This approach greatly reduces the number of
random sample to be drawn without introducing a big modeling error.




                                                                                      69
4 Eﬀects of functional synaptic plasticity
  on neural circuits
                                                     “Too much light often blinds
                                                     gentlemen of this sort. They cannot
                                                     see the forest for the trees.”
                                                           Christoph Martin Wieland,
                                                                            Musarion

In this chapter we describe the impact of introducing functional synaptic plasticity in
a digital reconstruction of neocortical tissue [Markram et al., 2015]. We simulated the
same neural network under two conditions: spontaneous activity and thalamocortical
stimulation. In the former case, the circuit was left free from external stimuli for the
whole duration of the simulation (1050 s). In the latter instead, we provided a external
stimulus through the innervating thalamocortical ﬁbers. The circuit was simulated at
in vivo-like calcium levels (1.2 mM).
   In the next sections we describe how the circuit and the stimuli were deﬁned (Sec.
4.1). Then we present several preliminary analyses, required to gather initial informa-
tion on network behavior during plastic simulations (Sec 4.2). Finally, we show how
the introduction of functional plasticity induced strengthening of connections between
neurons similarly innervated, as previously observed in vitro [Perin et al., 2011] (Sec.
4.3).


4.1 Fundamentals of microcircuit reconstruction and simulation
The process of reconstructing and simulating a neocortical microcircuit is described else-
where Markram et al. [2015]. In brief, reconstructed neurons were placed in space accord-
ing to data on dendrites, proportions and morphological features. Electrical types were
assigned by ﬁtting the reconstructed neurons to in vitro electrophysiological recordings
[Druckmann et al., 2007, Hay et al., 2011, Van Geit et al., 2016]. Neurons were con-
nected using a biophysically constrained probabilistic algorithm [Reimann et al., 2015].
Synaptic physiology was assigned to match the known features of each viable synaptic
pathway.
  The microcircuit used in this work is composed of 224000 neurons, organized in seven
cortical columns and six layers (Fig. 4.1A). The columns are arranged in a center-
surround fashion: a single central column, the main subject of all our analyses, is sur-
rounded by six other columns, providing a realistic neural environment without major




                                                                                       71
4 Eﬀects of functional synaptic plasticity on neural circuits

edge eﬀects. All the cortico-cortical excitatory synapses in the central column are plas-
tic, while those in the other columns, the thalamocortical synapses and all inhibitory
synapses are not. Thalamocortical ﬁbers innervate the seven columns matching known
constraints on target layers, synaptic physiology and number of synapses per connection
[Markram et al., 2015]. To provide spatially deﬁned stimuli to the network, we bundled
the ﬁbers in 100 clusters, based on the ﬁbers position (Fig. 4.1). Fibers belonging to
the same cluster are assumed to carry the same information into the cortex (see Sec.
4.1.1) The microcircuit was simulated using the NEURON simulator on the MIRA and

                $                              %




                                                                   ͑P


  Figure 4.1 – Circuit layout. A, the circuit used in our simulations is composed by seven
               columns organizes in a mosaic. Excitatory synapses in the central column
               are plastic. Thalamocortical ﬁbers (one exemplar showed as a green line)
               innervates the circuit matching several biological constraints [Markram et al.,
               2015]. B, neighbor thalamocortical ﬁbers are functionally clustered into 100
               groups. Fibers belonging to the same cluster are always stimulated together
               and represent a “pixel” of the network global receptive ﬁeld.

BBPIV supercomputers.

4.1.1 Thalamocortical stimulation
Thalamocortical stimuli were generated using a custom Python module. Each stimulus
is fed to the generator as a dotted pattern impinged on the microcircuit mosaic (see Fig.
4.2). The stimulus generator produces Poissonian spike trains for every ﬁber. The ﬁring
rate of every ﬁber is determined by the presence or absence of a “dot” in its receptive
ﬁeld.
   The whole set of stimuli used in this work is shown in Fig. 4.3. We generated
eight diﬀerent stimuli with identical statistical properties (i.e. number of dots, ﬁring
frequency). During simulation, the eight patterns were presented in a random order
every 200 ms. Every 30 s we left the network unstimulated for 1 s to test the eﬀects of
synaptic depression on the stimulus. That is, thalamocortical ﬁbers could recover from




72
                                    4.1 Fundamentals of microcircuit reconstruction and simulation


                               
        PHDQILULQJUDWH +]
                                


                                                                                     +]

                                
                                                                              PV


                                


                                                                                     +]

                                            ͑P                            PV



  Figure 4.2 – Exemplar stimulus. A randomly generated dotted pattern is converted into a
               Poissonian spike train for every thalamocortical ﬁber. Fibers having a dot in
               their receptive ﬁelds rapidly increase their ﬁring rate in the ﬁrst 75 ms after
               the presentation of the pattern (red ﬁring rate histogram). Unstimulated
               ﬁbers have instead a low and constant ﬁring rate (blue ﬁring rate histogram).


depression and provide a stronger stimulus to the cortex, emulating the abrupt arrival
of a stimulus.




                                                                                               73
4 Eﬀects of functional synaptic plasticity on neural circuits


                $                    %                   &                    '




                (                    )                   *                    +




                                                                                  ͑P




                                                            
                                     PHDQILULQJUDWH +]

  Figure 4.3 – Summary of all stimuli. Eight stimuli were randomly generated and labeled
               alphabetically. The stimuli share the same statistics in terms of intensity and
               stimulated area.


4.2 Preliminary analyses
Before evaluating the computational implications of having functional synaptic plasticity
in the network, we carried on several preliminary analyses. Our goal was to establish
whether the direction and magnitude of plastic changes have any global trends in the
network, to identify regions of interest and possible boundary conditions to take into
account.

4.2.1 Network balance
The main concern in any network simulation featuring Hebbian-like excitatory plasticity
is the tendency to runaway excitation, due to the disruption of excitation-inhibition
balance. To prevent this eﬀect, these models are usually coupled with some sort of
weight normalization mechanisms, resembling an homeostatic process. Although what
we presented in this work is not a traditional Hebbian model, we cannot rule out a priori
the possibility of runaway excitation. To verify the stability of excitation-inhibition
balance in out simulations, we analyzed the evolution of ﬁring rate over time in both the
spontaneous and stimulated condition.
   In the spontaneous case, activity levels remained constant throughout the whole du-



74
                                                                    4.2 Preliminary analyses

ration of the simulation (Fig. 4.4). The overall ﬁring rate almost never exceeded 2 Hz,
except for the initial network burst. This event was caused by the fully recovered state
of all the synapses at the beginning of every simulation. We did not try to prevent
network bursts as they are, in limited amount, part of the normal functioning of a neu-
ral network. Runaway excitation rather refers to the generation of these events at high
frequency, preventing the circuit from doing anything else.

       $




       %




  Figure 4.4 – Firing rate evolution during spontaneous network activity. A, Firing rate
               over time for all cells in the central column (bins size 20 ms); the y axis was
               truncated to 4 Hz to better resolve the low ﬁring rate; red square marks the
               time window in B; inset shows the initial network burst, axis units as in the
               main panel. B, close-up on network activity in the red window marked in A;
               blue line shows the ﬁring rate histogram (left y axis, truncated at 1.75 Hz);
               green dots show the spikes of individual neurons, sorted by depth in the tissue
               (right y axis).

  During thalamocortical stimulation, the circuit showed a more marked tendency to
burst (Fig. 4.5). However, the periodicity of these events did not reach alarming levels




                                                                                             75
4 Eﬀects of functional synaptic plasticity on neural circuits

(see network burst in Markram et al. [2015]). Vice versa, to our surprise, they acted
as an intrinsic homeostatic mechanism. As we later show in Sec. 4.2.2, network busts
bring a large fraction of layer 5 Thick-tufted pyramidal cells (TTPC) synapses down to
an intermediate eﬃcacy range. This eﬀect seems to prevent uncontrolled potentiation
and so the explosion of network activity. We then tested weather a speciﬁc stimulus

       $




       %




  Figure 4.5 – Firing rate evolution during stimulated network activity. A, Firing rate over
               time for all cells in the central column (bins size 20 ms); the y axis was
               truncated to 20 Hz to better resolve the low ﬁring rate; red square marks the
               time window in B; inset shows the initial network burst, axis units as in the
               main panel. B, close-up on network activity in the red window marked in
               A; blue line shows the ﬁring rate histogram (left y axis, truncated at 10 Hz);
               green dots show the spikes of individual neurons, sorted by depth in the tissue
               (right y axis); letters on top mark the stimulus applied at the corresponding
               time.

was more likely to generate a burst. Indeed half of them were generated by the stimulus
labeled as A (Fig. 4.6). We are currently investigating why this speciﬁc input pattern
was burst-prone, but no obvious reason was found yet.



76
                                                                            4.2 Preliminary analyses


                                        



                HYRNHGQHWZRUNEXUVWV   



                                        



                                        



                                        



                                         



                                         
                                             $   %   &   '     (    )   *       +
                                                         VWLPXOXV


  Figure 4.6 – Burst incidence per stimulus. Stimulus pattern “A” triggered half of the total
               network bursts.


   The careful reader will have noticed that our model of glutamatergic synapse lack an
explicit homeostatic component. We believe that the need of an homeostatic term to
balance most Spike Timing Dependent Plasticity (STDP) learning rules has created the
false expectation of a speciﬁc biological mechanism to counterbalance the emergence of
strong weights and runaway excitation. Based on our data, we suggest that the picture
is actually simpler than that. The network itself tend to stabilize its state, as network
bursts eﬀectively dump the “weights”. Obviously we do not question the existence
of homeostatic plasticity as a mechanism to rescue from catastrophic events, such as
brain lesions. We instead question the need of an homeostatic term to explain normal
network functioning. One might argue that our model does indeed have an homeostatic
mechanism, as the release probability cannot grow above the prescribed bounds (see
Sec. 3.1.6). However, we believe that assuming the existence of physiological bounds
for synaptic properties, such as release probability, is rather diﬀerent than the common
weight normalization.

4.2.2 Redistribution of synaptic eﬃcacy
We then examined how single synapses on the same postsynaptic neuron evolved over
time. For this analysis we considered a typical layer 5 TTPC. In the spontaneous case the



                                                                                                 77
4 Eﬀects of functional synaptic plasticity on neural circuits

vast majority of synapses converged to their ﬁxed points within the time of simulation
(Fig. 4.7). For this reason the initially unimodal distribution of release probabilities
became strongly bimodal. During thalamocortical stimulation, synapses showed more

       $                    %
                                                      




                                UHOHDVHSUREDELOLW\
                                                      

                                                      

                                                      

                                                      

                                                      


                                                                                          
                                                                                WLPH V


                            &                                     W VHF                 W VHF
                                                      
                                V\QDSVHFRXQW




                                                      


                                                      


                                                      


                                                         
                                                                                       
                                                             UHOHDVHSUREDELOLW\      UHOHDVHSUREDELOLW\


  Figure 4.7 – Redistribution of synaptic eﬃcacy during spontaneous activity. A, 2D mor-
               phological reconstruction of the cell used in this analysis. B, release probabil-
               ity evolution for a small sample of excitatory synapses on the cell (3 % shown).
               C, release probability distribution before (left) and after (right) spontaneous
               activity simulation.

interesting dynamics (Fig. 4.8). In particular, a third (apparent) equilibrium point
appeared for the release probability variable. This new point, approximatively located
at 75 % of the release probability range, attracted a large fraction of the available
synapses. It was interesting to notice that the existence of this point became clear
after the second network burst (see 4.5). To complement our observations, we analyzed
the spatial distribution of synapses on the cell morphology after simulation (Fig. 4.9).
We found that synapses having release probability within a standard deviation from
the mean were mostly located on the basal and oblique dendrites. This ﬁnding was



78
                                                                                            4.2 Preliminary analyses

       $                   %
                                                     




                               UHOHDVHSUREDELOLW\
                                                     

                                                     

                                                     

                                                     

                                                     


                                                                                         
                                                                               WLPH V


                           &                                     W VHF                 W VHF
                                                     
                               V\QDSVHFRXQW




                                                     


                                                     


                                                     


                                                        
                                                                                      
                                                            UHOHDVHSUREDELOLW\      UHOHDVHSUREDELOLW\


  Figure 4.8 – Redistribution of synaptic eﬃcacy during thalamocortical stimulation. A,
               2D morphological reconstruction of the cell used in this analysis. B, release
               probability evolution for a small sample of excitatory synapses on the cell (3
               % shown). C, release probability distribution before (left) and after (right)
               spontaneous activity simulation.


particularly interesting because it follows the characteristic innervation pattern of layer
5 TTPC to layer 5 TTPC [Markram et al., 1997a]. We are currently investigating what
functional implications this peculiar distribution might have.




                                                                                                                 79
4 Eﬀects of functional synaptic plasticity on neural circuits


                                  W VHF
                                                                   




                                                                   




                                                                             UHOHDVHSUREDELOLW\
                                                                   




                                                                   




  Figure 4.9 – Synaptic eﬃcacy distribution of a typical layer 5 TTPC. Synaptic locations
               are marked by circles and divided in three groups, depending on release proba-
               bility. Weak synapses in blue (below one standard deviation from the mean);
               average synapses in green (within one standard deviation from the mean);
               strong synapses in red (above one standard deviation from the mean). Most
               average release probability synapses are located on the basal and oblique den-
               drites.


4.2.3 Layer speciﬁc changes
We next examined how the origin layer of each cell aﬀects the average plastic change (see
Fig. 4.10). In the spontaneous case, the mean release probability did not change in any
layer, while the variance increased for all of them. This is a clear sign that these synapses
are simply converging to their stable ﬁxed points without major reorganization (see Sec.
4.2.2). In the stimulated case instead, only layer 2 and layer 3 remained “neutral”. The
result was largely expected, as ﬁring activity in layer 2 and 3 is very low even during
stimulation and so unlikely to cause large plastic changes. This phenomenon is partially
due to physiological features of neurons in these layers and partially to the missing
lateral spread of axons, not yet fully captured in our network model. The lower layers
are instead much more engaged, as reﬂected by the smaller variance change. Overall,
layer 5 cells showed a tendency to potentiate while layer 4 and 6 cells to depress. For
an in depth analysis of pathway-speciﬁc changes see Sec. 4.2.2.




80
                                                                             4.3 Common neighbors strength bias


       $                                                  %
                                                                        

                                                                        

                                                                        
           ̪PHDQ




                                                              ̪PHDQ
                                                                        

                         ಝ                                             ಝ

                         ಝ                                             ಝ

                         ಝ                                             ಝ

                                                                      

                                                                      

                                                                      
            ̪YDULDQFH




                                                               ̪YDULDQFH



                                                                      

                                                                      

                                                                      

                                                                      
                                                                                                
                                          OD\HU                                              OD\HU


  Figure 4.10 – Average release probability change after 10 minutes of simulation. A sponta-
                neous activity simulation. Mean (top) and variance (bottom) of the neuron
                average change, grouped by soma layer. B same of A but for the stimulated
                case.


4.3 Common neighbors strength bias
It has been shown that the Excitatory Post-Synaptic Potential (EPSP) amplitude be-
tween connected pyramidal cells in layer 5 is correlated with the number of common
neighbors Perin et al. [2011]. Given any pair of cells i and j, a common neighbor is
a cell k connected to both i and j (see Fig. 4.11A). That is, the average connection
strength increases with the number of common neighbors. Whether the common neigh-
bors strength bias is genetically prescribed or the consequence of plastic mechanisms at
work is currently unknown.
   The digital reconstruction of a neocortical circuit used in this work satisﬁes several
constraints on synaptic physiology (see Markram et al. [2015]), but the common neigh-
bors strength bias is not one of those. To test the hypothesis that this bias arises




                                                                                                             81
4 Eﬀects of functional synaptic plasticity on neural circuits

spontaneously in a plastic neural network, we analyzed the relationship between number
of common aﬀerent1 and mean release probability for a subset of TTPC connections in
layer 5. A pool of 300 connections was randomly selected from the central column.
Common aﬀerents instead were estimated considering all layer 5 TTPC, without any
restriction on the cortical column. To assess whether number of common aﬀerents and
mean release probability were correlated, we considered the Spearman’s rank correlation
coeﬃcient [Zwillinger and Kokoska, 2000]. This measure was chosen as it does not re-
quire any assumption on the marginal distributions of the paired variables to test (i.e.
normality). We found that in both the stimulated case (Fig. 4.11) and spontaneous ac-
tivity (Fig. 4.12) the bias emerged at the end of the simulation (t = 1050 s), particularly
strong in the former case.




 1
     A common aﬀerent is simply a common neighbor satisfying the property of being pre-synaptic to
     both the target neurons. Perin et al. [2011] showed that the strength bias is even more marked when
     common aﬀerents rather than common neighbors are considered.




82
                                                                                                    4.3 Common neighbors strength bias




$                                                         %                              

                        N                                                                 




                                                              FRQQHFWLRQFRXQW
                                              M                                           

                                                                                          

                                                                                          

                                L                                                          
                                                                                                                
                                                    K                                                  &$/773&


&                                                         '                              
                       
                                                              PHDQUHOHDVHSUREDELOLW\




                                    LQLW                                                                U 
                                 ILQDO
    FRQQHFWLRQFRXQW




                                                                                         
                       
                                                                                         
                       

                                                                                                     U 
                       

                                                                                       

                         
                                                                                            
                             PHDQUHOHDVHSUREDELOLW\                                                  &$/773&


     Figure 4.11 – Emergence of common neighbors strength bias (talamocortical stimulation).
                   A, example of common neighbor; neuron k is both a common neighbor (CN)
                   common aﬀerent (CA) for neurons i and j, while neuron h is not. B, common
                   aﬀerent (CA) distribution for the 300 connections considered in this analysis;
                   CA were computed considering the full population of layer 5 TTPC. C, mean
                   release probability distribution before (t=0) and after (t=1050) simulation.
                   D, common neighbors strength bias; after simulation a strong correlation
                   between common aﬀerents and mean release probability emerged.




                                                                                                                                   83
4 Eﬀects of functional synaptic plasticity on neural circuits




$                      
                                                           %                              
                                                               PHDQUHOHDVHSUREDELOLW\



                                                LQLW
                                             ILQDO                                     
    FRQQHFWLRQFRXQW




                                                                                       

                                                                                       
                       
                                                                                          
                                                                                                        U 
                                                                                          
                                                                                                          U 
                         
                                                                                                     
                             PHDQUHOHDVHSUREDELOLW\                                                &$/773&


     Figure 4.12 – Emergence of common neighbors strength bias (spontaneous activity). A,
                   mean release probability distribution before (t=0) and after (t=1050) sim-
                   ulation. B, common neighbors strength bias; after simulation a weak corre-
                   lation between common aﬀerents and mean release probability emerged.




84
5 Conclusion
                                                      “If you can meet with Triumph and
                                                      Disaster and treat those two
                                                      impostors just the same.”
                                                                     Rudyard Kipling, If

In this work we integrated state of the art data and theories on excitatory synaptic plas-
ticity to design a comprehensive model, suitable for large scale simulations of neocortical
circuits. Our main goal was to unify the fragmented knowledge on the subject and to
produce a centralized resource to study the biological foundation of learning and mem-
ory in a computer-based framework. We started by sketching a simpliﬁed description of
the synapse life-cycle on the basis of available experimental literature (see Chapter 2).
The model consists of a discrete set of states and transitions. Each state represent a
well deﬁned axo-dendritic conﬁguration, such as depressed synapse or bouton-ﬁlopodium
apposition. Transitions between states model the mechanisms of synaptic plasticity, the
main subject of this work. From the high-level description of the synapse life-cycle we
derived a mathematical model of plastic synapse and implemented it in NEURON (see
Chapter 3). Whenever possible, we relied on existing components and theories from
literature, limiting extensions to the minimum required. The Long Term Potentiation
(LTP)/Long Term Depression (LTD) component of the model was then parameterized
using a custom Genetic Algorithm (GA)-based strategy to reproduce in silico the ef-
fects of three standard STDP protocols on layer 5 TTPC connections [Markram et al.,
1997b]. The parameterization obtained for layer 5 TTPC connections was generalized
to all other excitatory-to-excitatory pathways. Small variations of the ﬁtted parameters
were necessary to adapt the model to the synaptic physiology of the other pathways.
Finally, we tested the eﬀects of functional and structural synaptic plasticity on neural
circuits by simulating a large-scale reconstruction of neocortical tissue [Markram et al.,
2015]. We considered two network conditions: spontaneous activity and thalamocortical
stimulation (see Chapter 4).
   The integrative nature of this thesis led us to challenge and revise many common as-
sumptions on synaptic plasticity. For example, accounting for stochastic vesicle release
and postsynaptic calcium dynamics had profound implications on the requirements for
LTP induction. The Graupner-Brunel (GB) model [Graupner and Brunel, 2012], revis-
ited in this work, explains the LTP’s frequency-dependence by postulating postsynaptic
calcium accumulation at high frequency-stimulation. This idea is an extension of the
popular “calcium control hypothesis” [Lisman, 1989]: prolonged but moderate elevation
of intracellular calcium triggers LTD, while larger transients result in LTP. However,
postsynaptic calcium accumulation is incompatible with the physiological features of



                                                                                          85
5 Conclusion

neocortical synapses. That is, fast calcium buﬀering dynamics, combined with low
presynaptic release probability and a small Readily-Releasable Pool (RRP), result in
rapid calcium clearance at the spine, even at high-frequency stimulation (see Sec. 3.1.6
for more details). We then concluded that the interplay of several calcium buﬀers is
a requirement for LTP’s frequency-dependence, in line with the results of other studies
[Holmes and Levy, 1990, Markram et al., 1998] (for a review see [Manninen et al., 2010]).
Our model captured LTP’s and LTD’s frequency dependence reported by Markram et al.
[1997b], Sjöström et al. [2001], without explicit ﬁtting on these properties. Further-
more, it reproduced the results on distance-dependent synaptic plasticity [Sjöström and
Häusser, 2006], even though this dataset was also never used during ﬁtting. In network
simulation, we observed the strengthening of connections between neurons similarly in-
nervated, as previously shown in vitro [Perin et al., 2011]. These results show how our
approach, based on the systematic enforcement of biological constraints, actually sim-
pliﬁes the problem of modeling synaptic dynamics: a small input dataset and a few new
constraints were suﬃcient to capture a broad range of long-term plastic dynamics for
several pathways.
   The main product of this work is a simulation-based framework to study synaptic
plasticity in morphologically detailed neural circuits. We hope that the availability of
such a resource will trigger more theoretical and experimental work on the subject. In
particular, our framework could be used to test speciﬁc hypotheses on network function-
ing during stimulation, other than the one we are already considering. Furthermore, the
modular structure of the synapse model allows the replacement of any component with
a new formulation. For this reason, GluSynapse could serve as a scaﬀolding environment
for the design of other models of speciﬁc synaptic mechanisms. Finally, experimental
work could be simpliﬁed by testing in the network model the expected outcome of a large
set of protocols. This approach could greatly reduce the burden of experiment design
and guide the intuition of researchers even during the experiments.
   Regardless of the speciﬁc implementation, choices and assumptions discussed in this
thesis, we believe that the underlying idea of this work has an important scientiﬁc value.
The brain is an organ and a dynamical system of exquisite complexity and, as such, vir-
tually any idea and theory about its functioning would ﬁnd some experimental support.
For this reason, the careful integration of data and theories from diﬀerent sources is a
crucial ingredient for well founded hypotheses on brain functioning. We attempted to
reduce the complexity of investigating synaptic plasticity by looking for mismatches and
limit cases. As previously mentioned, this approach forced us to re-evaluate the impact
of many well accepted assumptions on synaptic plasticity and, eventually, to reject some
of those. We have no doubt that also this model, like every other, is wrong at a cer-
tain level of detail. However, we think it represents a signiﬁcant step towards a better
understanding of synaptic plasticity and its implications on neural circuits.




86
                                                                         5.1 Limitations

5.1 Limitations
In this work we considered plasticity of excitatory-to-excitatory connections only. We
decided to focus on these connections because they are the most abundant and best char-
acterized, both theoretically and experimentally. However, plasticity of the other three
“families”, namely excitatory-to-inhibitory, inhibitory-to-excitatory and inhibitory-to-
inhibitory, certainly have a role in higher order brain functions such as learning and
memory.
   Although the excitatory-to-inhibitory case seems to be a minor extension, the other
two types have their own bulk of literature, data and theories and so represent an en-
tirely new subject. The structural plasticity component of the synapse model is still a
work in progress. For the mathematical description, we relied on the work of Fauth et al.
[2015] with minor modiﬁcations. This formulation was particularly convenient because
it guarantees the emergence of bimodal distributions for the number of synapses per
connection, under certain assumptions on the functional plasticity model. Assumptions
all met by the GB model, also used in this work. However, despite the mathematical
proof provided by Fauth et al. [2015], we would need to test this hypothesis in simula-
tion, as our model of functional changes departed quite substantially from the original
formulation by Graupner and Brunel [2012]. Furthermore, to improve the current model
parametrization we would need a larger experimental dataset, focusing on the speciﬁcity
of synaptogenesis and synapse elimination. That is, we do not have quantitative data
linking the emergence (disappearance) of a particular synapse with the activity of its
pre- and postsynaptic neurons.


5.2 Future work
The ﬁrst release of the GluSynapse model was centered around data integration. We
focused on capturing the most fundamental synaptic properties and orchestrating the
diﬀerent mechanisms in a coherent manner. In the second iteration will instead concen-
trate on parameter reﬁnement.
   The N-methyl-D-aspartate receptor (NMDAR) plays a crucial role in the induction
of functional synaptic plasticity. For this reason, accurately modeling its dynamics is of
primary importance. In this work we relied on the formulation proposed by Jahr and
Stevens [1990]. Despite popularity, this model has some important drawbacks. First,
the model is parametrized on hippocampal culture data, without any compensation for
liquid junction potential. Second, the model does not properly capture NMDA spike
generation. We plan to address these issues in the next iteration of our GluSynapse,
based on a set of recent studies on NMDAR dynamics. The model reﬁnement will be
carried on in two phases. Initially we will simply re-parametrize the current model using
neocortical data from in vitro experiments [Vargas-Caballero and Robinson, 2003]. In
the ﬁnal stage of the project, we will consider more detailed descriptions of the NMDAR
[Vargas-Caballero and Robinson, 2004, Moradi et al., 2013].
   The expression mechanisms of functional plasticity in the somatosensory cortex are




                                                                                       87
5 Conclusion

still debated. In particular, it is not completely clear whether STDP stimulation triggers
pre- or post-synaptic changes. We assumed the former case, based on a set of LTP exper-
iments employing an induction protocol similar to STDP [Markram and Tsodyks, 1996].
However, more recent experiments in the visual cortex seem to indicate both a pre-
and post-synaptic locus of expression [Sjöström et al., 2007]. To resolve this mismatch,
we are currently conducing an in depth study of STDP expression mechanisms in the
somatosensory cortex [Perin and Chindemi, in preparation]. The results of our in vitro
experiments will have important implications on the next reﬁnement of GluSynapse. In
particular, preliminary analysis of our data suggests the presence of a postsynaptic com-
ponent, complementing the (expected) presynaptic one. Adding a postsynaptic term to
our synapse model is then a necessary extension. The main eﬀect of this change will be
a reduction of the distance between the depressed and potentiated ﬁxed points for the
                        (d)        (p)
release probability, USE and USE respectively. That is, assuming a purely presynaptic
                                                                        (d)       (p)
expression mechanism forced us to allow a large diﬀerence between USE and USE to cap-
ture the full range of EPSP changes in our dataset. The introduction of a postsynaptic
component will obviously relax this constraint. Furthermore, the new data currently
acquired in our lab will provide a quantitative estimate of the relative contribution of
the pre- and post-synaptic contribution.
                                                                  (d)   (p)
   The ﬁrst release of GluSynapse featured ﬁxed values for USE , USE and NRRP . We
plan to transition from these mean values to more biologically plausible probability dis-
tributions. This reﬁnement is required to better capture several aspects of synaptic
reliability [Rahmon et al., in preparation]. Work on NRRP is already at a mature stage
and we only miss an analysis of the relationships with other synaptic parameters. That
is, we are currently investigating whether NRRP is independent from the other synaptic
variables (i.e. single vesicle release probability or postsynaptic area) or a tight coupling
exists. Furthermore, we need to establish its role in the expression of long-term func-
tional changes. The second release of GluSynapse will be parametrized using a new
generalization algorithm (see Sec. 3.2.5). The existing method does not allow to impose
soft constrains, such us penalizing very large or small parameters, and it is very sensitive
to initial conditions. The new algorithm is being designed to overcome these issues and
to accommodate all the changes described above.
   In this work we did not include any network simulation featuring rewiring. Simulating
rewiring in large networks is an extremely challenging task, mostly because of the large
number of (potential) synapses to consider. Although we already developed a prototype,
the computational resources required will became available only next year (INCITE 2018
grant). As previously discussed, we also want to improve our understanding of circuit
rewiring by producing a larger experimental dataset. Our main objective will be to
establish a clear link between structural and functional plasticity.
   Neuromodulation and Metaplasticity are the next logical steps down the road to a
unifying view of synaptic plasticity. The eﬀects of several neuromodulators on synaptic
activity are already under investigation at the Blue Brain Project (BBP). The synapse
model developed in this work could easily be adapted to support these forms of plasticity.
In particular, ρ is an ideal candidate for the extension. This parameter delimits the




88
                                                                       5.2 Future work

basin of attraction of the potentiated and depressed state. For this reason, it directly
aﬀects the likelihood to potentiate or depress a synapse. By adding the appropriate
dynamics to ρ , it is possible to alter the default behavior of one or more synapses to
implement neuron- and network-level forms of plasticity.




                                                                                     89
Bibliography
Jon I. Arellano, Ruth Benavides-Piccione, Javier DeFelipe, and Rafael Yuste. Ultra-
  structure of dendritic spines: correlation between synaptic and spine morphologies.
  Frontiers in Neuroscience, 1, 2007. ISSN 1662-453X. doi: 10.3389/neuro.01.1.1.010.
  2007. URL http://journal.frontiersin.org/article/10.3389/neuro.01.1.1.
  010.2007/full.

Ashish A. Bagal, Joseph P. Y. Kao, Cha-Min Tang, and Scott M. Thompson. Long-term
  potentiation of exogenous glutamate responses at single dendritic spines. Proceedings
  of the National Academy of Sciences of the United States of America, 102(40):14434–
  14439, October 2005. ISSN 0027-8424, 1091-6490. doi: 10.1073/pnas.0501956102.
  URL http://www.pnas.org/content/102/40/14434.

Natalia Bastrikova, Gregory A. Gardner, Jeﬀ M. Reece, Andreas Jeromin, and Ser-
 ena M. Dudek. Synapse elimination accompanies functional plasticity in hippocam-
 pal neurons. Proceedings of the National Academy of Sciences, 105(8):3123–3127,
 February 2008. ISSN 0027-8424, 1091-6490. doi: 10.1073/pnas.0800027105. URL
 http://www.pnas.org/content/105/8/3123.

Nadine Becker, Corette J. Wierenga, Rosalina Fonseca, Tobias Bonhoeﬀer, and
 U. Valentin Nägerl. LTD Induction Causes Morphological Changes of Presynaptic
 Boutons and Reduces Their Contacts with Spines. Neuron, 60(4):590–597, Novem-
 ber 2008. ISSN 0896-6273. doi: 10.1016/j.neuron.2008.09.018. URL http://www.
 sciencedirect.com/science/article/pii/S0896627308007666.

T. V. P. Bliss and G. L. Collingridge. A synaptic model of memory: long-term potentia-
  tion in the hippocampus. Nature, 361(6407):31–39, January 1993. ISSN 0028-0836. doi:
  10.1038/361031a0. URL https://www.nature.com/nature/journal/v361/n6407/
  abs/361031a0.html.

T. V. P. Bliss and T. Lømo. Long-lasting potentiation of synaptic transmission in
  the dentate area of the anaesthetized rabbit following stimulation of the perforant
  path. The Journal of Physiology, 232(2):331–356, July 1973. ISSN 1469-7793. doi:
  10.1113/jphysiol.1973.sp010273. URL http://onlinelibrary.wiley.com/doi/10.
  1113/jphysiol.1973.sp010273/abstract.

Brenda L Bloodgood and Bernardo L Sabatini. Ca2+ signaling in dendritic spines.
  Current Opinion in Neurobiology, 17(3):345–351, June 2007. ISSN 0959-4388. doi: 10.
  1016/j.conb.2007.04.003. URL http://www.sciencedirect.com/science/article/
  pii/S0959438807000578.



                                                                                    91
Bibliography

P. Brown, A. Hindmarsh, and L. Petzold. Using Krylov Methods in the Solution of
  Large-Scale Diﬀerential-Algebraic Systems. SIAM Journal on Scientiﬁc Computing,
  15(6):1467–1488, November 1994. ISSN 1064-8275. doi: 10.1137/0915088. URL
  http://epubs.siam.org/doi/abs/10.1137/0915088.

Nicholas T Carnevale and Michael L Hines. The NEURON book. Cambridge University
  Press, 2006.

Giuseppe Chindemi, Vincent Delattre, Adrien Devresse, Michael Doron, Michael Graup-
  ner, James King, Pramod Kumbhar, Max Nolte, Rodrigo Perin, Srikanth Ramaswamy,
  Michael Reimann, Christian Rossert, Werner Van Geit, Fabien Delalondre, Eilif
  Muller, Idan Segev, and Henry Markram. Biophysical modeling of calcium based
  synaptic plasticity. in preparation.

L. Niels Cornelisse, Ronald A. J. van Elburg, Rhiannon M. Meredith, Rafael Yuste,
  and Huibert D. Mansvelder. High Speed Two-Photon Imaging of Calcium Dynam-
  ics in Dendritic Spines: Consequences for Spine Calcium Kinetics and Buﬀer Ca-
  pacity. PLOS ONE, 2(10):e1073, October 2007. ISSN 1932-6203. doi: 10.1371/
  journal.pone.0001073. URL http://journals.plos.org/plosone/article?id=10.
  1371/journal.pone.0001073.

Vincenzo De Paola, Anthony Holtmaat, Graham Knott, Sen Song, Linda Wilbrecht,
  Pico Caroni, and Karel Svoboda. Cell Type-Speciﬁc Structural Plasticity of Ax-
  onal Branches and Boutons in the Adult Neocortex. Neuron, 49(6):861–875, March
  2006. ISSN 0896-6273. doi: 10.1016/j.neuron.2006.02.017. URL http://www.
  sciencedirect.com/science/article/pii/S0896627306001346.

Shaul Druckmann, Yoav Banitt, Albert Gidon, Felix Schürmann, Henry Markram, and
  Idan Segev. A Novel Multiple Objective Optimization Framework for Constraining
  Conductance-Based Neuron Models by Experimental Data. Frontiers in Neuroscience,
  1(1):7–18, October 2007. ISSN 1662-4548. doi: 10.3389/neuro.01.1.1.001.2007. URL
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2570085/.

Tarec Fares and Armen Stepanyants. Cooperative synapse formation in the neocortex.
  Proceedings of the National Academy of Sciences, 106(38):16463–16468, September
  2009. ISSN 0027-8424, 1091-6490. doi: 10.1073/pnas.0813265106. URL http://www.
  pnas.org/content/106/38/16463.

Michael Fauth, Florentin Wörgötter, and Christian Tetzlaﬀ. The Formation of Multi-
 synaptic Connections by the Interaction of Synaptic and Structural Plasticity and
 Their Functional Consequences. PLOS Computational Biology, 11(1):e1004031, Jan-
 uary 2015. ISSN 1553-7358. doi: 10.1371/journal.pcbi.1004031. URL http://
 journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004031.

Hagit Vardinon Friedman, Tal Bresler, Craig C Garner, and Noam E Ziv. Assembly of
 New Individual Excitatory Synapses: Time Course and Temporal Order of Synaptic
 Molecule Recruitment. Neuron, 27(1):57–69, July 2000. ISSN 0896-6273. doi: 10.1016/



92
                                                                            Bibliography

  S0896-6273(00)00009-X. URL http://www.sciencedirect.com/science/article/
  pii/S089662730000009X.

Robert C. Froemke, Mu-ming Poo, and Yang Dan. Spike-timing-dependent synaptic
  plasticity depends on dendritic location. Nature, 434(7030):221–225, March 2005.
  ISSN 0028-0836. doi: 10.1038/nature03366. URL http://www.nature.com/nature/
  journal/v434/n7030/abs/nature03366.html.

Galit Fuhrmann, Idan Segev, Henry Markram, and Misha Tsodyks. Coding of Temporal
 Information by Activity-Dependent Synapses. Journal of Neurophysiology, 87(1):140–
 148, January 2002. ISSN 0022-3077, 1522-1598. URL http://jn.physiology.org/
 content/87/1/140.

Michael Graupner and Nicolas Brunel. Calcium-based plasticity model explains sensi-
 tivity of synaptic changes to spike pattern, rate, and dendritic location. Proceedings of
 the National Academy of Sciences, 109(10):3991–3996, March 2012. ISSN 0027-8424,
 1091-6490. doi: 10.1073/pnas.1109359109. URL http://www.pnas.org/content/
 109/10/3991.

K. M. Harris and J. K. Stevens. Dendritic spines of CA 1 pyramidal cells in the rat
  hippocampus: serial electron microscopy with reference to their biophysical charac-
  teristics. Journal of Neuroscience, 9(8):2982–2997, August 1989. ISSN 0270-6474,
  1529-2401. URL http://www.jneurosci.org/content/9/8/2982.

Etay Hay, Sean Hill, Felix Schürmann, Henry Markram, and Idan Segev. Models of
  Neocortical Layer 5b Pyramidal Cells Capturing a Wide Range of Dendritic and Peri-
  somatic Active Properties. PLOS Computational Biology, 7(7):e1002107, July 2011.
  ISSN 1553-7358. doi: 10.1371/journal.pcbi.1002107. URL http://journals.plos.
  org/ploscompbiol/article?id=10.1371/journal.pcbi.1002107.

Bertil Hille et al. Ion channels of excitable membranes, volume 507. Sinauer Sunderland,
  MA, 2001.

Alan C. Hindmarsh, Peter N. Brown, Keith E. Grant, Steven L. Lee, Radu Serban,
  Dan E. Shumaker, and Carol S. Woodward. SUNDIALS: Suite of Nonlinear and
  Diﬀerential/Algebraic Equation Solvers. ACM Trans. Math. Softw., 31(3):363–396,
  September 2005. ISSN 0098-3500. doi: 10.1145/1089014.1089020. URL http://doi.
  acm.org/10.1145/1089014.1089020.

W. R. Holmes and W. B. Levy. Insights into associative long-term potentiation from com-
 putational models of NMDA receptor-mediated calcium inﬂux and intracellular cal-
 cium concentration changes. Journal of Neurophysiology, 63(5):1148–1168, May 1990.
 ISSN 0022-3077. doi: 10.1152/jn.1990.63.5.1148. URL http://www.physiology.org/
 doi/abs/10.1152/jn.1990.63.5.1148.

Anthony Holtmaat and Karel Svoboda. Experience-dependent structural synaptic plas-
 ticity in the mammalian brain. Nature Reviews Neuroscience, 10(9):647–658, Septem-



                                                                                       93
Bibliography

  ber 2009. ISSN 1471-003X. doi: 10.1038/nrn2699. URL http://www.nature.com/
  nrn/journal/v10/n9/abs/nrn2699.html.
Anthony J.G.D. Holtmaat, Joshua T. Trachtenberg, Linda Wilbrecht, Gordon M.
 Shepherd, Xiaoqun Zhang, Graham W. Knott, and Karel Svoboda. Transient
 and Persistent Dendritic Spines in the Neocortex In Vivo. Neuron, 45(2):279–
 291, January 2005. ISSN 0896-6273. doi: 10.1016/j.neuron.2005.01.003. URL
 http://www.sciencedirect.com/science/article/pii/S0896627305000048.
C. E. Jahr and C. F. Stevens. Voltage dependence of NMDA-activated macroscopic con-
  ductances predicted by single-channel kinetics. Journal of Neuroscience, 10(9):3178–
  3182, September 1990. ISSN 0270-6474, 1529-2401. URL http://www.jneurosci.
  org/content/10/9/3178.
C. E. Jahr and C. F. Stevens. Calcium permeability of the N-methyl-D-aspartate receptor
  channel in hippocampal neurons in culture. Proceedings of the National Academy of
  Sciences, 90(24):11573–11577, December 1993. ISSN 0027-8424, 1091-6490. URL
  http://www.pnas.org/content/90/24/11573.
Haruo Kasai, Tatsuya Hayama, Motoko Ishikawa, Satoshi Watanabe, Sho Yagishita, and
 Jun Noguchi. Learning rules and persistence of dendritic spines. The European journal
 of neuroscience, 32(2):241–249, July 2010. ISSN 1460-9568. doi: 10.1111/j.1460-9568.
 2010.07344.x.
Graham W. Knott, Anthony Holtmaat, Linda Wilbrecht, Egbert Welker, and Karel
  Svoboda. Spine growth precedes synapse formation in the adult neocortex in vivo.
  Nature Neuroscience, 9(9):1117–1124, September 2006. ISSN 1097-6256. doi: 10.1038/
  nn1747. URL http://www.nature.com/neuro/journal/v9/n9/full/nn1747.html.
H J Koester and B Sakmann. Calcium dynamics in single spines during coincident
  pre- and postsynaptic activity depend on relative timing of back-propagating action
  potentials and subthreshold excitatory postsynaptic potentials. Proceedings of the
  National Academy of Sciences of the United States of America, 95(16):9596–9601,
  August 1998. ISSN 0027-8424.
D. Kraft. A Software Package for Sequential Quadratic Programming. Deutsche
  Forschungs- und Versuchsanstalt für Luft- und Raumfahrt Köln: Forschungsbericht.
  Wiss. Berichtswesen d. DFVLR, 1988. URL https://books.google.ch/books?id=
  4rKaGwAACAAJ.
Dimitri Kullmann. Silent synapse. Scholarpedia, 6(8):10705, 2011. ISSN 1941-6016.
  doi: 10.4249/scholarpedia.10705. URL http://www.scholarpedia.org/article/
  Silent_Synapse.
Jean-Vincent Le Bé and Henry Markram. Spontaneous and evoked synaptic rewiring
  in the neonatal neocortex. Proceedings of the National Academy of Sciences, 103
  (35):13214–13219, August 2006. ISSN 0027-8424, 1091-6490. doi: 10.1073/pnas.
  0604691103. URL http://www.pnas.org/content/103/35/13214.



94
                                                                          Bibliography

Yann LeCun, Yoshua Bengio, and Geoﬀrey Hinton. Deep learning. Nature, 521(7553):
  436–444, May 2015. ISSN 0028-0836. doi: 10.1038/nature14539. URL http://www.
  nature.com/nature/journal/v521/n7553/abs/nature14539.html.

J. Lisman. A mechanism for the Hebb and the anti-Hebb processes underlying learning
   and memory. Proceedings of the National Academy of Sciences, 86(23):9574–9578,
   December 1989. ISSN 0027-8424, 1091-6490. URL http://www.pnas.org/content/
   86/23/9574.

Alex Loebel, Gilad Silberberg, Daniela Helbig, Henry Markram, Misha Tsodyks, and
  Magnus J. E. Richardson. Multiquantal Release Underlies the Distribution of Synaptic
  Eﬃcacies in the Neocortex. Frontiers in Computational Neuroscience, 3, November
  2009. ISSN 1662-5188. doi: 10.3389/neuro.10.027.2009. URL http://www.ncbi.nlm.
  nih.gov/pmc/articles/PMC2786302/.

J C Magee and D Johnston. Characterization of single voltage-gated Na+ and Ca2+
  channels in apical dendrites of rat CA1 pyramidal neurons. The Journal of Physiology,
  487(1):67–90, August 1995. ISSN 1469-7793. doi: 10.1113/jphysiol.1995.sp020862.
  URL http://onlinelibrary.wiley.com/doi/10.1113/jphysiol.1995.sp020862/
  abstract.

Robert C. Malenka and Mark F. Bear. LTP and LTD: An Embarrassment of
  Riches. Neuron, 44(1):5–21, September 2004. ISSN 0896-6273. doi: 10.1016/j.
  neuron.2004.09.012. URL http://www.sciencedirect.com/science/article/pii/
  S0896627304006087.

Tiina Manninen, Katri Hituri, Jeanette Hellgren Kotaleski, Kim T. Blackwell, and
  Marja-Leena Linne. Postsynaptic Signal Transduction Models for Long-Term Po-
  tentiation and Depression. Frontiers in Computational Neuroscience, 4, December
  2010. ISSN 1662-5188. doi: 10.3389/fncom.2010.00152. URL http://www.ncbi.nlm.
  nih.gov/pmc/articles/PMC3006457/.

H. Markram, J. Lübke, M. Frotscher, A. Roth, and B. Sakmann. Physiology and anatomy
  of synaptic connections between thick tufted pyramidal neurones in the developing rat
  neocortex. The Journal of Physiology, 500(Pt 2):409–440, April 1997a. ISSN 0022-
  3751, 1469-7793. URL http://jp.physoc.org/content/500/Pt_2/409.

Henry Markram and Misha Tsodyks. Redistribution of synaptic eﬃcacy between neocor-
  tical pyramidal neurons. , Published online: 29 August 1996; | doi:10.1038/382807a0,
  382(6594):807–810, August 1996. ISSN ${footerJournalISSN}. doi: 10.1038/382807a0.
  URL http://www.nature.com/nature/journal/v382/n6594/abs/382807a0.html.

Henry Markram, Joachim Lübke, Michael Frotscher, and Bert Sakmann. Regulation
  of Synaptic Eﬃcacy by Coincidence of Postsynaptic APs and EPSPs. Science, 275
  (5297):213–215, January 1997b. ISSN 0036-8075, 1095-9203. doi: 10.1126/science.
  275.5297.213. URL http://www.sciencemag.org/content/275/5297/213.




                                                                                    95
Bibliography

Henry Markram, Arnd Roth, and Fritjof Helmchen. Competitive Calcium Binding:
  Implications for Dendritic Calcium Signaling. Journal of Computational Neuroscience,
  5(3):331–348, July 1998. ISSN 0929-5313, 1573-6873. doi: 10.1023/A:1008891229546.
  URL http://link.springer.com/article/10.1023/A%3A1008891229546.
Henry Markram, Eilif Muller, Srikanth Ramaswamy, Michael W. Reimann, Marwan
  Abdellah, Carlos Aguado Sanchez, Anastasia Ailamaki, Lidia Alonso-Nanclares, Nico-
  las Antille, Selim Arsever, Guy Antoine Atenekeng Kahou, Thomas K. Berger, Ah-
  met Bilgili, Nenad Buncic, Athanassia Chalimourda, Giuseppe Chindemi, Jean-Denis
  Courcol, Fabien Delalondre, Vincent Delattre, Shaul Druckmann, Raphael Dumusc,
  James Dynes, Stefan Eilemann, Eyal Gal, Michael Emiel Gevaert, Jean-Pierre Gho-
  bril, Albert Gidon, Joe W. Graham, Anirudh Gupta, Valentin Haenel, Etay Hay,
  Thomas Heinis, Juan B. Hernando, Michael Hines, Lida Kanari, Daniel Keller, John
  Kenyon, Georges Khazen, Yihwa Kim, James G. King, Zoltan Kisvarday, Pramod
  Kumbhar, Sébastien Lasserre, Jean-Vincent Le Bé, Bruno R. C. Magalhães, An-
  gel Merchán-Pérez, Julie Meystre, Benjamin Roy Morrice, Jeﬀrey Muller, Alberto
  Muñoz-Céspedes, Shruti Muralidhar, Keerthan Muthurasa, Daniel Nachbaur, Tay-
  lor H. Newton, Max Nolte, Aleksandr Ovcharenko, Juan Palacios, Luis Pastor, Rodrigo
  Perin, Rajnish Ranjan, Imad Riachi, José-Rodrigo Rodrı́guez, Juan Luis Riquelme,
  Christian Rössert, Konstantinos Sfyrakis, Ying Shi, Julian C. Shillcock, Gilad Sil-
  berberg, Ricardo Silva, Farhan Tauheed, Martin Telefont, Maria Toledo-Rodriguez,
  Thomas Tränkler, Werner Van Geit, Jafet Villafranca Dı́az, Richard Walker, Yun
  Wang, Stefano M. Zaninetta, Javier DeFelipe, Sean L. Hill, Idan Segev, and Felix
  Schürmann. Reconstruction and Simulation of Neocortical Microcircuitry. Cell, 163
  (2):456–492, October 2015. ISSN 0092-8674. doi: 10.1016/j.cell.2015.09.029. URL
  http://www.sciencedirect.com/science/article/pii/S0092867415011915.
Gábor Molnár, Márton Rózsa, Judith Baka, Noémi Holderith, Pál Barzó, Zoltan Nusser,
  and Gábor Tamás. Human pyramidal to interneuron synapses are mediated by multi-
  vesicular release and multiple docked vesicles. eLife, 5:e18167, August 2016. ISSN
  2050-084X. doi: 10.7554/eLife.18167. URL https://elifesciences.org/articles/
  18167.
Keivan Moradi, Kamran Moradi, Mahin Ganjkhani, Mojtaba Hajihasani, Shahriar
  Gharibzadeh, and Gholamreza Kaka. A fast model of voltage-dependent NMDA recep-
  tors. Journal of Computational Neuroscience, 34(3):521–531, June 2013. ISSN 0929-
  5313, 1573-6873. doi: 10.1007/s10827-012-0434-4. URL https://link.springer.
  com/article/10.1007/s10827-012-0434-4.
Eilif Muller, Michael Reimann, and Srikanth Ramaswamy. Probabilistic AMPA and
  NMDA receptor with presynaptic short-term plasticity. Unpublished, 2011.
Thomas Nevian and Bert Sakmann. Spine Ca2+ Signaling in Spike-Timing-Dependent
 Plasticity. The Journal of Neuroscience, 26(43):11001–11013, October 2006. ISSN
  0270-6474, 1529-2401. doi: 10.1523/JNEUROSCI.1749-06.2006. URL http://www.
  jneurosci.org/content/26/43/11001.



96
                                                                           Bibliography

U.Valentin Nägerl, Nicola Eberhorn, Sidney B. Cambridge, and Tobias Bonhoef-
  fer. Bidirectional Activity-Dependent Morphological Plasticity in Hippocampal Neu-
  rons. Neuron, 44(5):759–767, December 2004. ISSN 0896-6273. doi: 10.1016/j.
  neuron.2004.11.016. URL http://www.sciencedirect.com/science/article/pii/
  S0896627304007299.

Daniel H. O’Connor, Gayle M. Wittenberg, and Samuel S.-H. Wang. Graded bidirec-
 tional synaptic plasticity is composed of switch-like unitary events. Proceedings of
 the National Academy of Sciences of the United States of America, 102(27):9679–
 9684, July 2005. ISSN 0027-8424, 1091-6490. doi: 10.1073/pnas.0502332102. URL
 http://www.pnas.org/content/102/27/9679.

Rodrigo Perin and Giuseppe Chindemi. Extracellular calcium concentration eﬀects on
 stdp. in preparation.

Rodrigo Perin, Thomas K. Berger, and Henry Markram. A synaptic organizing principle
 for cortical neuronal groups. Proceedings of the National Academy of Sciences, 108(13):
 5419–5424, March 2011. ISSN 0027-8424, 1091-6490. doi: 10.1073/pnas.1016051108.
 URL http://www.pnas.org/content/108/13/5419.

Carl C. H. Petersen, Robert C. Malenka, Roger A. Nicoll, and John J. Hopﬁeld. All-
  or-none potentiation at CA3-CA1 synapses. Proceedings of the National Academy
  of Sciences, 95(8):4732–4737, April 1998. ISSN 0027-8424, 1091-6490. URL http:
  //www.pnas.org/content/95/8/4732.

Dale Purves, George J Augustine, David Fitzpatrick, WC Hall, AS LaMantia, JO Mc-
 Namara, and LE White. Neuroscience, 2008. De Boeck, Sinauer, Sunderland, Mass,
 2014.

John Rahmon, Max Chindemi, Giuseppe Nolte, Rodrigo Perin, Srikanth Ramaswamy,
  Michael Reimann, Eilif Muller, and Henry Markram. Multi-vesicular release in the
  neocortex. in preparation.

Srikanth Ramaswamy, Sean L. Hill, James G. King, Felix Schürmann, Yun Wang, and
  Henry Markram. Intrinsic morphological diversity of thick-tufted layer 5 pyramidal
  neurons ensures robust and invariant properties of in silico synaptic connections. The
  Journal of Physiology, 590(4):737–752, February 2012. ISSN 1469-7793. doi: 10.
  1113/jphysiol.2011.219576. URL http://onlinelibrary.wiley.com/doi/10.1113/
  jphysiol.2011.219576/abstract.

Srikanth Ramaswamy, Jean-Denis Courcol, Marwan Abdellah, Stanislaw R. Adaszewski,
  Nicolas Antille, Selim Arsever, Guy Atenekeng, Ahmet Bilgili, Yury Brukau, Athanas-
  sia Chalimourda, Giuseppe Chindemi, Fabien Delalondre, Raphael Dumusc, Stefan
  Eilemann, Michael Emiel Gevaert, Padraig Gleeson, Joe W. Graham, Juan B. Her-
  nando, Lida Kanari, Yury Katkov, Daniel Keller, James G. King, Rajnish Ran-
  jan, Michael W. Reimann, Christian Rössert, Ying Shi, Julian C. Shillcock, Mar-
  tin Telefont, Werner Van Geit, Jafet Villafranca Diaz, Richard Walker, Yun Wang,



                                                                                     97
Bibliography

  Stefano M. Zaninetta, Javier DeFelipe, Sean L. Hill, Jeﬀrey Muller, Idan Segev,
  Felix Schürmann, Eilif B. Muller, and Henry Markram. The neocortical microcir-
  cuit collaboration portal: a resource for rat somatosensory cortex. Frontiers in
  Neural Circuits, 9, 2015. ISSN 1662-5110. doi: 10.3389/fncir.2015.00044. URL
  http://journal.frontiersin.org/article/10.3389/fncir.2015.00044/full.

Michael W. Reimann, James G. King, Eilif B. Muller, Srikanth Ramaswamy, and Henry
 Markram. An algorithm to predict the connectome of neural microcircuits. Frontiers in
 Computational Neuroscience, 9, October 2015. ISSN 1662-5188. doi: 10.3389/fncom.
 2015.00120. URL http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4597796/.

David A. Richards, José Maria Mateos, Sylvain Hugel, Vincenzo de Paola, Pico Caroni,
  Beat H. Gähwiler, and R. Anne McKinney. Glutamate induces the rapid formation
  of spine head protrusions in hippocampal slice cultures. Proceedings of the National
  Academy of Sciences, 102(17):6166–6171, April 2005. ISSN 0027-8424, 1091-6490. doi:
  10.1073/pnas.0501881102. URL http://www.pnas.org/content/102/17/6166.

Bernardo L. Sabatini and Karel Svoboda. Analysis of calcium channels in single spines
  using optical ﬂuctuation analysis. Nature, 408(6812):589–593, November 2000. ISSN
  0028-0836. doi: 10.1038/35046076. URL http://www.nature.com/nature/journal/
  v408/n6812/abs/408589a0.html.

Bernardo L. Sabatini, Thomas G. Oertner, and Karel Svoboda. The Life Cycle of
  Ca2+ Ions in Dendritic Spines. Neuron, 33(3):439–452, January 2002. ISSN 0896-
  6273. doi: 10.1016/S0896-6273(02)00573-1. URL http://www.sciencedirect.com/
  science/article/pii/S0896627302005731.

T. Schikorski and C. F. Stevens. Quantitative ﬁne-structural analysis of olfactory cortical
  synapses. Proceedings of the National Academy of Sciences, 96(7):4107–4112, March
  1999. ISSN 0027-8424, 1091-6490. doi: 10.1073/pnas.96.7.4107. URL http://www.
  pnas.org/content/96/7/4107.

Jackie Schiller, Guy Major, Helmut J. Koester, and Yitzhak Schiller. NMDA spikes in
  basal dendrites of cortical pyramidal neurons. Nature, 404(6775):285–289, March 2000.
  ISSN 0028-0836. doi: 10.1038/35005094. URL http://www.nature.com/nature/
  journal/v404/n6775/full/404285a0.html.

R. Schneggenburger, Z. Zhou, A. Konnerth, and E. Neher. Fractional contribution of
  calcium to the cation current through glutamate receptor channels. Neuron, 11(1):
  133–143, July 1993. ISSN 0896-6273. doi: 10.1016/0896-6273(93)90277-X. URL
  http://www.sciencedirect.com/science/article/pii/089662739390277X.

Menahem Segal. Dendritic spines and long-term plasticity. Nature Reviews Neuroscience,
 6(4):277–284, April 2005. ISSN 1471-003X. doi: 10.1038/nrn1649. URL http://www.
 nature.com/nrn/journal/v6/n4/abs/nrn1649.html.




98
                                                                           Bibliography

Per Jesper Sjöström and Michael Häusser. A Cooperative Switch Determines the Sign
  of Synaptic Plasticity in Distal Dendrites of Neocortical Pyramidal Neurons. Neuron,
  51(2):227–238, July 2006. ISSN 0896-6273. doi: 10.1016/j.neuron.2006.06.017. URL
  http://www.sciencedirect.com/science/article/pii/S0896627306004715.

Per Jesper Sjöström, Gina G Turrigiano, and Sacha B Nelson. Rate, Timing, and
  Cooperativity Jointly Determine Cortical Synaptic Plasticity. Neuron, 32(6):1149–
  1164, December 2001. ISSN 0896-6273. doi: 10.1016/S0896-6273(01)00542-6. URL
  http://www.sciencedirect.com/science/article/pii/S0896627301005426.

Per Jesper Sjöström, Gina G. Turrigiano, and Sacha B. Nelson. Multiple forms of long-
  term plasticity at unitary neocortical layer 5 synapses. Neuropharmacology, 52(1):176–
  184, January 2007. ISSN 0028-3908. doi: 10.1016/j.neuropharm.2006.07.021. URL
  http://www.sciencedirect.com/science/article/pii/S0028390806002310.

Armen Stepanyants, Patrick R. Hof, and Dmitri B. Chklovskii. Geometry and Structural
  Plasticity of Synaptic Connectivity. Neuron, 34(2):275–288, April 2002. ISSN 0896-
  6273. doi: 10.1016/S0896-6273(02)00652-9. URL http://www.sciencedirect.com/
  science/article/pii/S0896627302006529.

N. Toni, P.-A. Buchs, I. Nikonenko, C. R. Bron, and D. Muller. LTP promotes formation
  of multiple spine synapses between a single axon terminal and a dendrite. Nature,
  402(6760):421–425, November 1999. ISSN 0028-0836. doi: 10.1038/46574. URL
  http://www.nature.com/nature/journal/v402/n6760/abs/402421a0.html.

Nicolas Toni, E. Matthew Teng, Eric A. Bushong, James B. Aimone, Chunmei Zhao,
  Antonella Consiglio, Henriette van Praag, Maryann E. Martone, Mark H. Ellisman,
  and Fred H. Gage. Synapse formation on neurons born in the adult hippocampus. Na-
  ture Neuroscience, 10(6):727–734, June 2007. ISSN 1097-6256. doi: 10.1038/nn1908.
  URL https://www.nature.com/neuro/journal/v10/n6/full/nn1908.html.

Misha V. Tsodyks and Henry Markram. The neural code between neocortical pyramidal
 neurons depends on neurotransmitter release probability. Proceedings of the National
 Academy of Sciences, 94(2):719–723, January 1997. ISSN 0027-8424, 1091-6490. URL
 http://www.pnas.org/content/94/2/719.

Werner Van Geit, Michael Gevaert, Giuseppe Chindemi, Christian Rössert, Jean-Denis
 Courcol, Eilif B. Muller, Felix Schürmann, Idan Segev, and Henry Markram. BluePy-
 Opt: Leveraging Open Source Software and Cloud Infrastructure to Optimise Model
 Parameters in Neuroscience. Frontiers in Neuroinformatics, 10, 2016. ISSN 1662-5196.
 doi: 10.3389/fninf.2016.00017. URL http://journal.frontiersin.org/article/
 10.3389/fninf.2016.00017/abstract.

Mariana Vargas-Caballero and Hugh P. C. Robinson. A Slow Fraction of Mg2+ Un-
 block of NMDA Receptors Limits Their Contribution to Spike Generation in Cortical
 Pyramidal Neurons. Journal of Neurophysiology, 89(5):2778–2783, May 2003. ISSN




                                                                                     99
Bibliography

  0022-3077, 1522-1598. doi: 10.1152/jn.01038.2002. URL http://jn.physiology.
  org/content/89/5/2778.

Mariana Vargas-Caballero and Hugh P. C. Robinson. Fast and Slow Voltage-Dependent
 Dynamics of Magnesium Block in the NMDA Receptor: The Asymmetric Trapping
 Block Model. Journal of Neuroscience, 24(27):6171–6180, July 2004. ISSN 0270-6474,
 1529-2401. doi: 10.1523/JNEUROSCI.1380-04.2004. URL http://www.jneurosci.
 org/content/24/27/6171.

Noam E. Ziv and Craig C. Garner. Cellular and molecular mechanisms of presynap-
 tic assembly. Nature Reviews Neuroscience, 5(5):385–399, May 2004. ISSN 1471-
 003X. doi: 10.1038/nrn1370. URL http://www.nature.com/nrn/journal/v5/n5/
 abs/nrn1370.html.

Noam E. Ziv and Stephen J Smith. Evidence for a Role of Dendritic Filopodia in
 Synaptogenesis and Spine Formation. Neuron, 17(1):91–102, July 1996. ISSN 0896-
 6273. doi: 10.1016/S0896-6273(00)80283-4. URL http://www.sciencedirect.com/
 science/article/pii/S0896627300802834.

D. Zwillinger and S. Kokoska. CRC Standard Probability and Statistics Tables and
  Formulae. Chapman and Hall: New York, 2000.




100
Acronyms
AI Artiﬁcial Intelligence 15

AMPAR alpha-amino-3-hydroxy-5-methyl-4-isoxazole propionate receptor 20, 21, 23,
25, 26, 30–33, 38

ANN Artiﬁcial Neural Networks 15

AP action potential 19, 20, 36, 67

AZ Active Zone 20, 22

bAP backpropagating action potential 21

BBP Blue Brain Project 5, 7, 30, 31, 33, 67, 86

EPSP Excitatory Post-Synaptic Potential 23, 58, 61, 79, 86

GA Genetic Algorithm 17, 53, 55, 83

GB Graupner-Brunel 38–40, 43, 52, 85

GHK Goldman-Hodgkin-Katz 35, 37

HH Hodgkin-Huxley 34

HPC High Performance Computing 18

HVA High Voltage Activated 21, 31

LTD Long Term Depression 5, 7, 16, 17, 26, 31, 36, 40, 52, 53, 55–57, 60, 63, 83

LTP Long Term Potentiation 5, 7, 16, 17, 25, 26, 31, 33, 36, 37, 39, 40, 53, 55–57, 60,
63, 65, 83, 85

LVA Low Voltage Activated 21

MOO Multi-Objective Optimization 52

MSB Multi Synapse Bouton 21, 24

MVR Multi Vesicular Release 33, 34



                                                                                   101
Acronyms

 NMC Neocortical Microcircuit Collaboration 34

 NMDAR N-methyl-D-aspartate receptor 17, 20, 21, 23, 30–37, 39, 68, 85

 NMODL NEURON model description language 66, 67

 ODE Ordinary Diﬀerential Equation 32, 35, 67

 PC pyramidal cell 60–63

 PSD Postsynaptic Density 20, 22, 23, 56

 PSP Post-Synaptic Potential 53

 RRP Readily-Releasable Pool 40

 SLSQP Sequential Least SQuares Programming 59

 SNARE SNAP (Soluble NSF Attachment Protein) REceptor 20

 STDP Spike Timing Dependent Plasticity 5, 7, 36, 38, 40, 42, 47, 49–51, 53, 56, 60,
 63, 75, 83, 85

 TM Tsodyks-Markram 31–33

 TTPC Thick-tufted pyramidal cells 32, 40–42, 45, 47, 49–52, 58–63, 74, 75, 77, 78,
 80, 81, 83

 VDCC Voltage Dependent Calcium Channels 3, 20, 21, 30, 31, 34, 35, 39, 66




102
Glossary
 me-type is a morphological and electrical type 99

 pathway is the set of all connections between two given me-types 17, 18, 49, 59, 60,
 63, 65, 78, 83




                                                                                 103
                                 Giuseppe Chindemi
                        École Polytechnique Fédérale de Lausanne
             ENT-R BBP-CORE, Campus Biotech, Ch. des Mines 9, CH-1202 Genève
                                  Tel +41 21 69 37258
                               giuseppe.chindemi@epfl.ch


AREAS OF INTERESTS
Artificial Intelligence, Computer Science, Neuroscience, Synaptic Plasticity


EDUCATION
 2012 – Present         École Polytechnique Fédérale de Lausanne (EPFL), Switzerland
                        Blue Brain Project (BBP)
                        Ph.D. in Neuroscience
                        Dissertation: Towards a Unified Understanding of Synaptic Plasticity: Parsimonious
                        Modeling and Simulation of the Glutamatergic Synapse Life-Cycle (Advisor: Henry
                        Markram and Eilif B. Müller)

 2009 – 2012            Politecnico di Milano (Polimi), Italy
                        University of Illinois at Chicago (UIC), US
                        Double Degree MSc Program in Engineering of Computing Systems and
                        Computer Science
                        Master thesis: A Neural Network Study of Predicting Seizures in Epilepsy
                        (Advisor: Daniel Graupe and Simone Garatti)

 2006 – 2009            Politecnico di Milano (Polimi), Italy
                        BS degree in Engineering of Computing Systems
                        Bachelor thesis: Force Evaluation in Cloth Simulation Using Massive Parallelism
                        Hardware (Advisor: Giuseppina Gini)



TEACHING EXPERIENCE
Teaching assistant for the following graduate courses at EPFL:
   x “In Silico Neuroscience” (Prof. Felix Shürmann and Prof. Sean Hill) – Spring 2015
   x “Neurosciences II” (Prof. Carl Petersen, Prof. Ralf Schneggenburger and Prof. Felix
       Shürmann) – Fall 2013
   x “Pattern Classification and Machine Learning” (Prof. Matthias Seeger) – Spring 2013


PUBLICATIONS
  x Markram, Henry, Eilif Muller, Srikanth Ramaswamy, Michael W. Reimann, Marwan
     Abdellah, Carlos Aguado Sanchez, Anastasia Ailamaki et al. “Reconstruction and simulation
     of neocortical microcircuitry.” Cell 163, no. 2 (2015): 456-492.
  x Ramaswamy, Srikanth, Jean-Denis Courcol, Marwan Abdellah, Stanislaw R. Adaszewski,
     Nicolas Antille, Selim Arsever, Guy Atenekeng et al. “The neocortical microcircuit
     collaboration portal: a resource for rat somatosensory cortex.” Frontiers in neural circuits 9
     (2015).



                                                    1
    x   Van Geit, Werner, Michael Gevaert, Giuseppe Chindemi, Christian Rössert, Jean-Denis
        Courcol, Eilif B. Muller, Felix Schürmann, Idan Segev, and Henry Markram. “BluePyOpt:
        leveraging open source software and cloud infrastructure to optimise model parameters in
        neuroscience.” Frontiers in neuroinformatics 10 (2016).
    x   Rössert, Christian, Christian Pozzorini, Giuseppe Chindemi, Andrew P. Davison, Csaba
        Eroe, James King, Taylor H. Newton et al. “Automated point-neuron simplification of data-
        driven microcircuit models.” arXiv preprint arXiv:1604.00087 (2016).
    x   Reimann, Michael W., Max Nolte, Martina Scolamiero, Katharine Turner, Rodrigo Perin,
        *LXVHSSH&KLQGHPL3DZHâ'âRWNR5DQ/HYL.Dthryn Hess, and Henry Markram. “Cliques
        of Neurons Bound into Cavities Provide a Missing Link between Structure and Function.”
        Frontiers in Computational Neuroscience 11 (2017): 48.
    x   Doron, Michael, Giuseppe Chindemi, Eilif Muller, Henry Markram, and Idan Segev. “Timed
        Synaptic Inhibition Shapes NMDA Spikes, Influencing Local Dendritic Processing and
        Global I/O Properties of Cortical Neurons.” Cell reports 21, no. 6 (2017): 1550-1561.


PERSONAL INFORMATION
Date of birth: March 23th, 1988
Citizenship: Italian
Languages: Italian (native), English (fluent), French (intermediate)
Skills: Artificial Intelligence, Complex Systems Analysis, Genetic Algorithms, High Performance
Computing, Machine Learning, Modeling, Neuroscience, Synaptic Plasticity, Software Engineering




                                                2
