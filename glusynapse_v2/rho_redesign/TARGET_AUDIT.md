# Target audit: stage-3 cores and worst validation misses against the original papers (2026-10-02)

**Question.** Do the mean, SEM and n of the 20 stage-3 cores (STAGE3_DESIGN.md section 2) and the worst S1C validation
misses match the papers? Are any readouts or preparations off?
**Short answer.** One core is WRONG: `sjostrom_40hz_dt+10ms` takes 1.53 +- 0.14 (n=11) from the Fig 6B legend of
Sjostrom 2001, but that n=11 group pools +10 **and -10** ms at 40 Hz. The pure +10 value is Fig 1D: 1.54 +- 0.11 (n=6).
The other means and SEMs agree with the paper text, or with a figure read by eye to about +-0.03.
`sjostrom_40hz_dt0ms` SEM 0.038 is an SEM, not an SD. Every core comes from unitary paired recordings.

Sources were read page by page from the PDFs: Zilberter 2009, Letzkus 2006, Sjostrom 2003, Sjostrom 2007 and Sjostrom &
Hausser 2006 from plastyfire/papers/. Sjostrom 2001 came from the author's lab PDF, Markram 1997 from a course-site PDF.
PubMed gave only abstracts for Letzkus and Zilberter, so the PDFs were used. DOIs, via PubMed: Markram 1997
[10.1126/science.275.5297.213](https://doi.org/10.1126/science.275.5297.213); Sjostrom 2001
[10.1016/s0896-6273(01)00542-6](https://doi.org/10.1016/s0896-6273(01)00542-6); Sjostrom 2003
[10.1016/s0896-6273(03)00476-8](https://doi.org/10.1016/s0896-6273(03)00476-8); Sjostrom 2007
[10.1016/j.neuropharm.2006.07.021](https://doi.org/10.1016/j.neuropharm.2006.07.021); Sjostrom & Hausser 2006
[10.1016/j.neuron.2006.06.017](https://doi.org/10.1016/j.neuron.2006.06.017); Letzkus 2006
[10.1523/JNEUROSCI.2650-06.2006](https://doi.org/10.1523/JNEUROSCI.2650-06.2006); Zilberter 2009
[10.1093/cercor/bhn247](https://doi.org/10.1093/cercor/bhn247).
"fig" = read by eye from the figure (about +-0.03); "text" = printed in the paper.

## 1. The 20 cores

| # | target | current mean+-sem (n) | paper value | source location | status | fix |
|---|---|---|---|---|---|---|
| 1 | 10Hz_10ms | 1.201+-0.063 (6) | fig ~1.20-1.22 late points; n "six bidirectionally coupled neurons" | Markram 1997 Fig 3C open squares, p.214 | OK (fig) | none; n could be 3 if "neurons" is literal (UNVERIFIABLE n) |
| 2 | 10Hz_-10ms | 0.792+-0.026 (6) | fig ~0.78-0.80 at 30-40 min, bars ~+-0.02-0.03 | Markram 1997 Fig 3C closed squares | OK (fig) | as #1 |
| 3 | sj 0.1hz +10 | 0.97+-0.04 (7) | fig ~0.97, n=7; Chindemi csv 0.964+-0.051 | Sjostrom 2001 Fig 1D | OK (fig); SEM 0.04-0.05 | none (the 0.05 floor already covers it) |
| 4 | sj 0.1hz -10 | 0.69+-0.07 (7) | fig ~0.69, n=7 (csv 0.682+-0.077) | Sjostrom 2001 Fig 7B | OK (fig) | none |
| 5 | sj 20hz +10 | 1.31+-0.14 (3) | fig ~1.31, n=3 (csv 1.305+-0.139) | Sjostrom 2001 Fig 1D | OK (fig) | none |
| 6 | sj 40hz +10 | 1.53+-0.14 (11) | Fig 1D +10 only: n=6, fig ~1.54 (csv sjs01_04 1.540+-0.113). The 153+-14% n=11 is the "+-10 ms, 40 Hz" group | Sjostrom 2001 Fig 6B legend + symbol key; Fig 1D | **WRONG** (n, SEM) | 1.54+-0.113 (n=6), Fig 1D |
| 7 | sj 50hz -10 | 1.70+-0.19 (10) | fig ~1.70, n=10 (csv 1.692+-0.186) | Sjostrom 2001 Fig 7B | OK (fig) | none |
| 8 | sj 20hz -10 | 0.65+-0.09 (4) | fig ~0.65, n=4 (csv 0.639+-0.090) | Sjostrom 2001 Fig 7B | OK (fig) | none |
| 9 | sj 20hz -10 mglu_block | 1.02+-0.07 (7) | text "20 Hz (102 +- 7%, n = 7)", AM251 | Sjostrom 2003 p.647 col 1, Fig 8B | OK (text) | timing note below |
| 10 | sj07 step200 pair | 1.62+-0.07 (31) | text "162 +- 7%", Fig 1B "LTP, n=31" | Sjostrom 2007 p.179 sect. 3, Fig 1B | OK (text) | none |
| 11 | letzkus 1ap +10 | 0.72+-0.03 (15) | text "0.72 +- 0.03 (n = 15)", **pooled over locations** (rise 0.8-4 ms, Fig 2D) | Letzkus 2006 p.10423, Fig 2C,D | OK (text) | value OK; STAGE3 calls it "distal": it is pooled |
| 12 | letzkus 3ap +10 @prox | 1.30+-0.10 (10) | text "1.3 +- 0.1 (n = 10)" | Letzkus 2006 p.10425, Fig 5 | OK (text) | none |
| 13 | letzkus 3ap +10 @dist | 0.79+-0.06 (10) | text "0.79 +- 0.06 (n = 10)" | Letzkus 2006 p.10425, Fig 5 | OK (text) | none |
| 14 | letzkus 3ap -10 @prox | 0.89+-0.03 (16) | text "0.89 +- 0.03 (n = 16)" | Letzkus 2006 p.10425, Fig 5 | OK (text) | none |
| 15 | zilberter 1ap +10 | 0.64+-0.07 (6) | text "0.64 +- 0.07 of control, n = 6"; Fig 2 legend says n = 10 for 2B | Zilberter 2009 p.2311, Fig 2B, 3 | OK (text); n text vs legend | keep text n=6 |
| 16 | zilberter 1ap -10 | 0.56+-0.06 (4) | text "0.56 +- 0.06 of control, n = 4; data not shown" | Zilberter 2009 p.2311 | OK (text) | none |
| 17 | zilb train10 +4 _last | 1.49+-0.12 (11) | text "1.49 +- 0.12; n = 11", pre 3-5 ms before the 10th AP | Zilberter 2009 p.2312, Fig 4A,D | OK (text) | none |
| 18 | zilb train10 -10 _last | 0.72+-0.05 (13) | text "0.72 +- 0.05 of control, n = 13", pre 5-12 ms after the 10th AP | Zilberter 2009 p.2312, Fig 4C,D | OK (text) | none |
| 19 | zilb train10 -10 _last mglu_block | 0.73+-0.07 (7) | text "0.73 +- 0.07 of control, n = 7", AM251 2 uM | Zilberter 2009 p.2315, Fig 6C | OK (text) | none |
| 20 | zilberter 5ap 10hz +10 | 0.76+-0.07 (19) | text "0.76 +- 0.07 of control; n = 19" | Zilberter 2009 p.2311, Fig 2G (legend n=19) | OK (text) | none |
| 3N | sj07 pair no_block | 1.362+-0.074 (14) | L-NAME 130+-12% n=8, cPTIO 140+-9.3% n=6; pooled text "136%" | Sjostrom 2007 p.181 sect. 3.3, Fig 5A | OK (text, inverse-variance) | none |

## 2. Worst validation misses and the Sjostrom 2001 dt 0 / +-25 rows

| target | current | paper value | source location | status | fix |
|---|---|---|---|---|---|
| sj 40hz dt0 | 0.929+-0.038 (5) | author xlsx: freq 40, mean 92.89, column header **"SEM"** 3.84; legend "n's are between 3 and 5"; paper: "Means are reported as +- SEM" | Sjostrom 2001 Fig 7D; Exp. Proc. p.1162; lab DataPage Synch_firing.xlsx | OK (SEM, not SD); n=5 is the upper bound | none. With n 3-5 the SD is 0.067-0.086, which is normal. If it were an SD the SEM would be 0.017-0.022 |
| sj 50hz dt0 (dropped) | 0.923+-0.035 (5) | the same measurement as 40 Hz dt 0. The Fig 7C "40/50 Hz" dt 0 point is the 40 Hz synchronous set (text p.1156: "no LTP at 40 Hz (Figures 7C and 7D)") | Sjostrom 2001 Fig 7C,D | duplicate, drop confirmed | keep dropped |
| sj 20hz dt0 | 0.761+-0.099 (5) | xlsx 76.08 +- 9.94 SEM; n 3-5 | Fig 7D, xlsx | OK; n UNVERIFIABLE (3-5) | none |
| sj 20hz +25 | 0.652+-0.184 (5) | fig ~0.63-0.65, large bar; n 3-7 | Fig 7C 20 Hz panel | OK (fig); n UNVERIFIABLE | none |
| sj 20hz -25 | 0.705+-0.070 (5) | fig ~0.70; n 3-7 | Fig 7C 20 Hz panel | OK (fig); n UNVERIFIABLE | none |
| letzkus 3ap -10 @dist (DROPT) | 1.42+-0.09 (7) | text "1.42 +- 0.09 ... n = 7" | Letzkus p.10425 | OK | stays dropped |
| letzkus APV +10 / -10 @dist | 1.03+-0.07 / 0.97+-0.06 (n None) | 1.03+-0.07 n=5; 0.97+-0.06 n=4 | Letzkus p.10426, Fig 7A,B | OK; n missing | add n=5 / n=4 |
| letzkus nopost / burst +-500 | 0.98+-0.06 (6) / 0.99+-0.05 (5) | text, Methods p.10421 | Letzkus 2006 | OK | sign of the 500 ms offset not stated (already flagged) |
| S&H 50hz +10 @distal | 0.86+-0.09 (8) | legend "distal unitary connections typically depressed (open blue circles, 86 +- 9%, n = 8)"; distal = rise > **3 ms** | Sjostrom & Hausser 2006 Fig 4C legend | OK (legend) | none |
| L2/3->L5 50hz +10 | 1.06+-0.09 (19) | text "106% +- 9%, n = 19" | S&H 2006 Results p.228 | OK (text) | none |
| zilb train10 +5 (first AP) | 0.97+-0.06 (4) | text "0.97 +- 0.06 of control, n = 4" | Zilberter p.2312, Fig 4B | OK | none |
| zilb train10 post_only | 1.03+-0.04 (4) | text "1.03 +- 0.04 of control, n = 4" | Zilberter p.2312 | OK | none |
| zilb train4 +4 _last | 0.76+-0.07 (4) | text "0.76 +- 0.07 of control, n = 4, dt = 4 ms" | Zilberter p.2315, Fig 5D | OK | none |
| zilb train8 / train-LTP AM251 | 1.15+-0.08 (8) / 1.73+-0.24 (4) | text, same values | Zilberter p.2315 | OK | none |
| 10Hz_-10 AM251 (L5, Zilberter) | 1.07+-0.08 (3) | text "1.07 +- 0.08 of control, n = 3" | Zilberter p.2315 | OK | none |
| sj 0.1hz -10 mglu / nmdar | 1.01+-0.04 (4) / 1.07+-0.04 (7) | text "0.1 (101 +- 4%, n = 4)"; "107 +- 4%, n = 7" | Sjostrom 2003 p.647 | OK (text) | none |
| sj 20hz -10 nmdar | 1.02+-0.07 (8) | fig ~1.02, n=8 | Sjostrom 2003 Fig 8B squares | OK (fig) | none |
| sj 50hz +10 ifenprodil / APV | 1.59+-0.19 (7) / 1.04+-0.05 (4) | fig ~1.58 n=7 / ~1.04 n=4 (ctrl ~1.64 n=16) | Sjostrom 2003 Fig 7D | OK (fig) | mapping note below |
| sj 0.1hz -25; burst r50 | 0.65+-0.09 (6); 0.79+-0.06 (7) | fig, same | Sjostrom 2003 Fig 9C | OK (fig) | none |
| sj07 step pair AM251 | 2.13+-0.22 (10) | text "213 +- 22%", Fig 5A n=10 | Sjostrom 2007 p.181 | OK | none |
| sj07 pre_only / post_only | 1.00+-0.071 (12) | Fig 1B "Control, n=12"; "error bars are too small to show up" | Sjostrom 2007 Fig 1B | UNVERIFIABLE SEM (assumed 0.05) | keep; label as assumed |

## 3. Readout, protocol and preparation flags (no value change)

- **Readout window.** The model scores one post-induction steady state, the first-pulse EPSP ratio. The papers average
  different windows:
  - Sjostrom 2001 and S&H 2006: 10 min after induction to the end, >= 40 min (Sj01 p.1162; S&H Methods).
  - Sjostrom 2007: from 15 min, on the first response of a 30 Hz train (p.177). This matches a first-pulse readout.
  - Letzkus: 20-30 min over the 10 min baseline (p.10421).
  - Zilberter: 5 min to the end (p.2309).
  - Markram 1997: the "maximum deviation from the baseline at any time during a 10- to 50-min period" (note 9, p.215).
    The 1.20 / 0.79 csv values sit on the late Fig 3C points and are within 0.02 of that maximum.
  - So no target uses a readout incompatible with ours. The `_last` suffix in the Zilberter ids is the dt reference
    (dt to the last AP of the train), which matches the paper's protocol. It is not a readout.
- **Sjostrom 2003 drug arms (core 9).** The AM251 and NMDAR-block tLTD experiments used "-10 or -25 ms" at 0.1-20 Hz
  (Fig 3A legend; Fig 7D "-5 to -25 ms"). The 20 Hz AM251 group may mix -10 and -25. Both control timings give LTD
  (0.65 / 0.70), so mapping it to dt -10 is acceptable. This is an inference.
- **Ifenprodil mapped to mglu_block** (50 Hz +10). The paper shows that ifenprodil blocks presynaptic NR2B NMDARs in the
  tLTD chain and spares postsynaptic NMDARs at this age (p.644). Using the eCB-off arm is therefore defensible, but it is
  a modelling mapping.
- **Letzkus location.** Letzkus splits pairs at the median rise time of 2.7 ms (50/50). Our split puts 80 of 120 pairs in
  the distal group. Letzkus 1AP +10 (core 11) is pooled over 15 pairs with rise times of 0.8-4 ms, so scoring it on all
  120 model pairs weights distal pairs 2:1. S&H "distal" uses > 3 ms, not 2.7.
- **Zilberter selection.** "connections with low release probability were discarded" (p.2309). Chindemi matched
  Zilberter only after a similar exclusion, but our L2/3->L2/3 pairs are unselected. This affects every Zilberter core.
- **n uncertainties.** Sjostrom 2001 dt 0 / +-25 have n 3-5 (7D) or 3-7 (7C); the csv's 5 is assumed. Markram +-10 n=6
  rests on "six bidirectionally coupled neurons". Neither affects chi2 (SEMs come from the figures).
- **Paired recordings.** Every core and every row above comes from a unitary paired recording. Sjostrom 2003 Fig 8B
  control triangles pool 2001 pairs, which are still paired. No extracellular row is in the fit or validation csvs here
  (Egger stays validation; Nevian is not loaded).

## 4. Consequences and open points

1. Replace core 6 by Fig 1D: `("sjostrom_40hz_dt+10ms","control"): (1.54, 0.113, 6, "Sjostrom 2001 Fig 1D 40 Hz")`.
   This is the ebner_targets.csv row, which feeds targets.py through `_ebner_csv`. With S1C pred 1.368, z changes from
   -1.16 to -1.52 (z^2 1.3 -> 2.3). Acceptance "within 1.5 SEM" then becomes |d| < 0.17 (was 0.21). The edit belongs
   to the target-file owner. As it stands, the n=11 value double-counts the 40 Hz -10 cells, which are also the
   validation row 1.51 +- 0.31.
2. Correct the wording in STAGE3_DESIGN core 11 ("at distal inputs" -> "pooled over locations"). Decide whether to
   score it on a 50/50 rise-time subsample to match Letzkus' sampling. That is a user/orchestrator decision.
3. `sjostrom_40hz_dt0ms` stays a genuine miss. Its SEM is real, so a floor would be a scoring choice and not a data fix.
4. Open: Markram n (3 or 6); Zilberter 1AP n (text 6 vs legend 10); the Sj07 single-arm SEM (assumed). These are
   unresolvable from the papers.
