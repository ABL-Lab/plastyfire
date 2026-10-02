# Core fit targets and staged fitting (2026-10-02)

The user rejected the 65-target joint fit (R1D). It reached a low chi2 by pulling nearly every prediction towards an EPSP ratio of 1, and it missed the Markram 10 Hz +10 ms LTP (1.20 measured, 0.91 predicted). Two lessons follow. The rule must first reproduce the defining STDP features, and a few meaningful targets must be fitted instead of many weak ones. Every other target is validation.

## Core targets (18; everything else is validation)

L5->L5, STDP shape and frequency dependence (Markram 1997; Sjostrom 2001):
1. `10Hz_10ms|control` (LTP, 1.20)
2. `10Hz_-10ms|control` (LTD, 0.79)
3. `sjostrom_0.1hz_dt+10ms|control` (no LTP at low frequency, 0.97)
4. `sjostrom_0.1hz_dt-10ms|control` (LTD, 0.69)
5. `sjostrom_20hz_dt+10ms|control` (LTP, 1.31)
6. `sjostrom_40hz_dt+10ms|control` (LTP, 1.53)
7. `sjostrom_50hz_dt-10ms|control` (LTP at high frequency in either order, 1.70)

L5->L5, mechanisms (Sjostrom 2001/2003/2007):
8. `sjostrom_50hz_dt+10ms|nmdar_block` (LTP needs NMDAR, 1.04)
9. `sjostrom_0.1hz_dt-10ms|mglu_block` (timing LTD is eCB/mGluR dependent, 1.01)
10. `sjostrom07_step200ms_pair|control` (LTP with a depolarising step, 1.62)
11. `sjostrom07_step200ms_pair|no_block` (part of that LTP is NO dependent, 1.36)

L2/3->L5 (Letzkus 2006):
12. `letzkus_1ap_dt+10ms|control` (single-AP pairing at distal inputs, LTD, 0.72)
13. `letzkus_3ap_200hz_dt+10ms@proximal|control` (burst pairing, proximal, LTP, 1.30)
14. `letzkus_3ap_200hz_dt+10ms@distal|control` (burst pairing, distal, 0.79)

L2/3->L2/3 (Zilberter 2009):
15. `zilberter_1ap_dt+10ms|control` (single-AP pre-post, LTD, 0.64)
16. `zilberter_1ap_dt-10ms|control` (single-AP post-pre, LTD, 0.56)
17. `zilberter_train10_50hz_dt+4ms_last|control` (train pairing, LTP, 1.49)
18. `zilberter_train10_50hz_dt-10ms_last|control` (train post-pre, LTD, 0.72)

Egger 1999 is validation only (L4 spiny stellate cells, user decision 2026-10-02).

## Stages (each stage seeds the next)

1. **Stage 1: L5 STDP shape.** Targets 1-7. Free: a00, a01, a10, a11, gamma_d, gamma_p. Rules: Chindemi (no licence, no eCB, no NO), Chindemi + shaft licence, and shaft licence + eCB (k_E free). Pass: Markram +/-10 ms within 1 SEM, Sjostrom frequency points within 1.5 SEM. If Chindemi passes and the licence does not, the licence is what kills Markram LTP.
2. **Stage 2: L5 mechanisms.** Targets 1-11. Add k_E, A_NO and d_NO_max as needed. Seeded from stage 1.
3. **Stage 3: other pathways.** Targets 1-18 with the same uniform parameters. First rescore the stage 2 fit without refitting (a transfer test), then refit seeded from stage 2.
4. **Validation.** All other targets (47) are rescored, reported and never fitted.
