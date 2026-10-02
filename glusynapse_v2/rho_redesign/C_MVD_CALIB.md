# C-MVD calibration: absolute MVD threshold from the extracted c_VDCC (2026-10-02)

**Status: section 1 done (22293154, NO gap). Section 2, the cue search, is job 22296646 (pending).**

## Method (c_mvd_calib.py via run_c_mvd.sh)
- Input: the extracted split2 records, unweighted c_VDCC only. No plasticity fit.
- MVD arithmetic from gpu_v11, applied per synapse row:
  - pool: V_{k+1} = a_k V_k + tau_E1 (1 - a_k)/i_scale s_k, with tau_E1 100 ms and i_scale 1e-5;
  - own-glutamate state: Bg += arrivals, decaying with tau_d_NMDA 70 ms;
  - glutamate window: Bg > e^(-50/70) = 0.49. W = 50 ms (Marcaggi 2009).
- **Vg** = the peak of V inside the glutamate window. This is the quantity the lower MVD edge acts on.
- Vg is reported in two units:
  - absolute c_VDCC pool units, i.e. 1e5 x VDCC charge, the unit of P1_i;
  - P1_i units, with P1_i = 1e5 x cexp vdcc_q_post (NaN rows get the table median; the NaN counts are in the log).
- Rows:
  - **Should fire** (LTD): Zilberter train10 -10 last (control and mglu_block), train4 +4 last, 5ap 10 Hz +10.
  - **Optional**: Zilberter 1AP +10 and 1AP -10.
  - **Must not fire**: L5 Sj 20 Hz -10 and 0.1 Hz -10 (the AM251 lanes), Sj 50 Hz +10 (control, NMDAR block), Sj 40 Hz +10, Markram 10 Hz +/-10.
  - Zilberter train10 +4 last appears once: its APV lane should fire and its control lane must not.
- **Lanes share the control record's VDCC trace.** AM251 and APV are kernel lanes, not separate simulations. So the Zilberter APV train (fire) and its control (no fire) have the *same* pool. Only the potentiation gating (mvd_pot 0: MVD is silent while licensed pot runs) can separate them, never a pool threshold. The same holds for every LTP "no" row.
- **Threshold scan.** For each theta, the script computes the fraction of synapses per row with Vg > theta. Two comparisons are made:
  - the minimum fraction over the fire rows against the maximum over the strict no-fire rows (the L5 -10 rows);
  - the same against all no-fire rows.
  
  A margin above 0.5 counts as a usable gap.
- **t_a convention check.** The script also prints the edge delay and the lag from t_a to the next own VDCC current peak. Here t_a = pre spike + NetCon delay (model_v2.arrival_index), which is the glutamate onset at the synapse.

## 1. Results, pool magnitude (22293154: 1:29, MaxRSS 3.1 GB of 3G, 1 CPU)
- **P1 NaN.** L5L5: 19 of 191 cexp rows are NaN. L23L23: 0 of 489. The NaN rows get the table median: 107.6 for L5 and 47.3 for L2/3.
- **No gap in absolute units.** The L5 must-not-fire rows have the same or larger Vg at every threshold. Medians are 160-300 pool units for L5 against 57-180 for Zilberter.
  - Best scan against the strict rows: theta 4.7, min fire 0.98, max no-fire 0.99, margin -0.02.
  - Against all no-fire rows the margin is also -0.02.
- **P1 units: no usable gap.** The best margin is 0.10 against the strict rows and 0.06 against all rows, both at theta 0.97 P1.
  - A threshold near 1.7 P1 catches the trains (Z train10 -10 0.93, train4 0.93).
  - It misses Z 5ap 10 Hz and 1AP.
  - It also fires in 40/50 Hz +10 (0.72-0.74) and 20 Hz -10 AM251 (0.48).
- **Conclusion: pool magnitude is not a discriminating cue.** theta_MVD_abs has no valid value, and mvd_ref 2 stays a placeholder.
- **t_a check.** The edge delay median is 0.8 ms (L2/3) and 0.95 ms (L5). The median lag from t_a to the next own VDCC current peak is:
  - dt 0: **0.75 ms**;
  - Z +4: **4.0 ms**;
  - +10: 9.5-10 ms;
  - 40/50 Hz +10: 9.9-10.1 ms.
  
  So veto_t0 = 2 lies between the dt-0 and +4 peaks, and veto_Tv = 17 lies between the +15 and +20 bAPs at the synapse. Both are confirmed.
- **Sizing note.** MaxRSS was at the 3G limit. Next CPU runs use 4G.

## 2. Cue search (c_mvd_cues.py via `run_c_mvd.sh cues`, job 22296646, 4G 0:15:00, from 22293154)
- Records, rows and glutamate window are the same as in section 1. Fire sets:
  - fire3: Z train10 -10, train4 +4, 5ap 10 Hz;
  - fire5: fire3 plus 1AP +/-10.
- No-fire sets: strict (the L5 -10 rows) and all.
- Each cue is scanned in both directions (fire when the cue is above theta, or below it). A margin above 0.5 counts as a gap.
- Output: `/scratch/dhuruva/v11_test/c_mvd_cues/c_mvd_cues_scan.csv`. Log: grep `SCAN` -A60 in `logs/c_mvd_cues_<id>.out`.

| cue | definition (own quantities only) | anchor if it separates |
|---|---|---|
| a_lowd / a_lowp | max V in the window while own c* < theta_d,i (or theta_p,i) | Zilberter Fig 5E: BAPTA turns LTP into LTD at intermediate Ca. theta_d is a Chindemi parameter, so no new parameter. Qualitative only |
| b_int / b_dur1 | integral of V over the window / P1; time above 1 P1 | none direct. The 4/8/10 AP dose concerns LTP |
| c_m51 / c_m363 | peak of an own-arrival low-pass, tau 51.3 ms or 363 ms | Marcaggi 2009 (deactivation, sensitization; summation above 20 Hz). **Expected to fail**: the Zilberter LTD rows have one pre spike per pairing, the L5 rows have pre trains |
| d_ratio | max (V/P1)/(c*/c_post) in the window | none |
| e_sh / e_shint | own shaft dCa above rest in the window (uM, uM ms) | absolute uM, the v8 licence signal. The likely separation (Zilberter shafts reach 16 uM, L5 0.5-1 uM) reflects emodel dendrite diameter (thin L2/3 obliques), so it is partly a pathway proxy |
| e_varr | V at the own arrival / P1 | Nevian order: Ca before glutamate |

Results: pending 22296646. To fill in: the best margin per cue and whether any cue reaches 0.5.
If none does, the remaining Zilberter LTD misses are a uniformity limit (user decision).

## 2. Results of the cue search (job 22296646, 1:28, 4.2 GB)

No synapse-local cue separates the Zilberter LTD rows from the L5 rows that must not fire. A gap would need a margin above 0.5 (min fire fraction minus max no-fire fraction).

| cue | best margin | note |
|---|---|---|
| pool (absolute) | 0.12 | L5 pools are as large as Zilberter pools |
| pool in P1_i units | 0.10 | |
| (a) pool while own c* < theta_d / theta_p | 0.13 / 0.04 | BAPTA-style low-Ca gate does not help |
| (b) pool integral / time above 1 P1 | 0.05 / 0.21 | trains only; fails once 1AP rows are included |
| (c) own-arrival low-pass (51 / 363 ms) | 0 | L5 rows have pre trains too |
| (d) VDCC / spine Ca ratio | 0 | |
| (e) shaft dCa peak / integral | 0.11 / 0.01 | |

Conclusion: in the current emodels and extracted Ca, the Zilberter CB1-independent LTD has no synapse-local signature that the L5 no-change rows lack. MVD with any of these cues would depress L5 as well. MVD therefore stays off (mvd_mode 0). The remaining Zilberter LTD misses are a uniformity limit of the synapse-local rule under these inputs, and that decision goes to the user.
