# C-MVD calibration: absolute MVD threshold from the extracted c_VDCC (2026-10-02)

**Status: job 22293154 submitted (CPU pilot, 3G 0:30:00). The results section is filled in after the job.**

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

## Results
Pending 22293154. Log: `logs/c_mvd_22293154.out` (grep `P1 |SCAN|fraction|done`). Tables: `/scratch/dhuruva/v11_test/c_mvd/`.
To fill in after the job:
- per-row Vg quantiles (abs and P1 units);
- whether there is a gap, and theta_MVD_abs;
- the lag of the dt-0 and +10 bAP peaks against veto_t0 = 2;
- the seff line.

If there is no gap, MVD cannot be uniform with this emodel, and the user decides (ANCHORS_V11, open point 1).
