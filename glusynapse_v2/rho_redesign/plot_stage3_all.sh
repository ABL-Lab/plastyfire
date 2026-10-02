#!/bin/bash
# Plot every finished stage-3 style fit (val csv present) with plot_stage_fit.py; the fitted set is the rows of the fit csv.
# Usage (from plastyfire/, via sjob.sh, 1 CPU 300M 0:15:00): bash glusynapse_v2/rho_redesign/plot_stage3_all.sh
# MEASURED 22287347: 3 plots 241 MB, 15 s.
RS=/project/rrg-emuller/dhuruva/plastyfire/glusynapse_v2/rho_redesign
OUTD=/project/rrg-emuller/dhuruva/plastyfire/glusynapse_v2/figs/stage_fits/stage3
for spec in "s1C_w:/scratch/dhuruva/s1c_l23w/s1CW_w025:S1C rule refit, L2/3 weight x0.25" \
            "s3C_s::3C seeded" "s3C_u::3C unseeded" "s3N_s::3N seeded" "s3N_u::3N unseeded" \
            "s4N_s::4N seeded" "s4N_u::4N unseeded" "s5N_s::5N seeded (all targets)" \
            "s3V_s::3V seeded (v11 veto, peak)" "s3W_s::3W seeded (v11 veto, integral)" \
            "s4V_s::4V seeded (v11)" "s4V_u::4V unseeded (v11)" \
            "s5Nm8_s::5N, Markram +10 weight 8, seeded" "s5Nm8_u::5N, Markram +10 weight 8, unseeded" "s4Nm8_s::4N, Markram +10 weight 8, seeded"; do
  m=${spec%%:*}; rest=${spec#*:}; fc=${rest%%:*}; lbl=${rest#*:}; [ -n "$fc" ] || fc=$m
  [ -s $RS/${m}_val.csv ] || { echo "skip $m (no val csv)"; continue; }
  [ -s $OUTD/${m}_val.png ] && [ $OUTD/${m}_val.png -nt $RS/${m}_val.csv ] && { echo "up to date $m"; continue; }
  FIT=${m}_val FITCSV=$fc FITLBL="$lbl" OUT=$OUTD/${m}_val python $RS/plot_stage_fit.py
done
