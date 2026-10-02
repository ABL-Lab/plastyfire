#!/bin/bash
# MEASURED 22292756 core 1:49 7.9 GB; 22292757 L5-all (3 rescores) 7:31 37.6 GB, so l5all 47G 0:15:00.
# Tests of gpu_v11_rho.py (SPEC_V11.md): fixed-parameter rescores of s1C_s.json (--maxiter 0), L5 model only.
# Usage (from plastyfire/, via sjob.sh -g):  bash glusynapse_v2/rho_redesign/test_v11.sh core|l5all
#   core   the 7 L5 core targets, as run_stage_fit.sh C s:  a = v11 at defaults (must match s1C_s.csv to 1e-6: REPRO11 OK),
#          b = a + mvd_mode 1 (theta_MVD_lo 0.37 placeholder, no upper edge), c = a + veto_t0 3, veto_Tv 17, veto_peak 1.
#          SIZE: 11G 0:15:00 (L5-only 7-11-target fits / rescores MEASURED 22286005-17: 1:20-2:20, 7.1-8.5 GB; 3 runs in series).
#   l5all  the same a / b / c on all 40 L5 targets (only the DROPT rows dropped): a40 is compared with the L5 rows of
#          s1C_s_val.csv; b40 shows whether MVD breaks the L5 cores (20 Hz -10 AM251, 40 / 50 Hz LTP, Markram).
#          SIZE: PILOT 40G 0:20:00 (22289864/65, the same 40-target L5 rescores under v10, were OOM-killed at 10.0 GB).
# Outputs /scratch/dhuruva/v11_test/<tag>.{json,csv}. Log grep: "=== |REPRO11|v11 calib|v10 counts|NaN|chi2|Traceback".
set -euo pipefail
MODE=${1:-core}
R=/project/rrg-emuller/dhuruva/plastyfire; V2=glusynapse_v2; RS=$V2/rho_redesign; W=/scratch/dhuruva/s2g0321; X=$W/extracted; T=delta-split2-ljp25g0321-prefire
source glusynapse_v2/env_v3.sh
export CEXP_DIR=$W/cexp ANALYTICAL_BASIS_DIR=$W/basis_l5l5 L23_BASIS_DIR=$W/basis_l23l5 PYTHONPATH=$R
O=/scratch/dhuruva/v11_test; mkdir -p $O
DIRS=$X/ebner_$T-vca,$X/markram_$T-vca,$X/sj03_$T-vca,$X/sj03r50_$T-vca,$X/sj07_$T-vca,$X/l5extra_$T-vca
K5=$W/doublets/l5_keep_pairs.txt
PAIRS=$(comm -12 <(tr ',' '\n' < $V2/subset24_pairs.txt | sed '/^$/d' | sort -u) <(tr ',' '\n' < $K5 | sed '/^$/d' | sort -u) | paste -sd,)
G8X='"v5_mode": 2, "veto_T": 25.0, "ecb_ref": 2, "gate_src": 2, "gate_win": 100.0, "gate_theta": 0.15, "theta_V": 3.5751446275588403, "t_exact": 1, "v10_counts": 1'
FILTERS='{"pre_drive": 1, "i_scale": 1e-5, "t_drive": 4, "tau_E1": 100.0, "vamp_mode": 1, "v5_mode": 1, "rho_gamma": 1.0, "A_eCB": 1.0, "dpre_min": -0.29}'
DROPT="letzkus_3ap_200hz_dt-10ms@distal|control,letzkus_3ap_200hz_dt-10ms@distal|nmdar_block"
echo "=== test_v11 $MODE $(date)"; nvidia-smi -L
python -m py_compile $RS/gpu_v11_rho.py && echo "=== py_compile OK"

# P1_i source: NaN fraction of cexp vdcc_q_post per pathway table (rows without a value get the pooled median in fit_v6)
for pw in L5L5 L23L5 L23L23; do
  awk -F, -v pw=$pw 'NR==1{for(i=1;i<=NF;i++) if($i=="vdcc_q_post") c=i; next}
    {n++; v=tolower($c); if(v==""||v=="nan") m++} END{printf "=== P1 NaN %s: %d / %d rows (%.3f) of cexp vdcc_q_post\n", pw, m, n, (n?m/n:0)}' $W/cexp/$pw.csv
done

run(){ # tag set drop
  rm -f $O/$1.json $O/$1.csv
  python -u $RS/gpu_v11_rho.py --dirs $DIRS --pairs $PAIRS --groups paired_l5,sjostrom07,paired_l5_extra --filters "$FILTERS" \
    --free-filters theta_eCB --set "$2" --fit-gamma --drop-targets "$3" \
    --seed-fits $RS/s1C_s.json --seed-set '{}' --seed 5 --maxiter 0 --save $O/$1 || echo "exit $? ($1; REPRO DIFF expected)"
  if [ -s $O/$1.csv ]; then echo "=== $1 written"; else echo "=== $1 MISSING"; fi
}

cmp(){ # ref tol csv...   (REPRO11 line for the first csv vs ref, then a side-by-side table)
  python - "$@" <<'EOF'
import sys, pandas as pd
ref, tol, fs = sys.argv[1], float(sys.argv[2]), sys.argv[3:]
R = pd.read_csv(ref); K = ["target", "condition"]
import os
if not os.path.isfile(fs[0]):
    print(f"REPRO11 FAIL: {fs[0]} missing"); sys.exit(0)
fs = [f for f in fs if os.path.isfile(f)]
D = [pd.read_csv(f) for f in fs]
m = D[0].merge(R, on=K, suffixes=("", "_ref"))
d = float((m["pred"] - m["pred_ref"]).abs().max())
ok = len(m) == len(R) and d <= tol
print(f"REPRO11 {'OK' if ok else 'FAIL'}: {fs[0]} vs {ref}: {len(m)} of {len(R)} rows, max |pred - ref| {d:.3g} (tol {tol:g})")
T = R[K + ["target_mean", "target_sem"]].copy()
for f, df in zip(fs, D):
    tag = f.rsplit("/", 1)[-1][:-4]
    T = T.merge(df[K + ["pred"]].rename(columns={"pred": tag}), on=K, how="left")
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)
print(T.round(4).to_string(index=False))
for f, df in zip(fs, D):
    z = (df["pred"] - df["target_mean"]) / df["target_sem"]
    print(f"chi2 {f.rsplit('/', 1)[-1]}: {float((z ** 2).sum()):.3f} over {len(df)} targets")
EOF
}

SB="{$G8X, \"mvd_mode\": 1}"
SC="{$G8X, \"veto_t0\": 3.0, \"veto_Tv\": 17.0, \"veto_peak\": 1}"
case $MODE in
  core)
    KEEP='10Hz_10ms|control,10Hz_-10ms|control,sjostrom_0.1hz_dt+10ms|control,sjostrom_0.1hz_dt-10ms|control,sjostrom_20hz_dt+10ms|control,sjostrom_40hz_dt+10ms|control,sjostrom_50hz_dt-10ms|control'
    DROP=$(awk -F, -v keep="$KEEP" 'BEGIN{n=split(keep,k,","); for(i=1;i<=n;i++) K[k[i]]=1} NR>1{t=$1"|"$2; if(!(t in K)) printf "%s%s", (c++?",":""), t}' $RS/r1D_s.csv)
    run a "{$G8X}" "$DROPT,$DROP"
    run b "$SB" "$DROPT,$DROP"
    run c "$SC" "$DROPT,$DROP"
    cmp $RS/s1C_s.csv 1e-6 $O/a.csv $O/b.csv $O/c.csv ;;
  l5all)
    run a40 "{$G8X}" "$DROPT"
    run b40 "$SB" "$DROPT"
    run c40 "$SC" "$DROPT"
    cmp $RS/s1C_s_val.csv 1e-6 $O/a40.csv $O/b40.csv $O/c40.csv ;;
  *) echo "MODE core|l5all" >&2; exit 1 ;;
esac
echo "=== done $(date)"
