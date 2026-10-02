#!/bin/bash
# MEASURED 22306937: 15:07, 33.4 GB MaxRSS (REPRO + 7 L5 rescores), so a rerun is 42G 0:30:00 on -g.
# Tests of the eCB LP mode in gpu_v11_rho.py (ecb_lp_mode, SPEC_ECB_LP.md). L5 model only, fixed parameters (--maxiter 0).
# Usage (from plastyfire/, via sjob.sh -g):  bash glusynapse_v2/rho_redesign/test_ecb_lp.sh
#   1. repro  test_v11.sh core "a": v11 at defaults on the 7 L5 cores, seed s1C_s.json; must match s1C_s.csv to 1e-6 (REPRO11 OK).
#   2. l5all  all 40 L5 targets (only DROPT dropped), 7 rescores:
#        C0  s3C_u.json, its SET (G8X, v5_mode 2 = (1 - b)-weighted trigger, v7 25 ms veto)       mode off (reference)
#        C1  C0 + ecb_lp_mode 1 (10 ms LP trigger, no veto; Ebner)
#        C2  C0 + ecb_lp_mode 2 + veto window (2, 17] ms integral (veto_t0 2, veto_Tv 17, veto_peak 0, as s3W_s)
#        C2t C2 with theta_eCB 0.07 (seed json copy; probes the +-25 lag: one bAP 23.5 ms back gives L ~ e^-2.35 P1)
#        W0  s3W_s.json, its SET (G8V, v5_mode 1 = unweighted trigger, window (2, 17] integral)    mode off (reference)
#        W1  G8V + ecb_lp_mode 1
#        W2  W0 + ecb_lp_mode 2
# SIZE: 47G 0:30:00 on -g. Basis: test_v11.sh MEASURED 22292757 L5-all 3 rescores 7:31, 37.6 GB (so ~2.5 min each, same
#   memory per run since runs are serial); 22292756 core 3 runs 1:49. 7 x 2.5 + 0.6 = 18 min, + 50% = 27 -> 0:30:00.
# Outputs /scratch/dhuruva/ecb_lp_test/<tag>.{json,csv}. Log grep: "=== |REPRO11|v11|v10 counts|chi2|Traceback".
set -euo pipefail
R=/project/rrg-emuller/dhuruva/plastyfire; V2=glusynapse_v2; RS=$V2/rho_redesign; W=/scratch/dhuruva/s2g0321; X=$W/extracted; T=delta-split2-ljp25g0321-prefire
source glusynapse_v2/env_v3.sh
export CEXP_DIR=$W/cexp ANALYTICAL_BASIS_DIR=$W/basis_l5l5 L23_BASIS_DIR=$W/basis_l23l5 PYTHONPATH=$R
O=/scratch/dhuruva/ecb_lp_test; mkdir -p $O
DIRS=$X/ebner_$T-vca,$X/markram_$T-vca,$X/sj03_$T-vca,$X/sj03r50_$T-vca,$X/sj07_$T-vca,$X/l5extra_$T-vca
K5=$W/doublets/l5_keep_pairs.txt
PAIRS=$(comm -12 <(tr ',' '\n' < $V2/subset24_pairs.txt | sed '/^$/d' | sort -u) <(tr ',' '\n' < $K5 | sed '/^$/d' | sort -u) | paste -sd,)
G8X='"v5_mode": 2, "veto_T": 25.0, "ecb_ref": 2, "gate_src": 2, "gate_win": 100.0, "gate_theta": 0.15, "theta_V": 3.5751446275588403, "t_exact": 1, "v10_counts": 1'
G8V=${G8X/\"v5_mode\": 2/\"v5_mode\": 1}
WIN='"veto_t0": 2, "veto_Tv": 17, "veto_peak": 0, "veto_peak_k": 0.5'
FILTERS='{"pre_drive": 1, "i_scale": 1e-5, "t_drive": 4, "tau_E1": 100.0, "vamp_mode": 1, "v5_mode": 1, "rho_gamma": 1.0, "A_eCB": 1.0, "dpre_min": -0.29}'
DROPT="letzkus_3ap_200hz_dt-10ms@distal|control,letzkus_3ap_200hz_dt-10ms@distal|nmdar_block"
echo "=== test_ecb_lp $(date)"; nvidia-smi -L
python -m py_compile $RS/gpu_v11_rho.py && echo "=== py_compile OK"

run(){ # tag set drop seed
  rm -f $O/$1.json $O/$1.csv
  echo "=== run $1 seed $4 SET $2"
  python -u $RS/gpu_v11_rho.py --dirs $DIRS --pairs $PAIRS --groups paired_l5,sjostrom07,paired_l5_extra --filters "$FILTERS" \
    --free-filters theta_eCB --set "$2" --fit-gamma --drop-targets "$3" \
    --seed-fits $4 --seed-set '{}' --seed 5 --maxiter 0 --save $O/$1 || echo "exit $? ($1; REPRO DIFF expected)"
  if [ -s $O/$1.csv ]; then echo "=== $1 written"; else echo "=== $1 MISSING"; fi
}

# 1. REPRO (test_v11.sh core a)
KEEP='10Hz_10ms|control,10Hz_-10ms|control,sjostrom_0.1hz_dt+10ms|control,sjostrom_0.1hz_dt-10ms|control,sjostrom_20hz_dt+10ms|control,sjostrom_40hz_dt+10ms|control,sjostrom_50hz_dt-10ms|control'
DROP=$(awk -F, -v keep="$KEEP" 'BEGIN{n=split(keep,k,","); for(i=1;i<=n;i++) K[k[i]]=1} NR>1{t=$1"|"$2; if(!(t in K)) printf "%s%s", (c++?",":""), t}' $RS/r1D_s.csv)
run a "{$G8X}" "$DROPT,$DROP" $RS/s1C_s.json

# 2. L5-all rescores
sed -E 's/("theta_eCB": )[-0-9.eE+]+/\10.07/g' $RS/s3C_u.json > $O/s3C_u_th007.json
echo "=== theta_eCB 0.07 entries in seed copy: $(grep -c '"theta_eCB": 0.07' $O/s3C_u_th007.json || true)"
run C0  "{$G8X}"                               "$DROPT" $RS/s3C_u.json
run C1  "{$G8X, \"ecb_lp_mode\": 1}"           "$DROPT" $RS/s3C_u.json
run C2  "{$G8X, $WIN, \"ecb_lp_mode\": 2}"     "$DROPT" $RS/s3C_u.json
run C2t "{$G8X, $WIN, \"ecb_lp_mode\": 2}"     "$DROPT" $O/s3C_u_th007.json
run W0  "{$G8V, $WIN}"                         "$DROPT" $RS/s3W_s.json
run W1  "{$G8V, \"ecb_lp_mode\": 1}"           "$DROPT" $RS/s3W_s.json
run W2  "{$G8V, $WIN, \"ecb_lp_mode\": 2}"     "$DROPT" $RS/s3W_s.json

python - $RS/s1C_s.csv $O <<'EOF'
import sys, os, pandas as pd
ref, O = sys.argv[1], sys.argv[2]
K = ["target", "condition"]
f = f"{O}/a.csv"
if os.path.isfile(f):
    R = pd.read_csv(ref); m = pd.read_csv(f).merge(R, on=K, suffixes=("", "_ref"))
    d = float((m["pred"] - m["pred_ref"]).abs().max())
    ok = len(m) == len(R) and d <= 1e-6
    print(f"REPRO11 {'OK' if ok else 'FAIL'}: a.csv vs s1C_s.csv: {len(m)} of {len(R)} rows, max |pred - ref| {d:.3g} (tol 1e-6)")
else:
    print("REPRO11 FAIL: a.csv missing")
tags = ["C0", "C1", "C2", "C2t", "W0", "W1", "W2"]
D = {t: pd.read_csv(f"{O}/{t}.csv") for t in tags if os.path.isfile(f"{O}/{t}.csv")}
if not D:
    sys.exit("no l5all csv")
base = next(iter(D.values()))
Tb = base[K + ["target_mean", "target_sem"]].copy()
for t, df in D.items():
    Tb = Tb.merge(df[K + ["pred"]].rename(columns={"pred": t}), on=K, how="left")
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)
focus = ["sjostrom_20hz_dt0ms", "sjostrom_40hz_dt0ms", "sjostrom_20hz_dt+25ms", "sjostrom_20hz_dt-25ms", "10Hz_10ms",
         "10Hz_-10ms", "sjostrom_50hz_dt-10ms"]
F = Tb[Tb["target"].isin(focus) & (Tb["condition"] == "control") | Tb["target"].str.startswith("sjostrom07")]
print("=== focus rows (pred)"); print(F.round(3).to_string(index=False))
print("=== all L5 rows (pred)"); print(Tb.round(3).to_string(index=False))
for t, df in D.items():
    z = (df["pred"] - df["target_mean"]) / df["target_sem"]
    fz = z[df["target"].isin(focus) | df["target"].str.startswith("sjostrom07")]
    print(f"chi2 {t}: L5 total {float((z ** 2).sum()):.3f} over {len(df)} targets | focus rows {float((fz ** 2).sum()):.3f} "
          f"| |z| >= 2: {int((z.abs() >= 2).sum())}")
EOF
echo "=== done $(date)"
