#!/bin/bash
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/calib_anchors_%x_%j.out
# PARAM_ANCHORS.md cheap calibrations (no fit, no NEURON). MODE=CD or MODE=NO2 (env). Each run appends its result to
# PARAM_ANCHORS.md section 5 (flock).
#  CD : dpre_min readout check. rho frozen at rho0, d at every synapse, fit readout (batch_v2.PairBasis.ratio, the
#       Option B EPSP basis incl. the cv2 term), mean over the 24 L5 subset pairs as in the objective. Target: Sjostrom
#       2003 ACEA 0.71 +/- 0.03 (AEA 0.79). Reports the d that gives 0.71 / 0.79.
#  NO2: Padamsey 2017 NO-photolysis replay of the dpre ODE (single pre spike, NO uncaged 8.5 ms later, 30 / 60 pairings
#       at 5 Hz). Reports dpre for fitted theta_Z / tau_Z and the theta_Z = 0 variant, and the A_NO needed for
#       dPr = +0.29 at a grid of baseline Pr0 (Pr0 not yet read from their Fig 5G).
# Sizing: CD from calib_vdcc_kampa 22135045 (numpy over extracted npz: 0:53, 1.16 GB) -> --mem 1500M, 0:15.
#         NO2 is pure arithmetic (no npz) -> --mem 1G, 0:15.
set -euo pipefail
source .venv/bin/activate
export OPENBLAS_NUM_THREADS=1
DOC=glusynapse_v2/rho_redesign/PARAM_ANCHORS.md
OUT=$(MODE=$MODE python -u - <<'EOF'
import glob, os, sys, math
import numpy as np
V2 = "glusynapse_v2"; sys.path.insert(0, V2)
MODE = os.environ["MODE"]; JID = os.environ.get("SLURM_JOB_ID", "?")
if MODE == "CD":
    import batch_v2 as BV
    pairs = open(f"{V2}/subset24_pairs.txt").read().strip().split(",")
    fs = {os.path.basename(f).split("__")[0]: f for f in glob.glob(f"{V2}/extracted/ebner_delta-prefire-vca/*__sjostrom_0.1hz_dt-10ms.npz")}
    D = np.round(np.linspace(-0.40, 0.0, 401), 4); R = []; cv = []
    for p in pairs:
        if p not in fs:
            continue
        z = np.load(fs[p]); syn = np.asarray(z["syn"]).astype(int); rho0 = np.asarray(z["rho0"], float)
        pb = BV.PairBasis(p, syn)
        R.append([pb.ratio(rho0, rho0, np.full(len(syn), d)) for d in D])
        bm, bs = pb.epsp(rho0, np.zeros(len(syn))); cv.append(1.0 + min((bs / bm) ** 2, 0.25))
    R = np.array(R); M = np.nanmean(R, axis=0)
    at = lambda d: M[np.argmin(abs(D - d))]
    inv = lambda y: float(np.interp(y, M, D))       # M increases with d
    print(f"C-D (job {JID}): {len(R)} L5 pairs, mean(1+cv2) {np.mean(cv):.4f}. Mean ratio at d -0.29 / -0.245 / -0.21: "
          f"{at(-0.29):.3f} / {at(-0.245):.3f} / {at(-0.21):.3f}. d giving 0.71 (ACEA) {inv(0.71):.3f}, 0.79 (AEA) {inv(0.79):.3f}. "
          f"Lock rule: dpre_min = d(0.71).")
elif MODE == "NO2":
    sys.path.insert(0, "analytical_method"); from constants import TAU_IND_GB as TAU_IND   # 70 s
    A_NO, tauZ, thZ, tauNO = 4041.994931890931, 48.738611261385365, 1.7487344595586203, 6.7
    dt = 0.01; t = np.arange(0.0, 200.0, dt); delay, isi = 8.5, 200.0
    def I_per_pairing(thz, tz, N0):
        zres = math.exp(-isi / tz) / (1 - math.exp(-isi / tz))             # steady-state residual Z at 5 Hz
        Z = (1.0 + zres) * np.exp(-(t + delay) / tz)                         # t from the uncaging time
        N = N0 * np.exp(-t / tauNO)
        return float(np.sum(np.tanh(N) * np.tanh(np.clip(Z - thz, 0, None))) * dt)
    lines = []
    for lab, thz, tz in (("fit thZ 1.75, tauZ 48.7", thZ, tauZ), ("thZ 0, tauZ 48.7", 0.0, tauZ), ("thZ 0, tauZ 10", 0.0, 10.0)):
        for N0 in (10.0, 1000.0):
            I = I_per_pairing(thz, tz, N0); r = A_NO / (1e3 * TAU_IND)
            d30, d60 = 1 - math.exp(-r * 30 * I), 1 - math.exp(-r * 60 * I)
            need = {p0: (-math.log(1 - 0.29 / p0) / (60 * I) * 1e3 * TAU_IND if I > 0 and 0.29 / p0 < 1 else float("nan"))
                    for p0 in (0.3, 0.4, 0.5, 0.6)}
            lines.append(f"[{lab}, N0 {N0:g}] gate integral {I:.2f} ms/pairing; fitted A_NO gives d {d30:.3f} (30) / {d60:.3f} (60); "
                         "A_NO for dPr +0.29 after 60 at Pr0 0.3/0.4/0.5/0.6: " + "/".join(f"{v:.0f}" for v in need.values()))
    print(f"C-NO2 (job {JID}, dpre_max 1, tau_NO 6.7 ms, NO-before-spike control gives 0 by construction): " + " | ".join(lines))
EOF
)
echo "$OUT"
flock "$DOC" -c "printf -- '- %s\n' \"\$0\" >> $DOC" "$OUT"
echo "=== appended to $DOC"
