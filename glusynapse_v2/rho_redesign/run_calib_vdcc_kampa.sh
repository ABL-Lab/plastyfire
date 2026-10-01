#!/bin/bash
#SBATCH --job-name=calib_vdcc_kampa
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=1200M
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/calib_vdcc_kampa_%j.out
# PARAM_ANCHORS.md calibration C-V (theta_VDCC vs Kampa et al. 2006 basal-dendrite Ca imaging). No fit, no NEURON:
# reads existing extracted -ica_VDCC traces and computes the C1 gate signal c_VDCC (tau_E1 100 ms, i_scale 1e-5)
# after the 1st/3rd/5th postsynaptic AP of each burst, per synapse, for 1 AP, 3 AP @ 200 Hz (Letzkus -500 ms: EPSP
# 500 ms later, i.e. bAP-only) and 50 Hz trains (S&H +10 on L2/3->L5 synapses, Sjostrom -10/+10 on L5->L5).
# Lock rule: theta_VDCC* = the level that best separates 3 AP @ 200 Hz (supralinear, Ni-sensitive Ca, LTP) from
# 3 AP @ 50 Hz (linear Ca, no LTP) at basal synapses (Youden J). Also reports the supralinearity of c_VDCC.
# Sizing basis: diag_dltd 22132795 (numpy replay over extracted npz, 1 CPU): 0:12, 910 MB MaxRSS, 100% CPU
# -> 910 MB x 1.25 = 1.14 GB -> 1200M; time 0:12 x 1.5 -> rounded up to 0:15. Only vdcc/t/spikes are loaded
# (one npz member at a time, <= ~10 MB each), so memory should stay below the diag_dltd peak.
# MEASURED 22135045: 0:53, MaxRSS 1.16 GB (99% of 1.17 -> next run 1500M), CPU eff 60%. Result: null (no 200 Hz supralinearity).
set -euo pipefail
source .venv/bin/activate
export OPENBLAS_NUM_THREADS=1
mkdir -p /scratch/dhuruva/param_anchors
python -u - <<'EOF'
import glob, os, sys
import numpy as np, pandas as pd
V2 = "glusynapse_v2"; EX = f"{V2}/extracted"; OUT = "/scratch/dhuruva/param_anchors"
sys.path.insert(0, V2)
TAU, ISC = 100.0, 1e-5          # tau_E1 (ms), i_scale (nA): the fitted C1 gate constants
THETA_FIT = 5.500694980600818   # v4_C1Ajn_s5 theta_V
st = np.load(f"{V2}/syn_section_type.npz"); styp = dict(zip(st["syn"].tolist(), st["section_type"].tolist()))
dist5 = {}
for f in glob.glob(f"{V2}/local_t/out/*.npz"):
    z = np.load(f)
    if "og|ap1|sid" in z.files:
        for s, dd in zip(z["og|ap1|sid"], z["og|ap1|dist"]):
            dist5[int(s)] = float(dd)
g = pd.read_csv("ebner/pair_geometry_L23PC_L5TTPC.csv"); g["pair"] = g.pregid.astype(str) + "-" + g.postgid.astype(str)
ldist = dict(zip(g.pair, g.letzkus_distal.astype(bool)))
L23D, L5D = f"{EX}/ebner_l23l5_delta-prefire-vseg-rs", f"{EX}/ebner_delta-prefire-vca"
JOBS = [("L23", L23D, "letzkus_1ap_dt+10ms"), ("L23", L23D, "letzkus_3ap_200hz_dt-500ms"),
        ("L23", L23D, "letzkus_3ap_200hz_dt+10ms"), ("L23", L23D, "sjostrom_50hz_dt+10ms"),
        ("L5", L5D, "sjostrom_0.1hz_dt-10ms"), ("L5", L5D, "sjostrom_50hz_dt-10ms"), ("L5", L5D, "sjostrom_50hz_dt+10ms")]
rows, nneg = [], 0
for path, dd, proto in JOBS:
    fs = sorted(glob.glob(f"{dd}/*__{proto}.npz"))
    for f in fs:
        pair = os.path.basename(f).split("__")[0]
        d = np.load(f)
        t = np.asarray(d["t_ms"], float) if "t_ms" in d.files else None
        vd = np.asarray(d["vdcc"], float)
        if t is None:
            t = float(d["t0_ms"]) + round(float(d["dt_ms"]), 4) * np.arange(vd.shape[1])
        if vd.sum() < 0:
            vd = -vd; nneg += 1
        post = np.sort(np.asarray(d["postspikes"], float)); syn = np.asarray(d["syn"]).astype(int)
        if post.size == 0:
            continue
        cut = np.where(np.diff(post) > 60.0)[0] + 1
        bursts = np.split(post, cut)
        res = []                                   # per burst: (n_syn, 5) increments after AP k (nan if absent)
        for b in bursts:
            i0 = max(np.searchsorted(t, b[0] - 15.0) - 1, 0); i1 = min(np.searchsorted(t, b[-1] + 10.0) + 1, len(t) - 1)
            tw = t[i0:i1 + 1]; h = np.diff(tw); aV = np.exp(-h / TAU); bV = TAU * (1.0 - aV) / ISC
            V = np.zeros((vd.shape[0], len(tw)))
            for k in range(len(h)):
                V[:, k + 1] = aV[k] * V[:, k] + bV[k] * vd[:, i0 + k]
            r = np.full((vd.shape[0], 5), np.nan)
            for k in range(min(len(b), 5)):
                w = min(8.0, b[k + 1] - b[k]) if k + 1 < len(b) else 8.0
                m = (tw >= b[k]) & (tw <= b[k] + w)
                if m.any():
                    r[:, k] = V[:, m].max(axis=1)
            res.append(r)
        R = np.nanmedian(np.stack(res), axis=0)
        for i, s in enumerate(syn):
            rows.append(dict(path=path, proto=proto, pair=pair, syn=int(s), sec=styp.get(int(s), np.nan),
                             dist=dist5.get(int(s), np.nan), letzkus_distal=ldist.get(pair), n_bursts=len(bursts),
                             n_ap=len(bursts[0]), dV1=R[i, 0], dV3=R[i, 2], dV5=R[i, 4]))
        del d, vd, t
    print(f"{path} {proto}: {len(fs)} files", flush=True)
S = pd.DataFrame(rows)
miss = sorted(set(S.syn[S.sec.isna()]))
if miss:
    try:
        import h5py, model_v2 as MV
        with h5py.File(MV.EDGES, "r") as f:
            v = f[f"edges/{MV.EDGE_POP}/0/afferent_section_type"][np.array(miss)]
        S["sec"] = S.sec.fillna(S.syn.map(dict(zip(miss, v.astype(float)))))
    except Exception as e:
        print("section-type fallback failed:", e)
S.to_csv(f"{OUT}/calib_vdcc_kampa_syn.csv", index=False)
print("sign flips (vdcc stored as ica):", nneg, "; section codes:", S.sec.value_counts(dropna=False).to_dict())
pd.set_option("display.width", 250)
q = lambda x: pd.Series({"n": x.notna().sum(), "q10": x.quantile(.1), "q50": x.median(), "q90": x.quantile(.9)})
for col in ("dV1", "dV3", "dV5"):
    T = S.groupby(["path", "proto", "sec"])[col].apply(q).unstack()
    print(f"=== {col} (c_VDCC increment after AP k, units of theta_V; fit theta_V {THETA_FIT:.2f})"); print(T.round(3).to_string())
BAS = 2.0   # basal code (diag_l23_ltp convention: L5->L5 median 2, L2/3->L5 median 3)
L200 = 1 + np.exp(-5 / TAU) + np.exp(-10 / TAU); L50 = 1 + np.exp(-20 / TAU) + np.exp(-40 / TAU)
P = S[S.path == "L23"].pivot_table(index="syn", columns="proto", values=["dV1", "dV3"])
P.columns = [f"{a}|{b}" for a, b in P.columns]
P = P.join(S[S.path == "L23"].groupby("syn").sec.first())
for lab, sub in (("basal", P[P.sec == BAS]), ("non-basal", P[P.sec != BAS])):
    v1 = sub["dV1|letzkus_1ap_dt+10ms"]
    si200 = sub["dV3|letzkus_3ap_200hz_dt-500ms"] / (v1 * L200); si50 = sub["dV3|sjostrom_50hz_dt+10ms"] / (v1 * L50)
    print(f"=== L2/3->L5 {lab} (n {len(sub)}): supralinearity dV3/(dV1 x linear sum) median: 200 Hz {si200.median():.3f}"
          f" [q10 {si200.quantile(.1):.3f}, q90 {si200.quantile(.9):.3f}], 50 Hz {si50.median():.3f}")
A = P[P.sec == BAS]["dV3|letzkus_3ap_200hz_dt-500ms"].dropna().values
for lab, B in (("L2/3->L5 basal S&H 50 Hz (+10, same synapses)", P[P.sec == BAS]["dV3|sjostrom_50hz_dt+10ms"].dropna().values),
               ("L5->L5 basal Sjostrom 50 Hz -10", S[(S.path == "L5") & (S.proto == "sjostrom_50hz_dt-10ms") & (S.sec == BAS)].dV3.dropna().values)):
    if len(A) == 0 or len(B) == 0:
        print("no data for", lab); continue
    grid = np.geomspace(max(min(A.min(), B.min()), 1e-3), max(A.max(), B.max()), 400)
    J = np.array([(A > x).mean() - (B > x).mean() for x in grid]); k = int(np.argmax(J))
    ok = grid[J >= J[k] - 0.05]
    print(f"=== theta* vs {lab}: theta* {grid[k]:.3f} (J {J[k]:.3f}; within 0.05 of max: {ok.min():.3f}-{ok.max():.3f}); "
          f"at fit {THETA_FIT:.2f}: P(200 Hz open) {(A > THETA_FIT).mean():.3f}, P(50 Hz open) {(B > THETA_FIT).mean():.3f}")
print("units: theta x i_scale = VDCC charge in nA*ms (pC); 1 unit = 1e-5 pC = 31.2 Ca ions (z=2) = 0.60 uM total Ca in 0.087 um3")
EOF
