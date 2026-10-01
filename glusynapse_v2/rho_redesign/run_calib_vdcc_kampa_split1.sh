#!/bin/bash
#SBATCH --job-name=calib_kampa_split1
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=1500M
#SBATCH --time=00:15:00
#SBATCH --chdir=/lustre09/project/6070394/dhuruva/plastyfire
#SBATCH --output=/project/rrg-emuller/dhuruva/plastyfire/logs/calib_kampa_split1_%j.out
# Copy of run_calib_vdcc_kampa.sh (C-V, job 22135045) on the delta-split1 emodel (passes Kampa 200 Hz on 8/8 cells).
# MODE=V : C-V. theta_V of v5c (pool V, tau_E1 100 ms, i_scale 1e-5; same units as the C1 theta_VDCC) vs Kampa 2006:
#          3 bAPs @ 200 Hz (supralinear Ca, LTP) vs 3 bAPs @ 50 Hz (linear, no LTP) at basal synapses.
#          200 Hz arms: letzkus_3ap_200hz_dt-500ms (bAP-only; extracted by split1_ext_l23l5 step 2) and _dt+10ms
#          (EPSP-matched to S&H 50 Hz +10). Lock: Youden J (unpaired) and paired bracket F = P(dV3_50 < theta < dV3_200).
#          Also reports shaft-Ca (shaft_cai) 200/50 Hz ratio, Kampa's own observable.
# MODE=E : C-E. Maps a literature [Ca]i threshold for eCB release onto the v5c W pool (b = 0, i.e. no own glutamate,
#          as in Ca-clamp / photolysis experiments), with the model's own spine Ca (Chindemi constants: vol 0.087 um3,
#          eta 0.04, tauCa 12 ms). Steady mapping (Ca clamp) is analytic; transient mapping = V(t_bAP + 10 ms) / peak
#          free VDCC dCa of the same bAP, per synapse (L5 sjostrom_0.1hz_dt-10ms primary, L2/3 letzkus_1ap_dt+10ms).
# Inputs: glusynapse_v2/extracted/{ebner_l23l5_delta-split1-prefire-vseg-rs, ebner_delta-split1-prefire-vca}.
# Outputs: /scratch/dhuruva/param_anchors/calib_kampa_split1_{V,E}_syn.csv + this log.
# Sizing basis: 22135045 (same code path, 1 CPU): 0:53 elapsed, MaxRSS 1.16 GB (99 % of 1.17 G) -> 1.16 x 1.25 = 1.45 G
# -> 1500M; 0:53 x 1.5 -> 0:15 (minimum). Extra shaft_cai member read one at a time (<= 20 MB as float64).
# MEASURED: MODE=V 22140071: 0:46, MaxRSS 1.46 GB (97 % of 1.5 G -> next MODE=V run 1900M), CPU eff 61 %.
#           MODE=E 22140072: 0:19, MaxRSS 0.97 GB (65 %), CPU eff 74 %.
set -euo pipefail
source .venv/bin/activate
export OPENBLAS_NUM_THREADS=1
mkdir -p /scratch/dhuruva/param_anchors
MODE=${MODE:-V} python -u - <<'EOF'
import glob, os
import numpy as np, pandas as pd, h5py
MODE = os.environ["MODE"]; V2 = "glusynapse_v2"; EX = f"{V2}/extracted"; OUT = "/scratch/dhuruva/param_anchors"
TAU, ISC = 100.0, 1e-5
THETA_V5C, THETA_ECB, THETA_C1 = 0.245032531823422, 6.914872998545083, 5.500694980600818   # results/v5_V5c_s7.json, v4_C1Ajn_s5
F, VOL, ETA, TCA = 96485.33, 0.087e-15, 0.04, 12.0
KCA = ETA * 1e-9 / (2 * F * VOL) * 1e3          # uM/ms of free spine Ca per nA of VDCC current (2382)
L23D, L5D = f"{EX}/ebner_l23l5_delta-split1-prefire-vseg-rs", f"{EX}/ebner_delta-split1-prefire-vca"
EDG = {"L23": "/project/rrg-emuller/dhuruva/plastyfire/data/dhuruva_modified_edges_l23l5.h5",
       "L5": "/project/rrg-emuller/dhuruva/plastyfire/data/dhuruva_modified_edges.h5"}
POP = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"
BAS = 2.0   # SONATA afferent_section_type: 2 basal, 3 apical

def load(f):
    d = np.load(f)
    vd = np.asarray(d["vdcc"], float)
    t = np.asarray(d["t_ms"], float) if "t_ms" in d.files else float(d["t0_ms"]) + round(float(d["dt_ms"]), 4) * np.arange(vd.shape[1])
    if vd.sum() < 0:
        vd = -vd
    return d, t, vd, np.sort(np.asarray(d["postspikes"], float)), np.asarray(d["syn"]).astype(int)

def integ(x, tw, tau, gain):          # exact exponential integration, piecewise-constant input; x (n_syn, n_t)
    h = np.diff(tw); a = np.exp(-h / tau); b = tau * (1.0 - a) * gain
    Y = np.zeros((x.shape[0], len(tw)))
    for k in range(len(h)):
        Y[:, k + 1] = a[k] * Y[:, k] + b[k] * x[:, k]
    return Y

def sectypes(path, syns):
    u = np.unique(syns)
    with h5py.File(EDG[path], "r") as f:
        v = f[f"edges/{POP}/0/afferent_section_type"][u]
    return dict(zip(u.tolist(), v.astype(float).tolist()))

def files(dd, proto):
    return sorted(glob.glob(f"{dd}/*__{proto}.npz"))

q = lambda x: pd.Series({"n": x.notna().sum(), "q10": x.quantile(.1), "q50": x.median(), "q90": x.quantile(.9)})
pd.set_option("display.width", 250)

if MODE == "V":
    JOBS = [("L23", L23D, "letzkus_1ap_dt+10ms"), ("L23", L23D, "letzkus_3ap_200hz_dt-500ms"),
            ("L23", L23D, "letzkus_3ap_200hz_dt+10ms"), ("L23", L23D, "sjostrom_50hz_dt+10ms"),
            ("L5", L5D, "sjostrom_0.1hz_dt-10ms"), ("L5", L5D, "sjostrom_50hz_dt-10ms"), ("L5", L5D, "sjostrom_50hz_dt+10ms")]
    rows = []
    for path, dd, proto in JOBS:
        fs = files(dd, proto)
        for f in fs:
            pair = os.path.basename(f).split("__")[0]
            d, t, vd, post, syn = load(f)
            if post.size == 0:
                continue
            sh = np.asarray(d["shaft_cai"], float)
            bursts = np.split(post, np.where(np.diff(post) > 60.0)[0] + 1)
            res = []
            for b in bursts:
                i0 = max(np.searchsorted(t, b[0] - 15.0) - 1, 0); i1 = min(np.searchsorted(t, b[-1] + 10.0) + 1, len(t) - 1)
                tw = t[i0:i1 + 1]; V = integ(vd[:, i0:i1 + 1], tw, TAU, 1.0 / ISC)
                base = sh[:, max(np.searchsorted(t, b[0] - 1.0) - 1, 0)]
                r = np.full((vd.shape[0], 7), np.nan)
                for k in range(min(len(b), 5)):
                    w = min(8.0, b[k + 1] - b[k]) if k + 1 < len(b) else 8.0
                    m = (tw >= b[k]) & (tw <= b[k] + w)
                    if m.any():
                        r[:, k] = V[:, m].max(axis=1)
                        if k in (0, 2):
                            r[:, 5 + k // 2] = sh[:, i0:i1 + 1][:, m].max(axis=1) - base
                res.append(r)
            R = np.nanmedian(np.stack(res), axis=0)
            for i, s in enumerate(syn):
                rows.append(dict(path=path, proto=proto, pair=pair, syn=int(s), n_bursts=len(bursts), n_ap=len(bursts[0]),
                                 dV1=R[i, 0], dV3=R[i, 2], dV5=R[i, 4], dCa1=R[i, 5], dCa3=R[i, 6]))
            del d, vd, t, sh
        print(f"{path} {proto}: {len(fs)} files", flush=True)
    S = pd.DataFrame(rows)
    S["sec"] = np.nan
    for path in ("L23", "L5"):
        m = S.path == path
        if m.any():
            S.loc[m, "sec"] = S.loc[m, "syn"].map(sectypes(path, S.loc[m, "syn"].values))
    S.to_csv(f"{OUT}/calib_kampa_split1_V_syn.csv", index=False)
    print("section codes:", S.groupby("path").sec.value_counts().to_dict())
    for col in ("dV1", "dV3", "dCa1", "dCa3"):
        T = S.groupby(["path", "proto", "sec"])[col].apply(q).unstack()
        print(f"=== {col} (dV: pool V after AP k, v5c theta_V units; dCa: shaft_cai rise, mM)"); print(T.to_string(float_format=lambda v: f"{v:.4g}"))
    L200 = 1 + np.exp(-5 / TAU) + np.exp(-10 / TAU); L50 = 1 + np.exp(-20 / TAU) + np.exp(-40 / TAU)
    P = S[S.path == "L23"].pivot_table(index="syn", columns="proto", values=["dV1", "dV3", "dCa3"])
    P.columns = [f"{a}|{b}" for a, b in P.columns]
    P = P.join(S[S.path == "L23"].groupby("syn").sec.first())
    A200 = [p for p in ("letzkus_3ap_200hz_dt-500ms", "letzkus_3ap_200hz_dt+10ms") if f"dV3|{p}" in P.columns]
    for lab, sub in (("basal", P[P.sec == BAS]), ("non-basal", P[P.sec != BAS])):
        v1 = sub["dV1|letzkus_1ap_dt+10ms"]
        si50 = sub["dV3|sjostrom_50hz_dt+10ms"] / (v1 * L50)
        for p2 in A200:
            si = sub[f"dV3|{p2}"] / (v1 * L200)
            cr = sub[f"dCa3|{p2}"] / sub["dCa3|sjostrom_50hz_dt+10ms"]
            print(f"=== L2/3->L5 {lab} {p2} (n {si.notna().sum()}): pool SI dV3/(dV1 x linear) median 200 Hz {si.median():.3f} "
                  f"[q10 {si.quantile(.1):.3f}, q90 {si.quantile(.9):.3f}] vs 50 Hz {si50.median():.3f}; "
                  f"shaft dCa3 200/50 Hz median {cr.median():.3f} [q10 {cr.quantile(.1):.3f}, q90 {cr.quantile(.9):.3f}]")
    B = P[P.sec == BAS]
    for p2 in A200:
        for lab, Bv in (("L2/3->L5 basal S&H 50 Hz +10 (same synapses)", B["dV3|sjostrom_50hz_dt+10ms"]),
                        ("L5->L5 basal Sjostrom 50 Hz -10", S[(S.path == "L5") & (S.proto == "sjostrom_50hz_dt-10ms") & (S.sec == BAS)].dV3)):
            A = B[f"dV3|{p2}"].dropna().values; Bx = Bv.dropna().values
            if len(A) == 0 or len(Bx) == 0:
                print("no data for", p2, lab); continue
            grid = np.geomspace(max(min(A.min(), Bx.min()), 1e-3), max(A.max(), Bx.max()), 400)
            J = np.array([(A > x).mean() - (Bx > x).mean() for x in grid]); k = int(np.argmax(J)); ok = grid[J >= J[k] - 0.05]
            msg = (f"=== Youden {p2} vs {lab}: theta* {grid[k]:.3f} (J {J[k]:.3f}; J>=max-0.05: {ok.min():.3f}-{ok.max():.3f}); "
                   f"v5c {THETA_V5C:.3f}: P200 {(A > THETA_V5C).mean():.3f} P50 {(Bx > THETA_V5C).mean():.3f}; "
                   f"C1 {THETA_C1:.2f}: P200 {(A > THETA_C1).mean():.3f} P50 {(Bx > THETA_C1).mean():.3f}")
            if lab.startswith("L2/3"):
                pr = B[[f"dV3|{p2}", "dV3|sjostrom_50hz_dt+10ms"]].dropna().values
                Fb = np.array([((pr[:, 1] < x) & (pr[:, 0] > x)).mean() for x in grid]); kf = int(np.argmax(Fb))
                okf = grid[Fb >= Fb[kf] - 0.05]
                msg += (f"; paired n {len(pr)}, P(200>50) {(pr[:, 0] > pr[:, 1]).mean():.3f}, F* {Fb[kf]:.3f} at {grid[kf]:.3f} "
                        f"(F>=max-0.05: {okf.min():.3f}-{okf.max():.3f}); F at v5c {((pr[:, 1] < THETA_V5C) & (pr[:, 0] > THETA_V5C)).mean():.3f}")
            print(msg)

if MODE == "E":
    print(f"KCA {KCA:.1f} uM/ms per nA; steady mapping W_ss/dCa_ss = tauE1/(i_scale KCA tauCa) = {TAU / (ISC * KCA * TCA):.1f} units per uM")
    rows = []
    for path, dd, proto in (("L5", L5D, "sjostrom_0.1hz_dt-10ms"), ("L23", L23D, "letzkus_1ap_dt+10ms")):
        fs = files(dd, proto)
        for f in fs:
            pair = os.path.basename(f).split("__")[0]
            d, t, vd, post, syn = load(f)
            out = []
            for ap in post:
                i0 = max(np.searchsorted(t, ap - 15.0) - 1, 0); i1 = min(np.searchsorted(t, ap + 12.0) + 1, len(t) - 1)
                tw = t[i0:i1 + 1]; x = vd[:, i0:i1 + 1]
                V = integ(x, tw, TAU, 1.0 / ISC); C = integ(x, tw, TCA, KCA)
                m = (tw >= ap) & (tw <= ap + 8.0); j10 = min(np.searchsorted(tw, ap + 10.0), len(tw) - 1)
                out.append(np.stack([C[:, m].max(axis=1), V[:, j10], V[:, m].max(axis=1)], 1))
            R = np.nanmedian(np.stack(out), axis=0)
            for i, s in enumerate(syn):
                rows.append(dict(path=path, proto=proto, pair=pair, syn=int(s), dCa_peak_uM=R[i, 0], V_at10=R[i, 1], V_peak=R[i, 2]))
            del d, vd, t
        print(f"{path} {proto}: {len(fs)} files", flush=True)
    S = pd.DataFrame(rows); S["sec"] = np.nan
    for path in ("L23", "L5"):
        m = S.path == path
        if m.any():
            S.loc[m, "sec"] = S.loc[m, "syn"].map(sectypes(path, S.loc[m, "syn"].values))
    S["ratio"] = S.V_at10 / S.dCa_peak_uM
    S.to_csv(f"{OUT}/calib_kampa_split1_E_syn.csv", index=False)
    for col in ("dCa_peak_uM", "V_at10", "ratio"):
        print(f"=== {col}"); print(S.groupby(["path", "sec"])[col].apply(q).unstack().to_string(float_format=lambda v: f"{v:.4g}"))
    rT = S[S.path == "L5"].ratio.median(); rS = TAU / (ISC * KCA * TCA)
    print(f"transient mapping (L5 single bAP, median) {rT:.2f} units per uM peak free dCa; steady {rS:.1f}")
    print(f"fitted theta_eCB {THETA_ECB:.2f} <-> transient peak dCa {THETA_ECB / rT:.3f} uM, steady dCa {THETA_ECB / rS:.4f} uM")
    for lab, ca in (("0.1 uM", 0.1), ("0.2 uM", 0.2), ("0.3 uM", 0.3), ("0.5 uM", 0.5), ("1 uM", 1.0),
                    ("Wang&Zucker 2001 DSI half-max 3.6 uM", 3.6), ("3.9 uM", 3.9), ("Brenowitz&Regehr 2003 half-max 15 uM", 15.0)):
        tT, tS = ca * rT, ca * rS
        pT = (S[S.path == "L5"].V_at10 > tT).mean(); pS = (S[S.path == "L5"].V_at10 > tS).mean()
        print(f"Ca_thr {lab}: theta_eCB transient {tT:.2f} (P single-bAP W>theta at +10 ms, L5: {pT:.3f}); steady {tS:.1f} (P {pS:.3f})")
    print(f"at fitted {THETA_ECB:.2f}: P(single-bAP V_at10 > theta_eCB) L5 {(S[S.path == 'L5'].V_at10 > THETA_ECB).mean():.3f}, "
          f"L2/3->L5 {(S[S.path == 'L23'].V_at10 > THETA_ECB).mean():.3f}")
EOF
