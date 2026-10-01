"""LTD diagnosis of the v5c rule (LTD_DIAG.md). Read-only.

Replays the gpu_v5_rho kernel (v5_mode 2 = v5c) on the CPU, one record at a time (one BatchV2 per npz, as
diag_l23l23_v5.py), for every record of ONE pathway model (--model l5 | l23 | l23l23) at the parameters of ONE v5 fit
json, plus fixed-parameter counterfactual variants (no refit). Every variant is uniform and reads only own spine Ca:

  base     v5c as fitted (must reproduce the fit csv pred)
  noE      eCB step off (rho-only LTD)
  allE     eCB step at every own arrival (no trigger): the eCB bound with rho dynamics
  thE.3/.1/.01   theta_eCB x 0.3 / 0.1 / 0.01 (same W pool)
  Cs<f>    trigger pool Wc = (1 - b)-weighted copy of c* (own spine Ca influx, tau* = c*'s time constant);
           trigger at own arrival if Wc > f c_post (f new: fraction of the own bAP c*)
  Cf<f>    as Cs with the pool time constant tau_E1 (the V/W pool constant): fast copy of the same Ca input
  dm.40/.50      d_min -0.40 / -0.50 (base trigger)
  td.9/.8  theta_d x 0.9 / 0.8 (a00, a01 direction; Chindemi parameters, unrefit)
  Cf.2dm4  Cf0.2 with d_min -0.40
  round 2 (LTD_LIT.md point 5), bAP veto: a triggered own arrival is held for Tv ms and takes no eCB step if in that
  window (vV) own unweighted VDCC influx > theta_eCB, (vC) own c* influx > c_pre + kappa c_post, (vP) an own pot
  event occurs. Names: <trigger><veto><Tv>, trigger t01 = thE.01, t1 = thE.1, a = allE, Cf.5 = Cf0.5.

Per synapse (base): sec_type (edges afferent_section_type: 2 basal, 3 apical), dist (L5: local_t path distance um;
else edges distance_soma when present), rho0, rho_f, d0, ecb (got >= 1 eCB step), down (rho0 >= .5 -> rho_f < .5),
dep_any (c* > theta_d without pot at any step), Vpk (max V / theta_V), Wpk / W_arr (max W, max W at own arrivals,
in theta_eCB), vdcc_pk (max -ica_VDCC), pk_td / pk_tp (max c* / theta_d, theta_p), td_tp, Cs_arr / Cf_arr (max Wc at own
arrivals / c_post), n_arr, n_trig.
Outputs (--out, scratch): <tag>_<model>_ltd_{syn,targets,mech}.csv; summary in the job log.

    python glusynapse_v2/rho_redesign/diag_ltd.py --fit <v5 json> --model l5 --out /scratch/dhuruva/ltd_diag
"""
import argparse, glob, json, math, os, sys, time
import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__)); V2 = os.path.dirname(HERE); ROOT = os.path.dirname(V2)
sys.path.insert(0, V2); sys.path.insert(0, HERE)
import batch_v2                                   # noqa: E402
from batch_v2 import BatchV2, RHO_STAR_GB, TAU_IND_GB, TAU_EFFCA   # noqa: E402
import model_v2 as MV                             # noqa: E402
import rho_v4                                     # noqa: E402
from targets import load_targets                  # noqa: E402

NAN = float("nan")
# (name, trigger kind 0 W / 1 always / 2 Wc slow / 3 Wc fast / -1 off, theta_eCB scale, f (x c_post), d_min (nan = fit),
#  theta_d scale, veto kind 0 none / 1 own VDCC influx / 2 own c* influx / 3 own pot event, veto window Tv ms, kappa)
_B = (NAN, 1.0, 0, 0.0, 0.0)
VARIANTS = [("base", 0, 1.0, 0.0) + _B, ("noE", -1, 1.0, 0.0) + _B, ("allE", 1, 1.0, 0.0) + _B,
            ("thE.3", 0, 0.3, 0.0) + _B, ("thE.1", 0, 0.1, 0.0) + _B, ("thE.01", 0, 0.01, 0.0) + _B,
            ("Cs0.2", 2, 1.0, 0.2) + _B, ("Cs0.4", 2, 1.0, 0.4) + _B, ("Cs0.6", 2, 1.0, 0.6) + _B,
            ("Cs0.8", 2, 1.0, 0.8) + _B, ("Cf0.1", 3, 1.0, 0.1) + _B, ("Cf0.2", 3, 1.0, 0.2) + _B,
            ("Cf0.4", 3, 1.0, 0.4) + _B, ("Cf0.6", 3, 1.0, 0.6) + _B,
            ("dm.40", 0, 1.0, 0.0, -0.40, 1.0, 0, 0.0, 0.0), ("dm.50", 0, 1.0, 0.0, -0.50, 1.0, 0, 0.0, 0.0),
            ("td.9", 0, 1.0, 0.0, NAN, 0.9, 0, 0.0, 0.0), ("td.8", 0, 1.0, 0.0, NAN, 0.8, 0, 0.0, 0.0),
            ("Cf.2dm4", 3, 1.0, 0.2, -0.40, 1.0, 0, 0.0, 0.0),
            # round 2 (LTD_LIT.md 5): bAP veto. An own arrival whose trigger fired takes no eCB step if, within Tv ms
            # after it, (1) own unweighted VDCC influx sum bV s > theta_eCB (fitted), (2) own c* influx > c_pre +
            # kappa c_post, (3) an own pot event (c* > theta_p with the V gate) occurs
            ("vV25", 0, 1.0, 0.0, NAN, 1.0, 1, 25.0, 0.0), ("vP25", 0, 1.0, 0.0, NAN, 1.0, 3, 25.0, 0.0),
            ("t1vV25", 0, 0.1, 0.0, NAN, 1.0, 1, 25.0, 0.0),
            ("t01vV25", 0, 0.01, 0.0, NAN, 1.0, 1, 25.0, 0.0), ("t01vC25", 0, 0.01, 0.0, NAN, 1.0, 2, 25.0, 0.5),
            ("t01vP25", 0, 0.01, 0.0, NAN, 1.0, 3, 25.0, 0.0),
            ("aV15", 1, 1.0, 0.0, NAN, 1.0, 1, 15.0, 0.0), ("aV25", 1, 1.0, 0.0, NAN, 1.0, 1, 25.0, 0.0),
            ("aV35", 1, 1.0, 0.0, NAN, 1.0, 1, 35.0, 0.0), ("aC25", 1, 1.0, 0.0, NAN, 1.0, 2, 25.0, 0.5),
            ("aC25k1", 1, 1.0, 0.0, NAN, 1.0, 2, 25.0, 1.0), ("aP25", 1, 1.0, 0.0, NAN, 1.0, 3, 25.0, 0.0),
            ("Cf.5", 3, 1.0, 0.5) + _B, ("Cf.5vV25", 3, 1.0, 0.5, NAN, 1.0, 1, 25.0, 0.0),
            ("Cf.5vC25", 3, 1.0, 0.5, NAN, 1.0, 2, 25.0, 0.5), ("Cf.5vP25", 3, 1.0, 0.5, NAN, 1.0, 3, 25.0, 0.0)]
VN = [v[0] for v in VARIANTS]
VTK = np.array([v[1] for v in VARIANTS], np.int64)
VTE = np.array([v[2] for v in VARIANTS], np.float64)
VF = np.array([v[3] for v in VARIANTS], np.float64)
VDM = np.array([v[4] for v in VARIANTS], np.float64)
VTS = np.array([v[5] for v in VARIANTS], np.float64)
VVK = np.array([v[6] for v in VARIANTS], np.int64)
VTV = np.array([v[7] for v in VARIANTS], np.float64)
VKA = np.array([v[8] for v in VARIANTS], np.float64)
NPD = 32                                                     # pending-veto slots per variant
KEYV = ("thE.01", "allE", "Cf.5", "vV25", "t01vV25", "t01vP25", "aV25", "aC25", "aP25", "Cf.5vV25", "Cf.5vC25")


@njit(cache=False)
def replay(E, S, H, cnt, td, tp, r0, cpost, cpre, GD, GP, thV, iscV, tauE1, dmin, Ae, dp0, thE, tauD, RS, KT, tauS,
           TK, TE, FF, DM, TS, VK, TV, KA):
    """gpu_v5_rho kern mode 2 for one synapse, same step order; all variants at once (variant 0 = base).
    Veto variants: a triggered arrival is held pending for TV ms and its eCB step is applied unless vetoed.
    Returns rho_f, d0 per variant; base flags and per-synapse maxima."""
    n = E.shape[0]; nv = TK.shape[0]
    if thE <= 0.0:
        thE = thV
    r = np.full(nv, r0); d0 = np.full(nv, dp0)
    dmv = np.empty(nv)
    for v in range(nv):
        dmv[v] = dmin if DM[v] != DM[v] else DM[v]
    pdT = np.zeros((nv, NPD)); pdA = np.zeros((nv, NPD)); pdC = np.zeros((nv, NPD)); pdN = np.zeros(nv, np.int64)
    nveto = np.zeros(nv); nstep = np.zeros(nv)
    V = 0.0; W = 0.0; Bg = 0.0; Cs = 0.0; Cf = 0.0; tc = 0.0
    chV = -1.0; aV = 1.0; bV = 0.0; chD = -1.0; aD = 1.0; chS = -1.0; aS = 1.0; rCf = 1.0
    Vpk = 0.0; Wpk = 0.0; Warr = 0.0; spk = 0.0; Csa = 0.0; Cfa = 0.0; epk = 0.0
    dep_any = 0; n_arr = 0.0; n_trig = 0.0
    e = 0.0; s = 0.0; h = 0.0; hm = 0.0; c = 0.0; inc = 0.0
    for k in range(n + 1):
        last = k == n
        if not last:
            e = E[k]; s = S[k]; h = H[k]; hm = h * KT; c = cnt[k]
            if hm != chV:
                chV = hm
                aV = math.exp(-hm / tauE1); bV = tauE1 * (1.0 - aV) / iscV
            if hm != chS:
                chS = hm
                aS = math.exp(-hm / tauS)
                rCf = (tauE1 * (1.0 - aV)) / (tauS * (1.0 - aS)) if aS < 1.0 else 1.0
            inc = 0.0
            if k + 1 < n:
                inc = E[k + 1] - aS * e                              # c* influx of step k (cacr_from_effcai x b_k)
                if inc < 0.0:
                    inc = 0.0
        # resolve pending vetoes whose window closed (all at the end of the trace)
        for v in range(nv):
            if pdN[v] == 0:
                continue
            m = 0
            for j in range(pdN[v]):
                if last or tc > pdT[v, j]:
                    if VK[v] == 1:
                        vet = pdA[v, j] > thE
                    elif VK[v] == 2:
                        vet = pdA[v, j] > cpre + KA[v] * cpost
                    else:
                        vet = pdA[v, j] > 0.0
                    if vet:
                        nveto[v] += 1.0
                    else:
                        d0[v] = dmv[v] + (d0[v] - dmv[v]) * math.pow(1.0 - Ae, pdC[v, j])
                        nstep[v] += 1.0
                else:
                    pdT[v, m] = pdT[v, j]; pdA[v, m] = pdA[v, j]; pdC[v, m] = pdC[v, j]; m += 1
            pdN[v] = m
        if last:
            break
        if V > Vpk:
            Vpk = V
        if W > Wpk:
            Wpk = W
        if s > spk:
            spk = s
        if e > epk:
            epk = e
        potg = 1.0 if e > tp else 0.0
        if potg > 0.0 and thV > 0.0 and not (V > thV):
            potg = 0.0
        for v in range(nv):
            pot = potg
            dep = 1.0 if e > td * TS[v] else 0.0
            rv = r[v]
            rv = rv + h * (-rv * (1 - rv) * (RS - rv) + pot * GP * (1 - rv) - dep * (1 - pot) * GD * rv)
            r[v] = min(max(rv, 0.0), 1.0)
            if v == 0 and dep * (1 - pot) > 0.0:
                dep_any = 1
        if c > 0.0:
            Bg = Bg + c
            n_arr += c
            if W > Warr:
                Warr = W
            if Cs > Csa:
                Csa = Cs
            if Cf > Cfa:
                Cfa = Cf
            if W > thE:
                n_trig += c
            for v in range(nv):
                tk = TK[v]
                if tk == 0:
                    trig = W > thE * TE[v]
                elif tk == 1:
                    trig = True
                elif tk == 2:
                    trig = Cs > FF[v] * cpost
                elif tk == 3:
                    trig = Cf > FF[v] * cpost
                else:
                    trig = False
                if trig and Ae > 0.0:
                    if VK[v] == 0:
                        d0[v] = dmv[v] + (d0[v] - dmv[v]) * math.pow(1.0 - Ae, c)
                    elif pdN[v] < NPD:
                        j = pdN[v]; pdT[v, j] = tc + TV[v]; pdA[v, j] = 0.0; pdC[v, j] = c; pdN[v] = j + 1
        # window accumulation (the arrival step's own influx counts)
        for v in range(nv):
            for j in range(pdN[v]):
                if VK[v] == 1:
                    pdA[v, j] += bV * s
                elif VK[v] == 2:
                    pdA[v, j] += inc
                elif potg > 0.0:
                    pdA[v, j] = 1.0
        wb = 1.0 - (Bg if Bg < 1.0 else 1.0)
        W = aV * W + bV * (wb * s)
        Cs = aS * Cs + wb * inc
        Cf = aV * Cf + wb * inc * rCf
        if hm != chD:
            chD = hm
            aD = math.exp(-hm / tauD)
        Bg = Bg * aD
        V = aV * V + bV * s
        tc += hm
    return r, d0, Vpk, Wpk, Warr, spk, Csa, Cfa, epk, dep_any, n_arr, n_trig, nveto, nstep




def run_record(f, a, Q, prm, styp, dist):
    pair, proto = os.path.basename(f)[:-4].split("__")
    B = BatchV2([os.path.dirname(f)], protocols=[proto], pairs={pair}, fast=False, signals=("vdcc",), verbose=False)
    if not B.recs:
        return None
    B._stack(); r = B.recs[0]; H = B._H[0]
    t = r["t"]; n_t = len(t)
    td, tp = rho_v4.thetas(r["c_pre"], r["c_post"], a, rho_v4.opts(Q)[0])
    vd = r["vdcc"]
    thV = Q["theta_V"]; thE = Q["theta_eCB"] if Q["theta_eCB"] > 0 else thV
    rows = []; RF = []; D0 = []
    for i in range(len(r["syn"])):
        ai = np.asarray(r["arr"][i]).ravel().astype(np.int64)
        cnt = np.zeros(n_t); np.add.at(cnt, ai, 1.0)
        E = np.asarray(r["effcai"][i], np.float64)
        rf, d0, Vpk, Wpk, Warr, spk, Csa, Cfa, epk, dep_any, n_arr, n_trig, nveto, nstep = replay(
            E, np.asarray(vd[i], np.float64), H, cnt, float(td[i]), float(tp[i]), float(r["rho0"][i]),
            float(r["c_post"][i]), float(r["c_pre"][i]), *prm, VTK, VTE, VF, VDM, VTS, VVK, VTV, VKA)
        RF.append(rf.copy()); D0.append(d0.copy())
        s = int(r["syn"][i]); cp = float(r["c_post"][i])
        rows.append(dict(pair=pair, proto=proto, syn=s, sec_type=styp.get(s, np.nan), dist=dist.get(s, np.nan),
                         c_pre=float(r["c_pre"][i]), c_post=cp, td=float(td[i]), tp=float(tp[i]),
                         td_tp=float(td[i] / tp[i]), rho0=float(r["rho0"][i]), rho_f=float(rf[0]), d0=float(d0[0]),
                         ecb=bool(d0[0] < Q["dpre0"] - 1e-9), down=bool(r["rho0"][i] >= 0.5 and rf[0] < 0.5),
                         up=bool(r["rho0"][i] < 0.5 and rf[0] >= 0.5), dep_any=bool(dep_any),
                         Vpk=Vpk / thV, Wpk=Wpk / thE, W_arr=Warr / thE, vdcc_pk=spk, pk_td=epk / td[i],
                         pk_tp=epk / tp[i], Cs_arr=Csa / cp, Cf_arr=Cfa / cp, n_arr=n_arr, n_trig=n_trig,
                         **{f"ecb_{nm}": bool(d0[k] < Q["dpre0"] - 1e-9) for k, nm in enumerate(VN)},
                         **{f"nveto_{nm}": float(nveto[k]) for k, nm in enumerate(VN) if VVK[k] > 0}))
    b = B.basis(r); r0 = r["rho0"]; n = len(r0); z = np.zeros(n); dm = np.full(n, Q["dpre_min"])
    dp = np.full(n, Q["dpre0"]); RF = np.array(RF); D0 = np.array(D0)
    rec = dict(pair=pair, proto=proto, n_syn=n, rho_only=b.ratio(r0, RF[:, 0], dp), pre_only=b.ratio(r0, r0, D0[:, 0]),
               nmdar_block=b.ratio(r0, r0, dp), post_min=b.ratio(r0, z, dp), all_min=b.ratio(r0, z, dm),
               ecb_max=b.ratio(r0, r0, dm), ltp_max=b.ratio(r0, np.ones(n), dp))
    for k, nm in enumerate(VN):
        rec[f"control|{nm}"] = b.ratio(r0, RF[:, k], D0[:, k])
        rec[f"no_block|{nm}"] = rec[f"control|{nm}"]
        rec[f"mglu_block|{nm}"] = b.ratio(r0, RF[:, k], dp)
        rec[f"post_nmdar|{nm}"] = b.ratio(r0, r0, D0[:, k])
    return rows, rec


def locations(syns):
    """syn id -> section type (syn_section_type.npz, else edges afferent_section_type) and distance (L5 local_t path
    distance, else edges distance_soma if the edges file has it)."""
    import h5py
    st = np.load(os.path.join(V2, "syn_section_type.npz"))
    styp = dict(zip(st["syn"].tolist(), [float(x) for x in st["section_type"].tolist()]))
    dist = {}
    for f in glob.glob(os.path.join(V2, "local_t", "out", "*.npz")):
        z = np.load(f)
        for s, d in zip(z["og|ap1|sid"], z["og|ap1|dist"]):
            dist[int(s)] = float(d)
    miss_t = sorted(s for s in syns if s not in styp); miss_d = sorted(s for s in syns if s not in dist)
    with h5py.File(MV.EDGES, "r") as f:
        g = f[f"edges/{MV.EDGE_POP}/0"]
        print("edges fields:", sorted(g.keys()), flush=True)
        if miss_t:
            styp.update(zip(miss_t, g["afferent_section_type"][np.array(miss_t)].astype(float)))
        if miss_d and "distance_soma" in g:
            dist.update(zip(miss_d, g["distance_soma"][np.array(miss_d)].astype(float)))
    print(f"locations: {len(miss_t)} section types and {len(miss_d)} distances from edges", flush=True)
    return styp, dist


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True)
    ap.add_argument("--model", required=True, help="l5 | l23 | l23l23")
    ap.add_argument("--out", required=True)
    ap.add_argument("--models-from", default=os.path.join(HERE, "results", "v5_S1_V5cz_s5_39.json"))
    args = ap.parse_args()
    import pandas as pd
    t0 = time.time()
    fj = json.load(open(args.fit)); fa = fj["args"]
    tag = os.path.basename(args.fit)[:-5]
    batch_v2.BAP_GATE = fj.get("bap_gate")
    a = {k: fj["a"][k] for k in ("a00", "a01", "a10", "a11")}
    P = {**json.loads(fa["filters"]), **json.loads(fa["set"]), **fj["pre"]}
    Q = {**MV.DEFAULTS, "theta_eCB": 0.0, "A_eCB": 0.0, **P}
    assert int(Q["v5_mode"]) == 2, Q["v5_mode"]
    models = {m["name"]: m for m in fj.get("models", [])}
    if args.model not in models:
        models.update({m["name"]: m for m in json.load(open(args.models_from))["models"] if m["name"] not in models})
    M = models[args.model]
    batch_v2.BASIS_DIR = M["basis"] if os.path.isabs(M["basis"]) else os.path.join(ROOT, M["basis"])
    conds = set(fa["conditions"].split(","))
    drops = [d for d in (fa.get("drop_targets") or "").split(",") if d]
    drop = {d.split("/", 1)[1] if "/" in d else d for d in drops if "/" not in d or d.split("/", 1)[0] == args.model}
    T = {k: v for k, v in load_targets(tuple(M["groups"].split(","))).items()
         if k[1] in conds and f"{k[0]}|{k[1]}" not in drop}
    protos = sorted({k[0].split("@")[0] for k in T})
    loc = None; pairs = None
    if args.model == "l5":
        pairs = set(fa["pairs"].split(","))
    if M.get("geom"):
        g = pd.read_csv(M["geom"]); g["pair"] = g.pregid.astype(str) + "-" + g.postgid.astype(str); g = g[g.all_protocols]
        lc = [c for c in g.columns if c.endswith("_distal")]
        loc = {p: {c: bool(v) for c, v in zip(lc, row)} for p, row in zip(g.pair, g[lc].itertuples(index=False))}
        pairs = set(g.pair)
    gd, gp = rho_v4.rates(Q)
    prm = (gd, gp, float(Q["theta_V"]), float(Q["i_scale"]), float(Q["tau_E1"]), float(Q["dpre_min"]),
           float(Q["A_eCB"]), float(Q["dpre0"]), float(Q["theta_eCB"]), float(Q["tau_d_NMDA"]), float(RHO_STAR_GB),
           1e3 * TAU_IND_GB, float(TAU_EFFCA))
    print(f"{tag} {args.model}: {len(T)} targets, theta_V {Q['theta_V']:.4g} theta_eCB {Q['theta_eCB']:.4g} "
          f"d_min {Q['dpre_min']} dpre0 {Q['dpre0']} A_eCB {Q['A_eCB']} tau_E1 {Q['tau_E1']} tau* {TAU_EFFCA} "
          f"gd {gd:.1f} gp {gp:.1f} a {a}", flush=True)
    files = []
    for d in M["dirs"].split(","):
        d = d if os.path.isabs(d) else os.path.join(ROOT, d)
        for f in sorted(glob.glob(os.path.join(d, "*.npz"))):
            pr, pt = os.path.basename(f)[:-4].split("__")
            if pt in protos and (pairs is None or pr in pairs):
                files.append(f)
    print(f"{len(files)} record files", flush=True)
    syns = set()
    for f in files:
        with np.load(f) as z:
            syns.update(int(s) for s in z["syn"])
    styp, dist = locations(syns)
    rows, recs = [], []
    for n, f in enumerate(files):
        res = run_record(f, a, Q, prm, styp, dist)
        if res is None:
            continue
        rows += res[0]; recs.append(res[1])
        if n % 200 == 0:
            print(f"  {n}/{len(files)} records, {time.time() - t0:.0f} s", flush=True)
    os.makedirs(args.out, exist_ok=True)
    pre = os.path.join(args.out, f"{tag}_{args.model}_ltd")
    S = pd.DataFrame(rows); R = pd.DataFrame(recs)
    S.to_csv(pre + "_syn.csv", index=False)
    csv = args.fit[:-5] + ("" if args.model == "l5" else f"_{args.model}") + ".csv"
    fitp = {}
    if os.path.isfile(csv):
        c = pd.read_csv(csv); fitp = {(x.target, x.condition): x.pred for x in c.itertuples()}
    out, mech = [], []
    for (pid, cond), v in T.items():
        proto, _, where = pid.partition("@")
        x = R[R.proto == proto]
        if loc is not None:
            keep = []
            for p in x.pair:
                L = loc.get(p); ok = L is not None
                if ok and where and proto.startswith("letzkus") and L["letzkus_distal"] != (where == "distal"):
                    ok = False
                if ok and where == "distal" and not proto.startswith("letzkus") and not L["sh_distal"]:
                    ok = False
                keep.append(ok)
            x = x[np.array(keep, bool)] if len(x) else x
        elif where == "distal":
            continue
        q = dict(target=pid, condition=cond, data=v[0], sem=v[1], n_rec=len(x), fit_pred=fitp.get((pid, cond), np.nan))
        for nm in VN:
            col = "nmdar_block" if cond == "nmdar_block" else f"{cond}|{nm}"
            q[nm] = float(np.nanmean(x[col])) if len(x) else np.nan
        for k in ("rho_only", "pre_only", "post_min", "all_min", "ecb_max", "ltp_max"):
            q[k] = float(np.nanmean(x[k])) if len(x) else np.nan
        out.append(q)
        if cond == "control" and (v[0] < 1.0 or abs(v[0] - 1.0) <= 0.1):
            y = S[S.proto == proto]
            y = y[y.pair.isin(set(x.pair))]
            for grp, yy in [("all", y)] + [(f"sec{int(s)}", yy2) for s, yy2 in y.groupby("sec_type")] + \
                    ([(f"d{lo}-{hi}", y[(y.dist >= lo) & (y.dist < hi)]) for lo, hi in
                      ((0, 75), (75, 150), (150, 250), (250, 2000))] if y.dist.notna().any() else []):
                if not len(yy):
                    continue
                w = yy.W_arr.values
                mech.append(dict(target=pid, group=grp, n_syn=len(yy), f_rho0_1=float(np.mean(yy.rho0 >= .5)),
                                 f_ecb=float(yy.ecb.mean()), f_down=float(yy.down.mean()),
                                 f_down_of1=float(yy.down.sum() / max((yy.rho0 >= .5).sum(), 1)),
                                 f_dep_any=float(yy.dep_any.mean()), f_up=float(yy.up.mean()),
                                 f_none=float(np.mean(~yy.ecb & ~yy.down & ~yy.up)),
                                 W_lt001=float(np.mean(w < .01)), W_01_1=float(np.mean((w >= .01) & (w < 1))),
                                 W_1_10=float(np.mean((w >= 1) & (w < 10))), W_ge10=float(np.mean(w >= 10)),
                                 W_q50=float(np.median(w)), vdcc_q50=float(yy.vdcc_pk.median()),
                                 pk_td_q50=float(yy.pk_td.median()), f_pk_td_gt1=float(np.mean(yy.pk_td > 1)),
                                 td_tp_q50=float(yy.td_tp.median()), Cs_arr_q50=float(yy.Cs_arr.median()),
                                 Cf_arr_q50=float(yy.Cf_arr.median()), dist_q50=float(yy.dist.median()),
                                 **{f"fe_{nm}": float(yy[f"ecb_{nm}"].mean()) for nm in KEYV}))
    D = pd.DataFrame(out); MM = pd.DataFrame(mech)
    D.to_csv(pre + "_targets.csv", index=False); MM.to_csv(pre + "_mech.csv", index=False)
    pd.set_option("display.width", 320); pd.set_option("display.max_columns", 80); pd.set_option("display.max_rows", 500)
    print(f"\n=== {tag} {args.model}: targets x variants (record means)")
    print(D.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    if np.isfinite(D.fit_pred).any():
        print(f"repro max |base - fit_pred| = {np.nanmax(np.abs(D.base - D.fit_pred)):.2e}")
    ltd = D.data < 1.0
    print("chi2 per variant (all | data<1 | data>=1):")
    for nm in VN:
        z2 = ((D[nm] - D["data"]) / D["sem"]) ** 2
        print(f"  {nm:8s} {np.nansum(z2):8.2f} {np.nansum(z2[ltd]):8.2f} {np.nansum(z2[~ltd]):8.2f}")
    print(f"\n=== {tag} {args.model}: per-synapse mechanism, control LTD / no-change targets (base)")
    print(MM.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    if S.dist.notna().any():
        ok = S.dist.notna() & (S.vdcc_pk > 0)
        print(f"Spearman(dist, log vdcc_pk) = {pd.Series(S.dist[ok]).corr(np.log(S.vdcc_pk[ok]), method='spearman'):.3f}"
              f" over {ok.sum()} synapse-records")
    print(f"saved {pre}_{{syn,targets,mech}}.csv ({len(S)} synapse-records, {len(R)} records), {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
