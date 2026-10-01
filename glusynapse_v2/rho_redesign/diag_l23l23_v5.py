"""v5c diagnosis of the L2/3 -> L2/3 (Zilberter 2009) misfit (L23L23_DIAG.md section 7). Read-only.

Replays the gpu_v5_rho kernel (v5_mode 2 = v5b/v5c; v5_mode 3 = E2) on the CPU one record at a time (one BatchV2 per
npz, as diag_l23l23.py), for every record of ONE pathway model (--model l23l23 | l5 | l23) at the parameters of ONE v5
fit json (--fit), and in the same pass counterfactual rule variants with the SAME fixed parameters (no refit; every
variant is uniform, synapse-local and reads only own spine Ca quantities: c*, the own VDCC pool V, its (1 - b)-weighted
copy W, own arrivals):

  base   v5c as fitted (v5_mode of the json; must reproduce the fit csv pred)
  E2     v5_mode 3: an own VDCC excursion feeds W only if it began with the pool unprimed (V <= theta_V)
  noE    eCB step off (A_eCB 0): the post-only bound
  DW     post LTD: dep also while W > theta_eCB (the eCB trigger state drives Chindemi depression; AM251-insensitive)
  DV     as DW with V > theta_eCB
  PW     post LTD jump: at an own arrival with W > theta_eCB, rho -> 0 (and the fitted eCB step)
  PWo    PW without the eCB step (all trigger-LTD postsynaptic)
  GW     potentiation gate reads W > theta_V instead of V > theta_V (glutamate-unweighted VDCC Ca gates LTP)
  GWE    potentiation gate reads W > theta_eCB
  td.8, td.6, td.5   theta_d scaled x0.8 / 0.6 / 0.5 (a00, a01 direction; a Chindemi-parameter move, scored unrefit)
  E2+DW  combination

Per synapse: c_pre, c_post, td, tp, rho0, rho_f; pk_tp, pk_mM (pairing peak c*); t_tp, t_band, t_pot, t_conv, t_dep
(ms per pairing); Vpk, V_arr (theta_V units), VE_arr (V / theta_eCB at own arrivals), Wpk, W_arr (theta_eCB units),
trig (fraction of own arrivals with W > theta_eCB), Bg_arr; d0, d1; absolute spine Ca above rest (cai_CR - min, uM, from
batch_v2.cacr_from_effcai): bAP_uM (isolated post AP: no post AP in the previous 100 ms, no own arrival in [-150, +30] ms;
peak in [0, 15] ms), EPSP_uM / EPSP_max (isolated own arrival: no post AP in [-150, +50] ms, no own arrival in the
previous 100 ms; peak in [0, 50] ms; median / max over events), pair_uM (pairing-window peak, median).
Per record: ratios per variant x condition, decomposition (rho_only, pre_only), and ceilings: ltp_max (rho -> 1, d 0),
post_min (rho -> 0, d 0), all_min (rho -> 0, d_min), ecb_max (rho0, d_min).
Outputs (--out, scratch): <tag>_<model>_{syn,rec,targets,proto}.csv; log = the job log.

    python glusynapse_v2/rho_redesign/diag_l23l23_v5.py --fit <v5 json> --model l23l23 --out /scratch/dhuruva/l23l23_diag_v5
"""
import argparse, glob, json, math, os, sys, time
import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__)); V2 = os.path.dirname(HERE); ROOT = os.path.dirname(V2)
sys.path.insert(0, V2); sys.path.insert(0, HERE)
import batch_v2                                   # noqa: E402
from batch_v2 import BatchV2, RHO_STAR_GB, TAU_IND_GB, cacr_from_effcai   # noqa: E402
import model_v2 as MV                             # noqa: E402
import rho_v4                                     # noqa: E402
from targets import load_targets                  # noqa: E402

# (name, E2 trigger, eCB on, dep kind 0 / 1 W>thE / 2 V>thE, post jump 0 / 1 PW / 2 PWo, pot gate 0 V / 1 W>thV /
#  2 W>thE, theta_d scale)
VARIANTS = [("base", 0, 1, 0, 0, 0, 1.0), ("E2", 1, 1, 0, 0, 0, 1.0), ("noE", 0, 0, 0, 0, 0, 1.0),
            ("DW", 0, 1, 1, 0, 0, 1.0), ("DV", 0, 1, 2, 0, 0, 1.0), ("PW", 0, 1, 0, 1, 0, 1.0),
            ("PWo", 0, 0, 0, 2, 0, 1.0), ("GW", 0, 1, 0, 0, 1, 1.0), ("GWE", 0, 1, 0, 0, 2, 1.0),
            ("td.8", 0, 1, 0, 0, 0, 0.8), ("td.6", 0, 1, 0, 0, 0, 0.6), ("td.5", 0, 1, 0, 0, 0, 0.5),
            ("E2+DW", 1, 1, 1, 0, 0, 1.0)]
VN = [v[0] for v in VARIANTS]
VE2 = np.array([v[1] for v in VARIANTS], np.int64)
VAE = np.array([v[2] for v in VARIANTS], np.int64)
VDK = np.array([v[3] for v in VARIANTS], np.int64)
VPK = np.array([v[4] for v in VARIANTS], np.int64)
VGK = np.array([v[5] for v in VARIANTS], np.int64)
VTS = np.array([v[6] for v in VARIANTS], np.float64)


@njit(cache=False)
def replay(E, S, H, cnt, evk, td, tp, r0, GD, GP, thV, iscV, tauE1, dmin, Ae, dp0, thE, tauD, RS, KT,
           E2v, AEv, DKv, PKv, GKv, TSv):
    """gpu_v5_rho kern (mode 2 / 3) for one synapse, same step order, all variants at once (variant 0 = base).
    evk: 1 at the grid step of an own VDCC event onset (cev, unique steps). Returns rho_f, d0 per variant and the
    shared V / Bg traces, base W trace and base pot / dep / conversion flags."""
    n = E.shape[0]; nv = E2v.shape[0]
    if thE <= 0.0:
        thE = thV
    Vt = np.zeros(n); Wt = np.zeros(n); Bt = np.zeros(n)
    Pt = np.zeros(n, np.int8); Dt = np.zeros(n, np.int8); Ct = np.zeros(n, np.int8)
    r = np.full(nv, r0); d0 = np.full(nv, dp0); W = np.zeros(nv)
    V = 0.0; Bg = 0.0; chV = -1.0; aV = 1.0; bV = 0.0; chD = -1.0; aD = 1.0; gE = 1.0
    for k in range(n):
        e = E[k]; s = S[k]; h = H[k]; hm = h * KT; c = cnt[k]
        Vt[k] = V; Wt[k] = W[0]; Bt[k] = Bg
        for v in range(nv):
            pot = 1.0 if e > tp else 0.0
            dep = 1.0 if e > td * TSv[v] else 0.0
            raw = pot
            if pot > 0.0 and thV > 0.0:
                if GKv[v] == 0 and not (V > thV):
                    pot = 0.0
                elif GKv[v] == 1 and not (W[v] > thV):
                    pot = 0.0
                elif GKv[v] == 2 and not (W[v] > thE):
                    pot = 0.0
            if DKv[v] == 1 and W[v] > thE:
                dep = 1.0
            elif DKv[v] == 2 and V > thE:
                dep = 1.0
            rv = r[v]
            rv = rv + h * (-rv * (1 - rv) * (RS - rv) + pot * GP * (1 - rv) - dep * (1 - pot) * GD * rv)
            r[v] = min(max(rv, 0.0), 1.0)
            if v == 0:
                Pt[k] = int(pot); Dt[k] = int(dep * (1 - pot)); Ct[k] = int(raw > 0.0 and pot == 0.0)
        if c > 0.0:
            Bg = Bg + c
        for v in range(nv):
            if c > 0.0 and W[v] > thE:
                if AEv[v] == 1 and Ae > 0.0:
                    d0[v] = dmin + (d0[v] - dmin) * math.pow(1.0 - Ae, c)
                if PKv[v] >= 1:
                    r[v] = 0.0
        if hm != chV:
            chV = hm
            aV = math.exp(-hm / tauE1); bV = tauE1 * (1.0 - aV) / iscV
        wb = 1.0 - (Bg if Bg < 1.0 else 1.0)
        if evk[k] > 0.0:
            gE = 0.0 if V > thV else 1.0
        for v in range(nv):
            g = gE if E2v[v] == 1 else 1.0
            W[v] = aV * W[v] + bV * (wb * g * s)
        if hm != chD:
            chD = hm
            aD = math.exp(-hm / tauD)
        Bg = Bg * aD
        V = aV * V + bV * s
    return r, d0, Vt, Wt, Bt, Pt, Dt, Ct


def windows(r):
    ev = np.sort(np.concatenate([np.asarray(r["prespikes"], float).ravel(), np.asarray(r["postspikes"], float).ravel()]))
    if not len(ev):
        return []
    gs = np.split(ev, np.where(np.diff(ev) > 300.0)[0] + 1)
    out = []
    for j, g in enumerate(gs):
        hi = g[-1] + 500.0
        if j + 1 < len(gs):
            hi = min(hi, gs[j + 1][0] - 30.0)
        out.append((g[0] - 30.0, hi))
    return out


def med(x):
    x = np.asarray(x, float)
    return float(np.nanmedian(x)) if len(x) and np.any(np.isfinite(x)) else np.nan


def peaks(u, t, t0s, lo, hi):
    out = []
    for t0 in t0s:
        a, b = np.searchsorted(t, [t0 + lo, t0 + hi])
        if b > a:
            out.append(float(u[a:b].max()))
    return np.array(out)


def run_record(f, a, Q, prm):
    pair, proto = os.path.basename(f)[:-4].split("__")
    B = BatchV2([os.path.dirname(f)], protocols=[proto], pairs={pair}, fast=False, signals=("vdcc",), verbose=False)
    if not B.recs:
        return None
    B._stack(); r = B.recs[0]; H = B._H[0]
    t = r["t"]; n_t = len(t); h = np.diff(t)
    td, tp = rho_v4.thetas(r["c_pre"], r["c_post"], a, rho_v4.opts(Q)[0])
    vd = r["vdcc"]; W = windows(r); npair = max(len(W), 1)
    wsl = [slice(*np.searchsorted(t, [lo, hi])) for lo, hi in W]
    wsl = [w for w in wsl if w.stop > w.start]
    thV = Q["theta_V"]; thE = Q["theta_eCB"] if Q["theta_eCB"] > 0 else thV
    post = np.sort(np.asarray(r["postspikes"], float).ravel())
    rows = []; RF = []; D0 = []
    for i in range(len(r["syn"])):
        ai = np.asarray(r["arr"][i]).ravel().astype(np.int64)
        cnt = np.zeros(n_t); np.add.at(cnt, ai, 1.0)
        ev = np.asarray(r["cev"][i], np.float64); ev = ev[np.isfinite(ev)]
        evk = np.zeros(n_t)
        if len(ev):
            evk[np.unique(np.clip(np.searchsorted(t, ev - 1e-9), 0, n_t - 1))] = 1.0
        E = np.asarray(r["effcai"][i], np.float64)
        rf, d0, Vt, Wt, Bt, Pt, Dt, Ct = replay(E, np.asarray(vd[i], np.float64), H, cnt, evk, float(td[i]),
                                                float(tp[i]), float(r["rho0"][i]), *prm, VE2, VAE, VDK, VPK, VGK, VTS)
        RF.append(rf.copy()); D0.append(d0.copy())
        ta = t[ai] if len(ai) else np.zeros(0)
        # isolated bAPs / EPSPs for the absolute spine Ca
        iso_b = [x for x in post if not np.any((post < x) & (post > x - 100.0))
                 and not np.any((ta > x - 150.0) & (ta < x + 30.0))]
        iso_e = [x for x in ta if not np.any((post > x - 150.0) & (post < x + 50.0))
                 and not np.any((ta < x) & (ta > x - 100.0))]
        u = cacr_from_effcai(E[None, :], t)[0].astype(np.float64) * 1e3     # spine cai_CR above rest, uM
        pb = peaks(u, t, iso_b, 0.0, 15.0); pe = peaks(u, t, iso_e, 0.0, 50.0)
        hp = lambda m: float(np.sum(h[m[:-1]])) / npair
        q = dict(pair=pair, proto=proto, syn=int(r["syn"][i]), npair=len(W),
                 c_pre=float(r["c_pre"][i]), c_post=float(r["c_post"][i]), td=float(td[i]), tp=float(tp[i]),
                 rho0=float(r["rho0"][i]), rho_f=float(rf[0]),
                 pk_tp=med([E[w].max() / tp[i] for w in wsl]), pk_mM=med([E[w].max() for w in wsl]),
                 t_tp=hp(E > tp[i]), t_band=hp((E > td[i]) & (E <= tp[i])),
                 t_pot=float(np.sum(h * Pt[:-1])) / npair, t_conv=float(np.sum(h * Ct[:-1])) / npair,
                 t_dep=float(np.sum(h * Dt[:-1])) / npair,
                 Vpk=med([Vt[w].max() / thV for w in wsl]), V_arr=med(Vt[ai] / thV) if len(ai) else np.nan,
                 VE_arr=med(Vt[ai] / thE) if len(ai) else np.nan,
                 Wpk=med([Wt[w].max() / thE for w in wsl]), W_arr=med(Wt[ai] / thE) if len(ai) else np.nan,
                 trig=float(np.mean(Wt[ai] > thE)) if len(ai) else np.nan,
                 Bg_arr=med(Bt[ai]) if len(ai) else np.nan, t_V=hp(Vt > thV), t_W=hp(Wt > thE),
                 n_ev=len(ev) / npair, d0=float(d0[0]), d1=float(Q["dpre0"]),
                 bAP_uM=med(pb), bAP_n=len(pb), EPSP_uM=med(pe), EPSP_max=float(pe.max()) if len(pe) else np.nan,
                 EPSP_n=len(pe), pair_uM=med([u[w].max() for w in wsl]))
        for k, nm in enumerate(VN[1:], 1):
            q[f"rf_{nm}"] = float(rf[k]); q[f"d0_{nm}"] = float(d0[k])
        rows.append(q)
    b = B.basis(r); r0 = r["rho0"]; n = len(r0); z = np.zeros(n); dm = np.full(n, Q["dpre_min"])
    dp = np.full(n, Q["dpre0"])
    RF = np.array(RF); D0 = np.array(D0)
    rec = dict(pair=pair, proto=proto, n_syn=n, rho_only=b.ratio(r0, RF[:, 0], z), pre_only=b.ratio(r0, r0, D0[:, 0]),
               nmdar_block=b.ratio(r0, r0, z), ltp_max=b.ratio(r0, np.ones(n), z), post_min=b.ratio(r0, z, z),
               all_min=b.ratio(r0, z, dm), ecb_max=b.ratio(r0, r0, dm))
    for k, nm in enumerate(VN):
        rec[f"control|{nm}"] = b.ratio(r0, RF[:, k], D0[:, k])
        rec[f"no_block|{nm}"] = rec[f"control|{nm}"]
        rec[f"mglu_block|{nm}"] = b.ratio(r0, RF[:, k], dp)
        rec[f"post_nmdar|{nm}"] = b.ratio(r0, r0, D0[:, k])
    return rows, rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True)
    ap.add_argument("--model", required=True, help="l5 | l23 | l23l23")
    ap.add_argument("--out", required=True)
    ap.add_argument("--models-from", default=os.path.join(HERE, "results", "v5_V5cz_s5.json"),
                    help="json whose 'models' supplies a model spec missing from --fit (V5c_s7 has no l23l23)")
    ap.add_argument("--csv", default=None)
    args = ap.parse_args()
    import pandas as pd
    t0 = time.time()
    fj = json.load(open(args.fit)); fa = fj["args"]
    tag = os.path.basename(args.fit)[:-5]
    batch_v2.BAP_GATE = fj.get("bap_gate")
    a = {k: fj["a"][k] for k in ("a00", "a01", "a10", "a11")}
    P = {**json.loads(fa["filters"]), **json.loads(fa["set"]), **fj["pre"]}
    Q = {**MV.DEFAULTS, "theta_eCB": 0.0, "A_eCB": 0.0, **P}
    mode = int(Q["v5_mode"])
    assert mode in (2, 3), mode
    if mode == 3:                                  # fitted with the E2 trigger: base = E2
        VE2[0] = 1
    models = {m["name"]: m for m in fj.get("models", [])}
    if args.model not in models:
        models.update({m["name"]: m for m in json.load(open(args.models_from))["models"] if m["name"] not in models})
        print(f"model spec {args.model} from {args.models_from}", flush=True)
    M = models[args.model]
    basis = M["basis"] if os.path.isabs(M["basis"]) else os.path.join(ROOT, M["basis"])
    batch_v2.BASIS_DIR = basis
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
           1e3 * TAU_IND_GB)
    print(f"{tag} model {args.model}: {len(T)} targets, {len(protos)} protocols, basis {basis}, v5_mode {mode}", flush=True)
    print("params:", {k: round(float(Q[k]), 5) for k in ("theta_V", "theta_eCB", "A_eCB", "dpre_min", "dpre0", "tau_E1",
                                                            "i_scale", "tau_d_NMDA", "rho_gamma")},
          f"gamma_d {gd:.2f} gamma_p {gp:.2f}", {k: round(v, 4) for k, v in a.items()}, flush=True)
    files = []
    for d in M["dirs"].split(","):
        d = d if os.path.isabs(d) else os.path.join(ROOT, d)
        for f in sorted(glob.glob(os.path.join(d, "*.npz"))):
            pr, pt = os.path.basename(f)[:-4].split("__")
            if pt in protos and (pairs is None or pr in pairs):
                files.append(f)
    print(f"{len(files)} record files", flush=True)
    rows, recs = [], []
    for n, f in enumerate(files):
        res = run_record(f, a, Q, prm)
        if res is None:
            continue
        rows += res[0]; recs.append(res[1])
        if n % 200 == 0:
            print(f"  {n}/{len(files)} records, {time.time() - t0:.0f} s", flush=True)
    os.makedirs(args.out, exist_ok=True)
    pre = os.path.join(args.out, f"{tag}_{args.model}")
    S = pd.DataFrame(rows); R = pd.DataFrame(recs)
    # supralinearity: pairing peak / (own isolated bAP + own isolated EPSP), references pooled per synapse over records
    ref = S.groupby(["pair", "syn"]).agg(bAP_ref=("bAP_uM", "median"), EPSP_ref=("EPSP_uM", "median"),
                                         EPSPmax_ref=("EPSP_max", "max")).reset_index()
    S = S.merge(ref, on=["pair", "syn"], how="left")
    S["supra"] = S.pair_uM / (S.bAP_ref + S.EPSP_ref)
    S["bAP_EPSP"] = S.bAP_ref / S.EPSPmax_ref
    S["up"] = (S.rho0 < 0.5) & (S.rho_f >= 0.5); S["down"] = (S.rho0 >= 0.5) & (S.rho_f < 0.5)
    S.to_csv(pre + "_syn.csv", index=False); R.to_csv(pre + "_rec.csv", index=False)
    csv = args.csv or (args.fit[:-5] + ("" if args.model == "l5" else f"_{args.model}") + ".csv")
    fitp = {}
    if os.path.isfile(csv):
        c = pd.read_csv(csv); fitp = {(x.target, x.condition): x.pred for x in c.itertuples()}
    out = []
    for (pid, cond), v in T.items():
        proto, _, where = pid.partition("@")
        x = R[R.proto == proto]
        if loc is not None:
            keep = []
            for p in x.pair:
                L = loc.get(p)
                ok = L is not None
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
            q[nm] = float(np.nanmean(x[col])) if len(x) and col in x else np.nan
        for k in ("rho_only", "pre_only", "ltp_max", "post_min", "all_min", "ecb_max"):
            q[k] = float(np.nanmean(x[k])) if len(x) else np.nan
        out.append(q)
    D = pd.DataFrame(out)
    D["z_base"] = (D.base - D["data"]) / D["sem"]
    D["reach"] = np.where(D["data"] < np.where(D.condition == "mglu_block", D.post_min, D.all_min), "below_floor",
                          np.where(D["data"] > D.ltp_max, "above_ceil", "ok"))
    D.to_csv(pre + "_targets.csv", index=False)
    pd.set_option("display.width", 300); pd.set_option("display.max_columns", 80); pd.set_option("display.max_rows", 500)
    print(f"\n=== {tag} {args.model}: targets (pred = record mean; base must equal fit_pred when a fit csv exists)")
    print(D.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    if np.isfinite(D.fit_pred).any():
        print(f"repro max |base - fit_pred| = {np.nanmax(np.abs(D.base - D.fit_pred)):.2e}")
    chi = {nm: float(np.nansum(((D[nm] - D["data"]) / D["sem"]) ** 2)) for nm in VN}
    print("chi2 per variant: " + ", ".join(f"{k} {v:.2f}" for k, v in chi.items())
          + f"  (fit csv {np.nansum(((D.fit_pred - D['data']) / D['sem']) ** 2):.2f})")
    cols = ["c_pre", "c_post", "td", "tp", "pk_tp", "pk_mM", "t_tp", "t_band", "t_pot", "t_conv", "t_dep", "Vpk", "V_arr",
            "VE_arr", "Wpk", "W_arr", "trig", "Bg_arr", "t_V", "t_W", "n_ev", "d0", "bAP_uM", "EPSP_uM", "EPSP_max",
            "pair_uM", "bAP_ref", "EPSP_ref", "supra", "bAP_EPSP"]
    G = S.groupby("proto")
    Pq = G[cols].median()
    Pq["n_syn"] = G.size(); Pq["frac_rho0_1"] = G.rho0.apply(lambda z: float(np.mean(z >= 0.5)))
    Pq["P_up"] = G.apply(lambda z: z.up.sum() / max((z.rho0 < 0.5).sum(), 1))
    Pq["P_down"] = G.apply(lambda z: z.down.sum() / max((z.rho0 >= 0.5).sum(), 1))
    for nm in VN[1:]:
        Pq[f"Pup_{nm}"] = G.apply(lambda z, nm=nm: float(np.sum((z.rho0 < .5) & (z[f"rf_{nm}"] >= .5))) / max((z.rho0 < .5).sum(), 1))
        Pq[f"Pdn_{nm}"] = G.apply(lambda z, nm=nm: float(np.sum((z.rho0 >= .5) & (z[f"rf_{nm}"] < .5))) / max((z.rho0 >= .5).sum(), 1))
    Pq = Pq.reset_index(); Pq.to_csv(pre + "_proto.csv", index=False)
    print(f"\n=== {tag} {args.model}: per-protocol medians (ms per pairing; Vpk, V_arr in theta_V; VE_arr, Wpk, W_arr in "
          f"theta_eCB; Ca in uM above rest)")
    print(Pq[["proto", "n_syn"] + cols + ["frac_rho0_1", "P_up", "P_down"]].to_string(
        index=False, float_format=lambda v: f"{v:.3g}"))
    print("\nP_up / P_down per variant:")
    print(Pq[["proto"] + [c for c in Pq.columns if c.startswith(("Pup_", "Pdn_"))]].to_string(
        index=False, float_format=lambda v: f"{v:.2f}"))
    print(f"\nsaved {pre}_{{syn,rec,targets,proto}}.csv ({len(S)} synapses, {len(R)} records), {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
