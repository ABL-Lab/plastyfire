"""L2/3 -> L2/3 (Zilberter 2009) diagnosis of the 3-pathway score of C1Ajn_s5 (L23L23_DIAG.md). Read-only.

Replays the gpu_v4_rho kernel (vamp 1 = C1, t_drive 4, pre_drive 1) on the CPU one record at a time, exactly as
diag_dltd.py (which reproduced its fit csv to 4 decimals), for every record of ONE pathway model of the fit json
(--model l23l23 | l5 | l23), and in the same pass a few counterfactual rule variants with the SAME fixed parameters
(no refit, no new parameter: they reuse theta_V, tau_E1, tau_d_NMDA and the 0.5 count threshold of the (D) vgate):

  base       C1 as fitted (must reproduce the fit csv pred)
  C2         a blocked potentiation step is neutral (vamp_mode 2, existing kernel switch)
  Gb         pot also needs Gb > 0.5, Gb = own VDCC events weighted by b (own glutamate-bound NMDA state at the event,
             the complement of the eCB weight 1 - b), low-passed with tau_E1: "bAP while own glutamate is bound"
  Gbp        as Gb, but an event only counts if the own c_VDCC pool is already above theta_V (primed pool, Zilberter:
             the pre spike must come late in the train)
  E1         eCB step skipped at an own arrival while c_VDCC > theta_V (high-VDCC pool -> no eCB)
  E2         an own VDCC event only drives the eCB pool S1 while c_VDCC <= theta_V (eCB from unprimed events only)
  GbL        latched Gb: opens at Gb > 0.5 while c* > theta_p, stays open until c* <= theta_p (no C1 conversion
             of the c* tail after Gb decays)
  Gbp+E2, C2+Gbp, C2+Gbp+E2, GbL+E2   combinations

Per synapse (base rule; times per pairing in ms, pairing = pre arrivals + post APs closer than 300 ms):
  c_pre, c_post, td, tp, rho0, rho_f; pk_tp = pairing peak c*/theta_p; t_tp (c* > theta_p), t_band (theta_d < c* <=
  theta_p), t_pot (pot applied), t_conv (c* > theta_p but c_VDCC gate closed -> depresses under C1), t_dep (dep applied);
  Vpk = pairing peak c_VDCC / theta_V, V_arr = c_VDCC at own arrival / theta_V, t_V (c_VDCC > theta_V);
  n_ev (own VDCC events), n_ev_bef (before the pairing's first own arrival), w_ecb = sum (1 - b), w_b = sum b;
  Spk, Tpk (eCB pool / gate integrator peaks), xg_arr (T - theta_Tg at arrivals, median), open_arr (frac > 0);
  d0 (control), d1 (mglu_block), d2 (no_block) dpre.
Per record: ratio per variant x condition (control, mglu_block, no_block, nmdar_block), base decomposition rho_only
(rho_f, d 0), pre_only (rho0, d0), ecb_only (rho0, d2). Per target: pred (record mean, as gpu_v3/gpu_v4 chi2),
z, chi2 per variant. Outputs (--out, scratch): <model>_{syn,rec,targets,proto}.csv; log = the job log.

    python glusynapse_v2/rho_redesign/diag_l23l23.py --fit <json> --model l23l23 --out /scratch/dhuruva/l23l23_diag
"""
import argparse, glob, json, math, os, sys, time
import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__)); V2 = os.path.dirname(HERE); ROOT = os.path.dirname(V2)
sys.path.insert(0, V2); sys.path.insert(0, HERE)
import batch_v2                                   # noqa: E402
from batch_v2 import BatchV2, ca_impulses, RHO_STAR_GB, TAU_IND_GB   # noqa: E402
import model_v2 as MV                             # noqa: E402
import rho_v4                                     # noqa: E402
from targets import load_targets                  # noqa: E402

# (name, vamp_mode, gate kind 0 none / 1 Gb / 2 Gbp / 3 GbL latch, eCB kind 0 fitted / 1 E1 / 2 E2)
VARIANTS = [("base", 1, 0, 0), ("C2", 2, 0, 0), ("Gb", 1, 1, 0), ("Gbp", 1, 2, 0), ("E1", 1, 0, 1), ("E2", 1, 0, 2),
            ("Gbp+E2", 1, 2, 2), ("C2+Gbp", 2, 2, 0), ("C2+Gbp+E2", 2, 2, 2), ("GbL", 1, 3, 0), ("GbL+E2", 1, 3, 2)]
VN = [v[0] for v in VARIANTS]
VM = np.array([v[1] for v in VARIANTS], np.int64)
GK = np.array([v[2] for v in VARIANTS], np.int64)
EK = np.array([v[3] for v in VARIANTS], np.int64)
CONDS = ("control", "mglu_block", "no_block", "nmdar_block")


@njit(cache=False)
def replay(E, S, H, cnt, jmp, jev, td, tp, r0, dmin, dmax, am, r_no, thTg, thN, thZ, offT, scT, offN, sc, tauT, tauN,
           tauZ, tauE1, GD, GP, thV, iscV, RS, KT, VMv, GKv, EKv):
    """gpu_v4_rho kern for one synapse (VAMP 1 / 2, T_IMP, no_mode 0, use_s), same step order as diag_dltd.replay,
    for all variants at once (variant 0 = base). jmp: (1 - b) weighted own VDCC events per step (ca_impulses),
    jev: unweighted own VDCC events per step. Returns per-variant rho_f, d0, d1, d2 and base traces."""
    n = E.shape[0]; nv = VMv.shape[0]
    Vt = np.zeros(n); Tt = np.zeros(n); St = np.zeros(n)
    Pt = np.zeros(n, np.int8); Dt = np.zeros(n, np.int8); Ct = np.zeros(n, np.int8)
    xgA = np.full(n, np.nan)
    r = np.full(nv, r0); T = np.zeros(nv); S1 = np.zeros(nv); G = np.zeros(nv); La = np.zeros(nv)
    d0 = np.zeros(nv); d1 = np.zeros(nv); d2 = np.zeros(nv)
    N = 0.0; Z = 0.0; V = 0.0; hmp = 0.0
    ch = -1.0; aT = 1.0; bT = 0.0; aN = 1.0; bN = 0.0; aZ = 1.0; aV = 1.0
    chE = -1.0; aE = 1.0
    for k in range(n):
        e = E[k]; s = S[k]; h = H[k]; hm = h * KT; c = cnt[k]
        jb = jev[k] - jmp[k]
        for v in range(nv):
            if GKv[v] == 1 or GKv[v] == 3:
                G[v] += jb
            elif GKv[v] == 2 and V > thV:
                G[v] += jb
            if GKv[v] == 3:                    # latch: opens at G > 0.5 while c* > theta_p, closes when c* <= theta_p
                if e > tp:
                    if G[v] > 0.5:
                        La[v] = 1.0
                else:
                    La[v] = 0.0
            pot = 1.0 if e > tp else 0.0
            dep = 1.0 if e > td else 0.0
            raw = pot
            if pot > 0.0 and thV > 0.0 and not (V > thV):
                pot = 0.0
            if pot > 0.0 and (GKv[v] == 1 or GKv[v] == 2) and G[v] <= 0.5:
                pot = 0.0
            if pot > 0.0 and GKv[v] == 3 and La[v] < 0.5:
                pot = 0.0
            if VMv[v] == 2 and e > tp:
                dep = 0.0
            rv = r[v]
            rv = rv + h * (-rv * (1 - rv) * (RS - rv) + pot * GP * (1 - rv) - dep * (1 - pot) * GD * rv)
            r[v] = min(max(rv, 0.0), 1.0)
            if v == 0:
                Pt[k] = int(pot); Dt[k] = int(dep * (1 - pot)); Ct[k] = int(raw > 0.0 and pot == 0.0)
        if hm != ch:
            ch = hm
            aV = math.exp(-hm / tauE1)
            x = hm / tauT
            if x < 0.1:
                aT = 1.0 - x; bT = hm
            else:
                aT = math.exp(-x); bT = tauT * (1.0 - aT)
            x = hm / tauN
            if x < 0.1:
                aN = 1.0 - x; bN = hm
            else:
                aN = math.exp(-x); bN = tauN * (1.0 - aN)
            x = hm / tauZ
            aZ = 1.0 - x if x < 0.1 else math.exp(-x)
        V = aV * V + tauE1 * (1.0 - aV) / iscV * s
        Z = Z + c
        if c > 0.0:
            for v in range(nv):
                xg = T[v] - thTg
                if v == 0:
                    xgA[k] = xg
                if xg < 0.0:
                    xg = 0.0
                if EKv[v] == 1 and V > thV:
                    xg = 0.0
                fm = math.pow(1.0 - am * math.tanh(xg), c)
                d0[v] = dmin + (d0[v] - dmin) * fm
                d2[v] = dmin + (d2[v] - dmin) * fm
        xN = N - thN; xZ = Z - thZ
        if xN > 0.0 and xZ > 0.0 and hm > 0.0:
            g = math.tanh(xN) * math.tanh(xZ); ex = math.exp(-r_no * g * hm)
            for v in range(nv):
                d0[v] = dmax - (dmax - d0[v]) * ex
                d1[v] = dmax - (dmax - d1[v]) * ex
        if hmp != chE:
            chE = hmp
            x = hmp / tauE1
            aE = 1.0 - x if x < 0.1 else math.exp(-x)
        for v in range(nv):
            jj = jmp[k]
            if EKv[v] == 2 and V > thV:
                jj = 0.0
            S1[v] = aE * S1[v] + jj
            if S1[v] > offT:
                T[v] = aT * T[v] + bT * (S1[v] - offT) / scT
            else:
                T[v] = aT * T[v]
            G[v] = G[v] * aV
        if s > offN:
            N = aN * N + bN * (s - offN) / sc
        else:
            N = aN * N
        Z = aZ * Z
        hmp = hm
        Vt[k] = V; Tt[k] = T[0]; St[k] = S1[0]
    return r, d0, d1, d2, Vt, Tt, St, Pt, Dt, Ct, xgA


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


def run_record(f, a, P, Q, prm):
    pair, proto = os.path.basename(f)[:-4].split("__")
    B = BatchV2([os.path.dirname(f)], protocols=[proto], pairs={pair}, fast=False, signals=("vdcc",), verbose=False)
    if not B.recs:
        return None
    B._stack(); r = B.recs[0]; H = B._H[0]
    t = r["t"]; n_t = len(t); h = np.diff(t)
    imp = ca_impulses(r, P)
    td, tp = rho_v4.thetas(r["c_pre"], r["c_post"], a, rho_v4.opts(P)[0])
    delay = MV.edge_params(r["syn"])["delay"]
    vd = r["vdcc"]; W = windows(r); npair = max(len(W), 1)
    wsl = [slice(*np.searchsorted(t, [lo, hi])) for lo, hi in W]
    wsl = [w for w in wsl if w.stop > w.start]
    thV = Q["theta_V"]
    rows = []; RF = []; D0 = []; D1 = []; D2 = []
    for i in range(len(r["syn"])):
        cnt = np.zeros(n_t); np.add.at(cnt, r["arr"][i], 1.0)
        ev = np.asarray(r["cev"][i], np.float64); ev = ev[np.isfinite(ev)]
        jev = np.zeros(n_t)
        if len(ev):
            np.add.at(jev, np.clip(np.searchsorted(t, ev - 1e-9), 0, n_t - 1), 1.0)
        E = np.asarray(r["effcai"][i], np.float64)
        out = replay(E, np.asarray(vd[i], np.float64), H, cnt, np.asarray(imp[i], np.float64), jev,
                     float(td[i]), float(tp[i]), float(r["rho0"][i]), *prm, VM, GK, EK)
        rf, d0, d1, d2, Vt, Tt, St, Pt, Dt, Ct, xgA = out
        RF.append(rf.copy()); D0.append(d0.copy()); D1.append(d1.copy()); D2.append(d2.copy())
        ta = np.asarray(r["prespikes"], float) + delay[i]
        _, w = MV.glu_weight(ev, ta, Q["tau_d_NMDA"]) if len(ta) else (ev, np.ones(len(ev)))
        nb = 0
        for lo, hi in W:
            a_in = ta[(ta >= lo) & (ta < hi)]
            if len(a_in):
                nb += int(np.sum((ev >= lo) & (ev < a_in.min())))
        ai = np.asarray(r["arr"][i]).ravel()
        xg = xgA[ai] if len(ai) else np.array([np.nan])
        hp = lambda m: float(np.sum(h[m[:-1]])) / npair
        q = dict(pair=pair, proto=proto, syn=int(r["syn"][i]), npair=len(W),
                 c_pre=float(r["c_pre"][i]), c_post=float(r["c_post"][i]), td=float(td[i]), tp=float(tp[i]),
                 rho0=float(r["rho0"][i]), rho_f=float(rf[0]),
                 pk_tp=med([E[w].max() / tp[i] for w in wsl]), pk_mM=med([E[w].max() for w in wsl]),
                 t_tp=hp(E > tp[i]), t_band=hp((E > td[i]) & (E <= tp[i])),
                 t_pot=float(np.sum(h * Pt[:-1])) / npair, t_conv=float(np.sum(h * Ct[:-1])) / npair,
                 t_dep=float(np.sum(h * Dt[:-1])) / npair,
                 Vpk=med([Vt[w].max() / thV for w in wsl]), V_arr=med(Vt[ai] / thV) if len(ai) else np.nan,
                 t_V=hp(Vt > thV),
                 n_ev=len(ev) / npair, n_ev_bef=nb / npair, w_ecb=float(np.sum(w)) / npair,
                 w_b=float(np.sum(1.0 - w)) / npair,
                 Spk=med([St[w].max() for w in wsl]), Tpk=med([Tt[w].max() for w in wsl]),
                 xg_arr=med(xg), open_arr=float(np.mean(xg > 0)) if len(ai) else np.nan,
                 d0=float(d0[0]), d1=float(d1[0]), d2=float(d2[0]))
        for k, nm in enumerate(VN[1:], 1):
            q[f"rf_{nm}"] = float(rf[k]); q[f"d0_{nm}"] = float(d0[k]); q[f"d1_{nm}"] = float(d1[k])
        rows.append(q)
    b = B.basis(r); r0 = r["rho0"]
    RF = np.array(RF); D0 = np.array(D0); D1 = np.array(D1); D2 = np.array(D2)
    rec = dict(pair=pair, proto=proto, n_syn=len(r0), rho_only=b.ratio(r0, RF[:, 0], np.zeros(len(r0))),
               pre_only=b.ratio(r0, r0, D0[:, 0]), ecb_only=b.ratio(r0, r0, D2[:, 0]),
               nmdar_block=b.ratio(r0, r0, np.zeros(len(r0))))
    for k, nm in enumerate(VN):
        rec[f"control|{nm}"] = b.ratio(r0, RF[:, k], D0[:, k])
        rec[f"mglu_block|{nm}"] = b.ratio(r0, RF[:, k], D1[:, k])
        rec[f"no_block|{nm}"] = b.ratio(r0, RF[:, k], D2[:, k])
    return rows, rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True)
    ap.add_argument("--model", required=True, help="l5 | l23 | l23l23 (a model of the fit json)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--csv", default=None, help="fit csv to compare (default <fit>{,_l23,_l23l23}.csv)")
    args = ap.parse_args()
    import pandas as pd
    t0 = time.time()
    fj = json.load(open(args.fit)); fa = fj["args"]
    batch_v2.BAP_GATE = fj.get("bap_gate")
    a = {k: fj["a"][k] for k in ("a00", "a01", "a10", "a11")}
    P = {**json.loads(fa["filters"]), **json.loads(fa["set"]), **fj["pre"]}
    P["theta_Ti"] = P["theta_NOi"]
    Q = {**MV.DEFAULTS, **P}
    M = {m["name"]: m for m in fj["models"]}[args.model]
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
    prm = (Q["dpre_min"], Q["dpre_max"], Q["A_mglu"], Q["A_NO"] / (1e3 * TAU_IND_GB), Q["theta_Tg"], Q["theta_N"],
           Q["theta_Z"], Q["theta_Te"], Q["Te_scale"], Q["theta_NOi"], Q["i_scale"], Q["tau_T"], Q["tau_NO"], Q["tau_Z"],
           Q["tau_E1"], Q["gamma_d"], Q["gamma_p"], Q["theta_V"], Q["i_scale"], float(RHO_STAR_GB), 1e3 * TAU_IND_GB)
    print(f"model {args.model}: {len(T)} targets, {len(protos)} protocols, basis {basis}", flush=True)
    print("params:", {k: round(float(Q[k]), 5) for k in ("theta_Te", "tau_T", "theta_Tg", "A_mglu", "dpre_min", "A_NO",
                                                            "theta_Z", "tau_Z", "tau_E1", "theta_V", "rho_gamma",
                                                            "gamma_d", "gamma_p", "tau_d_NMDA")}, flush=True)
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
        res = run_record(f, a, P, Q, prm)
        if res is None:
            continue
        rows += res[0]; recs.append(res[1])
        if n % 200 == 0:
            print(f"  {n}/{len(files)} records, {time.time() - t0:.0f} s", flush=True)
    os.makedirs(args.out, exist_ok=True)
    pre = os.path.join(args.out, args.model)
    S = pd.DataFrame(rows); R = pd.DataFrame(recs)
    S["up"] = (S.rho0 < 0.5) & (S.rho_f >= 0.5); S["down"] = (S.rho0 >= 0.5) & (S.rho_f < 0.5)
    S.to_csv(pre + "_syn.csv", index=False); R.to_csv(pre + "_rec.csv", index=False)
    # ---- targets: record means as gpu_v3 / gpu_v4 (loc selections for l23)
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
        q["rho_only"] = float(np.nanmean(x.rho_only)); q["pre_only"] = float(np.nanmean(x.pre_only))
        q["ecb_only"] = float(np.nanmean(x.ecb_only))
        out.append(q)
    D = pd.DataFrame(out); D.to_csv(pre + "_targets.csv", index=False)
    pd.set_option("display.width", 260); pd.set_option("display.max_columns", 60); pd.set_option("display.max_rows", 500)
    print(f"\n=== {args.model}: targets (pred = record mean; base must equal fit_pred)")
    print(D.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print(f"repro max |base - fit_pred| = {np.nanmax(np.abs(D.base - D.fit_pred)):.2e}")
    chi = {nm: float(np.nansum(((D[nm] - D["data"]) / D["sem"]) ** 2)) for nm in VN}
    print("chi2 per variant: " + ", ".join(f"{k} {v:.2f}" for k, v in chi.items()) + f"  (fit csv {np.nansum(((D.fit_pred - D['data']) / D['sem']) ** 2):.2f})")
    # ---- per protocol local quantities (medians over synapses; P_up / P_down as fractions)
    cols = ["c_pre", "c_post", "td", "tp", "pk_tp", "pk_mM", "t_tp", "t_band", "t_pot", "t_conv", "t_dep", "Vpk", "V_arr",
            "t_V", "n_ev", "n_ev_bef", "w_ecb", "w_b", "Spk", "Tpk", "xg_arr", "open_arr", "d0", "d1", "d2"]
    G = S.groupby("proto")
    Pq = G[cols].median()
    Pq["n_syn"] = G.size(); Pq["frac_rho0_1"] = G.rho0.apply(lambda z: float(np.mean(z >= 0.5)))
    Pq["P_up"] = G.apply(lambda z: z.up.sum() / max((z.rho0 < 0.5).sum(), 1))
    Pq["P_down"] = G.apply(lambda z: z.down.sum() / max((z.rho0 >= 0.5).sum(), 1))
    for nm in VN[1:]:
        Pq[f"Pup_{nm}"] = G.apply(lambda z, nm=nm: float(np.sum((z.rho0 < .5) & (z[f"rf_{nm}"] >= .5))) / max((z.rho0 < .5).sum(), 1))
        Pq[f"Pdn_{nm}"] = G.apply(lambda z, nm=nm: float(np.sum((z.rho0 >= .5) & (z[f"rf_{nm}"] < .5))) / max((z.rho0 >= .5).sum(), 1))
    # rho0 = 0 synapses (can only potentiate) and rho0 = 1 (can only depress): the post-rule drive per state
    for st, sel in (("r0", S.rho0 < 0.5), ("r1", S.rho0 >= 0.5)):
        Gs = S[sel].groupby("proto")
        for k in ("t_tp", "t_pot", "t_conv", "t_dep", "pk_tp", "Vpk"):
            Pq[f"{k}_{st}"] = Gs[k].median()
            Pq[f"{k}_{st}_mean"] = Gs[k].mean()
    Pq = Pq.reset_index(); Pq.to_csv(pre + "_proto.csv", index=False)
    print(f"\n=== {args.model}: per-protocol medians (times = ms per pairing; Vpk, V_arr in theta_V units; pk_tp in theta_p)")
    print(Pq.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print(f"\nsaved {pre}_{{syn,rec,targets,proto}}.csv ({len(S)} synapses, {len(R)} records), {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
