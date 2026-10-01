"""Sjostrom 2004 dLTD diagnosis (DLTD_DIAG.md). Read-only: replays the gpu_v4_rho vamp-1 / t_drive-4 kernel on the CPU
(numba, same step order) for one fit, one record at a time, and prints per-synapse event, T-gate and c* numbers for
sjostrom04_dltd_step250ms (sj04 dir) and sjostrom_0.1hz_dt-10ms (ebner dir, the working eCB protocol).
Validation: the mean per-pair ratio must reproduce the fit csv pred (0.887 for sj04 control and mglu_block).
Run: sbatch --wrap (1 CPU). MEASURED 22132795 (v4_C1Ajs_s5, sj04 + ebner dt-10, 48 records): 910 MB MaxRSS, 0:12 elapsed,
100% CPU (requested 4G, 0:20) -> next run --mem=2G --time=0:15. Output: logs/diag_dltd_<id>.out; report DLTD_DIAG.md."""
import argparse, glob, json, math, os, sys
import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__)); V2 = os.path.dirname(HERE)
sys.path.insert(0, V2); sys.path.insert(0, HERE)
import batch_v2                                   # noqa: E402
from batch_v2 import BatchV2, ca_impulses, RHO_STAR_GB, TAU_IND_GB   # noqa: E402
import model_v2 as MV                             # noqa: E402
import rho_v4                                     # noqa: E402


@njit(cache=False)
def replay(E, S, H, cnt, jmp, td, tp, r0, dmin, dmax, am, r_no, thTg, thN, thZ, offT, scT, offN, sc, tauT, tauN, tauZ,
           tauE1, GD, GP, thV, iscV, RS, KT):
    """gpu_v4_rho kern for one synapse, VAMP 1, T_IMP, no_mode 0, use_s. Returns rho_f, d0, d1, d2 and traces
    S1 (eCB impulse state), T, the pot / dep flags actually applied, and the rho trace."""
    n = E.shape[0]
    Tt = np.zeros(n); St = np.zeros(n); Pt = np.zeros(n, np.int8); Dt = np.zeros(n, np.int8); Rt = np.zeros(n)
    xgA = np.full(n, np.nan)
    r = r0; T = 0.0; N = 0.0; Z = 0.0; V = 0.0
    d0 = 0.0; d1 = 0.0; d2 = 0.0; S1 = 0.0; hmp = 0.0
    for k in range(n):
        e = E[k]; s = S[k]; h = H[k]; hm = h * KT; c = cnt[k]
        pot = 1.0 if e > tp else 0.0
        dep = 1.0 if e > td else 0.0
        if pot > 0.0 and thV > 0.0 and not (V > thV):
            pot = 0.0
        r = r + h * (-r * (1 - r) * (RS - r) + pot * GP * (1 - r) - dep * (1 - pot) * GD * r)
        r = min(max(r, 0.0), 1.0)
        Pt[k] = int(pot); Dt[k] = int(dep * (1 - pot)); Rt[k] = r
        aV = math.exp(-hm / tauE1); V = aV * V + tauE1 * (1.0 - aV) / iscV * s
        Z = Z + c
        if c > 0.0:
            xg = max(T - thTg, 0.0); xgA[k] = T - thTg
            fm = math.pow(1.0 - am * math.tanh(xg), c)
            d0 = dmin + (d0 - dmin) * fm
            d2 = dmin + (d2 - dmin) * fm
        xN = N - thN; xZ = Z - thZ
        if xN > 0.0 and xZ > 0.0 and hm > 0.0:
            g = math.tanh(xN) * math.tanh(xZ); ex = math.exp(-r_no * g * hm)
            d0 = dmax - (dmax - d0) * ex
            d1 = dmax - (dmax - d1) * ex
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
        x = hmp / tauE1
        aE = 1.0 - x if x < 0.1 else math.exp(-x)
        S1 = aE * S1 + jmp[k]
        if S1 > offT:
            T = aT * T + bT * (S1 - offT) / scT
        else:
            T = aT * T
        if s > offN:
            N = aN * N + bN * (s - offN) / sc
        else:
            N = aN * N
        Z = aZ * Z
        hmp = hm
        Tt[k] = T; St[k] = S1
    return r, d0, d1, d2, St, Tt, Pt, Dt, Rt, xgA


def run_record(f, a, P):
    pair, proto = os.path.basename(f)[:-4].split("__")
    B = BatchV2([os.path.dirname(f)], protocols=[proto], pairs={pair}, fast=False, signals=("vdcc",), verbose=False)
    if not B.recs:
        return None
    B._stack(); r = B.recs[0]; H = B._H[0]
    t = r["t"]; n_t = len(t); h = np.diff(t)
    d = np.load(f); K = float(d["K_ca"]) if "K_ca" in d.files else MV.DEFAULTS["K_ca"]
    imp = ca_impulses(r, P)
    g, _, _ = rho_v4.opts(P)
    td, tp = rho_v4.thetas(r["c_pre"], r["c_post"], a, g)
    delay = MV.edge_params(r["syn"])["delay"]
    Q = {**MV.DEFAULTS, **P}
    vd = r["vdcc"]
    rows = []; d0s = []; d1s = []; rfs = []
    for i in range(len(r["syn"])):
        cnt = np.zeros(n_t); np.add.at(cnt, r["arr"][i], 1.0)
        out = replay(np.asarray(r["effcai"][i], np.float64), np.asarray(vd[i], np.float64), H, cnt,
                     np.asarray(imp[i], np.float64), float(td[i]), float(tp[i]), float(r["rho0"][i]),
                     Q["dpre_min"], Q["dpre_max"], Q["A_mglu"], Q["A_NO"] / (1e3 * TAU_IND_GB), Q["theta_Tg"],
                     Q["theta_N"], Q["theta_Z"], Q["theta_Te"], Q["Te_scale"], Q["theta_NOi"], Q["i_scale"], Q["tau_T"],
                     Q["tau_NO"], Q["tau_Z"], Q["tau_E1"], Q["gamma_d"], Q["gamma_p"], Q["theta_V"], Q["i_scale"],
                     float(RHO_STAR_GB), 1e3 * TAU_IND_GB)
        rf, d0, d1, d2, St, Tt, Pt, Dt, Rt, xgA = out
        rfs.append(rf); d0s.append(d0); d1s.append(d1)
        ta = np.asarray(r["prespikes"]) + delay[i]; npair = len(ta)
        ev = np.asarray(r["cev"][i], np.float64); ev = ev[np.isfinite(ev)]
        _, w = MV.glu_weight(ev, ta, Q["tau_d_NMDA"])
        # per pairing: events in [ta - 500, ta + 1000), before / after arrival
        nb = na = 0; ws = []; lat = []
        pk = []; plat = []; qv = []; Tpost = []; tdur = []
        sl = lambda lo, hi: slice(*np.searchsorted(t, [lo, hi]))
        for tj in ta:
            m = (ev >= tj - 500) & (ev < tj + 1000)
            nb += np.sum(m & (ev < tj)); na += np.sum(m & (ev >= tj)); ws += list(w[m]); lat += list(ev[m] - tj)
            pk.append(vd[i][sl(tj - 100, tj + 400)].max() / K)
            st = sl(tj + 50, tj + 250)
            plat.append(np.median(vd[i][st]) / K if st.stop > st.start else np.nan)
            q = sl(tj - 50, tj + 300); q = slice(q.start, min(q.stop, n_t - 1))
            qv.append(np.sum(vd[i][q] * h[q]))                            # nA ms = pC
            Tpost.append(Tt[sl(tj, tj + 400)].max() - Q["theta_Tg"])
            tdur.append(np.median(h[q]) if q.stop > q.start else np.nan)
        xg = xgA[np.isfinite(xgA)] if np.any(np.isfinite(xgA)) else np.array([np.nan])
        flip = np.sum(np.abs(np.diff((Rt >= 0.5).astype(int))))
        rows.append(dict(pair=pair, syn=int(r["syn"][i]), ev_before=nb / npair, ev_after=na / npair,
                         ev_lat_med=float(np.median(lat)) if lat else np.nan,
                         w_mean=float(np.mean(ws)) if ws else np.nan, pk_K=float(np.median(pk)), plat_K=float(np.nanmedian(plat)),
                         Q_pC=float(np.median(qv)), Smax=float(St.max()), Tmax=float(Tt.max()),
                         xg_arr_max=float(np.max(xg)), xgT_post400=float(np.max(Tpost)),
                         t_dep_ms=float(np.sum(h[(r["effcai"][i][:-1] > td[i])])),
                         t_pot_ms=float(np.sum(h[(r["effcai"][i][:-1] > tp[i])])),
                         t_dep_eff_ms=float(np.sum(h * Dt[:-1])), t_pot_eff_ms=float(np.sum(h * Pt[:-1])),
                         peak_cstar=float(r["effcai"][i].max()), td=float(td[i]), tp=float(tp[i]), c_post=float(r["c_post"][i]),
                         rho0=float(r["rho0"][i]), rho_f=float(rf), nflip=int(flip), d0=float(d0), d1=float(d1),
                         dt_grid=float(np.nanmedian(tdur))))
    b = B.basis(r)
    rho_f = np.array(rfs)
    rc = b.ratio(r["rho0"], rho_f, np.array(d0s)); rm = b.ratio(r["rho0"], rho_f, np.array(d1s))
    return rows, rc, rm, K, (float(vd.min()), float(vd.max())), len(r["prespikes"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True)
    ap.add_argument("--sets", default="sj04_delta-prefire-vca:sjostrom04_dltd_step250ms,"
                                      "ebner_delta-prefire-vca:sjostrom_0.1hz_dt-10ms")
    args = ap.parse_args()
    fj = json.load(open(args.fit))
    a = {k: fj["a"][k] for k in ("a00", "a01", "a10", "a11")}
    P = {**json.loads(fj["args"]["filters"]), **json.loads(fj["args"]["set"]), **fj["pre"]}
    P["theta_Ti"] = P["theta_NOi"]
    pairs = set(fj["args"]["pairs"].split(","))
    print("params:", {k: P[k] for k in ("t_drive", "theta_Te", "tau_T", "theta_Tg", "A_mglu", "tau_E1", "theta_V",
                                         "rho_gamma", "gamma_d", "gamma_p")}, flush=True)
    import pandas as pd
    for spec in args.sets.split(","):
        dname, proto = spec.split(":")
        rows, rcs, rms = [], [], []
        for f in sorted(glob.glob(os.path.join(V2, "extracted", dname, f"*__{proto}.npz"))):
            if os.path.basename(f).split("__")[0] not in pairs:
                continue
            res = run_record(f, a, P)
            if res is None:
                continue
            rr, rc, rm, K, vr, npre = res
            rows += rr; rcs.append(rc); rms.append(rm)
            print(f"{proto} {rr[0]['pair']}: ratio ctrl {rc:.4f} mglu {rm:.4f}  K {K:.3e}  vdcc range {vr}  n_pre {npre}",
                  flush=True)
        df = pd.DataFrame(rows)
        print(f"\n=== {proto} ({dname}): {len(set(df.pair))} pairs, {len(df)} synapses; pred ctrl {np.nanmean(rcs):.4f} "
              f"mglu {np.nanmean(rms):.4f}")
        pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40); pd.set_option("display.max_rows", 400)
        print(df.round(4).to_csv(index=False))
        num = df.drop(columns=["pair", "syn"])
        print("median:\n" + num.median().round(4).to_string())
        print("mean:\n" + num.mean().round(4).to_string())
        print(f"frac syn with any event after arrival {np.mean(df.ev_after > 0):.3f}, before {np.mean(df.ev_before > 0):.3f}; "
              f"frac xg_arr_max > 0 {np.mean(df.xg_arr_max > 0):.3f}; frac xgT_post400 > 0 {np.mean(df.xgT_post400 > 0):.3f}; "
              f"flips 1->0 {int(np.sum((df.rho0 >= .5) & (df.rho_f < .5)))}, 0->1 {int(np.sum((df.rho0 < .5) & (df.rho_f >= .5)))}"
              f" of {len(df)}; rho0=1 {int(np.sum(df.rho0 >= .5))}; t_dep_eff>0 {int(np.sum(df.t_dep_eff_ms > 0))}", flush=True)


if __name__ == "__main__":
    main()
