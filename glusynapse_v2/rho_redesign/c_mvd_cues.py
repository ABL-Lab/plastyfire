"""C-MVD cue search (C_MVD_CALIB.md section 2): synapse-local cues that might separate the Zilberter LTD rows from the L5
rows that must not get rho-LTD. Same records, rows, pool and glutamate window as c_mvd_calib.py (W 50 ms, b_c 0.49).
Per synapse row, over the steps k inside the own glutamate window (G_k: Bg > b_c):
  pool      max V                                   (the C-MVD baseline; abs and / P1)
  a_lowd    max V over G_k with c*_k < theta_d,i     (a: VDCC high while own spine Ca below theta_d; s1C_s a00/a01)
  a_lowp    max V over G_k with c*_k < theta_p,i
  b_int     sum V h over G_k / P1                    (b: time integral of the pool in the window, P1 ms)
  b_dur1    time in G_k with V > 1 P1 (ms)           (b: duration above threshold)
  c_m51     peak of own-arrival low-pass, tau 51.3 ms (c: mGluR1 deactivation, Marcaggi 2009); c_m363: tau 363 ms
  d_ratio   max over G_k of (V / P1) / (c*_k / c_post,i)  (d: VDCC pool vs total spine Ca, both in own single-bAP units)
  e_sh      max own shaft dCa above rest over G_k (uM); e_shint: sum over G_k of shaft dCa h (uM ms)
  e_varr    V at the own arrival / P1 (VDCC Ca before glutamate, max over arrivals)
Scan: for each cue and both signs (fire if cue > theta, or cue < theta), the best min(fire frac) - max(no-fire frac).
Run (1 CPU, from plastyfire/ via sjob.sh): bash glusynapse_v2/rho_redesign/run_c_mvd.sh cues
"""
import glob, json, os, sys, time
import numpy as np, pandas as pd
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import batch_v2
from c_mvd_calib import ROWS, W_, X, T, TAU_E1, ISC, TAU_D, BC, pairs, cexp_p1

A = json.load(open(os.path.join(HERE, "s1C_s.json")))["a"]
CUES = ["pool", "pool_P1", "a_lowd", "a_lowd_P1", "a_lowp_P1", "b_int", "b_dur1", "c_m51", "c_m363", "d_ratio",
        "e_sh", "e_shint", "e_varr"]


@njit(cache=False)
def _cues(a, bV, aD, h, s, e, x, c, bc, thd, thp, p1, cpost, am51, am363):
    V = 0.0; Bg = 0.0; m51 = 0.0; m363 = 0.0
    pool = 0.0; lowd = 0.0; lowp = 0.0; bint = 0.0; dur = 0.0; pk51 = 0.0; pk363 = 0.0; rat = 0.0
    sh = 0.0; shi = 0.0; varr = 0.0
    for k in range(s.shape[0]):
        Bg += c[k]
        m51 += c[k]; m363 += c[k]
        if m51 > pk51:
            pk51 = m51
        if m363 > pk363:
            pk363 = m363
        if c[k] > 0.0 and V > varr:
            varr = V
        if Bg > bc:
            if V > pool:
                pool = V
            if e[k] < thd and V > lowd:
                lowd = V
            if e[k] < thp and V > lowp:
                lowp = V
            bint += V * h[k]
            if V > p1:
                dur += h[k]
            if e[k] > 0.0:
                r = (V / p1) / (e[k] / cpost)
                if r > rat:
                    rat = r
            if x[k] > sh:
                sh = x[k]
            shi += x[k] * h[k]
        V = a[k] * V + bV[k] * s[k]
        Bg *= aD[k]; m51 *= am51[k]; m363 *= am363[k]
    return pool, lowd, lowp, bint, dur, pk51, pk363, rat, sh, shi, varr


def one(r, P1):
    t = r["t"]; h = np.append(np.diff(t), 0.0)
    a = np.exp(-h / TAU_E1); bV = TAU_E1 * (1.0 - a) / ISC; aD = np.exp(-h / TAU_D)
    am51 = np.exp(-h / 51.3); am363 = np.exp(-h / 363.0)
    S = np.asarray(r["vdcc"], np.float32).astype(np.float64)
    E = np.asarray(r["effcai"], np.float64)
    sh = np.asarray(r["shaft_cai"], np.float64)
    k0 = 0
    ev = np.concatenate([np.asarray(r["prespikes"], float).ravel(), np.asarray(r["postspikes"], float).ravel()])
    if ev.size:
        k0 = max(int(np.searchsorted(t, ev.min() - 1.0)) - 1, 0)
    dsh = np.maximum((sh - sh[:, k0:k0 + 1]) * 1e3, 0.0)           # uM above the pre-stimulus rest (v8 src 2)
    out = []
    for ii in range(S.shape[0]):
        arr = np.asarray(r["arr"][ii]).ravel().astype(np.int64)
        c = np.bincount(arr, minlength=len(t)).astype(float)
        cp, cq = float(r["c_pre"][ii]), float(r["c_post"][ii])
        thd = A["a00"] * cp + A["a01"] * cq; thp = A["a10"] * cp + A["a11"] * cq
        key = f"{r['pair']}:{int(r['syn'][ii])}"
        p1 = P1.get(key, np.nan)
        p1 = p1 if np.isfinite(p1) and p1 > 0 else np.nanmedian(list(P1.values()))
        v = _cues(a, bV, aD, h, S[ii], E[ii], dsh[ii], c, BC, thd, thp, p1, max(cq, 1e-12), am51, am363)
        pool, lowd, lowp, bint, dur, pk51, pk363, rat, shm, shi, varr = v
        out.append(dict(pair=r["pair"], syn=int(r["syn"][ii]), P1=p1, pool=pool, pool_P1=pool / p1, a_lowd=lowd,
                        a_lowd_P1=lowd / p1, a_lowp_P1=lowp / p1, b_int=bint / p1, b_dur1=dur, c_m51=pk51, c_m363=pk363,
                        d_ratio=rat, e_sh=shm, e_shint=shi, e_varr=varr / p1))
    return out


def scan(df, col, fire, nos):
    v = df[col].to_numpy(float); v = v[np.isfinite(v)]
    lo, hi = np.quantile(v, 0.01), np.quantile(v, 0.99)
    grid = np.geomspace(max(lo, 1e-9), max(hi, 2e-9), 80) if lo > 0 else np.linspace(lo, hi, 80)
    best = (-9.0, None, None, None, None)
    for sgn in (1, -1):
        F = {lab: np.array([((g[col] > th) if sgn > 0 else (g[col] < th)).mean() for th in grid])
             for lab, g in df.groupby("row", sort=False) if lab in fire + nos}
        m = np.min([F[l] for l in fire], axis=0) - np.max([F[l] for l in nos], axis=0)
        j = int(np.argmax(m))
        if m[j] > best[0]:
            best = (float(m[j]), sgn, float(grid[j]), float(np.min([F[l][j] for l in fire])),
                    float(np.max([F[l][j] for l in nos])))
    return best


def main():
    od = sys.argv[1]; os.makedirs(od, exist_ok=True); t0 = time.time()
    l5p = sorted(set(pairs(f"{os.path.dirname(HERE)}/subset24_pairs.txt")) & set(pairs(f"{W_}/doublets/l5_keep_pairs.txt")))
    zp = pairs(f"{W_}/doublets/l23l23_keep_pairs.txt")
    P1 = {"L5L5": cexp_p1("L5L5"), "L23L23": cexp_p1("L23L23")}
    rows = []
    for lab, grp, dd, proto in ROWS:
        if grp == "chk":
            continue
        z = dd.startswith("zilberter")
        batch_v2.BASIS_DIR = f"{W_}/basis_l23l23" if z else f"{W_}/basis_l5l5"
        d = X + "/" + dd % T; ok = set(zp if z else l5p)
        fs = [f for f in sorted(glob.glob(f"{d}/*__{proto}.npz")) if os.path.basename(f).split("__")[0] in ok]
        for f in fs:
            B = batch_v2.BatchV2(d, protocols=[proto], pairs=[os.path.basename(f).split("__")[0]], verbose=False,
                                 fast=False, signals=("vdcc", "shaft_cai"))
            for r in B.recs:
                for o in one(r, P1["L23L23" if z else "L5L5"]):
                    rows.append(dict(row=lab, group=grp, proto=proto, **o))
            del B
        print(f"{lab}: {len(fs)} files ({time.time() - t0:.0f} s)", flush=True)
    df = pd.DataFrame(rows); df.to_csv(f"{od}/c_mvd_cues_rows.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    print("per-row medians:\n" + df.groupby(["group", "row"], sort=False)[CUES].median().round(3).to_string(), flush=True)
    print("per-row q10 / q90:\n" + df.groupby(["group", "row"], sort=False)[CUES].quantile([.1, .9]).round(3).to_string(),
          flush=True)
    fire3 = [l for l, g, *_ in ROWS if g == "fire"]
    fire5 = fire3 + [l for l, g, *_ in ROWS if g == "opt"]
    strict = [l for l, g, *_ in ROWS if g == "no" and ("AM251" in l or "-10" in l)]
    allno = [l for l, g, *_ in ROWS if g == "no"]
    res = []
    for col in CUES:
        for fn, fire in (("fire3 (trains, 5ap)", fire3), ("fire5 (+ 1AP +/-10)", fire5)):
            for nn, nos in (("strict (L5 -10)", strict), ("all no-fire", allno)):
                m, sgn, th, ff, nf = scan(df, col, fire, nos)
                res.append(dict(cue=col, fire=fn, nofire=nn, margin=m, fires_if=">" if sgn > 0 else "<", theta=th,
                                min_fire_frac=ff, max_nofire_frac=nf))
    R = pd.DataFrame(res); R.to_csv(f"{od}/c_mvd_cues_scan.csv", index=False)
    print("SCAN (margin = min fire frac - max no-fire frac; gap if > 0.5):\n" + R.round(4).to_string(index=False), flush=True)
    print(f"done {time.time() - t0:.0f} s -> {od}", flush=True)


if __name__ == "__main__":
    main()
