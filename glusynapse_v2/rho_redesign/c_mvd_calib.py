"""C-MVD calibration (ANCHORS_V11 section 1): per-synapse MVD pool from the extracted records, unweighted c_VDCC only.

For each synapse row of the listed protocols, on the kernel grid and with the gpu_v11 MVD arithmetic:
  V_{k+1} = a_k V_k + tau_E1 (1 - a_k) / i_scale s_k,  a_k = exp(-h_k / tau_E1)   (s = step-mean -ica_VDCC, npz "vdcc")
  Bg += own arrivals at step k, then Bg *= exp(-h_k / tau_d_NMDA);  glutamate window = Bg > b_c = exp(-W / tau_d_NMDA)
  Vg = max of V_k over steps in the window (the value theta_MVD must stay below for MVD to fire at that synapse).
P1_i = 1e5 x cexp vdcc_q_post (ecb_ref 2; NaN -> table median). Vg is reported in pool units and in P1_i units.
Also prints the t_a convention check: edge delay, and the lag from t_a to the next own VDCC current peak.
Run (1 CPU, from plastyfire/ via sjob.sh): python glusynapse_v2/rho_redesign/c_mvd_calib.py OUTDIR
"""
import glob, os, sys, time
import numpy as np, pandas as pd
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import batch_v2
import model_v2 as MV

W_ = "/scratch/dhuruva/s2g0321"; X = W_ + "/extracted"; T = "delta-split2-ljp25g0321-prefire"
TAU_E1, ISC, TAU_D, WGLU = 100.0, 1e-5, 70.0, 50.0
BC = float(np.exp(-WGLU / TAU_D))
# (label, group, dir, proto). group: fire = LTD rows MVD should give; no = rows that must not get rho-LTD; opt = Z 1AP +/-10;
# chk = t_a timing check only. Lanes (mglu_block = AM251, nmdar_block = APV) share the control record's VDCC trace.
ROWS = [
    ("Z train10 -10 last (ctrl, mglu_block)", "fire", "zilberter_l23l23_%s-vseg-rs", "zilberter_train10_50hz_dt-10ms_last"),
    ("Z train4 +4 last", "fire", "zilberter_l23l23_%s-vseg-rs", "zilberter_train4_50hz_dt+4ms_last"),
    ("Z 5ap 10 Hz +10", "fire", "zilberter_l23l23_%s-vseg-rs", "zilberter_5ap_10hz_dt+10ms"),
    ("Z train10 +4 last: APV lane (fire) = ctrl (no)", "fire/no", "zilberter_l23l23_%s-vseg-rs", "zilberter_train10_50hz_dt+4ms_last"),
    ("Z 1AP +10", "opt", "zilberter_l23l23_%s-vseg-rs", "zilberter_1ap_dt+10ms"),
    ("Z 1AP -10", "opt", "zilberter_l23l23_%s-vseg-rs", "zilberter_1ap_dt-10ms"),
    ("L5 Sj 20 Hz -10 (AM251)", "no", "ebner_%s-vca", "sjostrom_20hz_dt-10ms"),
    ("L5 Sj 0.1 Hz -10 (AM251)", "no", "ebner_%s-vca", "sjostrom_0.1hz_dt-10ms"),
    ("L5 Sj 50 Hz +10 (ctrl LTP, NMDAR block)", "no", "ebner_%s-vca", "sjostrom_50hz_dt+10ms"),
    ("L5 Sj 40 Hz +10 (LTP)", "no", "ebner_%s-vca", "sjostrom_40hz_dt+10ms"),
    ("L5 Markram 10 Hz +10 (LTP)", "no", "markram_%s-vca", "10Hz_10ms"),
    ("L5 Markram 10 Hz -10", "no", "markram_%s-vca", "10Hz_-10ms"),
    ("L5 Sj 0.1 Hz +10 (t_a check)", "chk", "ebner_%s-vca", "sjostrom_0.1hz_dt+10ms"),
    ("L5 Sj 0.1 Hz dt 0 (t_a check)", "chk", "l5extra_%s-vca", "sjostrom_0.1hz_dt0ms"),
]


def pairs(path):
    return sorted({p for p in open(path).read().replace("\n", ",").split(",") if p})


def cexp_p1(name):
    df = pd.read_csv(f"{W_}/cexp/{name}.csv")
    q = df["vdcc_q_post"].to_numpy(float) * 1e5
    key = df.pre_gid.astype(str) + "-" + df.post_gid.astype(str) + ":" + df.syn_id.astype(np.int64).astype(str)
    nan = int(np.sum(~np.isfinite(q)))
    q = np.where(np.isfinite(q), q, np.nanmedian(q))
    print(f"P1 {name}: {nan} / {len(q)} NaN rows -> table median {np.nanmedian(q):.4g}", flush=True)
    return dict(zip(key, q))


@njit(cache=False)
def _scan(a, bV, aD, s, c, bc):
    """Kernel order per step: Bg jump, read (V_k, Bg), then the V / Bg step. Returns (Vg, Vmax)."""
    V = 0.0; Bg = 0.0; vg = 0.0; vm = 0.0
    for k in range(s.shape[0]):
        Bg += c[k]
        if Bg > bc and V > vg:
            vg = V
        if V > vm:
            vm = V
        V = a[k] * V + bV[k] * s[k]
        Bg *= aD[k]
    return vg, vm


def one(r, P1):
    t = r["t"]; h = np.append(np.diff(t), 0.0)
    a = np.exp(-h / TAU_E1); bV = TAU_E1 * (1.0 - a) / ISC; aD = np.exp(-h / TAU_D)
    S = np.asarray(r["vdcc"], np.float32).astype(np.float64)
    delay = MV.edge_params(r["syn"])["delay"]
    out = []
    for ii in range(S.shape[0]):
        s = S[ii]; arr = np.asarray(r["arr"][ii]).ravel().astype(np.int64)
        c = np.bincount(arr, minlength=len(t)).astype(float)
        vg, vm = _scan(a, bV, aD, s, c, BC)
        ta = t[arr]; lag = np.nan                    # t_a -> next own VDCC current peak (first arrival)
        if ta.size:
            i0, i1 = np.searchsorted(t, [ta[0] - 5.0, ta[0] + 30.0])
            if i1 > i0:
                lag = t[i0 + int(np.argmax(s[i0:i1]))] - ta[0]
        post = np.asarray(r["postspikes"], float).ravel()
        pl = float(post[np.argmin(np.abs(post - ta[0]))] - ta[0]) if (post.size and ta.size) else np.nan
        key = f"{r['pair']}:{int(r['syn'][ii])}"
        out.append(dict(pair=r["pair"], syn=int(r["syn"][ii]), P1=P1.get(key, np.nan), Vg=vg, Vmax=vm,
                        delay=float(delay[ii]), lag_peak=lag, post_minus_ta=pl, n_arr=int(arr.size)))
    return out


def main():
    od = sys.argv[1]; os.makedirs(od, exist_ok=True); t0 = time.time()
    l5p = sorted(set(pairs(f"{os.path.dirname(HERE)}/subset24_pairs.txt")) & set(pairs(f"{W_}/doublets/l5_keep_pairs.txt")))
    zp = pairs(f"{W_}/doublets/l23l23_keep_pairs.txt")
    P1 = {"L5L5": cexp_p1("L5L5"), "L23L23": cexp_p1("L23L23")}
    rows = []
    for lab, grp, dd, proto in ROWS:
        z = dd.startswith("zilberter")
        batch_v2.BASIS_DIR = f"{W_}/basis_l23l23" if z else f"{W_}/basis_l5l5"
        d = X + "/" + dd % T
        fs = [f for f in sorted(glob.glob(f"{d}/*__{proto}.npz")) if os.path.basename(f).split("__")[0] in set(zp if z else l5p)]
        n0 = len(rows)
        for f in fs:
            B = batch_v2.BatchV2(d, protocols=[proto], pairs=[os.path.basename(f).split("__")[0]], verbose=False,
                                 fast=False, signals=("vdcc",))
            for r in B.recs:
                for o in one(r, P1["L23L23" if z else "L5L5"]):
                    rows.append(dict(row=lab, group=grp, proto=proto, **o))
            del B
        print(f"{lab}: {len(fs)} files, {len(rows) - n0} synapse rows ({time.time() - t0:.0f} s)", flush=True)
    df = pd.DataFrame(rows); df["Vg_P1"] = df.Vg / df.P1
    df.to_csv(f"{od}/c_mvd_rows.csv", index=False)
    qs = [.1, .25, .5, .75, .9]
    pd.set_option("display.width", 250)
    S = df.groupby(["group", "row"], sort=False).agg(
        n=("Vg", "size"), **{f"Vg_q{int(q * 100)}": ("Vg", lambda v, q=q: v.quantile(q)) for q in qs},
        **{f"P1u_q{int(q * 100)}": ("Vg_P1", lambda v, q=q: v.quantile(q)) for q in qs},
        P1_q50=("P1", "median"), delay_q50=("delay", "median"), lag_q50=("lag_peak", "median"),
        post_ta_q50=("post_minus_ta", "median"))
    print(S.round(3).to_string(), flush=True); S.to_csv(f"{od}/c_mvd_summary.csv")
    # threshold scan: fraction of synapses with Vg > theta per row; the best worst-row margin for fire vs strict no-fire rows
    for col, unit in (("Vg", "pool units (1e5 x charge)"), ("Vg_P1", "P1_i units")):
        v = df[col].to_numpy(float); v = v[np.isfinite(v) & (v > 0)]
        grid = np.geomspace(np.quantile(v, 0.01), np.quantile(v, 0.99), 60)
        F = {lab: np.array([(g[col] > th).mean() for th in grid]) for lab, g in df[df.group != "chk"].groupby("row", sort=False)}
        fire = [lab for lab, grp, *_ in ROWS if grp == "fire"]
        strict = [lab for lab, grp, *_ in ROWS if grp == "no" and ("AM251" in lab or "-10" in lab)]
        allno = [lab for lab, grp, *_ in ROWS if grp == "no"]
        for name, nos in (("strict no-fire (L5 -10 rows)", strict), ("all no-fire rows", allno)):
            m = np.min([F[l] for l in fire], axis=0) - np.max([F[l] for l in nos], axis=0)
            j = int(np.argmax(m))
            print(f"SCAN {unit} vs {name}: best theta {grid[j]:.4g}: min fire frac {np.min([F[l][j] for l in fire]):.3f}, "
                  f"max no-fire frac {np.max([F[l][j] for l in nos]):.3f} (margin {m[j]:.3f}; gap if > 0.5)", flush=True)
        T2 = pd.DataFrame(F, index=np.round(grid, 4)).iloc[::6]
        print(f"fraction of synapses with {col} > theta ({unit}):\n{T2.round(3).T.to_string()}", flush=True)
    print(f"done {time.time() - t0:.0f} s -> {od}", flush=True)


if __name__ == "__main__":
    main()
