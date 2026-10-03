"""s6mech CPU probe (S6MECH_DESIGN.md section 1 "Needs a run" and section 3 kernel test 2/3): the coincidence event C per own
pre arrival, from the extracted L5 records, with the gpu_v11 peak-veto arithmetic that gpu_v12_rho uses as C.

C_a = [max own spine VDCC current step value over the steps overlapping (t_a + 2, t_a + 17] > 0.5 x I1_i]
  t_a   = prespike + edge delay (exact arrival, t_exact 1)
  steps = the record's sample intervals [t_j, t_j+1] (npz "vdcc" is the step value), as the kernel's window scan
  I1_i  = median peak step VDCC current over [x - 1, x + 8] ms of the synapse's isolated bAPs (no other post spike in
          (x - 100, x + 20)); tier 1 = no own arrival in (x - 150, x + 8), else tier 2; NaN -> pooled median. Computed
          here from the probed records of that synapse only (the kernel pools every L5 record of the fit), so a
          synapse's I1_i can differ slightly; the kernel's exact per-protocol C fraction is the "v12 C l5 <proto>" line
          of a gpu_v12 run with an arm on (or c_report 1).
Per row: C fraction over all (synapse, arrival); C fraction by position in the train (trains split at gaps > 100 ms);
the mean C count per 5-spike train (design C/5 column); and the B3 / B4 annul-eligible fraction (arrival without C whose
previous own arrival within T_m 25 / 50 ms had C), with the share of those that are the last pre of their train.
Expected (design lag table): 20 Hz +10 5/5; 40 / 50 Hz -10 4/5 (last none); 20 Hz 0 / +25 / -25 / -10 0/5; Markram -10 0/5,
+10 5/5; S07 step pair "most" (decides B3: a high S07 C fraction means B3 removes S07's hidden eCB factor).
Run (1 CPU, from plastyfire/ via sjob.sh): python glusynapse_v2/rho_redesign/diag_s6_c.py OUTDIR
"""
import glob, os, sys, time
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import batch_v2
import model_v2 as MV

W_ = "/scratch/dhuruva/s2g0321"; X = W_ + "/extracted"; T = "delta-split2-ljp25g0321-prefire"
VT0, VTV, K = 2.0, 17.0, 0.5          # gpu_v11 3V / 4V veto window and veto_peak_k (ANCHORS_V11 section 2)
TMS = (25.0, 50.0)                    # B3 / B4 T_m
# (label, dir, proto): the design's lag table rows, plus the 0.1 Hz rows as extra single-bAP (I1_i) sources
ROWS = [
    ("20 Hz +10", "ebner", "sjostrom_20hz_dt+10ms"), ("20 Hz 0", "l5extra", "sjostrom_20hz_dt0ms"),
    ("20 Hz +25", "l5extra", "sjostrom_20hz_dt+25ms"), ("20 Hz -25", "l5extra", "sjostrom_20hz_dt-25ms"),
    ("20 Hz -10", "ebner", "sjostrom_20hz_dt-10ms"), ("40 Hz +10", "ebner", "sjostrom_40hz_dt+10ms"),
    ("40 Hz -10", "ebner", "sjostrom_40hz_dt-10ms"), ("40 Hz 0", "l5extra", "sjostrom_40hz_dt0ms"),
    ("50 Hz +10", "ebner", "sjostrom_50hz_dt+10ms"), ("50 Hz -10", "ebner", "sjostrom_50hz_dt-10ms"),
    ("50 Hz 0", "l5extra", "sjostrom_50hz_dt0ms"), ("Markram +10", "markram", "10Hz_10ms"),
    ("Markram -10", "markram", "10Hz_-10ms"), ("S07 step pair", "sj07", "sjostrom07_step200ms_pair"),
    ("0.1 Hz +10", "ebner", "sjostrom_0.1hz_dt+10ms"), ("0.1 Hz -10", "ebner", "sjostrom_0.1hz_dt-10ms"),
]


def pairs(path):
    return sorted({p for p in open(path).read().replace("\n", ",").split(",") if p})


def train_pos(pre):
    """1-based position of each pre spike in its train (new train after a gap > 100 ms) and the train length."""
    pos = np.ones(len(pre), int); tid = np.zeros(len(pre), int)
    for k in range(1, len(pre)):
        if pre[k] - pre[k - 1] > 100.0:
            tid[k] = tid[k - 1] + 1
        else:
            tid[k] = tid[k - 1]; pos[k] = pos[k - 1] + 1
    ln = np.bincount(tid)[tid]
    return pos, ln


def one(r, lab, i1src):
    t = np.asarray(r["t"], np.float64); T_ = len(t)
    S = np.asarray(r["vdcc"], np.float32).astype(np.float64)
    pre = np.sort(np.asarray(r["prespikes"], float).ravel())
    post = np.sort(np.asarray(r["postspikes"], float).ravel())
    delay = np.asarray(MV.edge_params(r["syn"])["delay"], float)
    iso = np.array([not np.any((post != x) & (post > x - 100.0) & (post < x + 20.0)) for x in post], bool)
    pos, ln = train_pos(pre)
    out = []
    for ii in range(S.shape[0]):
        key = (r["pair"], int(r["syn"][ii]))
        tav = t[np.asarray(r["arr"][ii]).ravel().astype(np.int64)]
        L = i1src.setdefault(key, ([], []))
        for x in post[iso]:                                       # gpu_v11 veto_peak I1_i source
            a, b = np.searchsorted(t, [x - 1.0, x + 8.0])
            if b > a:
                L[0 if not np.any((tav > x - 150.0) & (tav < x + 8.0)) else 1].append(float(S[ii, a:b].max()))
        ta = pre + delay[ii]
        j0 = np.searchsorted(t, ta + VT0, side="right") - 1     # first step with t_j+1 > t_a + 2
        j1 = np.minimum(np.searchsorted(t, ta + VTV, side="left"), T_ - 1)   # steps j < j1 (t_j < t_a + 17), j <= T - 2
        pk = np.array([S[ii, max(a, 0):b].max() if b > max(a, 0) else 0.0 for a, b in zip(j0, j1)])
        for k in range(len(ta)):
            out.append((lab, r["pair"], int(r["syn"][ii]), k, float(ta[k]), int(pos[k]), int(ln[k]), float(pk[k])))
    return out


def main():
    od = sys.argv[1]; os.makedirs(od, exist_ok=True); t0 = time.time()
    l5p = sorted(set(pairs(f"{os.path.dirname(HERE)}/subset24_pairs.txt")) & set(pairs(f"{W_}/doublets/l5_keep_pairs.txt")))
    batch_v2.BASIS_DIR = f"{W_}/basis_l5l5"
    rows, i1src = [], {}
    for lab, dd, proto in ROWS:
        d = f"{X}/{dd}_{T}-vca"
        fs = [f for f in sorted(glob.glob(f"{d}/*__{proto}.npz")) if os.path.basename(f).split("__")[0] in set(l5p)]
        n0 = len(rows)
        for f in fs:
            B = batch_v2.BatchV2(d, protocols=[proto], pairs=[os.path.basename(f).split("__")[0]], verbose=False,
                                 fast=False, signals=("vdcc",))
            for r in B.recs:
                rows += one(r, lab, i1src)
            del B
        print(f"{lab}: {len(fs)} files, {len(rows) - n0} arrivals ({time.time() - t0:.0f} s)", flush=True)
    df = pd.DataFrame(rows, columns=["row", "pair", "syn", "k", "ta", "pos", "train_len", "peak"])
    p1 = {k: np.median(v[0]) if v[0] else (np.median(v[1]) if v[1] else np.nan) for k, v in i1src.items()}
    tiers = pd.Series({k: 1 if v[0] else (2 if v[1] else 0) for k, v in i1src.items()}).value_counts().to_dict()
    pool = float(np.nanmedian(list(p1.values())))
    df["I1"] = [p1.get((p, s), np.nan) for p, s in zip(df.pair, df.syn)]
    nnan = int(df.I1.isna().sum()); df["I1"] = df.I1.fillna(pool)
    df["C"] = df.peak > K * df.I1
    print(f"I1_i: {len(p1)} synapses, tiers (1 / 2 / none) {tiers}, pooled median {pool:.4g}; {nnan} arrival rows filled",
          flush=True)
    df = df.sort_values(["row", "pair", "syn", "ta"]).reset_index(drop=True)
    for tm in TMS:                                                # B3 / B4 annul eligibility
        el = np.zeros(len(df), bool)
        for _, g in df.groupby(["row", "pair", "syn"], sort=False):
            ta = g.ta.to_numpy(); c = g.C.to_numpy(); idx = g.index.to_numpy()
            lastc = -1e30
            for m in range(len(ta)):
                if not c[m] and ta[m] - lastc <= tm + 1e-3:
                    el[idx[m]] = True
                if c[m]:
                    lastc = ta[m]
        df[f"annul{int(tm)}"] = el
    df.to_csv(f"{od}/s6_c_arrivals.csv.gz", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    res = []
    for lab, *_ in ROWS:
        g = df[df.row == lab]
        if g.empty:
            continue
        d = dict(row=lab, synapses=g.groupby(["pair", "syn"]).ngroups, arrivals=len(g), C_frac=g.C.mean())
        five = g[g.train_len == 5]
        if len(five):
            d["C_per5"] = five.C.mean() * 5
        for p in range(1, 11):
            gp = g[g.pos == p]
            if len(gp):
                d[f"pos{p}"] = gp.C.mean()
        for tm in TMS:
            a = g[f"annul{int(tm)}"]
            d[f"annul{int(tm)}"] = a.mean()
            d[f"annul{int(tm)}_last"] = (a & (g.pos == g.train_len)).sum() / max(a.sum(), 1)
        res.append(d)
    R = pd.DataFrame(res); R.to_csv(f"{od}/s6_c_rows.csv", index=False)
    print("C per row (C_frac over synapse x arrival; C_per5 = mean C count per 5-spike train; posN = C frac of the N-th pre "
          "of a train; annulT = annul-eligible frac at T_m, annulT_last = share of those that are the last pre):\n"
          + R.round(3).to_string(index=False), flush=True)
    s7 = R[R.row == "S07 step pair"]
    if len(s7):
        print(f"S07 PROBE: C fraction of S07 step-pair pre arrivals = {float(s7.C_frac.iloc[0]):.4f} "
              f"(annul-eligible T_m 25 {float(s7.annul25.iloc[0]):.4f}, T_m 50 {float(s7.annul50.iloc[0]):.4f})", flush=True)
    print(f"done {time.time() - t0:.0f} s -> {od}", flush=True)


if __name__ == "__main__":
    main()
