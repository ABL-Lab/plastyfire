"""Representative traces of the v5c rule (fit v5_S1_V5c_s6_39, emodel delta-split1): BCL live (live_v5.run_task, prefire
phase, cooker expression, GluSynapseV5) vs the offline kernel (v5_rec_traces = compare_prefire_v5.v5_rec with time
series) on the same records. Plots: plot_traces_v5.py.

Groups (GROUPS): one pair per group, one record per member. Fixed pair 209736-199745 for the L5 Sjostrom groups: on the
partial run (results/S1_V5c_partial_records.csv) it is within 0.04 of the fit's prediction (fit csv pred) for 50 Hz +10,
0.1 Hz -10 (control and mglu_block), r50 -120 and the S07 pair (control and mglu_block). Markram and Letzkus (only the
pilot pairs in the partial run): auto, the candidate pair with the smallest sum |offline ratio - pred| / SEM over the
group's members (offline_v5), candidates = the fit's L5 pairs / the all_protocols letzkus_distal L2/3 pairs with records.

Per record (one npz in --out): offline kernel traces (fine window = one induction repetition on the record grid; rho and
dpre over the whole induction interpolated onto the live coarse grid; eCB step times; max effcai, max pool), live
NEURON vectors (soma v, rho_GB, dpre_GB, effcai_GB, cai_CR, Vg_GB = spine VDCC Ca pool, W_GB = eCB trigger pool,
bglu_GB, necb_GB, ica_VDCC at dt_fine in the window; rho_GB, dpre_GB, necb_GB at dt_coarse over the run), end states,
weight ratios through the pair's EPSP basis (offline_v5's b) and the fit csv target (mean, SEM, pred).
The repetition shown: the one with the first offline eCB step of the control record, else the middle one.

    python glusynapse_v2/bcl_validation/traces_v5.py --check            # task list, kernel-copy test, then exit
    python glusynapse_v2/bcl_validation/traces_v5.py --workers 11       # env: run_traces_S1_V5c.sh
"""
import argparse, gc, json, os, subprocess, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import live_v5 as L                                       # noqa: E402  (light: numpy only at import)

FIT = os.path.join(L.V2, "rho_redesign/results/v5_S1_V5c_s6_39.json")
OUT = "/scratch/dhuruva/bcl_traces_S1_V5c"
DT_FINE, N_COARSE, PAD0, PAD1, WMAX = 0.1, 30000, 250.0, 400.0, 2000.0
FINE = ("rho_GB", "dpre_GB", "effcai_GB", "cai_CR", "Vg_GB", "W_GB", "bglu_GB", "necb_GB", "ica_VDCC")
COARSE = ("rho_GB", "dpre_GB", "necb_GB")
L5FIX = "209736-199745"
# member = (proto, cond, fit-csv target, label)
GROUPS = [
    dict(name="markram10hz", path="L5", pair=None, fallback="180351-198084",
         title="L5$\\to$L5, Markram et al. 1997: 10 Hz pairing",
         members=[("10Hz_10ms", "control", "10Hz_10ms", "$\\Delta t$ +10 ms (LTP)"),
                  ("10Hz_-10ms", "control", "10Hz_-10ms", "$\\Delta t$ $-$10 ms (LTD)")]),
    dict(name="sj01freq", path="L5", pair=L5FIX, title="L5$\\to$L5, Sjöström et al. 2001: frequency dependence",
         members=[("sjostrom_50hz_dt+10ms", "control", "sjostrom_50hz_dt+10ms", "50 Hz, $\\Delta t$ +10 ms (LTP)"),
                  ("sjostrom_0.1hz_dt-10ms", "control", "sjostrom_0.1hz_dt-10ms", "0.1 Hz, $\\Delta t$ $-$10 ms (eCB tLTD)"),
                  ("sjostrom_0.1hz_dt-10ms", "mglu_block", "sjostrom_0.1hz_dt-10ms",
                   "0.1 Hz, $\\Delta t$ $-$10 ms, mGluR block")]),
    dict(name="sj03r50", path="L5", pair=L5FIX, title="L5$\\to$L5, Sjöström et al. 2003: r50 burst pairing",
         members=[("sjostrom_burst5x20hz_r50_dt-120ms", "control", "sjostrom_burst5x20hz_r50_dt-120ms",
                   "5$\\times$20 Hz burst, $\\Delta t$ $-$120 ms")]),
    dict(name="sj07step", path="L5", pair=L5FIX, title="L5$\\to$L5, Sjöström et al. 2007: step pairing",
         members=[("sjostrom07_step200ms_pair", "control", "sjostrom07_step200ms_pair", "step pairing, control"),
                  ("sjostrom07_step200ms_pair", "mglu_block", "sjostrom07_step200ms_pair", "step pairing, mGluR block")]),
    dict(name="letzkus", path="L23", pair=None, fallback="10149-186264",
         title="L2/3$\\to$L5, Letzkus et al. 2006: distal inputs",
         members=[("letzkus_1ap_dt+10ms", "control", "letzkus_1ap_dt+10ms", "1 AP, $\\Delta t$ +10 ms"),
                  ("letzkus_3ap_200hz_dt+10ms", "control", "letzkus_3ap_200hz_dt+10ms@distal",
                   "3 AP 200 Hz, $\\Delta t$ +10 ms (distal)"),
                  ("letzkus_3ap_200hz_dt+10ms", "nmdar_block", "letzkus_3ap_200hz_dt+10ms@distal",
                   "3 AP 200 Hz, $\\Delta t$ +10 ms, NMDAR block")]),
]


def _njit():
    from numba import njit

    @njit(cache=False)
    def v5_rec_traces(E, S, hs, cnt, td, tp, rho0, gd, gp, rs, kt, thV, isc, tauE1, dmin, Ae, dp0, thE, tauD, weighted,
                      o_rho, o_d, o_V, o_W, o_ev):
        """compare_prefire_v5.v5_rec, same steps and arithmetic, with the state after each grid step k stored in
        o_*[i, k] (rho, dpre, pool V, eCB trigger pool W; o_ev = own arrivals whose eCB step was applied at k)."""
        n, T = E.shape
        if thE <= 0.0:
            thE = thV
        for i in range(n):
            a_td = td[i]; a_tp = tp[i]
            r = rho0[i]; d0 = dp0
            V = 0.0; chV = -1.0; aV = 1.0; bV = 0.0
            W = 0.0; Bg = 0.0; chD = -1.0; aD = 1.0
            for k in range(T):
                e = E[i, k] * 1.0
                s = S[i, k] * 1.0
                h = hs[k]
                hm = h * kt
                c = cnt[i, k]
                pot = 1.0 if e > a_tp else 0.0
                dep = 1.0 if e > a_td else 0.0
                if pot > 0.0 and thV > 0.0 and not (V > thV):
                    pot = 0.0
                r = r + h * (-r * (1 - r) * (rs - r) + pot * gp * (1 - r) - dep * (1 - pot) * gd * r)
                if r < 0.0:
                    r = 0.0
                elif r > 1.0:
                    r = 1.0
                if weighted:
                    if c > 0.0:
                        Bg = Bg + c
                    trig = W > thE
                else:
                    trig = V > thE
                ev = 0.0
                if c > 0.0 and Ae > 0.0 and trig:
                    d0 = dmin + (d0 - dmin) * (1.0 - Ae) ** c
                    ev = c
                if hm != chV:
                    chV = hm
                    aV = np.exp(-hm / tauE1); bV = tauE1 * (1.0 - aV) / isc
                if weighted:
                    wb = 1.0 - (Bg if Bg < 1.0 else 1.0)
                    W = aV * W + bV * (wb * s)
                    if hm != chD:
                        chD = hm
                        aD = np.exp(-hm / tauD)
                    Bg = Bg * aD
                V = aV * V + bV * s
                o_rho[i, k] = r; o_d[i, k] = d0; o_V[i, k] = V; o_W[i, k] = W; o_ev[i, k] = ev
    return v5_rec_traces


# ---------------------------------------------------------------- worker (one BCL run, own process)
def make_hook():
    def hook(task, sim, cell, syns, gids, meta, h):
        tf = np.arange(task["win"][0], min(task["win"][1], meta["tstop"]), task["dt_fine"])
        tc = np.arange(0.0, meta["tstop"], task["dt_coarse"])
        hv = dict(f=h.Vector(tf), c=h.Vector(tc))
        rec = {}

        def vrec(ref, which):
            v = h.Vector(); v.record(ref, hv[which]); return v
        rec["v_soma_f"] = [vrec(cell.somatic[0](0.5)._ref_v, "f")]   # cell.somatic as live_v5 (fixhp)
        for k in FINE:
            rec[k + "_f"] = [vrec(getattr(s.hsynapse, "_ref_" + k), "f") for _, s in syns]
        for k in COARSE:
            rec[k + "_c"] = [vrec(getattr(s.hsynapse, "_ref_" + k), "c") for _, s in syns]

        def fin(out, _keep=hv):            # _keep: the h.Vector time bases must outlive the hook, else NEURON records 0 samples
            mf = min(len(v) for k, vs in rec.items() if k.endswith("_f") for v in vs)
            mc = min(len(v) for k, vs in rec.items() if k.endswith("_c") for v in vs)
            assert mf > 0 and mc > 0, f"live recording empty (n_tf {mf}, n_tc {mc})"
            arr = {k: np.array([np.asarray(v.to_python()[:(mf if k.endswith("_f") else mc)], np.float64) for v in vs], np.float32)
                   for k, vs in rec.items()}
            arr["v_soma_f"] = arr["v_soma_f"][0]
            np.savez_compressed(task["live_npz"], tf=tf[:mf], tc=tc[:mc], syn=np.asarray(gids, np.int64), **arr)
            out.update(live_npz=task["live_npz"], n_tf=int(mf), n_tc=int(mc), n_tf_req=len(tf), n_tc_req=len(tc))
        return fin
    return hook


def run_one(task, fit):
    res = L.run_task((task, fit), hook=make_hook())
    res.pop("globals", None)
    print("RESULT " + json.dumps(res), flush=True)


def _sub(task, fit, timeout):
    try:
        p = subprocess.run([sys.executable, os.path.abspath(__file__), "--fit", fit, "--one", json.dumps(task)],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=timeout)
        res = [l[len("RESULT "):] for l in p.stdout.splitlines() if l.startswith("RESULT ")]
        if res:
            return json.loads(res[-1])
        return dict(task, ok=False, error=f"exit {p.returncode}, no result", log_tail=p.stdout[-3000:])
    except subprocess.TimeoutExpired:
        return dict(task, ok=False, error=f"timeout {timeout} s")


# ---------------------------------------------------------------- main process (offline, selection, merge)
def params(CP, fit):
    fa = fit["args"]
    CP.batch_v2.BAP_GATE = fit.get("bap_gate")
    return {**CP.MV.DEFAULTS, **CP.V5_DEFAULTS, **json.loads(fa["filters"]), **json.loads(fa.get("set", "{}")),
            **fit["pre"]}


def fit_targets(fit):
    """{(path, target, cond): (mean, sem, pred)} from the fit csvs (pred = offline over all the fit's pairs)."""
    import pandas as pd
    T = {}
    for suf, path in (("", "L5"), ("_l23", "L23")):
        f = fit["args"]["save"] + suf + ".csv"
        f = f if os.path.isabs(f) else os.path.join(L.ROOT, f)
        for _, r in pd.read_csv(f).iterrows():
            T[(path, r.target, r.condition)] = (float(r.target_mean), float(r.target_sem), float(r.pred))
    return T


def tstop_of(path, pair, proto):
    f = os.path.join(L.PATHS[path]["sims"], pair, proto, "prefire_simulation_config.json")
    return float(json.load(open(f))["run"]["tstop"])


def has_all(path, pair, protos):
    return all(L.find_record(path, pair, p) and os.path.isfile(
        os.path.join(L.PATHS[path]["sims"], pair, p, "prefire_simulation_config.json")) for p in protos)


def candidates(fit, grp):
    import pandas as pd
    protos = sorted({m[0] for m in grp["members"]})
    if grp["path"] == "L5":
        ps = fit["args"]["pairs"].split(",")
    else:
        g = pd.read_csv(L.GEOM); g = g[g.all_protocols & g.letzkus_distal.astype(bool)]
        ps = list(g.pregid.astype(str) + "-" + g.postgid.astype(str))
    return [p for p in ps if has_all(grp["path"], p, protos)]


def select_pair(CP, fit, P, grp, T):
    path = grp["path"]; protos = sorted({m[0] for m in grp["members"]}); conds = sorted({m[1] for m in grp["members"]})
    cands = candidates(fit, grp); score = {}
    for i in range(0, len(cands), 8):
        chunk = cands[i:i + 8]
        off = CP.offline_v5(path, fit, P, set(chunk), set(protos), conds)
        for p in chunk:
            s = 0.0
            for proto, cond, tgt, _ in grp["members"]:
                o = off.get((p, proto, cond)); m, sem, pred = T[(path, tgt, cond)]
                s += abs(float(np.mean(o["ratio"])) - pred) / sem if o is not None else np.inf
            score[p] = s
        del off; gc.collect()
    best = min(score, key=score.get)
    print(f"[traces_v5] {grp['name']}: {len(cands)} candidates, best {best} (score {score[best]:.3f}; "
          f"median {np.median(list(score.values())):.3f})", flush=True)
    return best, len(cands), score[best]


def load_record(CP, path, pair, proto):
    bv = CP.batch_v2; bv.BASIS_DIR = CP.C.BASIS[path]
    B = bv.BatchV2(L.PATHS[path]["dirs"], protocols=[proto], pairs={pair}, fast=False, signals=("vdcc",))
    rs = [r for r in B.recs if r["pair"] == pair and r["proto"] == proto]
    assert len(rs) == 1, f"{path} {pair} {proto}: {len(rs)} records"
    return rs[0]


def kernel(kfun, CP, fit, P, r, Ae):
    rho_v4, bv = CP.rho_v4, CP.batch_v2
    n, T = r["effcai"].shape; t = np.asarray(r["t"], float)
    hs = np.zeros(T); hs[:T - 1] = np.diff(t) / 1000.0 / bv.TAU_IND_GB
    cnt = np.zeros((n, T)); np.add.at(cnt, (np.broadcast_to(np.arange(n)[:, None], r["arr"].shape), r["arr"]), 1.0)
    td, tp = rho_v4.thetas(r["c_pre"], r["c_post"], fit["a"], 1.0)
    gd, gp = rho_v4.rates(P); thV = rho_v4.vamp(P)[1]
    o = [np.empty((n, T), np.float32) for _ in range(5)]
    kfun(np.ascontiguousarray(r["effcai"], dtype=np.float32), np.ascontiguousarray(r["vdcc"], dtype=np.float32),
         hs, cnt, np.asarray(td, np.float64), np.asarray(tp, np.float64), r["rho0"].astype(np.float64), gd, gp,
         float(bv.RHO_STAR_GB), 1e3 * float(bv.TAU_IND_GB), thV, float(P["i_scale"]), float(P["tau_E1"]),
         float(P["dpre_min"]), float(Ae), 0.0 + float(P["dpre0"]), float(P["theta_eCB"]), float(P["tau_d_NMDA"]),
         int(P["v5_mode"]) == 2, *o)
    return o, np.asarray(td, float), np.asarray(tp, float)


def rep_window(r, ev_t, tstop):
    t = np.asarray(r["t"], float); T = len(t)
    a = np.asarray(r["arr"])[0]; a = a[(a >= 0) & (a < T - 1)]
    ta = np.unique(t[a])
    if len(ta) == 0:
        return 0.0, min(tstop, 1000.0), 0.0, [], 0, 0
    gaps = np.diff(ta)
    grp = np.split(ta, np.flatnonzero(gaps > 0.5 * gaps.max()) + 1) if len(gaps) else [ta]
    k = len(grp) // 2
    if len(ev_t):
        e0 = float(np.min(ev_t))
        hit = [j for j, g in enumerate(grp) if g[0] - 1.0 <= e0 <= g[-1] + 1.0]
        k = hit[0] if hit else k
    g = grp[k]
    w0 = max(0.0, g[0] - PAD0); w1 = min(tstop, g[-1] + PAD1, w0 + WMAX)
    return float(w0), float(w1), float(g[0]), [float(x) for x in g], len(grp), k


def offline_part(kfun, CP, fit, P, grp, pair, proto, conds):
    """-> {cond: (arrays, meta)} for the members of (pair, proto); window from the control kernel's eCB steps."""
    path = grp["path"]; r = load_record(CP, path, pair, proto)
    t = np.asarray(r["t"], float); n, T = r["effcai"].shape
    t_off = np.append(t[1:], t[-1])
    tstop = tstop_of(path, pair, proto)
    dc = float(max(1.0, np.ceil(tstop / N_COARSE))); tc = np.arange(0.0, tstop, dc)
    ctrl, td, tp = kernel(kfun, CP, fit, P, r, P["A_eCB"])
    ii, kk = np.nonzero(ctrl[4] > 0)
    w0, w1, trep, rep_arr, n_rep, k_rep = rep_window(r, t[kk], tstop)
    iw = np.flatnonzero((t >= w0) & (t <= w1))
    E = np.asarray(r["effcai"], np.float32); S = np.asarray(r["vdcc"], np.float32)
    offc = CP.offline_v5(path, fit, P, {pair}, {proto}, sorted(set(conds)))
    out = {}
    for cond in conds:
        if cond == "nmdar_block":     # gpu_v5_rho lane: rho frozen at rho0, d 0 (no kernel traces)
            o = None
        elif cond == "mglu_block":    # no eCB step (A_eCB 0); rho as control
            o = kernel(kfun, CP, fit, P, r, 0.0)[0]
        else:
            o = ctrl
        a = dict(syn=np.asarray(r["syn"], np.int64), rho0=np.asarray(r["rho0"], float), td=td, tp=tp,
                 c_pre=np.asarray(r["c_pre"], float), c_post=np.asarray(r["c_post"], float),
                 off_maxE=E.max(axis=1), off_maxV=ctrl[2].max(axis=1), off_maxW=ctrl[3].max(axis=1),
                 off_t_in=t[iw], off_E_f=E[:, iw], off_S_f=S[:, iw], tc_off=tc)
        if o is None:
            a.update(off_rho_c=np.tile(a["rho0"][:, None], (1, len(tc))).astype(np.float32),
                     off_d_c=np.zeros((n, len(tc)), np.float32), off_ev_syn=np.zeros(0, np.int64),
                     off_ev_t=np.zeros(0), off_rho_end=a["rho0"].copy(), off_d_end=np.zeros(n))
        else:
            ei, ek = np.nonzero(o[4] > 0)
            a.update(off_t_f=t_off[iw], off_rho_f=o[0][:, iw], off_d_f=o[1][:, iw], off_V_f=o[2][:, iw],
                     off_W_f=o[3][:, iw],
                     off_rho_c=np.stack([np.interp(tc, t_off, o[0][i]) for i in range(n)]).astype(np.float32),
                     off_d_c=np.stack([np.interp(tc, t_off, o[1][i]) for i in range(n)]).astype(np.float32),
                     off_ev_syn=ei.astype(np.int64), off_ev_t=t[ek], off_rho_end=o[0][:, -1].astype(float),
                     off_d_end=o[1][:, -1].astype(float))
        oc = offc[(pair, proto, cond)]
        m = dict(w0=w0, w1=w1, t_rep0=trep, rep_arrivals=rep_arr, n_rep=n_rep, k_rep=k_rep, dt_coarse=dc,
                 tstop=tstop, ratio_off=float(np.mean(oc["ratio"])),
                 copy_maxdiff_rho=float(np.max(np.abs(a["off_rho_end"] - np.asarray(oc["rho"], float)))),
                 copy_maxdiff_d=float(np.max(np.abs(a["off_d_end"] - np.asarray(oc["dpre"], float)))),
                 off_frozen=o is None, n_off_ecb=int(len(a["off_ev_t"])), n_grid=int(T))
        out[cond] = (a, m, oc)
    del ctrl, r, offc; gc.collect()
    return out


def merge(CP, P, a, m, oc, res):
    """Live npz (reordered to the record's synapses) + offline arrays -> final arrays; live ratio via the basis."""
    rho_v4 = CP.rho_v4; sigma = rho_v4.opts(P)[2]
    z = np.load(res["live_npz"]); ls = list(np.asarray(z["syn"]))
    idx = np.array([ls.index(g) if g in ls else -1 for g in a["syn"]])

    def ro(x):
        y = np.asarray(x, np.float32)[np.clip(idx, 0, None)]
        y[idx < 0] = np.nan
        return y
    for k in z.files:
        if k in ("tf", "tc", "v_soma_f", "syn"):
            a["live_" + k] = np.asarray(z[k])
        else:
            a["live_" + k] = ro(z[k])
    rl = ro(res["end"]["rho_GB"]).astype(float); dl = ro(res["end"]["dpre_GB"]).astype(float)
    a["live_rho_end"], a["live_d_end"] = rl, dl
    m["ratio_live"] = float(np.mean(rho_v4.ratio(oc["b"], oc["rho0"], rl, dl, sigma)))
    m["rho_end_maxdiff_live_off"] = float(np.nanmax(np.abs(rl - a["off_rho_end"])))
    m["d_end_maxdiff_live_off"] = float(np.nanmax(np.abs(dl - a["off_d_end"])))
    m["n_live_ecb"] = int(np.nansum(a["live_necb_GB_c"][:, -1])) if a["live_necb_GB_c"].size else 0
    m["n_post_ind"] = res.get("n_post_ind"); m["live_wall_s"] = res.get("wall_s")
    m["n_missing_live"] = int((idx < 0).sum())
    # time-base check: peak time of the effective Ca of the synapse with the largest peak, live vs record
    if "off_t_in" in a and a["off_E_f"].size and a["live_effcai_GB_f"].size:
        i = int(np.argmax(a["off_maxE"]))
        m["t_peak_offset_ms"] = float(a["live_tf"][np.nanargmax(a["live_effcai_GB_f"][i])]
                                      - a["off_t_in"][np.argmax(a["off_E_f"][i])])
    os.remove(res["live_npz"])
    return a, m


def kernel_copy_test(kfun, CP):
    rng = np.random.default_rng(0); n, T = 3, 3000
    E = rng.random((n, T)).astype(np.float32); S = (rng.random((n, T)) * 2e-5).astype(np.float32)
    hs = np.full(T, 0.25 / 1000 / 70.0); hs[-1] = 0.0
    cnt = (rng.random((n, T)) < 0.01).astype(float)
    td = np.full(n, 0.3); tp = np.full(n, 0.7); r0 = np.array([0.0, 1.0, 0.0])
    worst = 0.0
    for w in (True, False):
        for thE in (0.0, 2.0):
            ro_, do_ = np.empty(n), np.empty(n)
            CP.v5_rec(E, S, hs, cnt, td, tp, r0, 35.0, 507.0, 0.5, 7e4, 1.0, 1e-5, 100.0, -0.29, 1.0, 0.0, thE, 70.0,
                      w, ro_, do_)
            o = [np.empty((n, T), np.float32) for _ in range(5)]
            kfun(E, S, hs, cnt, td, tp, r0, 35.0, 507.0, 0.5, 7e4, 1.0, 1e-5, 100.0, -0.29, 1.0, 0.0, thE, 70.0, w, *o)
            worst = max(worst, float(np.abs(o[0][:, -1] - ro_).max()), float(np.abs(o[1][:, -1] - do_).max()))
    assert worst < 1e-5, f"v5_rec_traces != v5_rec (max {worst})"
    print(f"[traces_v5] check: v5_rec_traces end state = v5_rec (max diff {worst:.2e})", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", default=FIT)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--timeout", type=float, default=1500.0)
    ap.add_argument("--check", action="store_true", help="build the task list, test the kernel copy, exit")
    ap.add_argument("--one", default=None, help="internal: one BCL run (json task)")
    a = ap.parse_args()
    if a.one:
        run_one(json.loads(a.one), os.path.abspath(a.fit)); return
    t0 = time.time()
    import compare_prefire_v5 as CP                       # env ANALYTICAL_BASIS_DIR / L23_BASIS_DIR set by the job
    fit = json.load(open(a.fit)); P = params(CP, fit); T = fit_targets(fit)
    kfun = _njit()
    missing = [(g["path"], m[2], m[1]) for g in GROUPS for m in g["members"] if (g["path"], m[2], m[1]) not in T]
    assert not missing, f"targets not in the fit csvs: {missing}"
    if a.check:
        kernel_copy_test(kfun, CP)
        n = 0
        for g in GROUPS:
            protos = sorted({m[0] for m in g["members"]})
            if g["pair"]:
                assert has_all(g["path"], g["pair"], protos), f"{g['name']}: {g['pair']} lacks a record / workdir"
                pp = g["pair"]
            else:
                c = candidates(fit, g); assert c, f"{g['name']}: no candidate pair"
                pp = f"auto ({len(c)} candidates, fallback {g['fallback']})"
            for m in g["members"]:
                n += 1
                print(f"[traces_v5] check {g['name']}: {g['path']} {pp} {m[0]} {m[1]} target {m[2]} "
                      f"(data {T[(g['path'], m[2], m[1])][0]:.2f}, pred {T[(g['path'], m[2], m[1])][2]:.3f})", flush=True)
        import plot_traces_v5                              # noqa: F401  (import test of the plot script)
        print(f"[traces_v5] check OK: {n} records, out {a.out}", flush=True)
        return
    os.makedirs(os.path.join(a.out, "live"), exist_ok=True)
    plan = []
    for g in GROUPS:
        pair, nc, sc = g["pair"], None, None
        if pair is None:
            try:
                pair, nc, sc = select_pair(CP, fit, P, g, T)
            except Exception as e:
                pair = g["fallback"]; print(f"[traces_v5] {g['name']}: selection failed ({e}); fallback {pair}", flush=True)
        plan.append((g, pair, nc, sc))
    from concurrent.futures import ThreadPoolExecutor
    ex = ThreadPoolExecutor(a.workers); jobs = []
    for g, pair, nc, sc in plan:
        for proto in dict.fromkeys(m[0] for m in g["members"]):
            mem = [(j, m) for j, m in enumerate(g["members"]) if m[0] == proto]
            off = offline_part(kfun, CP, fit, P, g, pair, proto, [m[1] for _, m in mem])
            for j, (_, cond, tgt, lab) in mem:
                arr, meta, oc = off[cond]
                ex_, sem, pred = T[(g["path"], tgt, cond)]
                meta.update(group=g["name"], group_title=g["title"], member=j, label=lab, path=g["path"], pair=pair,
                            proto=proto, cond=cond, target=tgt, exp=ex_, sem=sem, pred=pred, n_candidates=nc,
                            select_score=sc, theta_V=float(CP.rho_v4.vamp(P)[1]),
                            theta_eCB=float(P["theta_eCB"]) if float(P["theta_eCB"]) > 0 else float(CP.rho_v4.vamp(P)[1]),
                            rho_star=float(CP.batch_v2.RHO_STAR_GB), dpre_min=float(P["dpre_min"]))
                key = f"{g['name']}__{j}_{pair}__{proto}__{cond}"
                task = dict(path=g["path"], pair=pair, proto=proto, cond=cond, phase="prefire", express="cooker",
                            win=[meta["w0"], meta["w1"]], dt_fine=DT_FINE, dt_coarse=meta["dt_coarse"],
                            live_npz=os.path.join(a.out, "live", key + ".npz"))
                print(f"[traces_v5] offline {key}: ratio_off {meta['ratio_off']:.3f} (pred {pred:.3f}, data {ex_:.2f}) "
                      f"window {meta['w0']:.0f}-{meta['w1']:.0f} ms (rep {meta['k_rep']}/{meta['n_rep']}), "
                      f"copy diff rho {meta['copy_maxdiff_rho']:.1e} d {meta['copy_maxdiff_d']:.1e}, "
                      f"{meta['n_off_ecb']} eCB steps; t {time.time() - t0:.0f}s", flush=True)
                jobs.append((key, arr, meta, oc, ex.submit(_sub, task, os.path.abspath(a.fit), a.timeout)))
            del off; gc.collect()
    summ = []
    for key, arr, meta, oc, fut in jobs:
        res = fut.result()
        if not res.get("ok"):
            print(f"[traces_v5] FAILED {key}: {res.get('error')}\n{res.get('traceback', res.get('log_tail', ''))}",
                  flush=True)
            meta.update(ok=False, error=res.get("error")); summ.append(meta); continue
        arr, meta = merge(CP, P, arr, meta, oc, res); meta["ok"] = True
        np.savez_compressed(os.path.join(a.out, key + ".npz"), meta=np.array(json.dumps(meta)), **arr)
        summ.append(meta)
        print(f"[traces_v5] done {key}: ratio data {meta['exp']:.2f}+/-{meta['sem']:.2f} pred {meta['pred']:.3f} "
              f"off {meta['ratio_off']:.3f} BCL {meta['ratio_live']:.3f}; |rho_live - rho_off| "
              f"{meta['rho_end_maxdiff_live_off']:.3f}, |d| {meta['d_end_maxdiff_live_off']:.3f}, eCB live "
              f"{meta['n_live_ecb']} off {meta['n_off_ecb']}, t_peak offset {meta.get('t_peak_offset_ms')} ms, "
              f"BCL {meta['live_wall_s']:.0f}s", flush=True)
    ex.shutdown()
    try:
        os.rmdir(os.path.join(a.out, "live"))
    except OSError:
        pass
    json.dump(summ, open(os.path.join(a.out, "summary.json"), "w"), indent=1)
    print(f"[traces_v5] {sum(m.get('ok', False) for m in summ)}/{len(summ)} records ok, {time.time() - t0:.0f}s; "
          f"wrote {a.out}", flush=True)
    if not any(m.get("ok", False) for m in summ):
        sys.exit("[traces_v5] no record ok")


if __name__ == "__main__":
    main()
