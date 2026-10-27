"""Per-step state traces of the 4V rule (gpu_v12_rho) for the trace panels of figures LB..LF (FIG_SPEC.md).

v12_tr is a copy of compare_prefire_v12.v12_rec (the CPU port of gpu_v12_rho._make_kernel_v10.kern; NEURON GluSynapseV12 =
port on 280/280 synapses) with recording statements added and nothing else: same arithmetic, same step order. For the
synapses with slot[i] >= 0 it records per step (TR_FIELDS):
  c        effcai as used (c*), compared with theta_d / theta_p (per synapse, rho_v4.thetas of c_pre / c_post);
  vdcc     own VDCC current (record step mean, S), the peak-veto statistic;
  shaft_dCa  shaft dCa as the kernel reads it (X, uM), compared with theta_G;
  L_G      shaft licence (t_k - t_last(X > theta_G) <= gate_win);   L_C  C-licence incl. the hold (c_lic > 0, else NaN);
  pot_raw / pot / dep  potentiation flag before the licences / after them / depression flag;
  pool     the eCB trigger pool at the start of the step (Vg for v5_mode 1, W for v5_mode 2); threshold thE = k_E P1 x uE;
  N        NO drive used in the step's d_NO integral;  d_NO, d_eCB, rho  after the step;  clk  kernel clock (ms) at the step.
and per own arrival (EV_FIELDS): sample k, clock, t_a (= clk - off), count, off, first-over-threshold offset
tq, thP = veto_peak_k spk_i, C (qa > thP), trigger value (pool - f x up) and its threshold, trig, veto, annulled, eCB step
taken, d_eCB after, L_G / L_C / pot at the step, tAC before the step, scheduled C time t_a + tq (NaN if no C).
The lanes are formed from the traces as compare_prefire_v12.offline_v12: control / post_nmdar d = clip(d_eCB + d_NO,
dpre_min, d_NO_max), mglu_block d = clip(dpre0 + d_NO), no_block d = d_eCB, nmdar_block d = 0 (no NO arm: d = d_eCB / dpre0).

Inputs are prepared exactly as offline_v12 (records of the s2g0321 tree through batch_v2.BatchV2, shaft rest and X
quantisation, arrival fractions, uE map, thetas); spk_i = spk_v12.spk_from_recs over every record of the pair in the
model's target protocols (load_targets of the fit model's groups; the GPU fit's batch for that pair) + the traced ones,
NaN -> median of that pair batch (kernel: median over all models; flagged per synapse as spk_tier 0).

Built-in check: every traced record is also run through compare_prefire_v12.v12_rec on the same inputs; the end state of
every synapse (rho, d, necb, nveto, nsh, tsh, dn, nno, nar, nC, nann, ncb) must be identical (max |diff| <= 1e-12) and the
last traced rho / d_eCB / d_NO must equal the end state. Printed per record; any failure -> exit 3 after all records.

Choice by contrast (selection_lat.json: choice dict per trace, population numbers under "stats"; pop_<path>.csv = end
states of every synapse of every scanned record): every synapse of every fit pair of the protocol (perpair csv) is run
through v12_rec under the mechanism-on fit and its knock-out; the contrast is the end-state difference; among synapses
where the mechanism acts (contrast > 0.01 and the figure's mask) the one closest to the MEDIAN contrast is traced.
  LB L2/3->L5 sjostrom_50hz_dt+10ms, V2B4 vs KO_lic (population from pop_L23.csv when present, no rescan): two synapses
     traced, preferably of one pair and rho0 = 0: 'closed' (KO_lic rho_end > 0.8, V2B4 rho_end < 0.3) and 'open' (both
     > 0.8, shaft dCa crosses theta_G), chosen by score (fig_LB doc), not by the median rule; train shown = the one with
     the most blocked / licensed pot time. Population: classes KO pot & B4 dep / both pot / both dep / other with the mean
     fraction of pot > 0 steps blocked; the earlier per-row lines (50 Hz all / @distal, 3AP @proximal / @distal) kept;
  LC L5 sjostrom_0.1hz_dt-10ms vs dt+10ms, same synapse: d_eCB_end(+10) - d_eCB_end(-10), stepping at -10 only; arrival =
     first stepping arrival at -10 and the same arrival index at +10; lanes control / mglu_block;
  LD L5 sjostrom_20hz_dt+10ms: |d_eCB_end(KO_veto) - d_eCB_end(V2B4)|, vetoed under V2B4 (preferably also stepping at
     dt-10ms); arrival = first vetoed triggering arrival (same index under KO_veto) + first stepping arrival at dt-10ms;
  LE L5 sjostrom07_step200ms_pair: |d_end control(V2B4) - d_end(KO_no)| with NO events; V2B4, V2B4_019, KO_no traced, plus
     the pair EPSP-ratio readout applied along the induction (proxy; control and no_block lanes);
  LF L5 sjostrom_50hz_dt-10ms: |d_eCB_end(KO_b4) - d_eCB_end(V2B4)| with an annulment under V2B4; Zilberter L2/3->L2/3
     train10 -10 last and post_only pooled: |rho_end(V2B1B4) - rho_end(V2B1B4_clic0)| (variants/V2B1B4_clic0.json =
     V2B1B4 with c_lic 0).
Variants (KO_lic, KO_b4, ...) are found at run time as /scratch/dhuruva/latest_figs/variants/*<name>*.json; missing ones
are listed and skipped. Variant jsons carry no v8 / v10 / v12 blocks: theta_G = pre gate_theta, switches from args.set;
lic_src 0 (KO_lic, latest_figs/SWITCH_MAP.md) = licence off, run as theta_G 0 (the only place the base licence acts in
v12_rec); this is the one departure from offline_v12's asserts, which accept gate_src 2 only.

    <interpreter> /scratch/dhuruva/latest_figs/traces/traces_v12_lat.py [--out DIR] [--figs LB,LC,LD,LE,LF]
(env as run_valid_v12_cmp.sh: BCL_EMODEL, ANALYTICAL_BASIS_DIR, L23_BASIS_DIR, L23L23_BASIS_DIR, CEXP_DIR)
"""
import argparse, gc, glob, json, os, sys, time
import numpy as np

HERE = "/lustre09/project/6070394/dhuruva/plastyfire/glusynapse_v2/bcl_validation"                 # read-only imports from bcl_validation (this copy lives on scratch)
sys.path.insert(0, HERE)
import compare_prefire_v12 as C12                         # noqa: E402  (imports live_v12 -> s2g0321 record dirs)
from compare_prefire_v5 import V5_DEFAULTS                # noqa: E402
from numba import njit                                    # noqa: E402

C8, C7, C, MV, rho_v4, batch_v2, PATHS = C12.C8, C12.C7, C12.C, C12.MV, C12.rho_v4, C12.batch_v2, C12.PATHS
live_v12, spk_v12 = C12.live_v12, C12.spk_v12
V2 = os.path.dirname(HERE)
sys.path.insert(0, V2)
from targets import load_targets                          # noqa: E402

PP = "/scratch/dhuruva/perpair"
FITS = dict(V2B4=f"{PP}/w_z10b4n14_u47_fit.json", V2B1B4=f"{PP}/w_z10bbn14_u45_fit.json",
            V2B4_019=f"{PP}/w_z9b4dNO_s42_fit.json")
VARDIR = "/scratch/dhuruva/latest_figs/variants"
VARIANTS = ("KO_lic", "KO_ecb", "KO_veto", "KO_no", "KO_b4", "CHP", "V2B1B4_clic0")
OUT = "/scratch/dhuruva/latest_figs/traces"
CEXP = os.environ.get("CEXP_DIR", "/scratch/dhuruva/s2g0321/cexp")
MODEL = {"L5": "l5", "L23": "l23", "L23L23": "l23l23"}      # path -> fit model name / per-pair pathway
VALCSV = {"L5": "_val.csv", "L23": "_val_l23.csv", "L23L23": "_val_l23l23.csv"}
CEXPCSV = {"L5": "L5L5.csv", "L23": "L23L5.csv", "L23L23": "L23L23.csv"}
TR_FIELDS = ("c", "vdcc", "shaft_dCa", "L_G", "L_C", "pot_raw", "pot", "dep", "pool", "N", "d_NO", "d_eCB", "rho", "clk")
EV_FIELDS = ("k", "clk", "t_a", "c", "off", "tq", "thP", "C", "trigval", "thE", "trig", "veto", "annul", "step",
             "d_eCB", "L_G", "L_C", "pot", "tAC_prev", "t_C")
OUTK = ("rho", "d", "ne", "nv", "nsh", "tsh", "dn", "nno", "nar", "nC", "nann", "ncb")
TOL = 1e-12
FAIL = []


@njit(cache=False)
def v12_tr(E, S, X, F, hs, cnt, td, tp, rho0, gd, gp, rs, kt, thG, win, tex, isc, tauE1, dmin, Ae, dp0, thE, tauD, Q,
           veto, trigw, vpk, vt0, vtv, thP, cl, cand, winc, ann, tm, nom, Ano, tauN, dnmax,
           slot, TR, EV, NEV,
           out_rho, out_d, out_necb, out_nveto, out_nsh, out_tsh, out_dn, out_nno, out_nar, out_nC, out_nann, out_ncb):
    """compare_prefire_v12.v12_rec + recording (lines marked # trace). slot (n,) int64: trace row of synapse i or -1;
    TR (nsel, len(TR_FIELDS), T), EV (nsel, maxA, len(EV_FIELDS)) float64 (NaN-filled), NEV (nsel,) arrivals recorded."""
    n, T = E.shape
    for i in range(n):
        si = slot[i]                                                      # trace
        ja = 0                                                            # trace
        a_td = td[i]; a_tp = tp[i]
        r = rho0[i]; d0 = dp0
        V = 0.0; chV = -1.0; aV = 1.0; bV = 0.0
        W = 0.0; Bg = 0.0; chD = -1.0; aD = 1.0
        ip = 0.0; up = 0.0
        clk = 0.0; tl = -1e30; lp = False; nsh = 0.0; tsh = 0.0
        ne = 0.0; nv = 0.0; na_ = 0.0
        Nn = 0.0; dn = 0.0; nno = 0.0; chN = -1.0; aN = 1.0; bN = 0.0
        nar = 0.0; nC = 0.0; ncb = 0.0
        tcl = -1e30; pC1 = 1e30; pC2 = 1e30; Lt = False; tAC = -1e30
        hpv = 0.0
        thEi = thE[i]; thPi = thP[i]
        for k in range(T):
            e = E[i, k] * 1.0
            s = S[i, k] * 1.0
            h = hs[k]
            hm = h * kt
            c = cnt[i, k]
            fx = F[i, k] if tex else 0.0
            qa = 0.0
            Cf = False; tan = clk
            r_off = 0.0; r_tq = 1e30; r_licC = np.nan; r_stp = False       # trace
            r_pool = W if trigw else V                                    # trace (pool at the step start)
            if c > 0.0:
                if vpk:                                      # 389-419: window (t_a + vt0, t_a + vtv], peak step mean
                    off = fx * hpv if tex else 0.0
                    tq = 1e30
                    if off > 0.0 and k > 0:
                        lo = vt0
                        hi = off if off < vtv else vtv
                        if hi > lo:
                            qa = S[i, k - 1] * 1.0
                            if qa > thPi:
                                tq = lo
                    cj = off; j = k
                    while j < T and cj < vtv:
                        hj = hs[j] * kt
                        lo = cj if cj > vt0 else vt0
                        hi = (cj + hj) if (cj + hj) < vtv else vtv
                        if hj > 0.0 and hi > lo:
                            sj = S[i, j] * 1.0
                            if sj > qa:
                                qa = sj
                            if sj > thPi and tq > 1e29:
                                tq = lo
                        cj += hj
                        j += 1
                    Cf = qa > thPi                           # 420-433
                    tan = clk - off
                    r_off = off; r_tq = tq                                # trace
                    if Cf and cl > 0:
                        tn = tan + tq
                        if pC1 > 1e29:
                            pC1 = tn
                        elif tn >= pC1:
                            if pC2 > 1e29 or tn > pC2:
                                pC2 = tn
                        else:
                            if pC2 > 1e29:
                                pC2 = pC1
                            pC1 = tn
                elif veto:
                    qa = Q[i, k]                             # v7 veto integral over [t_k, t_k + veto_T]
            if cl > 0:                                       # 446-450
                if pC1 <= clk:
                    tcl = pC1; pC1 = pC2; pC2 = 1e30
                    if pC1 <= clk:
                        tcl = pC1; pC1 = 1e30
            if c > 0.0:
                nar += c
                if Cf:
                    nC += c
            pot = 1.0 if e > a_tp else 0.0
            dep = 1.0 if e > a_td else 0.0
            lic = (clk - tl) <= win
            if lic and not lp:
                nsh += 1.0
            if lic:
                tsh += hm
            lp = lic
            pq = pot
            if pot > 0.0 and thG > 0.0 and not lic:
                pot = 0.0
            if cl > 0:                                       # 483-496
                licC = (clk - tcl) <= winc
                if cl == 2:
                    if licC:
                        Lt = True
                    elif Lt and e <= a_tp:
                        Lt = False
                    licC = licC or Lt
                r_licC = 1.0 if licC else 0.0                             # trace
                if not cand:
                    pot = pq
                if pot > 0.0 and not licC:
                    pot = 0.0
                    ncb += 1.0
            if trigw and c > 0.0:
                Bg = Bg + c
            r = r + h * (-r * (1 - r) * (rs - r) + pot * gp * (1 - r) - dep * (1 - pot) * gd * r)
            if r < 0.0:
                r = 0.0
            elif r > 1.0:
                r = 1.0
            if trigw:                                        # 528-531
                trig = (W - fx * ip) > thEi
                r_tv = W - fx * ip                                        # trace
            else:
                trig = (V - fx * up) > thEi
                r_tv = V - fx * up                                        # trace
            if tex and not vpk:
                qa = qa + fx * up
            if vpk:
                vet = qa > thPi                              # 539-540
            else:
                vet = veto and qa > thEi
            an = False
            r_tAC = tAC                                                   # trace
            if ann and c > 0.0 and (tan - tAC) <= tm:       # 547-549
                an = True
            if c > 0.0 and Ae > 0.0 and trig and not vet and not an:
                d0 = dmin + (d0 - dmin) * (1.0 - Ae) ** c
                ne += c
                r_stp = True                                              # trace
            elif c > 0.0 and Ae > 0.0 and trig:
                if vet:
                    nv += c
                else:
                    na_ += c
            if ann and c > 0.0 and Cf:
                tAC = tan
            if nom == 3 and c > 0.0 and pot > 0.0:          # 575-578
                Nn = Nn + c
                nno += 1.0
            if nom == 4 and c > 0.0 and Cf:                 # 579-582
                Nn = Nn + c
                nno += 1.0
            r_N = Nn                                                      # trace
            if hm != chV:
                chV = hm
                aV = np.exp(-hm / tauE1); bV = tauE1 * (1.0 - aV) / isc
            if trigw:
                wb = 1.0 - (Bg if Bg < 1.0 else 1.0)
                ip = bV * (wb * s)
                W = aV * W + ip
                if hm != chD:
                    chD = hm
                    aD = np.exp(-hm / tauD)
                Bg = Bg * aD
            up = bV * s
            V = aV * V + up
            if nom == 3 or nom == 4:                         # 631-638
                if hm != chN:
                    chN = hm
                    aN = np.exp(-hm / tauN); bN = tauN * (1.0 - aN) / kt
                dn = dn + Ano * Nn * bN
                if dn > dnmax:
                    dn = dnmax
                Nn = aN * Nn
            if si >= 0:                                                   # trace: per-step record
                TR[si, 0, k] = e; TR[si, 1, k] = s; TR[si, 2, k] = X[i, k]
                TR[si, 3, k] = 1.0 if lic else 0.0; TR[si, 4, k] = r_licC
                TR[si, 5, k] = pq; TR[si, 6, k] = pot; TR[si, 7, k] = dep; TR[si, 8, k] = r_pool
                TR[si, 9, k] = r_N; TR[si, 10, k] = dn; TR[si, 11, k] = d0; TR[si, 12, k] = r; TR[si, 13, k] = clk
                if c > 0.0 and ja < EV.shape[1]:                          # trace: per-arrival record
                    EV[si, ja, 0] = k; EV[si, ja, 1] = clk; EV[si, ja, 2] = tan; EV[si, ja, 3] = c
                    EV[si, ja, 4] = r_off; EV[si, ja, 5] = r_tq; EV[si, ja, 6] = thPi
                    EV[si, ja, 7] = 1.0 if Cf else 0.0; EV[si, ja, 8] = r_tv; EV[si, ja, 9] = thEi
                    EV[si, ja, 10] = 1.0 if trig else 0.0; EV[si, ja, 11] = 1.0 if vet else 0.0
                    EV[si, ja, 12] = 1.0 if an else 0.0; EV[si, ja, 13] = 1.0 if r_stp else 0.0
                    EV[si, ja, 14] = d0; EV[si, ja, 15] = 1.0 if lic else 0.0; EV[si, ja, 16] = r_licC
                    EV[si, ja, 17] = pot; EV[si, ja, 18] = r_tAC
                    EV[si, ja, 19] = (tan + r_tq) if Cf else np.nan
                    ja += 1
                    NEV[si] = ja
            if X[i, k] > thG:                                # 644-652 (gate_src 2, window)
                tl = clk
            clk += hm
            hpv = hm
        out_rho[i] = r
        out_d[i] = d0
        out_necb[i] = ne
        out_nveto[i] = nv
        out_nsh[i] = nsh
        out_tsh[i] = tsh
        out_dn[i] = dn
        out_nno[i] = nno
        out_nar[i] = nar
        out_nC[i] = nC
        out_nann[i] = na_
        out_ncb[i] = ncb


# ---------------------------------------------------------------- inputs (as compare_prefire_v12.offline_v12)
def fit_params(fit):
    fa = fit["args"]
    return {**MV.DEFAULTS, **V5_DEFAULTS, **json.loads(fa["filters"]), **json.loads(fa.get("set", "{}")), **fit["pre"]}


def set_basis(path):
    if path == "L23L23":
        C.BASIS["L23L23"] = os.environ.get("L23L23_BASIS_DIR") or C.BASIS.get("L23L23", "")
    batch_v2.BASIS_DIR = C.BASIS[path]


def model_groups(fit, path):
    for m in fit.get("models") or []:
        if m["name"] == MODEL[path]:
            return tuple(m["groups"].split(","))
    raise KeyError(f"fit has no model {MODEL[path]}")


def target_protos(fit, path):
    conds = set(fit["args"]["conditions"].split(","))
    gs = model_groups(fit, path)
    try:
        tg = load_targets(gs)
    except Exception:                                    # pathway groups that load_targets will not combine
        tg = {}
        for g in gs:
            tg.update(load_targets((g,)))
    return {pid.partition("@")[0] for (pid, cond) in tg if cond in conds}


_SPK = {}


def pair_spk(fit, path, pair, extra):
    """spk_i of every synapse of the pair: spk_v12.spk_from_recs over the pair's records of the model's target protocols
    (the GPU fit's batch for that pair; the traced protocols only if there are no targets), NaN -> median of this batch.
    One fixed set per (path, pair), so scans and traces use the same spk_i."""
    protos = frozenset(target_protos(fit, path) or set(extra))
    key = (path, pair, protos)
    if key not in _SPK:
        set_basis(path)
        B = batch_v2.BatchV2(PATHS[path]["dirs"], protocols=sorted(protos), pairs={pair}, fast=False, signals=("vdcc",),
                             verbose=False)
        tab = spk_v12.spk_from_recs(B.recs)
        spv, f, nn = spk_v12.fill(tab)
        _SPK[key] = (spv, {k: v[1] for k, v in tab.items()}, len(B.recs), f, nn)
        print(f"  spk {path} {pair}: {len(B.recs)} records, {len(tab)} synapses, NaN {nn} (filled {f:.4g} nA)", flush=True)
        del B
    return _SPK[key]


_REC = {}


def load_rec(path, pair, proto):
    key = (path, pair, proto)
    if key not in _REC:
        set_basis(path)
        B = batch_v2.BatchV2(PATHS[path]["dirs"], protocols=[proto], pairs={pair}, fast=False,
                             signals=("vdcc", "shaft_cai"), verbose=False)
        assert len(B.recs) == 1, f"{path} {pair} {proto}: {len(B.recs)} records"
        _REC[key] = B.recs[0]
    return _REC[key]


def drop_recs():
    _REC.clear()


def prep(fit, path, r, spv):
    """(args of v12_rec without the outputs, info) for one record, as offline_v12."""
    P = fit_params(fit)
    batch_v2.BAP_GATE = fit.get("bap_gate")
    o12 = live_v12.v12_opts(fit)
    nm = o12["no_mode"]
    v10 = fit.get("v10") or {}
    fa = fit["args"]
    ex = {**json.loads(fa["filters"]), **json.loads(fa.get("set", "{}")), **fit["pre"]}     # explicit switches only
    lic_src = int(ex.get("lic_src", v10.get("lic_src", 2)))
    assert int(ex.get("ecb_src", v10.get("ecb_src", 0))) == 0, "v12: ecb_src 0 only"
    win = float(P.get("gate_win", 100.0))
    if lic_src == 0:
        # variants/KO_lic (SWITCH_MAP: lic_src 0 -> kernel SRC 0, lic = True, pot never gated). In v12_rec the base
        # licence only acts through `pot > 0 and thG > 0 and not lic`, so theta_G 0 is the same rule; the traced L_G /
        # nsh / tsh are then only the shaft-window diagnostic (the kernel reports lic = True).
        thG = 0.0
    else:
        assert lic_src == 2 and int(P.get("gate_src", 0)) == 2 and win > 0, "v12: gate_src / lic_src 2 (shaft) or lic_src 0"
        thG = float(fit["v8"]["theta_G"]) if fit.get("v8") else float(P["gate_theta"])
    tex = bool(int(P.get("t_exact", 0)))
    v5m = int(P["v5_mode"]); assert v5m in (1, 2), f"v5_mode {v5m}"
    gamma, _, sigma = rho_v4.opts(P); assert gamma == 1.0
    gd, gp = rho_v4.rates(P); thV = rho_v4.vamp(P)[1]
    kt = 1e3 * float(batch_v2.TAU_IND_GB); rs = float(batch_v2.RHO_STAR_GB)
    Tv = float(P.get("veto_T", 0.0)); tauE1 = float(P["tau_E1"]); isc = float(P["i_scale"])
    vpk = bool(o12["veto_peak"])
    assert vpk or Tv > 0 or not tex, "t_exact needs veto_T > 0 or the peak veto"
    thE0 = float(P["theta_eCB"]); um = live_v12.L7.ue_map(fit, path)
    Ano = float(P.get("A_NO", 0.0)) if nm else 0.0
    tauN = float(P.get("tau_NO", 6.7)); dnmax = float(P.get("d_NO_max", 0.5)) if nm else 1e30
    dmin, dp0 = float(P["dpre_min"]), 0.0 + float(P["dpre0"])
    cand = o12["c_lic_and"] if o12["c_lic_and"] >= 0 else (1 if o12["c_lic"] == 1 else 0)
    n, T = r["effcai"].shape
    hs = np.zeros(T); hs[:T - 1] = np.diff(r["t"]) / 1000.0 / batch_v2.TAU_IND_GB
    cnt = np.zeros((n, T)); np.add.at(cnt, (np.broadcast_to(np.arange(n)[:, None], r["arr"].shape), r["arr"]), 1.0)
    td, tp = rho_v4.thetas(r["c_pre"], r["c_post"], fit["a"], 1.0)
    pre, post = map(int, r["pair"].split("-"))
    uE = np.ones(n) if um is None else np.array([um[(pre, post, int(s))] for s in r["syn"]])
    thE = uE * thE0 if thE0 > 0 else np.full(n, thV)
    Q = C7.veto_Q(r, cnt, Tv, tauE1, isc) if (Tv > 0 and not vpk) else np.zeros((n, T))
    sp = np.array([spv[(r["pair"], int(s))] for s in r["syn"]]) if vpk else np.ones(n)
    thP = o12["veto_peak_k"] * sp
    sh = np.asarray(r["shaft_cai"], np.float64)
    rest = C8.shaft_rest(r, sh)
    du = (sh - rest[:, None]) * 1e3
    X = np.clip(np.rint(du * 1e3), 0, 65535).astype(np.uint16).astype(np.float64) * 1e-3
    F = C8.arrival_F(r, n, T) if tex else np.zeros((1, 1))
    args = (np.ascontiguousarray(r["effcai"], dtype=np.float32), np.ascontiguousarray(r["vdcc"], dtype=np.float32),
            X, F, hs, cnt, np.asarray(td, np.float64), np.asarray(tp, np.float64), r["rho0"].astype(np.float64), gd,
            gp, rs, kt, thG, win, tex, isc, tauE1, dmin, float(P["A_eCB"]), dp0,
            np.ascontiguousarray(thE, np.float64), float(P["tau_d_NMDA"]), Q, Tv > 0, v5m == 2, vpk,
            o12["veto_t0"], o12["veto_Tv"], np.ascontiguousarray(thP, np.float64), o12["c_lic"], bool(cand),
            o12["c_hold"], bool(o12["c_annul"]), o12["c_Tm"] + C12.T_TOL, nm, Ano, tauN, dnmax)
    gth = float(P.get("gate_theta", thG))
    if lic_src != 0 and abs(gth - thG) > 1e-12:
        print(f"  WARNING: pre gate_theta {gth} != v8 theta_G {thG}; the port (and kernel) read v8 theta_G", flush=True)
    info = dict(td=np.asarray(td, float), tp=np.asarray(tp, float), thE=thE, thP=thP, spk=sp, uE=uE, rest=rest,
                thG=thG, win=win, gd=gd, gp=gp, dmin=dmin, dp0=dp0, dnmax=dnmax, nm=nm, Ano=Ano, tauN=tauN, Ae=float(P["A_eCB"]),
                vt0=o12["veto_t0"], vtv=o12["veto_Tv"], k_peak=o12["veto_peak_k"], c_lic=o12["c_lic"], cand=int(cand),
                c_hold=o12["c_hold"], c_annul=o12["c_annul"], c_Tm=o12["c_Tm"], v5_mode=v5m, t_exact=int(tex),
                veto_peak=int(vpk), theta_eCB=thE0, lic_src=lic_src, sigma=float(sigma), cnt=cnt)
    return args, info


def ref_run(args, n):
    o = {k: np.empty(n) for k in OUTK}
    C12.v12_rec(*args, *[o[k] for k in OUTK])
    return o


def trace_run(args, n, T, isel, cnt):
    isel = np.asarray(isel, np.int64)
    slot = np.full(n, -1, np.int64); slot[isel] = np.arange(len(isel))
    maxA = max(1, int(max(((cnt[i] > 0).sum() for i in isel), default=1)))
    TR = np.full((len(isel), len(TR_FIELDS), T), np.nan); EV = np.full((len(isel), maxA, len(EV_FIELDS)), np.nan)
    NEV = np.zeros(len(isel), np.int64)
    o = {k: np.empty(n) for k in OUTK}
    v12_tr(*args, slot, TR, EV, NEV, *[o[k] for k in OUTK])
    return TR, EV, NEV, o


def check(tag, ref, tr, TR, isel):
    d = {k: float(np.max(np.abs(ref[k] - tr[k]))) if len(ref[k]) else 0.0 for k in OUTK}
    last = dict(rho=TR[:, TR_FIELDS.index("rho"), -1], d=TR[:, TR_FIELDS.index("d_eCB"), -1],
                dn=TR[:, TR_FIELDS.index("d_NO"), -1])
    dl = {k: float(np.max(np.abs(v - ref[k][isel]))) if len(isel) else 0.0 for k, v in last.items()}
    mx = max(list(d.values()) + list(dl.values()))
    ok = mx <= TOL
    print(f"  PORT CHECK {tag}: max |traced - v12_rec| over {len(ref['rho'])} syn = {mx:.3e} "
          f"({', '.join(f'{k} {v:.1e}' for k, v in d.items() if v > 0) or 'all end states identical'}; "
          f"last traced rho/d/dn vs end {max(dl.values()):.1e}) -> {'OK' if ok else 'FAIL'}", flush=True)
    if not ok:
        FAIL.append(tag)
    return dict(maxdiff=mx, per_key=d, last=dl, ok=bool(ok))


# ---------------------------------------------------------------- population scans (v12_rec end states of every synapse)
_FJ = {}


def FJ(key):
    if key not in _FJ:
        _FJ[key] = json.load(open(FITS[key] if key in FITS else VFILES[key]))
    return _FJ[key]


def fpath(key):
    return FITS[key] if key in FITS else VFILES[key]


def pp_pairs(fitkey, path, target, cond="control"):
    """Pairs of the fit's per-pair table (<fit>_val_perpair.csv) for one target row."""
    import pandas as pd
    fp = FITS[fitkey]
    s = os.path.basename(fp)[:-len("_fit.json")]
    pp = pd.read_csv(f"{PP}/{s}_val_perpair.csv")
    pp = pp[(pp.pathway == MODEL[path]) & (pp.target == target) & (pp.condition == cond)]
    return sorted(set(pp.pair))


def lanes(info, d, dn):
    if info["nm"]:
        return np.clip(d + dn, info["dmin"], info["dnmax"]), np.clip(info["dp0"] + dn, info["dmin"], info["dnmax"])
    return np.asarray(d, float).copy(), np.full(np.shape(d), info["dp0"])


def scan(path, pairs, jobs, spkfit, hooks=None):
    """jobs {proto: (fit keys, pair set or None)} -> DataFrame, one row per (pair, proto, fit, synapse) with the v12_rec
    end state, the lanes (dctl = control d, dmg = mglu_block d; no_block d = d) and hook columns."""
    import pandas as pd
    rows, t0 = [], time.time()
    hooks = hooks or {}
    for ip, p in enumerate(pairs):
        have = [q for q in sorted(jobs) if (jobs[q][1] is None or p in jobs[q][1]) and has_rec(path, p, q)]
        if not have:
            continue
        spv, tiers, *_ = pair_spk(FJ(spkfit), path, p, have)
        for q in have:
            r = load_rec(path, p, q)
            n = r["effcai"].shape[0]
            for fk in jobs[q][0]:
                args, info = prep(FJ(fk), path, r, spv)
                o = ref_run(args, n)
                dc, dm = lanes(info, o["d"], o["dn"])
                ext = hooks[(q, fk)](r, args, info, o) if (q, fk) in hooks else {}
                for i in range(n):
                    rows.append(dict(pair=p, proto=q, fit=fk, i=i, syn=int(r["syn"][i]), rho0=float(r["rho0"][i]),
                                     dctl=float(dc[i]), dmg=float(dm[i]), spk_tier=int(tiers.get((p, int(r["syn"][i])), 0)),
                                     **{k: float(o[k][i]) for k in OUTK}, **{k: float(v[i]) for k, v in ext.items()}))
            _REC.pop((path, p, q), None); del r; gc.collect()     # one record in memory at a time
        drop_recs(); gc.collect()
        if ip % 10 == 0 or ip == len(pairs) - 1:
            print(f"  scan {path}: {ip + 1}/{len(pairs)} pairs, {len(rows)} rows ({time.time() - t0:.0f} s)", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(A.out, f"pop_{path}.csv"), index=False)
    return df


def wide(df, specs):
    """Join (pair, i) rows of several (proto, fit, tag) selections; columns get the suffix _tag."""
    out = None
    for proto, fit, tag in specs:
        x = df[(df.proto == proto) & (df.fit == fit)].set_index(["pair", "i"]).drop(columns=["proto", "fit"])
        x = x.add_suffix("_" + tag)
        out = x if out is None else out.join(x, how="inner")
    return out.reset_index()


def choose(w, key, mask, eps, what):
    """Among the synapses with mask and contrast w[key] > eps, the one closest to the median contrast."""
    m = mask & (w[key] > eps)
    c = w[m]
    info = dict(contrast=what, eps=eps, n_syn=int(len(w)), n_cand=int(len(c)),
                rule="closest to the median contrast among synapses where the mechanism acts (contrast > eps)")
    if c.empty:
        return None, info
    med = float(c[key].median())
    row = c.iloc[int(np.argmin(np.abs(c[key].to_numpy() - med)))]
    info.update(median=med, chosen=float(row[key]), max=float(c[key].max()), q25=float(c[key].quantile(.25)),
                q75=float(c[key].quantile(.75)), pair=row["pair"], i=int(row["i"]))
    return row, info


def lic_hook(r, args, info, o):
    """Per synapse: steps with pot before the licences > 0 (npq) and of those with the shaft licence closed (nblk);
    for the C-licence: pot > 0, shaft licence open, C-licence closed (nblkC). Traced kernel, 2 synapses per call."""
    n, T = r["effcai"].shape
    out = {k: np.zeros(n) for k in ("npq", "nblk", "nblkC")}
    for c0 in range(0, n, 2):
        isel = list(range(c0, min(n, c0 + 2)))
        TR, _, _, _ = trace_run(args, n, T, isel, info["cnt"])
        pq = TR[:, F_["pot_raw"]] > 0.5; lg = TR[:, F_["L_G"]] > 0.5; lc = TR[:, F_["L_C"]] > 0.5
        for s, i in enumerate(isel):
            out["npq"][i] = pq[s].sum(); out["nblk"][i] = (pq[s] & ~lg[s]).sum()
            out["nblkC"][i] = (pq[s] & lg[s] & ~lc[s]).sum() if info["c_lic"] else 0.0
        del TR
    return out


def ratio_hook(r, args, info, o):
    """EPSP-ratio readout of the pair (rho_v4.ratio on the pair basis) for the control and no_block lanes."""
    n = r["effcai"].shape[0]
    dc, _ = lanes(info, o["d"], o["dn"])
    try:
        b = batch_v2.PairBasis(r["pair"], r["syn"])
        rc = rho_v4.ratio(b, r["rho0"], o["rho"], dc, info["sigma"])
        rn = rho_v4.ratio(b, r["rho0"], o["rho"], o["d"], info["sigma"])
    except Exception as e:
        print(f"  ratio {r['pair']}: {type(e).__name__} {e}", flush=True)
        rc = rn = np.nan
    return dict(ratio_ctl=np.full(n, rc), ratio_nob=np.full(n, rn))


def ratio_course(fit, path, r, spv, npts=2000):
    """Pair EPSP-ratio time course (control / no_block lanes) from the traced rho, d_eCB, d_NO of every synapse."""
    args, info = prep(fit, path, r, spv)
    n, T = r["effcai"].shape
    idx = np.unique(np.linspace(0, T - 1, npts).astype(int))
    R = np.zeros((n, len(idx))); DE = np.zeros_like(R); DN = np.zeros_like(R)
    for c0 in range(0, n, 2):
        isel = list(range(c0, min(n, c0 + 2)))
        TR, _, _, _ = trace_run(args, n, T, isel, info["cnt"])
        R[isel] = TR[:, F_["rho"]][:, idx]; DE[isel] = TR[:, F_["d_eCB"]][:, idx]; DN[isel] = TR[:, F_["d_NO"]][:, idx]
        del TR
    b = batch_v2.PairBasis(r["pair"], r["syn"])
    rc = np.array([rho_v4.ratio(b, r["rho0"], R[:, j], lanes(info, DE[:, j], DN[:, j])[0], info["sigma"])
                   for j in range(len(idx))])
    rn = np.array([rho_v4.ratio(b, r["rho0"], R[:, j], DE[:, j], info["sigma"]) for j in range(len(idx))])
    return np.asarray(r["t"], float)[idx], rc, rn


def has_rec(path, pair, proto):
    set_basis(path)
    b = os.path.join(batch_v2.BASIS_DIR, "basis_{}_{}.csv".format(*pair.split("-")))
    return live_v12.L.find_record(path, pair, proto) is not None and os.path.isfile(b)


def train_start(r, ta, gap=150.0):
    """Start of the spike cluster (pre and post spikes, gaps > gap ms) that contains time ta (record ms)."""
    ev = np.sort(np.concatenate([np.asarray(r["prespikes"], float).ravel(), np.asarray(r["postspikes"], float).ravel()]))
    ev = ev[ev <= ta + 1e-9]
    if not ev.size:
        return float(ta)
    k = len(ev) - 1
    while k > 0 and ev[k] - ev[k - 1] <= gap:
        k -= 1
    return float(ev[k])


def first_step(r, TR, cond):
    k = np.flatnonzero(cond)
    return float(r["t"][k[0]]) if k.size else None


def first_ev(EV, NEV, s, *fields):
    for j in range(int(NEV[s])):
        if all(EV[s, j, EV_FIELDS.index(f)] > 0.5 for f in fields):
            return j
    return None


def ev_time(r, EV, s, j):
    return float(EV[s, j, EV_FIELDS.index("t_a")] + r["t"][0])


# ---------------------------------------------------------------- one saved trace
def save_trace(fig, vkey, path, r, spv, tiers, isel, roles, crit, event, win, extra=None, arrays=None):
    fit = FJ(vkey)
    args, info = prep(fit, path, r, spv)
    n, T = r["effcai"].shape
    tag = f"{fig} {r['pair']} {r['proto']} {vkey}"
    ref = ref_run(args, n)
    TR, EV, NEV, tr = trace_run(args, n, T, isel, info["cnt"])
    chk = check(tag, ref, tr, TR, np.asarray(isel))
    t = np.asarray(r["t"], np.float64)
    delay = np.asarray(MV.edge_params(r["syn"])["delay"], float)
    pre = np.asarray(r["prespikes"], float).ravel()
    i_ = np.asarray(isel)
    dctl, dmg = lanes(info, TR[:, F_["d_eCB"]], TR[:, F_["d_NO"]])
    meta = dict(fig=fig, variant=vkey, fit=fpath(vkey), path=path, pair=r["pair"], proto=r["proto"], roles=roles,
                syn=[int(r["syn"][i]) for i in isel], criterion=crit, t_event=event, window_ms=win,
                t0=float(t[0]), port_check=chk,
                params={k: (float(v) if np.isscalar(v) else None) for k, v in info.items()
                        if k not in ("td", "tp", "thE", "thP", "spk", "uE", "rest", "cnt")},
                spk_tier=[int(tiers.get((r["pair"], int(r["syn"][i])), 0)) for i in isel],
                tr_fields=TR_FIELDS, ev_fields=EV_FIELDS, **(extra or {}))
    pkt = dict(meta=json.dumps(meta), t=t, prespikes=pre, postspikes=np.asarray(r["postspikes"], float).ravel(),
               arrivals=pre[None, :] + delay[i_][:, None], syn_all=np.asarray(r["syn"]),
               EV=EV, NEV=NEV, d_control=dctl, d_mglu=dmg, d_noblock=TR[:, F_["d_eCB"]].copy(),
               theta_d=info["td"][i_], theta_p=info["tp"][i_], thE=info["thE"][i_], thP=info["thP"][i_],
               spk=info["spk"][i_], uE=info["uE"][i_], rest=info["rest"][i_],
               c_pre=np.asarray(r["c_pre"], float)[i_], c_post=np.asarray(r["c_post"], float)[i_], **(arrays or {}))
    for j, f in enumerate(TR_FIELDS):
        pkt[f] = TR[:, j] if f in ("rho", "d_eCB", "d_NO", "clk") else TR[:, j].astype(np.float32)
    for k in OUTK:
        pkt["end_" + k] = ref[k]
    fn = os.path.join(A.out, f"{fig}__{r['pair']}__{r['proto']}__{vkey}.npz")
    np.savez_compressed(fn, **pkt)
    print(f"  wrote {fn} (syn {meta['syn']} {roles}, event {event})", flush=True)
    SEL.append(dict(file=fn, **{k: meta[k] for k in ("fig", "variant", "pair", "proto", "syn", "roles", "criterion",
                                                      "t_event", "spk_tier")}, port_ok=chk["ok"], port_maxdiff=chk["maxdiff"]))
    return TR, EV, NEV


def stat(fig, msg, **kw):
    print(f"STATS {fig}: {msg}", flush=True)
    STATS.setdefault(fig, []).append(dict(msg=msg, **{k: (float(v) if isinstance(v, (int, float, np.number)) else v)
                                                      for k, v in kw.items()}))


def frac(a, b):
    return float(a) / float(b) if b else float("nan")


def pick_trace(fig, path, pair, protos, spkfit):
    spv, tiers, *_ = pair_spk(FJ(spkfit), path, pair, protos)
    return spv, tiers, {q: load_rec(path, pair, q) for q in protos}


# ---------------------------------------------------------------- figures
P50, P3A = "sjostrom_50hz_dt+10ms", "letzkus_3ap_200hz_dt+10ms"


LB_POT, LB_DEP = 0.8, 0.3                                 # rho_end > 0.8: potentiated; rho_end < 0.3: depressed


def lb_class(w):
    """KO pot & B4 dep / both pot / both dep / other (rho_end thresholds LB_POT, LB_DEP)."""
    kp, kd = w.rho_off > LB_POT, w.rho_off < LB_DEP
    bp, bd = w.rho_on > LB_POT, w.rho_on < LB_DEP
    return np.where(kp & bd, "KOpot_B4dep", np.where(kp & bp, "both_pot", np.where(kd & bd, "both_dep", "other")))


def train_starts(r, gap=150.0):
    """Starts of the spike clusters (pre + post spikes, gaps > gap ms), record ms."""
    ev = np.sort(np.concatenate([np.asarray(r["prespikes"], float).ravel(), np.asarray(r["postspikes"], float).ravel()]))
    return ev[np.r_[True, np.diff(ev) > gap]] if ev.size else np.zeros(1)


def best_train(r, TRs, mask):
    """Start of the train holding the most time (ms) of mask (per-step bool over the record)."""
    st = train_starts(r)
    t = np.asarray(r["t"], float); dt = np.r_[np.diff(t), 0.0]
    j = np.clip(np.searchsorted(st, t, side="right") - 1, 0, len(st) - 1)
    tt = np.bincount(j[mask], weights=dt[mask], minlength=len(st))
    return float(st[int(np.argmax(tt))]), float(tt.max()) if tt.size else 0.0


def lb_metrics(r, TRs, s, tp):
    """Per traced synapse s: pot-step time (c* > theta_p, ms), blocked fraction of that time, peak shaft dCa (uM),
    peak c* / theta_p, licensed pot time."""
    t = np.asarray(r["t"], float); dt = np.r_[np.diff(t), 0.0]
    pq = TRs[s, F_["pot_raw"]] > 0.5; lg = TRs[s, F_["L_G"]] > 0.5
    tpq = float(dt[pq].sum()); tbl = float(dt[pq & ~lg].sum())
    return dict(tpot=tpq, fblk_t=frac(tbl, tpq), pkX=float(np.nanmax(TRs[s, F_["shaft_dCa"]])),
                cpk=float(np.nanmax(TRs[s, F_["c"]]) / tp), tlic=float(dt[pq & lg].sum()))


def lb_rec_info(r, i):
    """Location fields of the record / edge table that are cheaply at hand (none built here)."""
    out = {}
    n = len(r["syn"])
    for k in r.keys() if hasattr(r, "keys") else []:
        if any(s in str(k).lower() for s in ("dist", "soma", "branch", "sec", "apic", "basal", "loc")):
            try:
                a = np.asarray(r[k])
                out[f"rec.{k}"] = a[i].item() if (a.ndim == 1 and a.shape[0] == n) else f"shape {a.shape}"
            except Exception:
                out[f"rec.{k}"] = "?"
    try:
        ep = MV.edge_params(r["syn"])
        for k in (ep.keys() if hasattr(ep, "keys") else []):
            if any(s in str(k).lower() for s in ("dist", "soma", "branch", "sec", "apic", "basal", "loc", "type")):
                a = np.asarray(ep[k])
                out[f"edge.{k}"] = a[i].item() if (a.ndim == 1 and a.shape[0] == n) else f"shape {a.shape}"
    except Exception as e:
        out["edge_err"] = f"{type(e).__name__}"
    return out


def fig_LB():
    """L2/3->L5 50 Hz +10 ms, V2B4 vs KO_lic. Population from pop_<path>.csv when present (the full run's scan; else
    rescanned). Two traced synapses, preferably of one pair, both rho0 = 0: col 1 'closed' (KO_lic rho_end > LB_POT,
    V2B4 rho_end < LB_DEP), col 2 'open' (both > LB_POT). Candidates of the best pairs are traced under V2B4 and scored:
    closed = blocked fraction of pot time + pot time / max + (1 - min(peak shaft dCa / theta_G, 1));
    open = (1 - blocked fraction) + pot time / max + min(peak shaft dCa / theta_G, 3) / 3."""
    import pandas as pd
    path = "L23"
    s50 = set(pp_pairs("V2B4", path, P50)); s50d = set(pp_pairs("V2B4", path, P50 + "@distal"))
    s3p = set(pp_pairs("V2B4", path, P3A + "@proximal")); s3d = set(pp_pairs("V2B4", path, P3A + "@distal"))
    pc = os.path.join(A.out, f"pop_{path}.csv")
    if os.path.isfile(pc):
        df = pd.read_csv(pc)
        print(f"  population from {pc} ({len(df)} rows, no rescan)", flush=True)
    else:
        jobs = {P50: (["V2B4", "KO_lic"], s50 | s50d), P3A: (["V2B4", "KO_lic"], s3p | s3d)}
        hooks = {(P50, "V2B4"): lic_hook, (P3A, "V2B4"): lic_hook}
        df = scan(path, sorted(s50 | s50d | s3p | s3d), jobs, "V2B4", hooks)
    w50 = wide(df, [(P50, "V2B4", "on"), (P50, "KO_lic", "off")]); w50["dr"] = w50.rho_off - w50.rho_on
    w3 = wide(df, [(P3A, "V2B4", "on"), (P3A, "KO_lic", "off")]); w3["dr"] = w3.rho_off - w3.rho_on
    w3["pot_on"] = w3.rho_on - w3.rho0_on
    for name, w, ps in (("50Hz +10 (all)", w50, s50), ("50Hz +10 @distal", w50, s50d),
                        ("3AP +10 @proximal", w3, s3p), ("3AP +10 @distal", w3, s3d)):
        x = w[w.pair.isin(ps)]
        stat("LB", f"{name}: licence closes a pot>0 interval on {int((x.nblk_on > 0).sum())}/{len(x)} synapses "
                   f"({frac((x.nblk_on > 0).sum(), len(x)):.3f}; pot>0 at all on {int((x.npq_on > 0).sum())}); "
                   f"blocked pot steps {frac(x.nblk_on.sum(), x.npq_on.sum()):.3f} of pot steps; mean rho_end V2B4 "
                   f"{x.rho_on.mean():.3f} vs KO_lic {x.rho_off.mean():.3f}; synapses with rho_KO - rho_V2B4 > 0.01: "
                   f"{int((x.dr > 0.01).sum())}", row=name, n_syn=len(x), n_block=int((x.nblk_on > 0).sum()),
             rho_on=x.rho_on.mean(), rho_off=x.rho_off.mean())
    # ---- population of the traced row (all 598 synapses of 50 Hz +10, every pair scanned): four classes
    wp = w50.copy()
    wp["cls"] = lb_class(wp)
    wp["fblk"] = np.where(wp.npq_on > 0, wp.nblk_on / wp.npq_on.where(wp.npq_on > 0, 1), np.nan)
    cn = {"KOpot_B4dep": "KO_lic pot & V2B4 dep", "both_pot": "both pot", "both_dep": "both dep", "other": "other"}
    parts = []
    for c, nm in cn.items():
        x = wp[wp.cls == c]; xf = x.fblk.dropna()
        parts.append(f"{nm} {len(x)} (mean fraction of pot>0 steps blocked {xf.mean():.3f} on {len(xf)} with pot>0)")
    stat("LB", f"50Hz +10 population ({len(wp)} synapses, {wp.pair.nunique()} pairs; pot = rho_end > {LB_POT}, "
               f"dep = rho_end < {LB_DEP}): " + "; ".join(parts),
         n_syn=len(wp), **{f"n_{c}": int((wp.cls == c).sum()) for c in cn},
         **{f"fblk_{c}": float(wp[wp.cls == c].fblk.mean()) for c in cn})
    x0 = wp[wp.rho0_on == 0]
    stat("LB", f"50Hz +10, rho0 = 0 synapses ({len(x0)}): KO_lic pot & V2B4 dep {int((x0.cls == 'KOpot_B4dep').sum())}, "
               f"both pot {int((x0.cls == 'both_pot').sum())}, both dep {int((x0.cls == 'both_dep').sum())}, "
               f"other {int((x0.cls == 'other').sum())}")
    # ---- trace candidates: closed (KOpot_B4dep) and open (both_pot), rho0 = 0 preferred, same pair preferred
    def cands(cls):
        c = wp[wp.cls == cls]
        return c[c.rho0_on == 0] if (c.rho0_on == 0).any() else c
    cc, co = cands("KOpot_B4dep"), cands("both_pot")
    same = sorted(set(cc.pair) & set(co.pair))
    pairs = same or sorted(set(cc.pair) | set(co.pair))
    print(f"  LB candidates: closed {len(cc)}, open {len(co)} (rho0 = 0 where available); pairs holding both: "
          f"{len(same)}", flush=True)
    M = []                                                # one dict per traced candidate
    thG = None
    for p in pairs:
        ic = list(cc[cc.pair == p].i.astype(int)); io = list(co[co.pair == p].i.astype(int))
        spv, tiers, recs = pick_trace("LB", path, p, [P50], "V2B4")
        r = recs[P50]
        args, info = prep(FJ("V2B4"), path, r, spv); thG = info["thG"]
        n, T = r["effcai"].shape
        for c0 in range(0, len(ic + io), 2):
            isel = (ic + io)[c0:c0 + 2]
            TRs, _, _, _ = trace_run(args, n, T, isel, info["cnt"])
            for s, i in enumerate(isel):
                rw = wp[(wp.pair == p) & (wp.i == i)].iloc[0]
                M.append(dict(pair=p, i=i, syn=int(r["syn"][i]), role="closed" if i in ic else "open",
                              rho0=float(rw.rho0_on), rho_B4=float(rw.rho_on), rho_KO=float(rw.rho_off),
                              fblk=float(rw.fblk), **lb_metrics(r, TRs, s, float(info["tp"][i]))))
            del TRs
        _REC.pop((path, p, P50), None); del r; gc.collect()
    M = pd.DataFrame(M)
    for role in ("closed", "open"):
        m = M.role == role
        tn = M.loc[m, "tpot"] / max(M.loc[m, "tpot"].max(), 1e-9)
        xg = M.loc[m, "pkX"] / thG
        M.loc[m, "score"] = (M.loc[m, "fblk_t"] + tn + (1 - xg.clip(upper=1.0)) if role == "closed" else
                             (1 - M.loc[m, "fblk_t"]) + tn + xg.clip(upper=3.0) / 3)
    M.loc[(M.role == "open") & (M.pkX <= thG), "score"] -= 10      # open column: shaft dCa must cross theta_G
    M["score"] = M.score.fillna(-1e9)
    with pd.option_context("display.width", 200, "display.max_columns", 20, "display.max_rows", 200):
        print(M.sort_values(["role", "score"], ascending=[True, False]).to_string(index=False, float_format="%.4g"),
              flush=True)
    if same:
        best = {p: M[(M.pair == p) & (M.role == "closed")].score.max() + M[(M.pair == p) & (M.role == "open")].score.max()
                for p in same}
        pb = max(best, key=best.get)
        pick = {ro: M[(M.pair == pb) & (M.role == ro)].sort_values("score").iloc[-1] for ro in ("closed", "open")}
    else:
        pick = {ro: M[M.role == ro].sort_values("score").iloc[-1] for ro in ("closed", "open")}
    CHOICE["LB"] = dict(rule="closed: KO_lic rho_end > 0.8 & V2B4 rho_end < 0.3; open: both > 0.8; rho0 = 0 preferred; "
                             "same pair preferred; max score (fig_LB doc)",
                        candidates=M.to_dict("records"), picked={k: v.to_dict() for k, v in pick.items()})
    for ro, rw in pick.items():
        stat("LB", f"traced '{ro}': pair {rw.pair} syn {int(rw.syn)} (@distal pair: {'yes' if rw.pair in s50d else 'no'}), "
                   f"rho0 {rw.rho0:.0f}, rho_end V2B4 {rw.rho_B4:.3f} / KO_lic {rw.rho_KO:.3f}, peak shaft dCa "
                   f"{rw.pkX:.3f} uM (theta_G {thG:.2f}), fraction of pot>0 steps blocked {rw.fblk:.3f} (of pot time "
                   f"{rw.fblk_t:.3f}; c* > theta_p for {rw.tpot:.0f} ms, peak c*/theta_p {rw.cpk:.2f})",
             role=ro, pair=rw.pair, syn=int(rw.syn), rho0=rw.rho0, rho_B4=rw.rho_B4, rho_KO=rw.rho_KO, pkX=rw.pkX,
             fblk=rw.fblk)
    # ---- save: one npz per (pair, variant) holding both traced synapses (closed first); per-column train in meta
    by = {}
    for ro in ("closed", "open"):
        by.setdefault(pick[ro].pair, []).append(ro)
    for p, roles in by.items():
        spv, tiers, recs = pick_trace("LB", path, p, [P50], "V2B4")
        r = recs[P50]
        isel = [int(pick[ro].i) for ro in roles]
        args, info = prep(FJ("V2B4"), path, r, spv)
        TRs, _, _, _ = trace_run(args, r["effcai"].shape[0], r["effcai"].shape[1], isel, info["cnt"])
        evs, locs = [], []
        for s, ro in enumerate(roles):
            pq = TRs[s, F_["pot_raw"]] > 0.5; lg = TRs[s, F_["L_G"]] > 0.5
            ev, tt = best_train(r, TRs, (pq & ~lg) if ro == "closed" else (pq & lg))
            evs.append(ev); locs.append(lb_rec_info(r, isel[s]))
            print(f"LB {ro}: {p} syn {int(r['syn'][isel[s]])}, train {ev:.1f} ms ({tt:.1f} ms of "
                  f"{'blocked' if ro == 'closed' else 'licensed'} pot in it); location fields {locs[-1] or 'none'}",
                  flush=True)
        for v in ("V2B4", "KO_lic"):
            save_trace("LB", v, path, r, spv, tiers, isel, roles, CHOICE["LB"]["rule"], evs[0], [-20.0, 160.0],
                       dict(layout="v2", t_event_cols=evs, locations=locs, theta_G=thG,
                            picked=[{k: (v_.item() if hasattr(v_, "item") else v_) for k, v_ in pick[ro].items()}
                                    for ro in roles]))
        drop_recs()


L5P = dict(m01="sjostrom_0.1hz_dt-10ms", p01="sjostrom_0.1hz_dt+10ms", p20="sjostrom_20hz_dt+10ms",
           z20="sjostrom_20hz_dt0ms", m20="sjostrom_20hz_dt-10ms", s07="sjostrom07_step200ms_pair",
           m50="sjostrom_50hz_dt-10ms")
DF5 = {}


def l5_scan():
    if "df" not in DF5:
        pairs = sorted(set().union(*[pp_pairs("V2B4", "L5", L5P[k]) for k in ("m01", "p20", "s07", "m50")]))
        jobs = {L5P["m01"]: (["V2B4"], None), L5P["p01"]: (["V2B4"], None), L5P["p20"]: (["V2B4", "KO_veto"], None),
                L5P["z20"]: (["V2B4"], None), L5P["m20"]: (["V2B4"], None),
                L5P["s07"]: (["V2B4", "V2B4_019", "KO_no"], None), L5P["m50"]: (["V2B4", "KO_b4"], None)}
        hooks = {(L5P["s07"], f): ratio_hook for f in ("V2B4", "V2B4_019", "KO_no")}
        DF5["df"] = scan("L5", pairs, jobs, "V2B4", hooks)
    return DF5["df"]


def fig_LC():
    df = l5_scan(); path = "L5"; m, p = L5P["m01"], L5P["p01"]
    for q in (m, p):
        x = df[(df.proto == q) & (df.fit == "V2B4")]
        stat("LC", f"{q}: arrivals stepping {int(x['ne'].sum())}/{int(x.nar.sum())} = {frac(x['ne'].sum(), x.nar.sum()):.3f}; "
                   f"synapses with a step {int((x['ne'] > 0).sum())}/{len(x)}; mean d_eCB end {x.d.mean():.3f}, "
                   f"mean control d {x.dctl.mean():.3f}, mglu_block d {x.dmg.mean():.3f}", proto=q,
             frac_step=frac(x['ne'].sum(), x.nar.sum()))
    w = wide(df, [(m, "V2B4", "m"), (p, "V2B4", "p")]); w["dd"] = w.d_p - w.d_m
    row, ci = choose(w, "dd", (w.ne_m > 0) & (w.ne_p == 0), 0.01, "d_eCB_end(+10) - d_eCB_end(-10), same synapse")
    if row is None:
        row, ci = choose(w, "dd", w.ne_m > 0, 0.01, "d_eCB_end(+10) - d_eCB_end(-10) (relaxed: steps at +10 allowed)")
    CHOICE["LC"] = ci
    if row is None:
        print("LC: no candidate", flush=True); return
    pr, i = row["pair"], int(row["i"])
    spv, tiers, recs = pick_trace("LC", path, pr, [m, p], "V2B4")
    assert np.array_equal(recs[m]["syn"], recs[p]["syn"]), "LC: synapse order differs"
    args, info = prep(FJ("V2B4"), path, recs[m], spv)
    _, EV, NEV, _ = trace_run(args, recs[m]["effcai"].shape[0], recs[m]["effcai"].shape[1], [i], info["cnt"])
    j = first_ev(EV, NEV, 0, "step"); j = 0 if j is None else j
    print(f"LC: {pr} syn {int(recs[m]['syn'][i])}, contrast {ci['chosen']:.3f}, arrival {j}", flush=True)
    for q in (m, p):
        r = recs[q]
        a2, i2 = prep(FJ("V2B4"), path, r, spv)
        _, E2, N2, _ = trace_run(a2, r["effcai"].shape[0], r["effcai"].shape[1], [i], i2["cnt"])
        ta = ev_time(r, E2, 0, min(j, int(N2[0]) - 1))
        save_trace("LC", "V2B4", path, r, spv, tiers, [i], ["eCB"], ci["rule"] + "; " + ci["contrast"], ta, [-30.0, 60.0],
                   dict(choice=ci, arrival_index=int(j), t_arrival=ta, lanes=["control", "mglu_block"]))
    drop_recs()


def fig_LD():
    df = l5_scan(); path = "L5"; p, z, m = L5P["p20"], L5P["z20"], L5P["m20"]
    for q in (p, z, m):
        x = df[(df.proto == q) & (df.fit == "V2B4")]
        tr = x['ne'].sum() + x.nv.sum() + x.nann.sum()
        stat("LD", f"{q}: triggering arrivals vetoed {int(x.nv.sum())}/{int(tr)} = {frac(x.nv.sum(), tr):.3f} "
                   f"(stepped {int(x['ne'].sum())}, annulled {int(x.nann.sum())})", proto=q, frac_veto=frac(x.nv.sum(), tr))
    w = wide(df, [(p, "V2B4", "on"), (p, "KO_veto", "off"), (m, "V2B4", "m")]); w["dd"] = (w.d_off - w.d_on).abs()
    stat("LD", f"{p}: mean d_eCB end V2B4 {w.d_on.mean():.3f} vs KO_veto {w.d_off.mean():.3f}; synapses with "
               f"|diff| > 0.01: {int((w.dd > 0.01).sum())}/{len(w)}")
    row, ci = choose(w, "dd", (w.nv_on > 0) & (w.ne_m > 0), 0.01, "|d_eCB_end(KO_veto) - d_eCB_end(V2B4)| at 20Hz +10 "
                     "(synapse also steps at 20Hz -10)")
    if row is None:
        row, ci = choose(w, "dd", w.nv_on > 0, 0.01, "|d_eCB_end(KO_veto) - d_eCB_end(V2B4)| at 20Hz +10")
    CHOICE["LD"] = ci
    if row is None:
        print("LD: no candidate", flush=True); return
    pr, i = row["pair"], int(row["i"])
    spv, tiers, recs = pick_trace("LD", path, pr, [p, m], "V2B4")
    print(f"LD: {pr} syn {int(recs[p]['syn'][i])}, contrast {ci['chosen']:.3f}", flush=True)
    for q, v, flds in ((p, "V2B4", ("veto", "trig")), (p, "KO_veto", None), (m, "V2B4", ("step",))):
        r = recs[q]
        if flds is not None:
            a2, i2 = prep(FJ(v), path, r, spv)
            _, E2, N2, _ = trace_run(a2, r["effcai"].shape[0], r["effcai"].shape[1], [i], i2["cnt"])
            j = first_ev(E2, N2, 0, *flds); j = 0 if j is None else j
            if q == p:
                jp = j
        else:
            j = jp
            a2, i2 = prep(FJ(v), path, r, spv)
            _, E2, N2, _ = trace_run(a2, r["effcai"].shape[0], r["effcai"].shape[1], [i], i2["cnt"])
        ta = ev_time(r, E2, 0, j)
        save_trace("LD", v, path, r, spv, tiers, [i], ["window"], ci["rule"] + "; " + ci["contrast"], ta, [-15.0, 30.0],
                   dict(choice=ci, arrival_index=int(j), t_arrival=ta))
    drop_recs()


def fig_LE():
    df = l5_scan(); path = "L5"; q = L5P["s07"]
    w = wide(df, [(q, "V2B4", "on"), (q, "V2B4_019", "o19"), (q, "KO_no", "ko")]); w["dd"] = (w.dctl_on - w.d_on).abs()
    pr_ = w.drop_duplicates("pair")
    stat("LE", f"{q}: mean end total d control {w.dctl_on.mean():.3f} vs no_block {w.d_on.mean():.3f} (V2B4), control "
               f"{w.dctl_o19.mean():.3f} (V2B4_019), KO_no {w.dctl_ko.mean():.3f}; mean pair EPSP ratio control "
               f"{pr_.ratio_ctl_on.mean():.3f} / no_block {pr_.ratio_nob_on.mean():.3f} (V2B4), V2B4_019 "
               f"{pr_.ratio_ctl_o19.mean():.3f}, KO_no {pr_.ratio_ctl_ko.mean():.3f}; synapses with NO events "
               f"{int((w.nno_on > 0).sum())}/{len(w)}", d_ctl=w.dctl_on.mean(), d_nob=w.d_on.mean())
    row, ci = choose(w, "dd", w.nno_on > 0, 0.01, "|control d - no_block d| at the end (d_NO effect), V2B4")
    CHOICE["LE"] = ci
    if row is None:
        print("LE: no candidate", flush=True); return
    pr, i = row["pair"], int(row["i"])
    spv, tiers, recs = pick_trace("LE", path, pr, [q], "V2B4")
    r = recs[q]
    print(f"LE: {pr} syn {int(r['syn'][i])}, contrast {ci['chosen']:.3f}", flush=True)
    ev = float(np.min(np.asarray(r["prespikes"], float)))
    for v in ("V2B4", "V2B4_019", "KO_no"):
        tt, rc, rn = ratio_course(FJ(v), path, r, spv)
        save_trace("LE", v, path, r, spv, tiers, [i], ["NO"], ci["rule"] + "; " + ci["contrast"], ev, None,
                   dict(choice=ci, lanes=["control", "no_block"]), dict(t_ratio=tt, ratio_ctl=rc, ratio_nob=rn))
    drop_recs()


def fig_LF():
    df = l5_scan(); path = "L5"; q = L5P["m50"]
    w = wide(df, [(q, "V2B4", "on"), (q, "KO_b4", "off")]); w["dd"] = (w.d_on - w.d_off).abs()
    x = df[(df.proto == q) & (df.fit == "V2B4")]
    stat("LF", f"{q}: arrivals annulled {int(x.nann.sum())}/{int(x.nar.sum())} = {frac(x.nann.sum(), x.nar.sum()):.3f} "
               f"(C events {int(x.nC.sum())}); synapses with |d_KO_b4 - d_V2B4| > 0.01: {int((w.dd > 0.01).sum())}/{len(w)}; "
               f"mean d_eCB end V2B4 {w.d_on.mean():.3f} vs KO_b4 {w.d_off.mean():.3f}; mean control d "
               f"{w.dctl_on.mean():.3f} vs {w.dctl_off.mean():.3f}", frac_annul=frac(x.nann.sum(), x.nar.sum()))
    row, ci = choose(w, "dd", w.nann_on > 0, 0.01, "|d_eCB_end(KO_b4) - d_eCB_end(V2B4)|, 50Hz -10")
    CHOICE["LF_B4"] = ci
    if row is not None:
        pr, i = row["pair"], int(row["i"])
        spv, tiers, recs = pick_trace("LF", path, pr, [q], "V2B4")
        r = recs[q]
        args, info = prep(FJ("V2B4"), path, r, spv)
        _, EV, NEV, _ = trace_run(args, r["effcai"].shape[0], r["effcai"].shape[1], [i], info["cnt"])
        j = first_ev(EV, NEV, 0, "annul"); j = 0 if j is None else j
        ev = train_start(r, ev_time(r, EV, 0, j))
        print(f"LF B4: {pr} syn {int(r['syn'][i])}, contrast {ci['chosen']:.3f}, train {ev:.1f} ms", flush=True)
        for v in ("V2B4", "KO_b4"):
            save_trace("LF", v, path, r, spv, tiers, [i], ["C"], ci["rule"] + "; " + ci["contrast"], ev, [-20.0, 140.0],
                       dict(choice=ci, t_annulled=ev_time(r, EV, 0, j)))
        drop_recs()
    # B1: Zilberter L2/3 -> L2/3, V2B1B4 vs the same parameters with c_lic 0
    path = "L23L23"; zs = ("zilberter_train10_50hz_dt-10ms_last", "zilberter_train10_50hz_post_only")
    pairs = sorted(set(pp_pairs("V2B1B4", path, zs[0])) | set(pp_pairs("V2B1B4", path, zs[1])))
    jobs = {z: (["V2B1B4", "V2B1B4_clic0"], None) for z in zs}
    dz = scan(path, pairs, jobs, "V2B1B4")
    ws = []
    for z in zs:
        wz = wide(dz, [(z, "V2B1B4", "on"), (z, "V2B1B4_clic0", "off")]); wz["dr"] = (wz.rho_off - wz.rho_on).abs()
        wz["proto"] = z; ws.append(wz)
        stat("LF", f"{z}: synapses where the C-licence changes rho_end by > 0.01: {int((wz.dr > 0.01).sum())}/{len(wz)} "
                   f"({frac((wz.dr > 0.01).sum(), len(wz)):.3f}); C-licence blocks pot on "
                   f"{int((wz.ncb_on > 0).sum())} synapses; mean rho_end V2B1B4 {wz.rho_on.mean():.3f} vs c_lic 0 "
                   f"{wz.rho_off.mean():.3f}", proto=z)
    import pandas as pd
    wa = pd.concat(ws, ignore_index=True)
    row, ci = choose(wa, "dr", wa.ncb_on > 0, 0.01, "|rho_end(c_lic 0) - rho_end(V2B1B4)|, both Zilberter records pooled")
    CHOICE["LF_B1"] = ci
    if row is None:
        print("LF B1: no candidate", flush=True); return
    pr, i, z = row["pair"], int(row["i"]), row["proto"]
    ci["proto"] = z
    spv, tiers, recs = pick_trace("LF", path, pr, [z], "V2B1B4")
    r = recs[z]
    args, info = prep(FJ("V2B1B4"), path, r, spv)
    TR, _, _, _ = trace_run(args, r["effcai"].shape[0], r["effcai"].shape[1], [i], info["cnt"])
    pq = TR[0, F_["pot_raw"]] > 0.5; lg = TR[0, F_["L_G"]] > 0.5; lc = TR[0, F_["L_C"]] > 0.5
    tk = first_step(r, TR, pq & lg & ~lc)
    if tk is None:
        tk = first_step(r, TR, pq & ~lc)
    ev = train_start(r, tk if tk is not None else float(r["t"][0]))
    print(f"LF B1: {pr} {z} syn {int(r['syn'][i])}, contrast {ci['chosen']:.3f}, train {ev:.1f} ms", flush=True)
    for v in ("V2B1B4", "V2B1B4_clic0"):
        save_trace("LF", v, path, r, spv, tiers, [i], ["C_licence"], ci["rule"] + "; " + ci["contrast"], ev,
                   [-20.0, 320.0], dict(choice=ci))
    drop_recs()


F_ = {f: j for j, f in enumerate(TR_FIELDS)}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--figs", default="LB,LC,LD,LE,LF")
    A = ap.parse_args()
    os.makedirs(A.out, exist_ok=True)
    SEL, MISSING, STATS, CHOICE = [], set(), {}, {}
    VFILES = {}
    for v in VARIANTS:
        c = sorted(glob.glob(os.path.join(VARDIR, f"*{v}*.json")))
        c = [x for x in c if not x.endswith(("_val.json", "_spk.json"))]
        if c:
            ex = [x for x in c if os.path.basename(x) in (f"{v}.json", f"{v}_fit.json")]
            VFILES[v] = (ex or c)[0]
    print(f"fits {FITS}\nvariants found {VFILES}\nmissing {[v for v in VARIANTS if v not in VFILES]}", flush=True)
    t0 = time.time()
    for fg in A.figs.split(","):
        print(f"\n=== {fg}", flush=True)
        try:
            globals()[f"fig_{fg}"]()
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f"  {fg}: ERROR {type(e).__name__}: {e}", flush=True)
            FAIL.append(f"{fg} error: {type(e).__name__}: {e}")
        drop_recs()
        print(f"  ({time.time() - t0:.0f} s)", flush=True)
    json.dump(dict(selection=SEL, choice=CHOICE, stats=STATS, variant_files=VFILES, port_failures=FAIL,
                   missing_variants=[v for v in VARIANTS if v not in VFILES]),
              open(os.path.join(A.out, "selection_lat.json"), "w"), indent=1, default=str)
    print(f"\nPORT CHECK SUMMARY: {sum(s['port_ok'] for s in SEL)}/{len(SEL)} traced records identical to v12_rec "
          f"(max diff {max([s['port_maxdiff'] for s in SEL] or [0]):.3e}); failures {FAIL or 'none'}", flush=True)
    if FAIL:
        sys.exit(3)
