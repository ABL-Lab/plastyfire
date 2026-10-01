"""GPU objective v8 (GATE_ALT.md, GATE_IMPL.md): v7 (gpu_v7_rho, not modified) with the LTP licence read from a chosen
synapse-local Ca quantity. Option gate_src in --set json:
  0  Vg (own spine VDCC pool, tau_E1): every call is the v7 code path unchanged (W2_N3 repro at MAXITER 0).
  1  spine [Ca] above rest (cai_CR - min_ca_CR, recovered exactly from effcai, batch_v2.cacr_from_effcai), uM.
  2  shaft [Ca] above the synapse's pre-stimulus rest (npz shaft_cai = cai of the segment carrying the synapse), uM.
  3  c* (effcai) itself, uM s (mM ms), licence = [c* > theta_G].
  gate_win > 0 (ms, src 1/2): peak-in-window, licence = [max of the signal over [t - gate_win, t] > theta_G] (GATE_ALT G1);
  gate_win 0 (src 1/2): 100 ms low-pass, dG/dt = (x - G) / tau_E1, licence = [G > theta_G] (GATE_ALT G2).
Threshold theta_G (natural unit, one value for every pathway): gate_theta > 0 fixes it (theta_V then unused); otherwise
theta_G = kappa x theta_V (DE slot theta_V), kappa = Q50 / gate_conv_thV, Q50 = median over synapses (pooled over all
models) of each synapse's median gate quantity at the moments the Vg gate opens (Vg crosses gate_conv_thV, default W2_N3
theta_V 3.575), so the W2_N3 seed starts at the converted threshold of the median synapse. gate_kappa > 0 pins kappa.
The saved json gets a "v8" block (gate_src, gate_win, kappa, Q50, theta_G and unit). For gate_src != 0 the fit_v6
repro check (seed rule spec) fails by design; that REPRO DIFF exit is caught here and the run exits 0.
src 2 calibration print (GATE_ALT checklist): single isolated bAP peak shaft Delta[Ca] per synapse (bap_ref selection,
peak over [x - 1, x + 30] ms) by cexp loc / proximal (<= 100 um), Spearman vs cexp vdcc_q_post; 3AP 200 Hz record peaks.
t_exact 1 (SET): GPUModelV8 inherits gpu_v7x_rho.GPUModelV7X; the eCB trigger and veto use the v7x exact-arrival corrections
(live convention); t_exact 0 = the v7 arithmetic unchanged (G8 results bit-for-bit). v5_mode 1/2 only (no E2 latch). Shaft is stored on the GPU as uint16 nM above rest (clip 0..65535 nM).
Run: python gpu_v8_rho.py <fit_v6 args>.
"""
import json, math, os, sys
import numpy as np
from numba import cuda, njit

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import gpu_v7_rho, gpu_v7x_rho
from gpu_v7_rho import KT
from gpu_v7x_rho import GPUModelV7X, arrival_frac
from batch_v2 import RHO_STAR_GB, TAU_EFFCA

UNIT = {1: "uM", 2: "uM", 3: "uM s"}
_INST = []
_STATE = dict(kappa=None)


def _make_kernel_v8(weighted, src, win, tex=False):
    """gpu_v7_rho._make_kernel_v7 (no E2) with the potentiation licence read from gate quantity src (see module doc)."""
    RS = float(RHO_STAR_GB)
    TRIG_W = bool(weighted)
    SRC = int(src)
    WIN = float(win)
    TS = float(TAU_EFFCA)
    TEX = bool(tex)                       # v7x arrival-time convention (gpu_v7x_rho, SET t_exact 1)

    @cuda.jit
    def kern(E, S, Hg, eoff, slen, hoff, aptr, astep, acnt, rho0, order, td, tp, prm, out_rho, out_d, uV, uE, vq,
             X, gsc, gfix, fa):
        p = cuda.blockIdx.x * cuda.blockDim.x + cuda.threadIdx.x
        q = cuda.blockIdx.y * cuda.blockDim.y + cuda.threadIdx.y
        if p >= td.shape[0] or q >= order.shape[0]:
            return
        i = order[q]
        GD = prm[p, 0]; GP = prm[p, 1]; thV0 = prm[p, 2]; iscV = prm[p, 3]; tauE1 = prm[p, 4]
        dmin = prm[p, 5]; Ae = prm[p, 6]; dp0 = prm[p, 7]
        thE0 = prm[p, 8]; tauD = prm[p, 9]
        thV = thV0 * uV[i]
        if thE0 <= 0.0:
            thE = thV
        else:
            thE = thE0 * uE[i]
        thG = gfix if gfix > 0.0 else thV0 * gsc
        a_td = td[p, i]; a_tp = tp[p, i]
        r = rho0[i]; d0 = dp0
        V = 0.0; chV = -1.0; aV = 1.0; bV = 0.0
        W = 0.0; Bg = 0.0; chD = -1.0; aD = 1.0
        Gq = 0.0; clk = 0.0; tl = -1e30; chS = -1.0; aS = 1.0; bS = 0.0
        ip = 0.0; up = 0.0                                      # v7x: W / V input of the previous step
        e0 = eoff[i]; n = slen[i]; h0 = hoff[i]
        ka = aptr[i]; kae = aptr[i + 1]
        na = astep[ka] if ka < kae else -1
        for k in range(n):
            e = E[e0 + k] * 1.0
            s = S[e0 + k] * 1.0
            h = Hg[h0 + k]
            hm = h * KT
            c = 0.0
            qa = 0.0
            fx = 0.0
            if k == na:
                c = acnt[ka]; qa = vq[ka]
                if TEX:
                    fx = fa[ka]
                ka += 1
                na = astep[ka] if ka < kae else -1
            pot = 1.0 if e > a_tp else 0.0
            dep = 1.0 if e > a_td else 0.0
            if SRC == 3:
                lic = e > thG
            elif WIN > 0.0:
                lic = (clk - tl) <= WIN
            else:
                lic = Gq > thG
            if pot > 0.0 and thG > 0.0 and not lic:
                pot = 0.0
            r = r + h * (-r * (1 - r) * (RS - r) + pot * GP * (1 - r) - dep * (1 - pot) * GD * r)
            if r < 0.0:
                r = 0.0
            elif r > 1.0:
                r = 1.0
            if TRIG_W:
                if c > 0.0:
                    Bg = Bg + c
                trig = (W - fx * ip) > thE if TEX else W > thE      # v7x: W at the exact arrival
            else:
                trig = (V - fx * up) > thE if TEX else V > thE
            if TEX:
                qa = qa + fx * up                                    # v7x: veto window from t_a
            if c > 0.0 and Ae > 0.0 and trig and not (qa > thE):
                d0 = dmin + (d0 - dmin) * math.pow(1.0 - Ae, c)
            if hm != chV:
                chV = hm
                aV = math.exp(-hm / tauE1); bV = tauE1 * (1.0 - aV) / iscV
            if TRIG_W:
                wb = 1.0 - (Bg if Bg < 1.0 else 1.0)
                if TEX:
                    ip = bV * (wb * s)
                    W = aV * W + ip
                else:
                    W = aV * W + bV * (wb * s)
                if hm != chD:
                    chD = hm
                    aD = math.exp(-hm / tauD)
                Bg = Bg * aD
            if TEX:
                up = bV * s
                V = aV * V + up
            else:
                V = aV * V + bV * s
            if SRC == 1 or SRC == 2:
                if SRC == 1:
                    x = 0.0
                    if hm > 0.0 and k + 1 < n:
                        if hm != chS:
                            chS = hm
                            aS = math.exp(-hm / TS); bS = 1e3 / (TS * (1.0 - aS))
                        x = (E[e0 + k + 1] * 1.0 - aS * e) * bS          # step mean of cai_CR - min_ca_CR, uM
                    tx = clk + hm
                else:
                    x = X[e0 + k] * 1e-3                                  # nM -> uM
                    tx = clk
                if WIN > 0.0:
                    if x > thG:
                        tl = tx
                else:
                    Gq = aV * Gq + (1.0 - aV) * x
            clk += hm
        out_rho[p, i] = r
        out_d[p, i, 0] = d0; out_d[p, i, 1] = dp0; out_d[p, i, 2] = d0

    return kern


@njit(cache=False)
def _open_q(t, s, e, x, src, win, tauE1, isc, thV, TS):
    """Gate quantity at each Vg gate opening (first step k with V_k > thV after V_{k-1} <= thV), as the kernel sees it."""
    n = t.shape[0]
    res = np.empty(n); m = 0
    xs = np.zeros(n); txs = np.zeros(n)
    V = 0.0; Gq = 0.0; clk = 0.0; prev = False
    for k in range(n):
        hm = (t[k + 1] - t[k]) if k + 1 < n else 0.0
        op = V > thV
        if op and not prev:
            if src == 3:
                qv = e[k]
            elif win > 0.0:
                qv = 0.0; j = k - 1
                while j >= 0 and txs[j] >= clk - win:
                    if xs[j] > qv:
                        qv = xs[j]
                    j -= 1
            else:
                qv = Gq
            res[m] = qv; m += 1
        prev = op
        aV = math.exp(-hm / tauE1); bV = tauE1 * (1.0 - aV) / isc
        V = aV * V + bV * s[k]
        xk = 0.0; tx = clk
        if src == 1:
            if hm > 0.0 and k + 1 < n:
                aS = math.exp(-hm / TS); xk = (e[k + 1] - aS * e[k]) * 1e3 / (TS * (1.0 - aS))
            tx = clk + hm
        elif src == 2:
            xk = x[k]
        xs[k] = xk; txs[k] = tx
        Gq = aV * Gq + (1.0 - aV) * xk
        clk += hm
    return res[:m]


def _shaft_rest(r, sh):
    """Per-synapse shaft cai (mM) at the last grid sample >= 1 ms before the first pre or post spike."""
    t = np.asarray(r["t"], np.float64)
    ev = np.concatenate([np.asarray(r["prespikes"], float).ravel(), np.asarray(r["postspikes"], float).ravel()])
    k = 0
    if ev.size:
        k = max(int(np.searchsorted(t, ev.min() - 1.0)) - 1, 0)
    return sh[:, k].astype(np.float64)


class GPUModelV8(GPUModelV7X):
    def __init__(self, B, targets, filters, **kw):
        super().__init__(B, targets, filters, **kw)
        P0 = self.P0
        self.gate_src = int(P0.get("gate_src", 0))
        if self.gate_src == 0:
            return
        if self.gate_src not in (1, 2, 3):
            raise ValueError(f"gate_src {self.gate_src}")
        if not self.v5 or self.v5_mode == 3:
            raise ValueError("v8 gate_src needs v5_mode 1 or 2")
        self.gate_win = float(P0.get("gate_win", 100.0)) if self.gate_src in (1, 2) else 0.0
        self.gate_theta = float(P0.get("gate_theta", 0.0))
        thc = float(P0.get("gate_conv_thV", 3.5751446275588403))
        tauE1 = float(P0.get("tau_E1", 100.0)); isc = float(P0.get("i_scale", 1e-5))
        if self.veto_T <= 0:
            na = int(np.asarray(self.d["aptr"].copy_to_host())[-1])
            self.vq = cuda.to_device(np.zeros(max(na, 1), np.float64))
        eoff = np.asarray(self.d["eoff"].copy_to_host())
        N = self.n_ctrl
        X = np.zeros(int(eoff[-1]) if self.gate_src == 2 else 1, np.uint16)
        oq = np.full(N, np.nan); self.burst_sh = np.full(N, np.nan)
        rests, nclip, ev, keys = [], 0, {}, []
        for r in B.recs:
            t = np.asarray(r["t"], np.float64); sl = r["sl"]; n, T = r["effcai"].shape
            S = np.asarray(r["vdcc"], np.float32).astype(np.float64)
            if self.gate_src == 2:
                if r.get("shaft_cai") is None:
                    raise ValueError(f"{r['pair']}/{r['proto']}: gate_src 2 needs shaft_cai (BatchV2 signals)")
                sh = np.asarray(r["shaft_cai"], np.float64)
                rest = _shaft_rest(r, sh); rests.append(rest)
                du = (sh - rest[:, None]) * 1e3                              # uM above rest
                nm = np.rint(du * 1e3); nclip += int(np.sum(nm > 65535))
                post = np.sort(np.asarray(r["postspikes"], float).ravel())
                iso = np.array([not np.any((post != x) & (post > x - 100.0) & (post < x + 20.0)) for x in post], bool)
                b3 = "3ap_200hz" in r["proto"]
            for ii in range(n):
                i = sl.start + ii
                e = np.asarray(r["effcai"][ii], np.float64)
                keys.append((r["pair"], int(r["syn"][ii])))
                if self.gate_src == 2:
                    X[eoff[i]:eoff[i] + T] = np.clip(nm[ii], 0, 65535).astype(np.uint16)
                    xr = np.clip(nm[ii], 0, 65535) * 1e-3                    # what the kernel reads, uM
                    ta = t[np.asarray(r["arr"][ii]).ravel().astype(np.int64)]
                    L = ev.setdefault(keys[-1], [])
                    for x in post[iso]:                                      # isolated bAPs as gpu_v7_rho.bap_ref
                        if np.any((ta > x - 150.0) & (ta < x + 8.0)):
                            continue
                        a, b = np.searchsorted(t, [x - 1.0, x + 30.0])
                        if b > a:
                            L.append(float(du[ii, a:b].max()))
                    if b3:
                        self.burst_sh[i] = float(du[ii].max())
                else:
                    xr = np.zeros(1)
                v = _open_q(t, S[ii], e, xr, self.gate_src, self.gate_win, tauE1, isc, thc, float(TAU_EFFCA))
                if v.size:
                    oq[i] = np.median(v)
            if self.gate_src == 2:
                r["shaft_cai"] = None
        self.bap_sh = np.array([np.median(ev[k]) if ev.get(k) else np.nan for k in keys]) if self.gate_src == 2 \
            else np.full(N, np.nan)                                          # per physical synapse, broadcast to rows
        self.open_q = oq
        self.X = cuda.to_device(X) if self.gate_src == 2 else cuda.to_device(np.zeros(1, np.float32))
        del X
        self.tex = bool(self.t_exact)
        if not self.tex:
            self.fa = cuda.to_device(np.zeros(1, np.float64))
        self._kern8 = _make_kernel_v8(weighted=self.v5_mode == 2, src=self.gate_src, win=self.gate_win, tex=self.tex)
        _INST.append(self)
        msg = (f"v8 gate_src {self.gate_src} win {self.gate_win} ms, gate_theta {self.gate_theta}: Vg-opening gate quantity "
               f"per synapse q10/50/90 {np.nanquantile(oq, [.1, .5, .9]).round(4) if np.isfinite(oq).any() else '-'} "
               f"{UNIT[self.gate_src]} ({int(np.isfinite(oq).sum())} of {N} synapses open at theta_V {thc:.4g})")
        if self.gate_src == 2:
            rr = np.concatenate(rests) * 1e6
            msg += (f"; shaft rest q10/50/90 {np.quantile(rr, [.1, .5, .9]).round(1)} nM, {nclip} samples clipped at 65.5 uM; "
                    f"uint16 shaft {self.X.nbytes / 1e9:.2f} GB on GPU")
        print(msg, flush=True)

    @staticmethod
    def _kappa(P0):
        if _STATE["kappa"] is None:
            q = np.concatenate([m.open_q for m in _INST]); q = q[np.isfinite(q)]
            thc = float(P0.get("gate_conv_thV", 3.5751446275588403))
            Q50 = float(np.median(q)) if q.size else float("nan")
            kap = float(P0.get("gate_kappa", 0.0)) or Q50 / thc
            _STATE.update(kappa=kap, Q50=Q50, theta_conv=thc, n_open=int(q.size))
            print(f"v8 kappa {kap:.6g} {UNIT.get(int(P0.get('gate_src', 0)), '')} per theta_V unit (Q50 {Q50:.6g} at theta_V "
                  f"{thc:.4g}, {q.size} synapses pooled over {len(_INST)} models)", flush=True)
            if int(P0.get("gate_src", 0)) == 2:
                _calib()
        return _STATE["kappa"]

    def rho_dpre_syn(self, td, tp, Ps):
        if self.gate_src == 0:
            return super().rho_dpre_syn(td, tp, Ps)
        kap = self._kappa(self.P0)
        P, N = td.shape
        out_rho = cuda.device_array((P, N), np.float64); out_d = cuda.device_array((P, N, 3), np.float64)
        bx, by = self.block
        g = self.d
        self._kern8[(-(-P // bx), -(-N // by)), (bx, by)](
            g["E"], g["S"], g["Hg"], g["eoff"], g["slen"], g["hoff"], g["aptr"], g["astep"], g["acnt"], g["rho0"],
            g["order"], cuda.to_device(np.ascontiguousarray(td, np.float64)),
            cuda.to_device(np.ascontiguousarray(tp, np.float64)), cuda.to_device(self._params5(Ps)), out_rho, out_d,
            *self.u6, self.vq, self.X, float(kap), float(self.gate_theta), self.fa)
        return out_rho.copy_to_host(), out_d.copy_to_host()


def _calib():
    """GATE_ALT checklist: single-bAP shaft Delta[Ca] by cexp location, proximal median, Spearman vs vdcc_q_post."""
    try:
        import fit_v6
        from scipy.stats import spearmanr
        cdir = os.environ.get("CEXP_DIR", "/scratch/dhuruva/split1/cexp")
        allb, allq = [], []
        for G, name in fit_v6._M:
            if G not in _INST:
                continue
            df = fit_v6.join_cexp(G, name, cdir)
            u = np.unique([f"{p}:{s}" for p, s in G._keys], return_index=True)[1]     # one row per physical synapse
            b = G.bap_sh[u]; loc = df["loc"].astype(str).to_numpy()[u]; dist = df["dist_um"].to_numpy(float)[u]
            qp = df["vdcc_q_post"].to_numpy(float)[u]
            q = lambda v: np.nanquantile(v, [.1, .5, .9]).round(4) if np.isfinite(v).any() else "-"
            parts = [f"{L} n{int(np.isfinite(b[loc == L]).sum())} {q(b[loc == L])}" for L in sorted(set(loc))]
            prox = b[dist <= 100.0]
            ok = np.isfinite(b) & np.isfinite(qp)
            rs = spearmanr(b[ok], qp[ok])[0] if ok.sum() > 2 else float("nan")
            bu = G.burst_sh
            print(f"v8 calib {name}: single-bAP peak shaft dCa uM q10/50/90 {q(b)} ({int(np.isfinite(b).sum())} of {len(b)} "
                  f"synapses) | by loc {'; '.join(parts)} | proximal <=100 um median "
                  f"{np.nanmedian(prox) if np.isfinite(prox).any() else float('nan'):.4f} (n {int(np.isfinite(prox).sum())}) | "
                  f"Spearman vs vdcc_q_post {rs:.3f} (n {int(ok.sum())}) | 3AP 200 Hz record peak q10/50/90 {q(bu)} "
                  f"(n {int(np.isfinite(bu).sum())})", flush=True)
            allb.append(b[ok]); allq.append(qp[ok])
        if allb:
            b = np.concatenate(allb); qp = np.concatenate(allq)
            print(f"v8 calib pooled: Spearman single-bAP shaft dCa vs vdcc_q_post {spearmanr(b, qp)[0]:.3f} (n {len(b)})",
                  flush=True)
    except Exception as ex:             # diagnostics only, never stops a fit
        print(f"v8 calib failed: {type(ex).__name__}: {ex}", flush=True)


if __name__ == "__main__":
    _s = json.loads(sys.argv[sys.argv.index("--set") + 1]) if "--set" in sys.argv else {}
    if int(_s.get("ecb_ref", 0)) == 2:                  # as gpu_v7_rho.__main__
        os.environ.setdefault("CEXP_DIR", "/scratch/dhuruva/split1/cexp")
        _e = os.environ.get("SCALE_E") or '{"vdcc_q_post": 100000.0}'
        assert json.loads(_e) == {"vdcc_q_post": 1e5}, f"ecb_ref 2 needs SCALE_E vdcc_q_post x 1e5, got {_e}"
        os.environ["SCALE_E"] = _e
    elif int(_s.get("ecb_ref", 0)) == 1:
        assert not os.environ.get("SCALE_E"), "ecb_ref 1 is not combined with SCALE_E"
    src = int(_s.get("gate_src", 0))
    import fit_v6
    fit_v6.gpu_v6_rho.GPUModelV6 = GPUModelV8
    if src == 2:                                        # keep the shaft trace in the records (popped after the build)
        _BV = fit_v6.BatchV2
        fit_v6.BatchV2 = lambda *a, signals=("vdcc",), **k: _BV(*a, signals=tuple(signals) + ("shaft_cai",), **k)
    save = sys.argv[sys.argv.index("--save") + 1]
    code = 0
    try:
        fit_v6.main()
    except SystemExit as ex:
        if src == 0 or str(ex) != "REPRO DIFF":
            raise
        print("REPRO DIFF expected (gate_src != 0: a different rule than the seed)", flush=True)
    if src and os.path.isfile(save + ".json"):
        out = json.load(open(save + ".json"))
        thV = float(out["v5"]["theta_V"]); gt = float(_s.get("gate_theta", 0.0))
        out["v8"] = dict(gate_src=src, t_exact=int(_s.get("t_exact", 0)), gate_win=float(_s.get("gate_win", 100.0)) if src in (1, 2) else 0.0, gate_theta=gt,
                         kappa=_STATE.get("kappa"), Q50=_STATE.get("Q50"), theta_conv=_STATE.get("theta_conv"),
                         n_open=_STATE.get("n_open"), theta_G=gt if gt > 0 else thV * float(_STATE.get("kappa") or np.nan),
                         unit=UNIT[src])
        json.dump(out, open(save + ".json", "w"), indent=1)
        print(f"v8: theta_G {out['v8']['theta_G']:.6g} {UNIT[src]} -> {save}.json", flush=True)
    sys.exit(code)
